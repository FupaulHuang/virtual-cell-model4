#!/usr/bin/env python3
"""Small raw-source response adapter; no H1 response input or label is permitted.

Two independent ridge heads predict explicit log2(CPM+1) differences and natural
log1p(50000*p) bulk differences. Source-context LOSO is used for regularization
selection, not claimed as an unbiased final performance estimate. Every outer
held-out context is absent from BOTH training labels and training features.

The returned effect/coverage tuples have the Atlas interface. Learned effects
already predict destination response magnitude: use amplitude=1 downstream,
not the incumbent's additional .6/.3 attenuation. Confidence is a diagnostic
probability of source-mean sign agreement, not a DE probability or a hard gate.
CD4 is intentionally omitted in v1; no correlated raw/DE duplicate is accepted.
"""
from __future__ import annotations

import os
for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_key, '1')
import argparse
from dataclasses import dataclass, field
import hashlib
import json
from pathlib import Path
import re

import numpy as np

FAMILIES = ('K562', 'HCT116', 'HEK293T')
SPACES = ('log2fc', 'bulk_delta')
UNITS = {
    'log2fc': 'log2(target arithmetic CPM + 1) - log2(matched control arithmetic CPM + 1)',
    'bulk_delta': 'ln(1 + 50000 * target native-panel count probability) - ln(1 + 50000 * matched control probability)',
}
DEFAULT_LAMBDAS = (0.001, 0.03, 0.3, 3.0)
MODIFIERS = ('constant', 'availability_fraction', 'mean_log_cell_count', 'sign_agreement',
             'log_effect_disagreement', 'destination_control_expression', 'source_control_gap',
             'arithmetic_log2_se', 'arithmetic_se_available')
FEATURE_NAMES = tuple(f'{effect}:{modifier}' for effect in ('specific', 'shared') for modifier in MODIFIERS)
CONFIDENCE_NAMES = ('intercept', 'log_abs_source_mean', *MODIFIERS[1:])


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def _family(source, path=''):
    # Require an allowed raw-source identity, even when the caller supplies a path alias.
    text = f'{source} {path}'
    if re.search(r'(?i)(?<![a-z0-9])h1(?![a-z0-9])', text):
        raise ValueError('H1 response banks are forbidden in adapter labels and inputs')
    found = [f for f in FAMILIES if re.search(rf'(?i)(?<![a-z0-9]){f}(?![a-z0-9])', str(source))]
    if len(found) != 1:
        raise ValueError('Only one identified raw K562/HCT116/HEK293T bank is allowed')
    return found[0]


@dataclass
class RawBank:
    source: str
    targets: np.ndarray
    genes: np.ndarray
    n_cells: np.ndarray
    counts: np.ndarray
    mean_cpm: np.ndarray
    control_probability: np.ndarray
    control_mean_cpm: np.ndarray
    measured: np.ndarray
    path: str = ''
    digest: str = ''
    mean_variance: np.ndarray | None = None
    control_mean_variance: np.ndarray | None = None
    _cache: dict = field(default_factory=dict, repr=False)

    def __post_init__(self):
        self.family = _family(self.source, self.path)
        self.targets = np.asarray(self.targets).astype(str)
        self.genes = np.asarray(self.genes).astype(str)
        original_n = np.asarray(self.n_cells)
        if not np.isfinite(original_n).all() or (original_n != np.floor(original_n)).any():
            raise ValueError('Source cell counts must be finite integers')
        self.n_cells = original_n.astype(np.int64)
        self.measured = np.asarray(self.measured, dtype=bool)
        shape = (len(self.targets), len(self.genes))
        if len(set(self.targets)) != shape[0] or len(set(self.genes)) != shape[1]:
            raise ValueError('Raw source target and gene axes must be unique')
        if self.n_cells.shape != (shape[0],) or (self.n_cells < 0).any() or self.measured.shape != (shape[1],):
            raise ValueError('Invalid source counts/measurement mask')
        for key in ('counts', 'mean_cpm', 'control_probability', 'control_mean_cpm'):
            value = np.asarray(getattr(self, key))
            if value.shape == (shape[1],) and key.startswith('control_'):
                value = np.broadcast_to(value, shape)
            if value.shape != shape or not np.isfinite(value).all() or (value < 0).any():
                raise ValueError(f'Invalid raw source array: {key}')
            setattr(self, key, value)
        if (self.counts != np.floor(self.counts)).any():
            raise ValueError('Raw counts must be actual nonnegative integers')
        if ((self.n_cells > 0) & (self.counts.sum(axis=1) <= 0)).any():
            raise ValueError('Measured source groups require positive raw count mass')
        if not np.allclose(self.control_probability.sum(axis=1), 1, rtol=1e-5, atol=1e-7):
            raise ValueError('Matched control probabilities must sum to one on the saved native axis')
        if (self.mean_variance is None) != (self.control_mean_variance is None):
            raise ValueError('Both arithmetic mean-variance arrays are required together')
        for key in ('mean_variance', 'control_mean_variance'):
            value = getattr(self, key)
            if value is not None:
                value = np.asarray(value)
                if value.shape != shape or not np.isfinite(value).all() or (value < 0).any():
                    raise ValueError(f'Invalid observed uncertainty: {key}')
                setattr(self, key, value)
        self._rows = {t: i for i, t in enumerate(self.targets)}
        self._cols = {g: i for i, g in enumerate(self.genes)}

    def effects(self, space, prior_counts=100000.0, minimum_cells=20):
        """Native-axis effects and full-source common component; own gene omitted in centering."""
        if space not in SPACES or not np.isfinite(prior_counts) or prior_counts < 0:
            raise ValueError('Invalid response space/prior')
        key = (space, float(prior_counts), int(minimum_cells))
        if key in self._cache:
            return self._cache[key]
        total = self.counts.sum(axis=1, dtype=np.float64)
        fraction = np.divide(total, total + prior_counts, out=np.zeros_like(total), where=total + prior_counts > 0)
        if space == 'log2fc':
            treated = fraction[:, None] * self.mean_cpm + (1 - fraction[:, None]) * self.control_mean_cpm
            effects = np.log2(treated + 1) - np.log2(self.control_mean_cpm + 1)
        else:
            probability = (self.counts.astype(np.float64) + prior_counts * self.control_probability)
            probability = np.divide(probability, (total + prior_counts)[:, None],
                                    out=np.zeros_like(probability), where=(total + prior_counts)[:, None] > 0)
            effects = np.log1p(50000 * probability) - np.log1p(50000 * self.control_probability)
        eligible = self.n_cells >= minimum_cells
        valid = eligible[:, None] & self.measured[None, :]
        centering = valid.copy()
        for i, t in enumerate(self.targets):
            if t in self._cols:
                centering[i, self._cols[t]] = False
        common = np.divide(np.where(centering, effects, 0).sum(axis=0), centering.sum(axis=0),
                           out=np.zeros(len(self.genes)), where=centering.sum(axis=0) > 0)
        result = (np.where(valid, effects, 0).astype(np.float32), common, eligible)
        self._cache[key] = result
        return result

    def pair_data(self, targets, genes, space, prior_counts, minimum_cells):
        ri = np.fromiter((self._rows.get(str(t), -1) for t in targets), dtype=int, count=len(targets))
        ci = np.fromiter((self._cols.get(str(g), -1) for g in genes), dtype=int, count=len(genes))
        valid = (ri >= 0) & (ci >= 0)
        r, c = np.maximum(ri, 0), np.maximum(ci, 0)
        valid &= (self.n_cells[r] >= minimum_cells) & self.measured[c]
        effect, common, _ = self.effects(space, prior_counts, minimum_cells)
        se = np.zeros(len(r)); known = np.zeros(len(r), dtype=bool)
        if self.mean_variance is not None:
            # These bank fields are variance of the arithmetic mean, not cell variance.
            se = np.sqrt(self.mean_variance[r, c] / (self.mean_cpm[r, c] + 1)**2
                         + self.control_mean_variance[r, c] / (self.control_mean_cpm[r, c] + 1)**2) / np.log(2)
            known[:] = True
        return {'effect': np.where(valid, effect[r, c], 0), 'common': np.where(valid, common[c], 0),
                'valid': valid, 'n': np.where(valid, self.n_cells[r], 0),
                'control_log2_cpm': np.where(valid, np.log2(self.control_mean_cpm[r, c] + 1), 0),
                'arithmetic_se': np.where(valid, se, 0), 'se_known': valid & known}

    def provenance(self):
        return {'family': self.family, 'source': self.source, 'path': self.path, 'sha256': self.digest,
                'native_genes': len(self.genes), 'targets': len(self.targets),
                'eligible_targets_20_cells': int((self.n_cells >= 20).sum()),
                'arithmetic_mean_uncertainty_available': self.mean_variance is not None}


def load_raw_bank(path):
    path = Path(path).resolve()
    with np.load(path, allow_pickle=False) as d:
        source = str(d['source'])
        _family(source, path)
        # Required count schema excludes DE-only and pseudo-count representations.
        return RawBank(source, d['targets'], d['genes'], d['n_cells'], d['target_count_sums'],
                       d['target_mean_cpm'], d['matched_control_probability'], d['matched_control_mean_cpm'],
                       d['measured_genes'], str(path), sha256(path),
                       d['target_mean_cpm_variance'] if 'target_mean_cpm_variance' in d else None,
                       d['matched_control_mean_cpm_variance'] if 'matched_control_mean_cpm_variance' in d else None)


def validate_families(banks, minimum=1):
    names = [_family(b.source, b.path) for b in banks]
    if len(names) < minimum or len(set(names)) != len(names):
        raise ValueError('Use distinct raw source-context families; duplicate representations are forbidden')
    return names


def build_features(banks, targets, genes, destination_control_cpm, space,
                   prior_counts=100000.0, minimum_cells=20):
    """Paired target/gene features; destination controls are allowed, destination responses are not."""
    validate_families(banks)
    targets, genes = np.asarray(targets).astype(str), np.asarray(genes).astype(str)
    control = np.asarray(destination_control_cpm, dtype=float)
    if targets.ndim != 1 or genes.shape != targets.shape or control.shape != targets.shape:
        raise ValueError('Features require equally sized target/gene/control vectors')
    if not np.isfinite(control).all() or (control < 0).any():
        raise ValueError('Destination controls must be finite nonnegative arithmetic CPM')
    data = [b.pair_data(targets, genes, space, prior_counts, minimum_cells) for b in banks]
    valid = np.stack([d['valid'] for d in data]); count = valid.sum(axis=0)
    def average(key):
        return np.divide(np.stack([d[key] for d in data]).sum(axis=0), count,
                         out=np.zeros(len(targets)), where=count > 0)
    effects = np.stack([d['effect'] for d in data])
    mean, shared = average('effect'), average('common')
    specific = mean - shared
    absolute = np.divide(np.abs(effects).sum(axis=0), count, out=np.zeros(len(targets)), where=count > 0)
    agreement = np.divide(np.abs(mean), absolute, out=np.zeros(len(targets)), where=absolute > 1e-12)
    disagreement = np.sqrt(np.divide(np.where(valid, (effects - mean)**2, 0).sum(axis=0), count,
                                    out=np.zeros(len(targets)), where=count > 0))
    logn = np.divide(np.stack([np.log1p(np.minimum(d['n'], 1000)) for d in data]).sum(axis=0),
                     count, out=np.zeros(len(targets)), where=count > 0) / np.log(1001)
    destination = np.log2(control + 1) / 16
    gap = (average('control_log2_cpm') / 16) - destination
    known = np.stack([d['se_known'] for d in data]).sum(axis=0)
    uncertainty = np.divide(np.stack([d['arithmetic_se'] for d in data]).sum(axis=0), known,
                            out=np.zeros(len(targets)), where=known > 0)
    modifiers = np.column_stack([np.ones(len(targets)), count / len(banks), logn, agreement,
                                np.log1p(disagreement), destination, gap, np.log1p(uncertainty), known / len(banks)])
    x = np.column_stack([specific[:, None] * modifiers, shared[:, None] * modifiers])
    confidence = np.column_stack([np.ones(len(targets)), np.log1p(np.abs(mean)), modifiers[:, 1:]])
    x[count == 0] = 0
    return x, confidence, count, mean


def _sample_pairs(destination, inputs, max_pairs, minimum_cells, seed):
    """Response-blind, equal-target quotas, uniformly sampled genes; own coordinates omitted."""
    if max_pairs < 1:
        raise ValueError('max_pairs must be positive')
    validate_families(inputs)
    if destination.family in {b.family for b in inputs}:
        raise ValueError('Destination responses cannot be their own input features')
    input_columns = [np.isin(destination.genes, b.genes[b.measured]) for b in inputs]
    eligible = []
    for t, n in zip(destination.targets, destination.n_cells):
        if n < minimum_cells:
            continue
        measured = np.zeros(len(destination.genes), dtype=bool)
        for b, columns in zip(inputs, input_columns):
            row = b._rows.get(t)
            if row is not None and b.n_cells[row] >= minimum_cells:
                measured |= columns
        columns = np.flatnonzero(destination.measured & (destination.genes != t) & measured)
        if len(columns):
            eligible.append((str(t), columns))
    if not eligible:
        raise ValueError('No eligible cross-context target/gene pairs')
    if max_pairs < len(eligible):
        raise ValueError('Sampling budget must include every eligible target')
    rng = np.random.default_rng(seed)
    labels, genes, weights = [], [], []
    quotient, remainder = divmod(max_pairs, len(eligible))
    for i, (target, columns) in enumerate(eligible):
        size = min(len(columns), quotient + (i < remainder))
        chosen = np.sort(rng.choice(columns, size=size, replace=False))
        labels.extend([target] * size); genes.extend(destination.genes[chosen])
        weights.extend([1 / (len(eligible) * size)] * size)
    return np.asarray(labels), np.asarray(genes), np.asarray(weights)


def examples(destination, inputs, space, *, max_pairs=30000, seed=20261021,
             prior_counts=100000.0, minimum_cells=20):
    target, genes, weights = _sample_pairs(destination, inputs, max_pairs, minimum_cells, seed)
    ri = np.asarray([destination._rows[t] for t in target]); ci = np.asarray([destination._cols[g] for g in genes])
    x, confidence, count, source_mean = build_features(inputs, target, genes,
        destination.control_mean_cpm[ri, ci], space, prior_counts, minimum_cells)
    # Labels use observed arithmetic means / raw count sums without the source input prior.
    y = destination.effects(space, prior_counts=0, minimum_cells=minimum_cells)[0][ri, ci].astype(float)
    if not (count > 0).all():
        raise ValueError('Sampled training pair has no measured input source')
    return {'x': x, 'confidence_x': confidence, 'y': y, 'weight': weights,
            'confidence_y': ((source_mean * y) > 0).astype(float), 'source_mean': source_mean,
            'targets': target, 'genes': genes, 'destination': destination.family,
            'input_families': [b.family for b in inputs]}


def fold_training_examples(banks, heldout_family, space, **kwargs):
    validate_families(banks, minimum=3)
    if heldout_family not in {b.family for b in banks}:
        raise ValueError('Unknown held-out source context')
    remaining = [b for b in banks if b.family != heldout_family]
    # Crucial: heldout is removed before constructing every training feature AND label task.
    return [examples(destination, [b for b in remaining if b.family != destination.family], space, **kwargs)
            for destination in remaining]


def fit_ridge(tasks, regularization, *, confidence=False):
    if not np.isfinite(regularization) or regularization <= 0:
        raise ValueError('Strictly positive ridge regularization required')
    key, label = ('confidence_x', 'confidence_y') if confidence else ('x', 'y')
    x = np.concatenate([t[key] for t in tasks]); y = np.concatenate([t[label] for t in tasks])
    weight = np.concatenate([t['weight'] / len(tasks) for t in tasks])
    if not np.isfinite(x).all() or not np.isfinite(y).all() or not np.isclose(weight.sum(), 1):
        raise ValueError('Invalid weighted training examples')
    scale = np.sqrt(np.sum(weight[:, None] * x*x, axis=0))
    scale = np.maximum(scale, 1e-8)
    z = x / scale
    # Uncentered features and no intercept preserve zero-effect identity of the response head.
    penalty = np.eye(x.shape[1]) * regularization
    if confidence:
        penalty[0, 0] = regularization * .01
    coefficient = np.linalg.solve((z * weight[:, None]).T @ z + penalty, (z * weight[:, None]).T @ y)
    return {'coefficient': coefficient.tolist(), 'scale': scale.tolist(),
            'regularization': float(regularization), 'confidence': confidence}


def apply_head(head, x):
    value = (np.asarray(x) / np.asarray(head['scale'])) @ np.asarray(head['coefficient'])
    return np.clip(value, 0, 1) if head['confidence'] else value


def evaluate(task, head, confidence_head):
    predicted = apply_head(head, task['x'])
    probability = apply_head(confidence_head, task['confidence_x'])
    y, w = task['y'], task['weight']
    return {'pairs': len(y), 'targets': len(set(task['targets'])),
            'mse': float(np.sum(w * (predicted-y)**2)),
            'mae': float(np.sum(w * np.abs(predicted-y))),
            'zero_response_mse': float(np.sum(w * y*y)),
            'unweighted_source_mean_mse': float(np.sum(w * (task['source_mean']-y)**2)),
            'direction_accuracy': float(np.sum(w * (np.sign(predicted) == np.sign(y)))),
            'sign_confidence_brier': float(np.sum(w * (probability-task['confidence_y'])**2))}


def fit_loso(banks, *, lambdas=DEFAULT_LAMBDAS, max_pairs=30000, seed=20261021,
             prior_counts=100000.0, minimum_cells=20):
    families = validate_families(banks, minimum=3)
    if set(families) != set(FAMILIES):
        raise ValueError('v1 requires exactly the three registered non-H1 raw contexts')
    banks = sorted(banks, key=lambda b: FAMILIES.index(b.family))
    lambdas = tuple(sorted(set(float(x) for x in lambdas)))
    if not lambdas or any(not np.isfinite(x) or x <= 0 for x in lambdas):
        raise ValueError('A nonempty positive fixed ridge grid is required')
    settings = dict(max_pairs=max_pairs, seed=seed, prior_counts=prior_counts, minimum_cells=minimum_cells)
    heads, report = {}, {}
    for space in SPACES:
        evaluations = {str(lam): [] for lam in lambdas}
        for heldout in banks:
            train = fold_training_examples(banks, heldout.family, space, **settings)
            test = examples(heldout, [b for b in banks if b.family != heldout.family], space, **settings)
            for lam in lambdas:
                head = fit_ridge(train, lam)
                confidence = fit_ridge(train, lam, confidence=True)
                evaluations[str(lam)].append({'heldout_context': heldout.family,
                    'training_label_contexts': [t['destination'] for t in train],
                    'training_input_contexts': sorted({b for t in train for b in t['input_families']}),
                    'test_input_contexts': test['input_families'],
                    **evaluate(test, head, confidence)})
        # Equal weight for each held-out context, then prefer stronger ridge on exact ties.
        selected = min(lambdas, key=lambda lam: (np.mean([x['mse'] for x in evaluations[str(lam)]]), -lam))
        final = [examples(destination, [b for b in banks if b.family != destination.family], space, **settings)
                 for destination in banks]
        heads[space] = {'response': fit_ridge(final, selected),
                        'sign_confidence': fit_ridge(final, selected, confidence=True)}
        report[space] = {'selected_regularization': selected,
            'selection_criterion': 'equal-context mean held-out response MSE; larger ridge for exact ties',
            'grid': evaluations, 'selected_mean_context_mse': float(np.mean([x['mse'] for x in evaluations[str(selected)]]))}
    model = {'schema': 'raw_response_adapter_v1', 'units': UNITS, 'heads': heads,
        'response_features': list(FEATURE_NAMES), 'confidence_features': list(CONFIDENCE_NAMES),
        'settings': {**settings, 'ridge_grid': list(lambdas)},
        'source_banks': [b.provenance() for b in banks], 'forbidden_response_contexts': ['H1'],
        'implementation_sha256': sha256(__file__), 'downstream_amplitude': 1.0,
        'confidence_semantics': 'diagnostic probability that pooled source-mean sign agrees with destination response; not DE significance; not multiplied into response',
        'limitations': [
            'Only three source contexts; LOSO selects regularization and is not an unbiased final estimate.',
            'LOSO training has one input context per label task; held-out tests have two, final inference may have three.',
            'Native-panel count probabilities differ across assays; arithmetic CPM retains each bank denominator.',
            'Observed mean uncertainty exists only for K562; unavailable uncertainty is masked, not inferred as zero noise.',
            'Shared components use all eligible non-H1 source targets excluding own-target coordinates.',
            'No target/gene identity features, H1 responses, CD4, DE-only duplicate banks or official scores are used.',
            'Context folds reuse perturbation identities; this is context transfer, not unseen-target validation.',
        ]}
    return model, {'loso': report, 'settings': model['settings'], 'source_banks': model['source_banks'],
                   'limitations': model['limitations'], 'h1_responses_used': False}


def predict_effects(model, banks, targets, genes, destination_control_mean_cpm, *, chunk_targets=16):
    """Return ({space: (effect, coverage)}, {space: sign_confidence}) on exact requested axes.

    No destination response object is accepted. The caller supplies only measured
    arithmetic control CPM. Destination bulk controls are used by desired_mean.
    """
    validate_families(banks)
    if model.get('schema') != 'raw_response_adapter_v1' or chunk_targets < 1:
        raise ValueError('Invalid adapter schema/chunk size')
    if (model['implementation_sha256'] != sha256(__file__) or model['units'] != UNITS
            or model['response_features'] != list(FEATURE_NAMES)
            or model['confidence_features'] != list(CONFIDENCE_NAMES)):
        raise ValueError('Adapter implementation, features or response units changed after fitting')
    expected = {p['family']: p for p in model['source_banks']}
    if set(expected) != {b.family for b in banks}:
        raise ValueError('Inference source families must match the trained adapter')
    for b in banks:
        if b.digest != expected[b.family]['sha256'] or b.source != expected[b.family]['source']:
            raise ValueError('Inference source bank identity differs from training')
    targets, genes = np.asarray(targets).astype(str), np.asarray(genes).astype(str)
    if targets.ndim != 1 or genes.ndim != 1 or len(set(targets)) != len(targets) or len(set(genes)) != len(genes):
        raise ValueError('Unique one-dimensional target/gene axes are required')
    control = np.asarray(destination_control_mean_cpm, dtype=float)
    if control.shape == (len(genes),):
        control = np.broadcast_to(control, (len(targets), len(genes)))
    if control.shape != (len(targets), len(genes)) or not np.isfinite(control).all() or (control < 0).any():
        raise ValueError('Destination arithmetic control CPM axes/values are invalid')
    settings = model['settings']; outputs, confidence_outputs = {}, {}
    for space in SPACES:
        effect = np.zeros(control.shape, dtype=np.float32)
        coverage = np.zeros(control.shape, dtype=np.float32)
        confidence = np.zeros(control.shape, dtype=np.float32)
        for start in range(0, len(targets), chunk_targets):
            end = min(start + chunk_targets, len(targets))
            pair_targets = np.repeat(targets[start:end], len(genes))
            pair_genes = np.tile(genes, end-start)
            x, cx, n, _ = build_features(banks, pair_targets, pair_genes, control[start:end].ravel(), space,
                                        settings['prior_counts'], settings['minimum_cells'])
            value = apply_head(model['heads'][space]['response'], x)
            certainty = apply_head(model['heads'][space]['sign_confidence'], cx)
            value[n == 0] = 0; certainty[n == 0] = 0
            effect[start:end] = value.reshape(end-start, len(genes))
            coverage[start:end] = n.reshape(end-start, len(genes))
            confidence[start:end] = certainty.reshape(end-start, len(genes))
        outputs[space] = effect, coverage
        confidence_outputs[space] = confidence
    return outputs, confidence_outputs


def effects_to_profiles(effects, control_mean_cpm, control_bulk_probability, *, clip=3.0):
    """Invert exact response units, then normalize on the destination gene axis.

    The arithmetic inverse includes its +1 pseudocount. This helper avoids an
    accidental extra incumbent amplitude or approximate CPM-ratio conversion.
    """
    if not np.isfinite(clip) or clip <= 0:
        raise ValueError('Positive finite response clipping required')
    result = []
    for space, control in (('log2fc', control_mean_cpm), ('bulk_delta', control_bulk_probability)):
        effect = np.asarray(effects[space][0], dtype=np.float64)
        control = np.asarray(control, dtype=np.float64)
        if effect.ndim != 2 or control.shape not in (effect.shape, (effect.shape[1],)):
            raise ValueError('Destination effect/control axes differ')
        if not np.isfinite(effect).all() or not np.isfinite(control).all() or (control < 0).any():
            raise ValueError('Invalid destination effects/controls')
        if space == 'log2fc':
            desired = np.maximum((control + 1) * np.exp2(np.clip(effect, -clip, clip)) - 1, 0)
        else:
            if not np.allclose(control.sum(axis=-1), 1, atol=1e-7, rtol=1e-5):
                raise ValueError('Bulk controls must sum to one')
            desired = np.expm1(np.maximum(np.log1p(50000 * control) + np.clip(effect, -clip, clip), 0))
        total = desired.sum(axis=1, keepdims=True)
        if (total <= 0).any() or not np.isfinite(desired).all():
            raise ValueError('Inverted response has no finite positive expression mass')
        result.append(desired / total)
    return tuple(result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bank', action='append', required=True, help='Raw non-H1 NPZ; exactly one per allowed context')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--report', type=Path, required=True)
    p.add_argument('--pairs-per-context', type=int, default=30000)
    p.add_argument('--seed', type=int, default=20261021)
    args = p.parse_args()
    if args.output.exists() or args.report.exists() or args.output.resolve() == args.report.resolve():
        raise FileExistsError('Use distinct new adapter and report paths')
    model, report = fit_loso([load_raw_bank(x) for x in args.bank], max_pairs=args.pairs_per_context, seed=args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True); args.report.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(model, indent=2, allow_nan=False) + '\n')
    report.update(model_path=str(args.output.resolve()), model_sha256=sha256(args.output))
    args.report.write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({s: report['loso'][s]['selected_regularization'] for s in SPACES}), flush=True)


if __name__ == '__main__':
    main()
