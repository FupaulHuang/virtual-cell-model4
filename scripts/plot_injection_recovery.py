#!/usr/bin/env python3
"""Descriptive plots of the two immutable oracle reports; no new selection gate."""
from pathlib import Path
import hashlib
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def main():
    dev_path = ROOT / "reports/injection_recovery_development.json"
    conf_path = ROOT / "reports/injection_recovery_confirmation.json"
    dev = json.loads(dev_path.read_text())
    conf = json.loads(conf_path.read_text())
    if conf["development_report"]["sha256"] != hashlib.sha256(dev_path.read_bytes()).hexdigest():
        raise ValueError("Development report changed after confirmation")
    if dev["analysis_freeze_sha256"] != conf["analysis_freeze_sha256"]:
        raise ValueError("Analysis changed between stages")
    rows = []
    for report in (dev, conf):
        rows += [dict(r, method="reference") for r in report["reference_results"]]
        rows += report["candidate_results"]
    if any(r.get("status", "complete") != "complete" for r in rows):
        raise ValueError("Incomplete arm: do not draw a selected successful subset")
    families = [("composition", "Real-cell composition"),
                ("dirichlet_multinomial", "Dirichlet–multinomial")]
    methods = [("reference", "Independent reference", "#374151"),
               ("current", "Current pooled emitter", "#c84f4f"),
               ("preserving_uniform", "Preserving emitter", "#1b8a77")]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8.5), constrained_layout=True)
    for col, (family, title) in enumerate(families):
        ax = axes[0, col]
        largest_null_count = 5
        for j, (method, label, color) in enumerate(methods):
            subset = [r for r in rows if r["family"] == family and r["method"] == method and r["is_null"]]
            if len(subset) != 4:
                raise ValueError("NULL panel incomplete")
            for offset, r in zip(np.linspace(-.12, .12, 4), subset):
                marker = "D" if r["partition"] == "development" else "o"
                largest_null_count = max(largest_null_count, r["metrics"]["significant_calls"])
                ax.scatter(j + offset, r["metrics"]["significant_calls"], color=color,
                           marker=marker, s=45, zorder=3)
        ax.axhline(5, color="#737373", linestyle="--", linewidth=1)
        ax.text(.02, .95, "Dashed line: reference NULL limit", transform=ax.transAxes,
                va="top", fontsize=9, color="#555555")
        ax.set_yscale("symlog", linthresh=1)
        ax.set_ylim(bottom=-.15, top=largest_null_count * 2.5)
        ax.set_xticks(range(3), ["Reference", "Current", "Preserving"])
        ax.set_title(title + ": zero-effect calls", weight="bold")
        ax.set_ylabel("Significant genes (q < 0.05)")
        ax.grid(axis="y", alpha=.18)
        ax = axes[1, col]
        for j, (method, label, color) in enumerate(methods):
            centers, means = [], []
            for x, regime in enumerate(("WEAK", "MODERATE", "STRONG")):
                subset = [r for r in rows if r["family"] == family and r["method"] == method and r["regime"] == regime]
                if len(subset) != 4:
                    raise ValueError("Signal panel incomplete")
                center = x + (j - 1) * .22
                vals = []
                for offset, r in zip(np.linspace(-.055, .055, 4), subset):
                    value = r["metrics"]["signed_recovery"]
                    if value is None:
                        continue
                    vals.append(value)
                    marker = "D" if r["partition"] == "development" else "o"
                    ax.scatter(center + offset, value, color=color, marker=marker,
                               s=30, alpha=.65, zorder=3)
                if vals:
                    centers.append(center)
                    means.append(float(np.mean(vals)))
            ax.plot(centers, means, color=color, linewidth=1.7, label=label)
        ax.set_xticks(range(3), ["Weak", "Moderate", "Strong"])
        ax.set_ylim(-.02, 1.02)
        ax.set_ylabel("Correct-sign significant calls / true-effect genes")
        ax.set_title(title + ": signed recovery", weight="bold")
        ax.grid(axis="y", alpha=.18)
    axes[1, 0].legend(loc="upper left", fontsize=8.5)
    fig.suptitle("Known-moment diagnostic · 400 cells per group\n"
                 "Dots: fixed template/seed cases; diamonds: development. Lines: descriptive means.\n"
                 "True-effect denominator: |log₂FC| ≥ 0.1 and comparison CPM > 5. Undefined ratios omitted. No official score.",
                 fontsize=12)
    prefix = ROOT / "reports/injection_recovery_20261021"
    fig.savefig(prefix.with_suffix(".png"), dpi=180)
    fig.savefig(prefix.with_suffix(".svg"))
    plt.close(fig)
    print(prefix)


if __name__ == "__main__":
    main()
