"""Control-only cube-style flight followed by fixed-size pivotal landing.

Independent implementation of the balancing-martingale construction; no
external implementation was copied. In exact arithmetic, every coordinate
retains inclusion probability n/N and sample size is n. Flight balances the
retained auxiliary span; pivotal landing relaxes non-intercept constraints.
Final balance is approximate and is audited, not optimized after sampling.

References: Deville & Tille (2004), doi:10.1093/biomet/91.4.893;
Chauvet (2012), https://arxiv.org/abs/1211.5442.

The caller must provide a verified control-only population. This module reads
no datasets, perturbation responses, labels, scores or file paths. It returns
indices of actual cells and never modifies their counts.
"""
from __future__ import annotations
import os
for _key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(_key,'1')
from dataclasses import dataclass
import hashlib
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import svds

LOCAL_SNAP_TOL=1e-12
TERMINAL_SNAP_TOL=1e-9
MAX_SNAP_L1=1e-8
NULL_TOL=1e-10
RANK_RTOL=1e-12
PCA_RANK_RTOL=1e-10
PCA_TOL=1e-8
PCA_MAXITER=1000
FEATURE_NAMES=('intercept','scaled_depth','PC1','PC2','PC3','PC4')


def array_sha256(value):
    x=np.ascontiguousarray(value)
    h=hashlib.sha256();h.update(str(x.dtype).encode());h.update(str(x.shape).encode());h.update(x.tobytes())
    return h.hexdigest()


def _blocks(matrix,block=256):
    for start in range(0,matrix.shape[0],block):
        part=matrix[start:start+block]
        yield start,(part.toarray() if sparse.issparse(part) else np.asarray(part)).astype(np.float64,copy=True)


def _validated_counts(raw):
    if sparse.issparse(raw):
        matrix=sparse.csr_matrix(raw,dtype=np.float64,copy=True)
        matrix.sum_duplicates();matrix.sort_indices()
    else:
        matrix=np.asarray(raw)
    if matrix.ndim!=2 or min(matrix.shape)<1: raise ValueError('Nonempty cell-by-gene counts required')
    depth=np.empty(matrix.shape[0]);h=hashlib.sha256();h.update(str(matrix.shape).encode());h.update(b'canonical_int64_counts')
    for start,part in _blocks(matrix):
        if not np.isfinite(part).all() or (part<0).any() or (part!=np.floor(part)).any():
            raise ValueError('Controls must contain finite nonnegative integer counts')
        total=part.sum(1)
        if (total<=0).any() or (total>2**53-1).any(): raise ValueError('Positive exactly representable row depths required')
        depth[start:start+len(part)]=total
        h.update(part.astype('<i8').tobytes())
    return matrix,depth,h.hexdigest()


def fit_features(rawcounts,seed=20261026):
    """Return N x 6 features and a JSON-serializable reconstruction/provenance state.

    Arithmetic CPM uses each cell's full supplied-axis library depth. Genes
    require population mean CPM>5; zero/numerically constant standard deviations
    are dropped. Four numerical leading PCs use seeded truncated ARPACK SVD
    (dense SVD only if min(shape)<=64). Unavailable PC directions are zero padded.
    Freeze these actual features/state before sampling outcomes are examined.
    """
    matrix,depth,counts_hash=_validated_counts(rawcounts);n,g=matrix.shape
    mean=np.zeros(g)
    for start,part in _blocks(matrix):
        mean+=(part/depth[start:start+len(part),None]*1e6).sum(0)
    mean/=n;gate=np.flatnonzero(mean>5.)
    variance=np.zeros(len(gate))
    for start,part in _blocks(matrix):
        cpm=part[:,gate]/depth[start:start+len(part),None]*1e6
        variance+=((cpm-mean[gate])**2).sum(0)
    variance/=n;std=np.sqrt(variance)
    variable=std>1e-12*np.maximum(1.,np.abs(mean[gate]))
    selected=gate[variable];selected_std=std[variable]
    z=np.empty((n,len(selected)),dtype=np.float64)
    for start,part in _blocks(matrix):
        z[start:start+len(part)]=(part[:,selected]/depth[start:start+len(part),None]*1e6-mean[selected])/selected_std
    loading=np.zeros((len(selected),4));singular=[];residual_left=[];residual_right=[];anchors=[]
    method='no_variable_control_genes';components=0
    if len(selected) and n>1:
        count=min(4,n-1,len(selected))
        if min(z.shape)<=64:
            u,s,vt=np.linalg.svd(z,full_matrices=False);u,s,vt=u[:,:count],s[:count],vt[:count]
            method='dense_svd_small_matrix'
        else:
            v0=np.random.default_rng(seed).normal(size=min(z.shape))
            u,s,vt=svds(z,k=count,which='LM',solver='arpack',tol=PCA_TOL,maxiter=PCA_MAXITER,v0=v0)
            order=np.argsort(-s,kind='stable');u,s,vt=u[:,order],s[order],vt[order]
            method='seeded_numerical_truncated_arpack_svd'
        cutoff=PCA_RANK_RTOL*max(float(s.max()),1.)
        for j,value in enumerate(s):
            if value<=cutoff: continue
            vector=vt[j].copy();left=u[:,j].copy();anchor=int(np.argmax(np.abs(vector)))
            if vector[anchor]<0:vector=-vector;left=-left
            loading[:,components]=vector
            singular.append(float(value));anchors.append(int(selected[anchor]))
            residual_left.append(float(np.linalg.norm(z@vector-value*left)/value))
            residual_right.append(float(np.linalg.norm(z.T@left-value*vector)/value))
            components+=1
    scores=z@loading
    score_mean=scores.mean(0);score_std=np.sqrt(((scores-score_mean)**2).mean(0))
    standardized=np.divide(scores-score_mean,score_std,out=np.zeros_like(scores),where=score_std>0)
    depth_mean=float(depth.mean());depth_std=float(np.sqrt(np.mean((depth-depth_mean)**2)))
    scaled_depth=(depth-depth_mean)/depth_std if depth_std>0 else np.zeros(n)
    features=np.column_stack([np.ones(n),scaled_depth,standardized])
    if not np.isfinite(features).all(): raise ValueError('Nonfinite control balancing features')
    state={'schema':'control_cube_features_v1','seed':int(seed),'rows':n,'genes':g,
        'input_counts_sha256':counts_hash,'feature_matrix_sha256':array_sha256(features),
        'feature_names':list(FEATURE_NAMES),'normalization':'Arithmetic CPM=1e6*raw_count/full supplied-axis cell depth',
        'gene_gate':'control population arithmetic mean CPM strictly >5',
        'standardization':'Population gene standard deviation; constant threshold1e-12*max(1,abs(mean))',
        'gate_gene_indices':gate.tolist(),'selected_gene_indices':selected.tolist(),
        'selected_gene_mean_cpm':mean[selected].tolist(),'selected_gene_std_cpm':selected_std.tolist(),
        'dropped_constant_gate_genes':int((~variable).sum()),'depth_mean':depth_mean,'depth_std':depth_std,
        'pca_method':method,'pca_components_requested':4,'pca_components_retained':components,
        'pca_tolerance':PCA_TOL,'pca_maxiter':PCA_MAXITER,'pca_rank_rtol':PCA_RANK_RTOL,
        'singular_values':singular,'singular_relative_residual_left':residual_left,
        'singular_relative_residual_right':residual_right,'loading_anchor_gene_indices':anchors,
        'loading_sign_rule':'Largest absolute loading positive; first supplied gene index breaks exact ties',
        'singular_tie_rule':'Seeded numerical solver basis; freeze actual loadings/scores, do not claim a unique basis at repeated singular values',
        'loadings':loading.tolist(),'pc_score_means':score_mean.tolist(),'pc_score_std':score_std.tolist(),
        'only_controls_permitted':True,'caller_must_verify_control_row_identity':True,
        'feature_fitting_uses_response_labels':False}
    return features,state


def _constraints(features):
    n=features.shape[0];one=np.ones(n)/np.sqrt(n)
    centered=features[:,1:]-features[:,1:].mean(0)
    if centered.shape[1]:
        u,s,_=np.linalg.svd(centered,full_matrices=False)
        cutoff=RANK_RTOL*max(float(s.max(initial=0)),1.)
        keep=s>cutoff;u=u[:,keep]
        # Reproject tiny centering roundoff away from the exact size constraint.
        u-=one[:,None]*(one@u)[None,:]
        if u.shape[1]:u,_=np.linalg.qr(u,mode='reduced')
        constraint=np.vstack([one,u.T])
        singular=s.tolist()
    else:constraint=one[None,:];singular=[]
    return constraint,{'rank':len(constraint),'rank_rtol':RANK_RTOL,'centered_feature_singular_values':singular,
                       'constraint_matrix_sha256':array_sha256(constraint)}


def _snap(probabilities,tolerance=LOCAL_SNAP_TOL):
    old=np.asarray(probabilities,dtype=float);new=old.copy()
    if not np.isfinite(old).all() or (old < -tolerance).any() or (old>1+tolerance).any():
        raise ArithmeticError('Probability left the unit cube beyond numerical tolerance')
    new[np.abs(new)<=tolerance]=0.;new[np.abs(new-1)<=tolerance]=1.
    return new,float(np.abs(new-old).sum())


def null_direction(constraint):
    """One canonical right null vector; full SVD retains the r+1 null row."""
    a=np.asarray(constraint,dtype=float)
    if a.shape[1]<=a.shape[0]:raise ValueError('Flight needs at least r+1 active units')
    _,_,vh=np.linalg.svd(a,full_matrices=True)
    v=vh[-1].copy()
    # The first constraint is a constant intercept; enforce its zero-sum form.
    v[-1]=-np.sum(v[:-1])
    largest=int(np.argmax(np.abs(v)))
    if v[largest]<0:v=-v
    v/=np.max(np.abs(v))
    residual=float(np.max(np.abs(a@v)));size_error=float(abs(v.sum()))
    if residual>NULL_TOL or size_error>NULL_TOL:raise ArithmeticError('SVD null direction failed its balance/size check')
    return v,{'null_residual':residual,'direction_sum_error':size_error}


def flight_branches(probabilities,constraint):
    """Return (p_plus, p_minus, P(plus), numerical diagnostics) for enumeration/tests."""
    p=np.asarray(probabilities,dtype=float)
    if ((p<=0)|(p>=1)).any():raise ValueError('Flight coordinates must be strictly fractional')
    v,diagnostic=null_direction(constraint)
    pos=v>0;neg=v<0
    plus=min(np.min((1-p[pos])/v[pos],initial=np.inf),np.min(-p[neg]/v[neg],initial=np.inf))
    minus=min(np.min(p[pos]/v[pos],initial=np.inf),np.min((p[neg]-1)/v[neg],initial=np.inf))
    if not np.isfinite([plus,minus]).all() or min(plus,minus)<=0:raise ArithmeticError('No finite feasible flight step')
    a,sa=_snap(p+plus*v);b,sb=_snap(p-minus*v)
    if np.all((a>0)&(a<1)) or np.all((b>0)&(b<1)):raise ArithmeticError('Flight failed to fix a coordinate')
    diagnostic.update(lambda_plus=float(plus),lambda_minus=float(minus),plus_snap_l1=sa,minus_snap_l1=sb)
    return a,b,float(minus/(plus+minus)),diagnostic


def pivotal_branches(a,b):
    """Fixed-sum two-coordinate martingale landing; no balance optimization."""
    if not 0<a<1 or not 0<b<1:raise ValueError('Pivotal coordinates must be strictly fractional')
    total=a+b
    if total<1:
        first=np.array([total,0.]);second=np.array([0.,total]);probability=a/total
    else:
        first=np.array([1.,total-1]);second=np.array([total-1,1.]);probability=(1-b)/(2-total)
    first,s1=_snap(first);second,s2=_snap(second)
    if not 0<=probability<=1 or not np.isfinite(probability):raise ArithmeticError('Invalid pivotal probability')
    return first,second,float(probability),{'first_snap_l1':s1,'second_snap_l1':s2}


@dataclass
class PreparedSampler:
    features:np.ndarray
    cells:int=400

    def __post_init__(self):
        self.features=np.array(self.features,dtype=np.float64,copy=True)
        if self.features.ndim!=2 or min(self.features.shape)<1 or not np.isfinite(self.features).all():
            raise ValueError('Finite nonempty population feature matrix required')
        if not np.array_equal(self.features[:,0],np.ones(len(self.features))):raise ValueError('First feature must be an exact intercept1')
        if isinstance(self.cells,(bool,np.bool_)) or not isinstance(self.cells,(int,np.integer)) or not 0<=self.cells<=len(self.features):
            raise ValueError('Sample size must be an integer between0 andN')
        self.cells=int(self.cells);self.population=len(self.features);self.pi=self.cells/self.population
        self.constraint,self.rank_info=_constraints(self.features);self.rank=len(self.constraint)
        self.features.flags.writeable=False;self.constraint.flags.writeable=False
        self.feature_sha256=array_sha256(self.features)

    def select(self,seed):
        rng=np.random.default_rng(seed);order=rng.permutation(self.population)
        p=np.full(self.population,self.pi);initial=p.copy();active=[];cursor=0
        flight_steps=0;landing_steps=0;snap_l1=0.;max_null=0.;max_sum=0.
        if 0<self.cells<self.population:
            while True:
                while len(active)<self.rank+1 and cursor<self.population:
                    i=int(order[cursor]);cursor+=1
                    if 0<p[i]<1:active.append(i)
                if len(active)<=self.rank:break
                indices=np.asarray(active,dtype=int)
                plus,minus,probability,d=flight_branches(p[indices],self.constraint[:,indices])
                choose_plus=rng.random()<probability
                p[indices]=plus if choose_plus else minus
                snap_l1+=d['plus_snap_l1'] if choose_plus else d['minus_snap_l1']
                max_null=max(max_null,d['null_residual']);max_sum=max(max_sum,d['direction_sum_error'])
                flight_steps+=1
                active=[i for i in active if 0<p[i]<1]
            flight_end=p.copy();landing_fractional=len(active)
            while len(active)>1:
                i,j=active[:2]
                a,b,probability,d=pivotal_branches(float(p[i]),float(p[j]))
                first=rng.random()<probability;p[[i,j]]=a if first else b
                snap_l1+=d['first_snap_l1'] if first else d['second_snap_l1']
                landing_steps+=1;active=[k for k in active if 0<p[k]<1]
            terminal_snap=0.
            if active:
                i=active[0];forced=self.cells-int(np.count_nonzero(p==1))
                if forced not in (0,1) or abs(p[i]-forced)>TERMINAL_SNAP_TOL:
                    raise ArithmeticError('Fixed-size intercept lost beyond the declared terminal roundoff allowance')
                terminal_snap=float(abs(p[i]-forced));p[i]=forced;snap_l1+=terminal_snap
        else:
            flight_end=p.copy();landing_fractional=0;terminal_snap=0.
        if snap_l1>MAX_SNAP_L1:raise ArithmeticError('Accumulated probability rounding exceeds the declared budget')
        if not np.all((p==0)|(p==1)):raise ArithmeticError('Landing did not reach an integer sample')
        indices=np.flatnonzero(p==1)
        if len(indices)!=self.cells:raise ArithmeticError('Incorrect cardinality; no subset repair is permitted')
        target=self.features.T@initial;flight_residual=self.features.T@flight_end-target
        final_residual=self.features.T@p-target
        constraint_residual=self.constraint@(flight_end-initial)
        if np.max(np.abs(constraint_residual))>MAX_SNAP_L1:
            raise ArithmeticError('Accumulated flight-constraint drift exceeds the numerical budget')
        audit={'schema':'cube_flight_pivotal_landing_v1','seed':int(seed),'population':self.population,
            'sample_size':self.cells,'inclusion_probability':self.pi,'first_order_probability_guarantee':'Exact-arithmetic coordinatewise martingale; finite-precision snap corrections explicitly bounded below',
            'fixed_size_verified':True,'without_replacement_verified':len(set(indices))==self.cells,
            'feature_sha256':self.feature_sha256,**self.rank_info,
            'flight_steps':flight_steps,'fractional_units_before_landing':landing_fractional,'landing_steps':landing_steps,
            'initial_probability_sum_minus_integer_size':float(initial.sum()-self.cells),
            'max_null_residual':max_null,'max_direction_sum_error':max_sum,
            'local_snap_tolerance':LOCAL_SNAP_TOL,'terminal_snap_tolerance':TERMINAL_SNAP_TOL,
            'terminal_snap_l1':terminal_snap,'accumulated_snap_l1':snap_l1,'maximum_snap_l1':MAX_SNAP_L1,
            'flight_feature_total_residual':flight_residual.tolist(),
            'flight_retained_constraint_residual':constraint_residual.tolist(),
            'final_feature_total_residual':final_residual.tolist(),
            'landing_feature_total_change':(final_residual-flight_residual).tolist(),
            'sample_mean_feature_residual':(final_residual/self.cells).tolist() if self.cells else None,
            'HT_population_feature_total_residual':(final_residual/self.pi).tolist() if self.pi else None,
            'auxiliary_balance_claim':'Approximate after pivotal landing; no exact balance claim or nearest-subset optimization',
            'unbiased_estimands':'Finite-population linear means in exact arithmetic; no unbiased pooled-ratio or log-profile claim',
            'raw_cells_modified':False,'selection_indices_sha256':array_sha256(indices)}
        return indices,audit


def prepare_sampler(features,cells=400):
    return PreparedSampler(features,cells)


def select_balanced(features,seed,cells=400):
    return prepare_sampler(features,cells).select(seed)
