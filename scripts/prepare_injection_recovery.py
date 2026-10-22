#!/usr/bin/env python3
"""Independent control-only composition and latent DM reference laws.

Intentionally imports no candidate emitter, donor sampler, count transport,
moment-repair kernel, or prediction builder.
"""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ.setdefault(key,'1')
from pathlib import Path
import argparse
from datetime import datetime,timezone
import hashlib
import importlib.metadata
import json
import time
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import gammaln,logsumexp
from sklearn.utils.extmath import randomized_svd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'evaluations/injection_recovery_400'
PROTOCOL=ROOT/'docs/INJECTION_RECOVERY_PROTOCOL_20261021.md'
SEEDS=(20261026,20261027)
REGIMES=('NULL','WEAK','MODERATE','STRONG')
DM_STRENGTH=(0.,.1,.3,.6)
TAUS=(0.,float(np.log(1.25)),float(np.log(2)),float(np.log(4)))
KAPPA=20000.
GENE_NULL_TOL=1e-14


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
    return h.hexdigest()


def named_seed(seed,name):
    return (seed+int.from_bytes(hashlib.sha256(name.encode()).digest()[:4],'little'))%2**32


def write_json(path,value):
    path=Path(path)
    if path.exists():raise FileExistsError(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def read(path):return json.loads(Path(path).read_text())


def checked_counts(value):
    x=sparse.csr_matrix(value)
    if not np.isfinite(x.data).all() or (x.data<0).any() or (x.data>1_000_000).any() or (x.data!=np.floor(x.data)).any():raise ValueError('Non-count or out-of-range input')
    x=x.astype(np.int32);x.sum_duplicates();x.eliminate_zeros()
    depth=np.asarray(x.sum(axis=1,dtype=np.int64)).ravel()
    if (depth<=0).any() or (depth>1_000_000).any():raise ValueError('Invalid integer library depth')
    return x,depth


def population_moments(counts):
    x,d=checked_counts(counts)
    p=np.asarray(x.astype(np.float64).multiply((1/d)[:,None]).mean(axis=0)).ravel()
    c=np.asarray(x.mean(axis=0,dtype=np.float64)).ravel()
    return dict(p=p,b=c/d.mean(),mean_counts=c,mean_depth=float(d.mean()))


def fisher_allocation(n_low,n_high,cells,tau):
    if any(type(v) not in (int,np.int64,np.int32) for v in (n_low,n_high,cells)) or min(n_low,n_high)<0 or not 1<=cells<=n_low+n_high or not np.isfinite(tau):
        raise ValueError('Invalid finite-population allocation')
    k=np.arange(max(0,cells-n_low),min(cells,n_high)+1,dtype=np.int64)
    def logchoose(n,r):return gammaln(n+1)-gammaln(r+1)-gammaln(n-r+1)
    logp=logchoose(n_high,k)+logchoose(n_low,cells-k)+tau*k
    prob=np.exp(logp-logsumexp(logp));prob/=prob.sum()
    return k,prob


def composition_moments(counts,high,cells,tau):
    high=np.asarray(high,bool)
    if len(high)!=counts.shape[0] or not high.any() or high.all():raise ValueError('Both reference states required')
    ks,prob=fisher_allocation(int((~high).sum()),int(high.sum()),cells,tau)
    weight=float(ks@prob/cells)
    lo=population_moments(counts[~high]);hi=population_moments(counts[high])
    if tau==0:
        result=population_moments(counts)
        weight=float(high.mean())
    else:
        c=(1-weight)*lo['mean_counts']+weight*hi['mean_counts'];depth=(1-weight)*lo['mean_depth']+weight*hi['mean_depth']
        result=dict(p=(1-weight)*lo['p']+weight*hi['p'],b=c/depth,mean_counts=c,mean_depth=float(depth))
    return result,dict(expected_high_fraction=weight,k=ks,probability=prob)


def draw_composition(high,cells,tau,seed):
    high=np.asarray(high,bool);ks,prob=fisher_allocation(int((~high).sum()),int(high.sum()),cells,tau)
    rng=np.random.default_rng(seed);k=int(rng.choice(ks,p=prob))
    chosen=np.r_[rng.choice(np.flatnonzero(high),k,replace=False),rng.choice(np.flatnonzero(~high),cells-k,replace=False)]
    return rng.permutation(chosen)


def dm_population_moments(state_means,atom_states,atom_depths):
    mu=np.asarray(state_means,dtype=np.float64);states=np.asarray(atom_states);depth=np.asarray(atom_depths,dtype=np.float64)
    if mu.ndim!=2 or not np.isfinite(mu).all() or (mu<0).any() or not np.allclose(mu.sum(1),1,rtol=0,atol=1e-12):raise ValueError('Invalid state probabilities')
    if len(states)!=len(depth) or len(depth)==0 or not np.issubdtype(states.dtype,np.integer) or (states<0).any() or (states>=len(mu)).any():raise ValueError('Invalid joint atoms')
    if not np.isfinite(depth).all() or (depth<=0).any() or (depth!=np.floor(depth)).any() or (depth>1e6).any():raise ValueError('Invalid atom depths')
    w=np.bincount(states,minlength=len(mu))/len(states);mass=np.bincount(states,weights=depth,minlength=len(mu))/len(states)
    c=mass@mu
    return dict(p=w@mu,b=c/depth.mean(),mean_counts=c,mean_depth=float(depth.mean()))


def transfer_program(state_means,up,down,fraction):
    mu=np.asarray(state_means,dtype=np.float64);up=np.asarray(up,int);down=np.asarray(down,int)
    if mu.ndim!=2 or not np.isfinite(mu).all() or (mu<0).any() or not np.allclose(mu.sum(1),1,rtol=0,atol=1e-12):raise ValueError('Invalid state means')
    if not np.isfinite(fraction) or not 0<=fraction<1 or len(up)==0 or len(down)==0 or len(set(up)&set(down)):
        raise ValueError('Invalid probability-transfer program')
    if len(set(up))!=len(up) or len(set(down))!=len(down):raise ValueError('Duplicate program indices')
    if min(up.min(),down.min())<0 or max(up.max(),down.max())>=mu.shape[1]:raise ValueError('Program gene index outside axis')
    q=mu.copy()
    if fraction==0:return q
    for state in range(len(mu)):
        u=float(mu[state,up].sum());d=float(mu[state,down].sum())
        if min(u,d)<=0:raise ValueError('Program requires expressed up/down support in every state')
        mass=fraction*min(u,d)
        q[state,up]*=1+mass/u;q[state,down]*=1-mass/d
    q/=q.sum(axis=1,keepdims=True)
    return q


def simulate_dm(state_means,atom_states,atom_depths,n_cells,kappa,seed):
    dm_population_moments(state_means,atom_states,atom_depths)
    if type(n_cells) is not int or n_cells<1 or not np.isfinite(kappa) or kappa<=0:raise ValueError('Invalid simulation size/concentration')
    mu=np.asarray(state_means,dtype=np.float64);depth=np.asarray(atom_depths,dtype=np.int64);states=np.asarray(atom_states)
    rng=np.random.default_rng(seed);atoms=rng.integers(0,len(states),size=n_cells)
    result=np.zeros((n_cells,mu.shape[1]),dtype=np.int32)
    supports=[np.flatnonzero(row>0) for row in mu]
    for i,atom in enumerate(atoms):
        s=int(states[atom]);support=supports[s]
        theta=rng.dirichlet(kappa*mu[s,support])
        result[i,support]=rng.multinomial(int(depth[atom]),theta)
    if not np.array_equal(result.sum(1,dtype=np.int64),depth[atoms]):raise AssertionError('Simulator lost library depths')
    return result,atoms


def control_roles():
    prior_path=ROOT/'evaluations/emission_null_400/manifest.json';prior=read(prior_path)
    inference=Path(prior['inference_controls']);reference=ROOT/'evaluations/emission_null_400/reference.h5ad'
    if sha(inference)!=prior['inference_sha256'] or sha(reference)!=prior['reference_sha256']:raise ValueError('Frozen control sources changed')
    d=ad.read_h5ad(inference);r=ad.read_h5ad(reference)
    if set(d.obs.target_gene.astype(str))!={'non-targeting'} or set(r.obs.source_target_gene.astype(str))!={'non-targeting'}:raise ValueError('Only true non-targeting controls are permitted')
    comp=np.flatnonzero(r.obs.target_gene.astype(str).to_numpy()=='non-targeting');ref=np.flatnonzero(r.obs.target_gene.astype(str).to_numpy()!='non-targeting')
    if d.n_obs!=4000 or len(comp)!=4000 or len(ref)!=4000 or not d.var_names.equals(r.var_names):raise ValueError('Unexpected control counts or axes')
    sets=[set(d.obs_names),set(r.obs_names[comp]),set(r.obs_names[ref])]
    if any(len(s)!=4000 for s in sets) or any(sets[i]&sets[j] for i in range(3) for j in range(i)):raise ValueError('Inference/comparison/reference source IDs overlap')
    for x in (d.X,r.X):checked_counts(x)
    provenance=dict(input_manifest_path=str(prior_path),input_manifest_sha256=sha(prior_path),
        inference_path=str(inference),inference_sha256=sha(inference),old_null_reference_path=str(reference),old_null_reference_sha256=sha(reference),
        comparison_rows=comp.tolist(),reference_pool_rows=ref.tolist(),all_source_cells_non_targeting=True,all_three_roles_disjoint=True,
        source_ids_sha256=[hashlib.sha256('\0'.join(sorted(s)).encode()).hexdigest() for s in sets])
    return d,r,comp,ref,provenance


def fit_control_states(counts):
    x,depth=checked_counts(counts);linear=x.astype(np.float64).multiply((1e6/depth)[:,None]).tocsr()
    means=np.asarray(linear.mean(axis=0)).ravel();gate=means>5
    matrix=linear[:,gate].toarray();np.log1p(matrix,out=matrix)
    center=matrix.mean(0);scale=matrix.std(0,ddof=0);keep=scale>1e-8
    selected=np.flatnonzero(gate)[keep];matrix=matrix[:,keep];center=center[keep];scale=scale[keep]
    matrix-=center;matrix/=scale
    _,_,components=randomized_svd(matrix,n_components=2,n_iter=5,random_state=20261026)
    for c in components:
        if c[np.argmax(np.abs(c))]<0:c*=-1
    scores=matrix@components.T
    return dict(selected_gene_indices=selected,feature_center=center,feature_scale=scale,components=components,
        split_medians=np.median(scores,axis=0),dm_quartiles=np.quantile(scores[:,0],[.25,.5,.75],method='linear')),scores


def project_states(counts,parameters):
    x,depth=checked_counts(counts)
    q=x[:,parameters['selected_gene_indices']].astype(np.float64).multiply((1e6/depth)[:,None]).toarray()
    q=np.log1p(q);q-=parameters['feature_center'];q/=parameters['feature_scale']
    return q@parameters['components'].T


def select_program(genes,mean_cpm,detection,name,excluded=()):
    eligible=(mean_cpm>=10)&(mean_cpm<= (1000 if name=='BROAD' else 100))
    if name=='ZERO_RICH':eligible&=detection<=.5
    eligible[np.asarray(list(excluded),dtype=int)]=False
    rows=np.flatnonzero(eligible);rows=rows[np.lexsort((genes[rows],mean_cpm[rows]))]
    pairs=list(zip(rows[::2],rows[1::2]))
    def key(pair):return hashlib.sha256(f'20261026:{name}:{genes[pair[0]]}|{genes[pair[1]]}'.encode()).digest()
    pairs=sorted(pairs,key=key)
    if len(pairs)<32:raise ValueError('Insufficient predeclared gene-program support: '+name)
    up=[];down=[]
    for pair in pairs[:32]:
        a,b=pair if key(pair)[0]%2==0 else pair[::-1];up.append(int(a));down.append(int(b))
    return np.asarray(up),np.asarray(down)


def profile_bank(aliases,templates,regimes,strengths,baseline,moments,direct,marginal_null,genes):
    p=np.stack([m['p'] for m in moments]);b=np.stack([m['b'] for m in moments]);null=np.asarray(regimes)=='NULL'
    p[null]=baseline['p'];b[null]=baseline['b']
    dp=p-baseline['p'];db=b-baseline['b'];dp[null]=0;db[null]=0
    lfc=np.log2((1e6*p+1e-9)/(1e6*baseline['p']+1e-9));bd=np.log1p(50000*b)-np.log1p(50000*baseline['b']);lfc[null]=0;bd[null]=0
    return dict(aliases=np.asarray(aliases,dtype=str),targets=np.asarray(aliases,dtype=str),genes=np.asarray(genes,dtype=str),oracle_p=p,oracle_b=b,
        baseline_p=baseline['p'],baseline_b=baseline['b'],delta_p=dp,delta_b=db,truth_lfc_cpm_eps1e9=lfc,
        truth_bulk_delta_log1p50k=bd,expected_mean_counts=np.stack([m['mean_counts'] for m in moments]),
        expected_depth=np.asarray([m['mean_depth'] for m in moments]),baseline_expected_mean_counts=baseline['mean_counts'],
        baseline_expected_depth=np.asarray(baseline['mean_depth']),is_null=null,template_id=np.asarray(templates),regime=np.asarray(regimes),
        strength=np.asarray(strengths),direct_modified_gene_mask=np.asarray(direct,bool),marginal_null_gene_mask=np.asarray(marginal_null,bool))


def freeze():
    if OUT.exists():raise FileExistsError('Injection truth directory exists; no overwrite')
    started=time.time();d,r,comp,ref,source=control_roles();genes=np.asarray(d.var_names,dtype=str);G=len(genes)
    x,depth=checked_counts(d.X);xr,_=checked_counts(r.X[ref]);parameters,scores=fit_control_states(x);rscores=project_states(xr,parameters)
    comp_base=population_moments(xr);aliases=[];templates=[];regimes=[];strengths=[];moment_rows=[];direct=[];nullm=[];highs=[];weights=[]
    for template in range(2):
        high=rscores[:,template]>=parameters['split_medians'][template]
        if min(int(high.sum()),int((~high).sum()))<400:raise ValueError('Unsupported frozen reference state bin')
        highs.append(high)
        for regime,tau in zip(REGIMES,TAUS):
            tau=tau if template==0 else -tau
            m,a=composition_moments(xr,high,400,tau)
            aliases.append(f'COMP_T{template+1}_{regime}');templates.append(template+1);regimes.append(regime);strengths.append(tau)
            moment_rows.append(m);weights.append(a['expected_high_fraction']);direct.append(np.zeros(G,bool));nullm.append(np.full(G,regime=='NULL'))
    comp_profiles=profile_bank(aliases,templates,regimes,strengths,comp_base,moment_rows,direct,nullm,genes)
    state=np.searchsorted(parameters['dm_quartiles'],scores[:,0],side='right').astype(np.int64)
    global_p=population_moments(x)['p'];mus=[]
    for s in range(4):
        if not np.any(state==s):raise ValueError('Empty DM state')
        mu=.95*population_moments(x[state==s])['p']+.05*global_p;mu/=mu.sum();mus.append(mu)
    mu0=np.stack(mus);dm_base=dm_population_moments(mu0,state,depth);detection=np.asarray((x>0).mean(axis=0)).ravel()
    program1=select_program(genes,1e6*global_p,detection,'BROAD')
    program2=select_program(genes,1e6*global_p,detection,'ZERO_RICH',excluded=np.r_[*program1])
    aliases=[];templates=[];regimes=[];strengths=[];moment_rows=[];direct=[];nullm=[];program_means=[]
    for template,(up,down) in enumerate((program1,program2)):
        for regime,fraction in zip(REGIMES,DM_STRENGTH):
            mu=transfer_program(mu0,up,down,fraction);m=dm_population_moments(mu,state,depth)
            aliases.append(f'DM_T{template+1}_{regime}');templates.append(template+1);regimes.append(regime);strengths.append(fraction);moment_rows.append(m);program_means.append(mu)
            mask=np.zeros(G,bool);mask[np.r_[up,down]]=fraction!=0;direct.append(mask)
            # Conditional means AND concentration/joint depth-state law are unchanged.
            nullm.append(np.all(np.abs(mu-mu0)<=GENE_NULL_TOL,axis=0))
    dm_profiles=profile_bank(aliases,templates,regimes,strengths,dm_base,moment_rows,direct,nullm,genes)
    OUT.mkdir(parents=True)
    np.savez_compressed(OUT/'population_parameters.npz',genes=genes,**parameters,composition_reference_high=np.stack(highs),
        composition_expected_high_weights=np.asarray(weights),dm_atom_states=state,dm_atom_depths=depth,dm_baseline_state_means=mu0,
        dm_alias_state_means=np.stack(program_means),dm_kappa=np.asarray(KAPPA),program1_up=program1[0],program1_down=program1[1],
        program2_up=program2[0],program2_down=program2[1])
    families={}
    for name,bank,context in [('composition',comp_profiles,'COMPOSITION'),('dirichlet_multinomial',dm_profiles,'DM')]:
        folder=OUT/name;folder.mkdir();np.savez_compressed(folder/'profiles.npz',**bank)
        pd.DataFrame({'gene_name':genes}).to_csv(folder/'gene_names.csv',index=False)
        pd.DataFrame({'target_gene':bank['aliases'],'n_cells':400}).to_csv(folder/'pert_counts.csv',index=False)
        families[name]=dict(context=context,profiles_path=str(folder/'profiles.npz'),profiles_sha256=sha(folder/'profiles.npz'),aliases=bank['aliases'].tolist(),
            development_aliases=bank['aliases'][bank['template_id']==1].tolist(),heldout_template_aliases=bank['aliases'][bank['template_id']==2].tolist(),
            genes_path=str(folder/'gene_names.csv'),targets_path=str(folder/'pert_counts.csv'),seeds={})
    source['inference_empirical_p']=population_moments(x)['p'].tolist();source['inference_empirical_b']=population_moments(x)['b'].tolist()
    files=[Path(__file__),PROTOCOL,ROOT/'tests/test_injection_recovery.py']
    freeze_record=dict(created_utc=datetime.now(timezone.utc).isoformat(),source_provenance=source,
        files={str(p):sha(p) for p in files},seeds=list(SEEDS),seed_roles={'20261026':'development','20261027':'heldout_seed'},
        state_learning='Inference-only standardized log1pCPM,2PCs,fixed sign,n_iter5,seed20261026',
        dm_kappa=KAPPA,dm_variance_multiplier_at_depth50000=(50000+KAPPA)/(1+KAPPA),marginal_null_conditional_probability_tolerance=GENE_NULL_TOL,
        population_parameters_path=str(OUT/'population_parameters.npz'),population_parameters_sha256=sha(OUT/'population_parameters.npz'),
        families=families,request_semantics='Anchored probability deltas; explicitNULL identity baseline; root candidate builder owns clipping and audits',
        sign_threshold=.1,effect_bins=[.1,.25,.5],reference_null_max_significant_calls=5,
        candidate_or_reference_recovery_outcomes_read_before_freeze=False,
        versions={n:importlib.metadata.version(n) for n in ('numpy','scipy','scikit-learn','anndata')},seconds=time.time()-started)
    write_json(OUT/'freeze.json',freeze_record);write_json(ROOT/'reports/injection_recovery_truth_freeze.json',freeze_record)
    print(json.dumps({'frozen':True,'families':2,'aliases_each':8,'seconds':time.time()-started}),flush=True)


def verify_freeze():
    f=read(OUT/'freeze.json')
    for p,h in f['files'].items():
        if sha(p)!=h:raise ValueError('Frozen truth code/protocol changed: '+p)
    source=f['source_provenance']
    for path_key,hash_key in [('input_manifest_path','input_manifest_sha256'),('inference_path','inference_sha256'),
                             ('old_null_reference_path','old_null_reference_sha256')]:
        if sha(source[path_key])!=source[hash_key]:raise ValueError('Frozen control provenance changed: '+path_key)
    if sha(f['population_parameters_path'])!=f['population_parameters_sha256']:raise ValueError('Latent parameters changed')
    for family in f['families'].values():
        if sha(family['profiles_path'])!=family['profiles_sha256']:raise ValueError('Oracle profiles changed')
    return f


def save_adata(path,counts,genes,context,labels,ids,metadata=None):
    path=Path(path)
    if path.exists():raise FileExistsError(path)
    x,depth=checked_counts(counts);obs=pd.DataFrame({'target_gene':pd.Categorical(labels),'context':context},index=ids)
    if not obs.index.is_unique:raise ValueError('Output IDs must be unique')
    if metadata:
        for key,values in metadata.items():obs[key]=values
    obj=ad.AnnData(x,obs=obs,var=pd.DataFrame(index=genes));obj.uns['purpose']='Independent known-signal reference; never submit'
    path.parent.mkdir(parents=True,exist_ok=True);obj.write_h5ad(path,compression='gzip',compression_opts=1)
    return dict(path=str(path),sha256=sha(path),shape=list(obj.shape),stored_entries=int(x.nnz),max_library=int(depth.max()),counts_validated=True)


def generate(family,seed):
    f=verify_freeze();spec=f['families'][family];folder=OUT/family/f'seed_{seed}'
    if folder.exists():raise FileExistsError(folder)
    folder.mkdir();parameters=np.load(f['population_parameters_path'],allow_pickle=False);bank=np.load(spec['profiles_path'],allow_pickle=False)
    genes=bank['genes'].astype(str);aliases=bank['aliases'].astype(str);context=spec['context'];start=time.time();blocks=[];labels=[];ids=[];source_rows=[];role_audits={}
    if family=='composition':
        d,r,comp,ref,source=control_roles();xin,dep=checked_counts(d.X);xc,_=checked_counts(r.X[comp]);xr,_=checked_counts(r.X[ref])
        inference_ids=d.obs_names.astype(str).tolist();comparison_ids=r.obs_names[comp].astype(str).tolist()
        role_audits['inference']=save_adata(folder/'inference.h5ad',xin,genes,context,['non-targeting']*4000,
            [f'{context}_{seed}_infer_{i}' for i in range(4000)],{'original_source_cell_id':inference_ids})
        role_audits['comparison']=save_adata(folder/'comparison.h5ad',xc,genes,context,['non-targeting']*4000,
            [f'{context}_{seed}_comp_{i}' for i in range(4000)],{'original_source_cell_id':comparison_ids})
        for i,alias in enumerate(aliases):
            high=parameters['composition_reference_high'][int(bank['template_id'][i])-1]
            chosen=draw_composition(high,400,float(bank['strength'][i]),named_seed(seed,f'composition:{alias}:reference'))
            block=xr[chosen];blocks.append(block);labels.extend([alias]*400);ids.extend([f'{context}_{seed}_{alias}_{j}' for j in range(400)])
            source_rows.extend(dict(alias=alias,row=j,reference_pool_row=int(v),original_source_cell_id=str(r.obs_names[ref[v]]),state_high=bool(high[v])) for j,v in enumerate(chosen))
        pd.DataFrame(source_rows).to_csv(folder/'reference_source_rows.csv',index=False)
        inference_m=population_moments(xin);comparison_m=population_moments(xc)
    else:
        state=parameters['dm_atom_states'];depth=parameters['dm_atom_depths'];mu0=parameters['dm_baseline_state_means']
        xin,atoms_in=simulate_dm(mu0,state,depth,4000,KAPPA,named_seed(seed,'dm:inference'))
        xc,atoms_comp=simulate_dm(mu0,state,depth,4000,KAPPA,named_seed(seed,'dm:comparison'))
        role_audits['inference']=save_adata(folder/'inference.h5ad',xin,genes,context,['non-targeting']*4000,
            [f'{context}_{seed}_infer_{i}' for i in range(4000)],{'latent_depth_atom':atoms_in,'latent_state':state[atoms_in]})
        role_audits['comparison']=save_adata(folder/'comparison.h5ad',xc,genes,context,['non-targeting']*4000,
            [f'{context}_{seed}_comp_{i}' for i in range(4000)],{'latent_depth_atom':atoms_comp,'latent_state':state[atoms_comp]})
        for i,alias in enumerate(aliases):
            block,atoms=simulate_dm(parameters['dm_alias_state_means'][i],state,depth,400,KAPPA,named_seed(seed,f'dm:{alias}:reference'))
            blocks.append(sparse.csr_matrix(block));labels.extend([alias]*400);ids.extend([f'{context}_{seed}_{alias}_{j}' for j in range(400)])
            source_rows.extend(dict(alias=alias,row=j,latent_depth_atom=int(v),latent_state=int(state[v]),independent_synthetic_draw=True) for j,v in enumerate(atoms))
        pd.DataFrame(source_rows).to_csv(folder/'reference_latent_rows.csv',index=False)
        inference_m=population_moments(xin);comparison_m=population_moments(xc)
    empirical=[]
    for i,(alias,block) in enumerate(zip(aliases,blocks)):
        m=population_moments(block);empirical.append(m)
    blocks.append(sparse.csr_matrix(xc));labels.extend(['non-targeting']*4000);ids.extend([f'{context}_{seed}_comp_{i}' for i in range(4000)])
    role_audits['reference']=save_adata(folder/'reference.h5ad',sparse.vstack(blocks,format='csr'),genes,context,labels,ids)
    np.savez_compressed(folder/'empirical_profiles.npz',aliases=aliases,genes=genes,
        inference_p=inference_m['p'],inference_b=inference_m['b'],comparison_p=comparison_m['p'],comparison_b=comparison_m['b'],
        reference_p=np.stack([m['p'] for m in empirical]),reference_b=np.stack([m['b'] for m in empirical]),
        inference_baseline_p_offset=inference_m['p']-bank['baseline_p'],inference_baseline_b_offset=inference_m['b']-bank['baseline_b'],
        comparison_baseline_p_offset=comparison_m['p']-bank['baseline_p'],comparison_baseline_b_offset=comparison_m['b']-bank['baseline_b'])
    receipt=dict(family=family,seed=seed,context=context,aliases=aliases.tolist(),roles=role_audits,
        inference_path=role_audits['inference']['path'],comparison_path=role_audits['comparison']['path'],reference_path=role_audits['reference']['path'],
        inference_sha256=role_audits['inference']['sha256'],comparison_sha256=role_audits['comparison']['sha256'],reference_sha256=role_audits['reference']['sha256'],
        profiles_path=spec['profiles_path'],profiles_sha256=spec['profiles_sha256'],population_parameters_path=f['population_parameters_path'],
        population_parameters_sha256=f['population_parameters_sha256'],freeze_sha256=sha(OUT/'freeze.json'),
        empirical_profiles_path=str(folder/'empirical_profiles.npz'),empirical_profiles_sha256=sha(folder/'empirical_profiles.npz'),
        source_rows_path=str(folder/('reference_source_rows.csv' if family=='composition' else 'reference_latent_rows.csv')),
        seconds=time.time()-start,no_candidate_emitter_used=True,no_actual_perturbed_response_used=True,
        heldout_seed=seed==SEEDS[1],heldout_template_aliases=aliases[bank['template_id']==2].tolist())
    receipt['source_rows_sha256']=sha(receipt['source_rows_path']);write_json(folder/'receipt.json',receipt)
    print(json.dumps({'family':family,'seed':seed,'complete':True,'seconds':time.time()-start}),flush=True)


def finalize():
    f=verify_freeze();manifest=dict(freeze_sha256=sha(OUT/'freeze.json'),population_parameters_path=f['population_parameters_path'],
        population_parameters_sha256=f['population_parameters_sha256'],families=f['families'],seeds=list(SEEDS),
        candidate_request_semantics=f['request_semantics'],sign_threshold=.1,effect_bins=[.1,.25,.5],reference_null_max_significant_calls=5,
        no_candidate_or_DE_outcomes_used=True,notes='Absolute truth and ownbaseline remain separate from empirical control offsets; this is not biological validation.')
    for family in manifest['families']:
        for seed in SEEDS:
            rp=OUT/family/f'seed_{seed}'/'receipt.json';r=read(rp)
            for role in r['roles'].values():
                if sha(role['path'])!=role['sha256']:raise ValueError('Truth output hash changed')
            manifest['families'][family]['seeds'][str(seed)]=dict(r,receipt_path=str(rp),receipt_sha256=sha(rp))
    write_json(OUT/'manifest.json',manifest);write_json(ROOT/'reports/injection_recovery_truth_ready.json',manifest)
    print(json.dumps({'truth_ready':True,'families':2,'seeds_each':2}),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('stage',choices=['freeze','generate','finalize'])
    p.add_argument('--family',choices=['composition','dirichlet_multinomial']);p.add_argument('--seed',type=int,choices=SEEDS);a=p.parse_args()
    if a.stage=='generate' and (a.family is None or a.seed is None):p.error('--family and --seed required')
    if a.stage=='freeze':freeze()
    elif a.stage=='generate':generate(a.family,a.seed)
    else:finalize()


if __name__=='__main__':main()
