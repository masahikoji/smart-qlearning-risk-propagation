"""Finite-sample SMART verification with estimated geometry.

All procedures use the same generated data in each replicate. Stage-1 risk is
integrated exactly using its known population Gram matrix. No test-sample
noise is added. Monte Carlo standard errors are computed on paired gaps.
"""
from __future__ import annotations
import argparse,csv,json,time
from pathlib import Path
import numpy as np
from feedback_tuning import (fit_smart, simulate_smart, stage1_risk, aggregate,
                             population_geometry, default_rules)
from verify_gaussian import quadrature


def main():
    p=argparse.ArgumentParser();p.add_argument('--replicates',type=int,default=5000)
    p.add_argument('--sizes',default='250,1000,4000');p.add_argument('--seed',type=int,default=20260909)
    p.add_argument('--output',default='results');p.add_argument('--criterion',choices=('wald','logrss'),default='wald');p.add_argument('--case',type=int,default=-1);a=p.parse_args()
    ns=[int(x) for x in a.sizes.split(',')]
    outdir=Path(a.output);outdir.mkdir(parents=True,exist_ok=True)
    geom=population_geometry(.7,.3)
    deltas=[(0.,0.),(0.,1.),(1.,1.),(2.,3.)]
    columns=['risk_fa','risk_akaike_limit','risk_wide','risk_ablation',
             'fa_minus_akaike','fa_minus_ablation','reported_scaled_sure',
             'w1','w2','chi','kappa','common_risk_hat']
    rows=[]; clock=time.time()
    # Each configuration has an invocation-order-independent random stream.
    for s,delta in enumerate(deltas):
        if a.case >= 0 and s != a.case:
            continue
        lim=quadrature(np.array(delta),geom['chi'],geom['kappa'],order=1024)
        for k,n in enumerate(ns):
            rng=np.random.default_rng(np.random.SeedSequence([a.seed,s,n]))
            rec=np.full((a.replicates,len(columns)),np.nan)
            fail=0;switch=0; clipped=0;maxidentity=0.;mean_alpha=np.zeros(len(default_rules()))
            for r in range(a.replicates):
                data,truth=simulate_smart(n,*delta,rng)
                fit=fit_smart(data,criterion=a.criterion)
                if fit['stabilized']:
                    fail+=1
                    # Prespecified common zero fallback: retain its risk.
                    rr=float(n*stage1_risk(fit['coefficients1'],truth))
                    rec[r,:6]=[rr,rr,rr,rr,0,0]
                    continue
                cc=fit['candidate_coefficients1'].T
                risks=n*stage1_risk(cc,truth)
                rr=float(n*stage1_risk(fit['coefficients1'],truth))
                ab=aggregate(fit['w'],fit['chi'],fit['kappa'],omit_feedback_for_weights=True)
                rab=float(n*stage1_risk(ab['alpha']@cc,truth))
                rec[r]=[rr,risks[4],risks[-1],rab,rr-risks[4],rr-rab,
                        fit['scaled_risk_estimate'],*fit['w'],fit['chi'],fit['kappa'],
                        fit['common_risk_estimate']]
                switch+=int(fit['actions_differ_from_reference']>0)
                clipped+=int(fit['geometry_clipped'])
                maxidentity=max(maxidentity,fit['response_score_identity_error'])
                mean_alpha+=fit['alpha']
            result=dict(delta1=delta[0],delta2=delta[1],n=n,criterion=a.criterion,replicates=a.replicates,
                        seed=a.seed,stabilization_count=fail,switch_replicates=switch,
                        clipped_geometry_replicates=clipped,max_score_identity_error=maxidentity,
                        limit_fa=geom['common_risk']+lim['fa_risk'],
                        limit_akaike=geom['common_risk']+lim['akaike_limit_risk'],
                        limit_gap=lim['fa_risk']-lim['akaike_limit_risk'],
                        limit_ablation_gap=lim['fa_risk']-lim['ablation_risk'],
                        mean_alpha=(mean_alpha/a.replicates).tolist())
            for i,col in enumerate(columns):
                valid=rec[np.isfinite(rec[:,i]),i]
                result[col]=float(valid.mean()) if len(valid) else None
                result[col+'_mcse']=float(valid.std(ddof=1)/np.sqrt(len(valid))) if len(valid)>1 else None
            result['score_covariance']=np.cov(rec[:,7:9],rowvar=False).tolist() if not fail else None
            rows.append(result)
            print(f"delta={delta} n={n} R={a.replicates} FA={result['risk_fa']:.5f} "
                  f"(SE {result['risk_fa_mcse']:.5f}) limit={result['limit_fa']:.5f} "
                  f"SURE={result['reported_scaled_sure']:.5f} gap={result['fa_minus_akaike']:.5f} "
                  f"(SE {result['fa_minus_akaike_mcse']:.5f}), elapsed={time.time()-clock:.1f}s",flush=True)
            (outdir/'smart_results.json').write_text(json.dumps(rows,indent=2))
    keys=[k for k in rows[0] if k not in ('mean_alpha','score_covariance')]
    with (outdir/'smart_results.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=keys);wr.writeheader()
        wr.writerows([{k:r[k] for k in keys} for r in rows])
    (outdir/'smart_metadata.json').write_text(json.dumps(dict(geometry=geom,replicates=a.replicates,
        sizes=ns,seed=a.seed,elapsed_seconds=time.time()-clock,
        evaluation='Exact stage-1 population coefficient-Gram MSE',
        normalisation='n-scaled full Q1 risk; barsigma1=1',
        temperature=2.,rules=[r.__dict__ for r in default_rules()],
        finite_sample_weights=a.criterion),indent=2))

if __name__=='__main__':main()
