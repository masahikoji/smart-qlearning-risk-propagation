"""Reconstruct supplemental Tables S3/S4 with the audited v17 numerical core.

Added during the 2026-09-11 audit; this is not a recovered historical script.
The filename fills the reproduction entry point named by the manuscript.
"""
from pathlib import Path
import sys,argparse,json,csv
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
CORE=ROOT if (ROOT/'feedback_tuning.py').exists() else ROOT/'input_1/feedback_aware_simulation_v17'
sys.path.insert(0,str(CORE))
from feedback_tuning import population_geometry
from verify_gaussian import quadrature

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/supplement_gaussian')
    p.add_argument('--low-order',type=int,default=256)
    p.add_argument('--high-order',type=int,default=512)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    geom=population_geometry(.7,.3);rows=[];convergence=[]
    for T in (1.,2.,4.):
        for delta in ((0.,0.),(0.,1.),(1.,1.),(2.,3.)):
            lo=quadrature(np.array(delta),geom['chi'],geom['kappa'],a.low_order,T)
            hi=quadrature(np.array(delta),geom['chi'],geom['kappa'],a.high_order,T)
            hi['expected_dispersion']=hi['naive_weighted_candidate_sure']-hi['frozen_sure']
            hi['expected_net_correction']=hi['correction']-hi['expected_dispersion']
            keys=('fa_risk','fa_sure','correction','gradient_cost','frozen_sure','naive_weighted_candidate_sure','soft_oracle_bound','upper_risk_score')
            change=max(abs(lo[k]-hi[k])for k in keys)
            discrepancy=abs(hi['fa_risk']-hi['fa_sure'])
            hi['low_high_max_difference']=change
            convergence.append({'T':T,'delta':delta,'max_change':change,'risk_identity_error':discrepancy})
            rows.append(hi)
    (a.output/'supplement_gaussian_results.json').write_text(json.dumps(rows,indent=2))
    with(a.output/'table_s3_adaptation.csv').open('w',newline='') as f:
        keys=['delta1','delta2','correction','expected_dispersion','expected_net_correction','gradient_cost']
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:r[k]for k in keys}for r in rows if r['temperature']==2)
    with(a.output/'table_s4_tuning.csv').open('w',newline='') as f:
        keys=['delta1','delta2','temperature','fa_risk','akaike_limit_risk'];w=csv.DictWriter(f,fieldnames=keys)
        w.writeheader();w.writerows({k:r[k]for k in keys}for r in rows)
    result={'added_by_audit':True,'purpose':'Supplement Tables S3 and S4','orders':[a.low_order,a.high_order],
        'max_integral_change':max(r['max_change']for r in convergence),
        'max_risk_identity_error':max(r['risk_identity_error']for r in convergence),'cases':convergence}
    (a.output/'convergence.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
