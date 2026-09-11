"""Independent Gaussian quadrature checks and a preregistered diagnostic grid.

No local signal is passed to the fitting function. It appears only when
integrating its sampling risk against N(delta,I). Numerical convergence is
reported explicitly; quadrature is not substituted for a proof.
"""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np
from scipy.special import roots_hermitenorm, logsumexp, ndtr
from scipy.integrate import quad
from feedback_tuning import aggregate, target_matrix, default_rules, population_geometry


def hard_gaussian_risk(delta, chi, kappa):
    """External hard-AIC risk: analytic inner moments, split adaptive outer integral."""
    threshold=np.sqrt(2.)
    theta=float(delta[0]+chi*delta[1])
    normconst=1/np.sqrt(2*np.pi)
    def integrand(z):
        h2=z if abs(z)>threshold else 0.
        mu=delta[0]+chi*h2
        a=threshold-mu;b=-threshold-mu
        pa=normconst*np.exp(-a*a/2);pb=normconst*np.exp(-b*b/2)
        pout=ndtr(b)+ndtr(-a)
        m1=mu*pout+pa-pb
        m2=(mu*mu+1)*pout+(a+2*mu)*pa-(b+2*mu)*pb
        risk=m2-2*theta*m1+theta*theta+kappa*kappa*(h2-delta[1])**2
        return risk*normconst*np.exp(-(z-delta[1])**2/2)
    bounds=(-np.inf,-threshold,threshold,np.inf)
    vals=[quad(integrand,bounds[i],bounds[i+1],epsabs=1e-10,epsrel=1e-10,limit=200) for i in range(3)]
    return sum(t[0] for t in vals),sum(t[1] for t in vals)


def quadrature(delta, chi, kappa, order=640, temperature=2.):
    nodes, weights=roots_hermitenorm(order)
    weights=weights/np.sqrt(2*np.pi)
    keep=weights>1e-25
    nodes,weights=nodes[keep],weights[keep]
    x,y=np.meshgrid(nodes+delta[0],nodes+delta[1],indexing='ij')
    w=np.column_stack((x.ravel(),y.ravel()))
    iw=np.outer(weights,weights).ravel()
    A=target_matrix(chi,kappa); truth=A@delta
    out=aggregate(w,chi,kappa,temperature)
    ab=aggregate(w,chi,kappa,temperature,omit_feedback_for_weights=True)
    risk=lambda x: float(iw@np.sum((x-truth)**2,axis=-1))
    cand_risks=np.einsum('i,ij->j',iw,np.sum((out['f']-truth[None,None,:])**2,axis=-1))
    cand_score=iw@out['score']
    hardrisk,harderr=hard_gaussian_risk(delta,chi,kappa)
    return dict(delta1=float(delta[0]),delta2=float(delta[1]),chi=float(chi),kappa=float(kappa),
                order=order,temperature=temperature,fa_risk=risk(out['fbar']),
                fa_sure=float(iw@out['sure']),frozen_sure=float(iw@out['frozen_sure']),
                correction=float(iw@out['tuning_correction']),
                naive_weighted_candidate_sure=float(iw@np.sum(out['alpha']*out['score'],axis=-1)),
                ablation_risk=risk(ab['fbar']),ablation_sure=float(iw@ab['sure']),
                oracle_risk=float(np.min(cand_risks)),
                soft_oracle=float(-temperature*logsumexp(-cand_risks/temperature)+temperature*np.log(len(cand_risks))),
                soft_oracle_bound=float(-temperature*logsumexp(-cand_risks/temperature)+temperature*np.log(len(cand_risks))+iw@out['gradient_variance']/temperature**2),
                akaike_limit_risk=float(cand_risks[4]),wide_risk=float(cand_risks[-1]),
                hard_aic_risk=float(hardrisk),hard_integration_error_estimate=float(harderr),
                candidate_unbiasedness_error=float(np.max(np.abs(cand_risks-cand_score))),
                candidate_risks=cand_risks.tolist(),
                gradient_cost=float(iw@out['gradient_variance']/temperature**2),
                oracle_bound=float(np.min(cand_risks)+temperature*np.log(len(default_rules()))+
                                   iw@out['gradient_variance']/temperature**2),
                upper_risk_score=float(iw@out['risk_upper_score']),
                mean_alpha=(iw@out['alpha']).tolist())


def main():
    p=argparse.ArgumentParser(); p.add_argument('--order',type=int,default=1024)
    p.add_argument('--output',default='results'); args=p.parse_args()
    outdir=Path(args.output);outdir.mkdir(parents=True,exist_ok=True)
    geom=population_geometry(.7,.3)
    cases=[((0.,1.),0.,.5),((0.,0.),.6,.5),((0.,1.),.6,.5),
           ((1.,1.),.6,.5),((2.,3.),.6,.5),((0.,1.),1.2,.5)]
    cases += [(d,geom['chi'],geom['kappa']) for d in ((0.,0.),(0.,1.),(1.,1.),(2.,3.))]
    rows=[]; convergence=[]
    for delta,chi,kappa in cases:
        low=quadrature(np.array(delta),chi,kappa,args.order//2)
        hi=quadrature(np.array(delta),chi,kappa,args.order)
        error=max(abs(hi[k]-low[k]) for k in ('fa_risk','fa_sure','ablation_risk','ablation_sure','gradient_cost'))
        hi['smooth_quadrature_change']=error
        rows.append(hi)
        print(f"delta={delta} chi={chi:.4f} kappa={kappa:.4f} FA={hi['fa_risk']:.9f} "
              f"SURE={hi['fa_sure']:.9f} frozen={hi['frozen_sure']:.9f} "
              f"abl={hi['ablation_risk']:.9f} change={error:.2g}",flush=True)
    (outdir/'gaussian_results.json').write_text(json.dumps(rows,indent=2))
    scalar_keys=[k for k,v in rows[0].items() if not isinstance(v,list)]
    with (outdir/'gaussian_results.csv').open('w',newline='') as f:
        dw=csv.DictWriter(f,fieldnames=scalar_keys);dw.writeheader()
        dw.writerows([{k:r[k] for k in scalar_keys} for r in rows])
    print('Hard-AIC is an external comparator, integrated separately across its thresholds.')

if __name__=='__main__': main()
