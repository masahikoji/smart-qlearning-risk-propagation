"""Run one example or read observed E3-prime arrays from a CSV.

The CSV entry point does not validate the scientific model assumptions.
"""
import argparse
import json
import numpy as np
from feedback_tuning import fit_smart, simulate_smart, stage1_risk, default_rules

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--csv', default=None)
    parser.add_argument('--criterion', choices=('wald','logrss'), default='logrss')
    args=parser.parse_args()
    truth=None
    if args.csv:
        records=np.genfromtxt(args.csv,delimiter=',',names=True,dtype=float)
        required=('X1','A1','X2','A2','Y1','Y2')
        if records.dtype.names is None or not set(required).issubset(records.dtype.names):
            parser.error('CSV requires X1,A1,X2,A2,Y1,Y2 headers.')
        data={key:np.atleast_1d(records[key]) for key in required}
    else:
        data,truth=simulate_smart(1000,1.,1.,np.random.default_rng(20260909))
    fit=fit_smart(data,criterion=args.criterion)
    if fit['stabilized']:
        print(json.dumps({'stabilized':True,'reason':fit['reason']})); return
    result={
        'criterion':args.criterion,'n':fit['n'],
        'estimated_primitive_scores':fit['w'].tolist(),
        'estimated_chi':fit['chi'],'estimated_kappa':fit['kappa'],
        'temperature_fixed_in_advance':2.,
        'weights':{r.name:float(a) for r,a in zip(default_rules(),fit['alpha'])},
        'Q1_coefficients_1_x_a_ax_x2':fit['coefficients1'].tolist(),
        'asymptotic_scaled_risk_score':fit['scaled_risk_estimate'],
        'warning':'Risk score is first-order, not a confidence bound or an exactly unbiased finite-sample SMART score.'}
    if truth is not None:
        result['simulation_only_realized_scaled_prediction_loss']=float(fit['n']*stage1_risk(fit['coefficients1'],truth))
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
