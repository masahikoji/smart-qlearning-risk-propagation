"""Deterministic derivative, risk-identity, and implementation checks."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from feedback_tuning import (smooth_map, candidate_maps, aggregate, target_matrix,
    default_rules, fit_smart, simulate_smart, population_geometry)

class Tests(unittest.TestCase):
    def test_scalar_derivatives(self):
        x=np.linspace(-7,7,141); eps=1e-5
        for p in (.5,1.,2.,None):
            f,d,dd=smooth_map(x,p)
            fp,dp,_=smooth_map(x+eps,p); fm,dm,_=smooth_map(x-eps,p)
            self.assertLess(np.max(np.abs((fp-fm)/(2*eps)-d)),2e-8)
            self.assertLess(np.max(np.abs((dp-dm)/(2*eps)-dd)),2e-8)

    def test_candidate_derivatives(self):
        rng=np.random.default_rng(412)
        w=rng.normal(size=(100,2))*2
        eps=1e-5
        for chi,kappa in ((0,0),(.4,.7),(-.8,.2),(1.5,1.)):
            v=candidate_maps(w,chi,kappa)
            A=target_matrix(chi,kappa)
            direct=np.sum((v['f']-(w@A.T)[:,None,:])**2,axis=-1)+2*np.einsum('ab,...jab->...j',A,v['jac'])-np.sum(A*A)
            self.assertLess(np.max(np.abs(direct-v['score'])),1e-12)
            for k in range(2):
                e=np.zeros_like(w); e[:,k]=eps
                vp,vm=candidate_maps(w+e,chi,kappa),candidate_maps(w-e,chi,kappa)
                self.assertLess(np.max(np.abs((vp['f']-vm['f'])/(2*eps)-v['jac'][...,k])),2e-7)
                self.assertLess(np.max(np.abs((vp['score']-vm['score'])/(2*eps)-v['grad_score'][...,k])),2e-6)
                self.assertLess(np.max(np.abs((vp['feedback']-vm['feedback'])/(2*eps)-v['grad_feedback'][...,k])),2e-6)

    def test_aggregate_identity_and_oracle_bound(self):
        rng=np.random.default_rng(821)
        w=rng.normal(size=(500,2))*2
        eps=1e-5
        for chi,kappa in ((0,0),(.4,.7),(-.8,.2),(1.5,1.)):
            A=target_matrix(chi,kappa)
            for temp in (.5,2.,8.):
                for ablation in (False,True):
                    v=aggregate(w,chi,kappa,temp,omit_feedback_for_weights=ablation)
                    direct=np.sum((v['fbar']-w@A.T)**2,axis=-1)+2*np.einsum('ab,...ab->...',A,v['jacbar'])-np.sum(A*A)
                    self.assertLess(np.max(np.abs(v['sure']-direct)),2e-12)
                    if not ablation:
                        self.assertTrue(np.all(v['sure']<=v['risk_upper_score']+1e-11))
                    for k in range(2):
                        e=np.zeros_like(w); e[:,k]=eps
                        vp=aggregate(w+e,chi,kappa,temp,omit_feedback_for_weights=ablation)
                        vm=aggregate(w-e,chi,kappa,temp,omit_feedback_for_weights=ablation)
                        self.assertLess(np.max(np.abs((vp['fbar']-vm['fbar'])/(2*eps)-v['jacbar'][...,k])),4e-6)

    def test_no_feedback_ablation(self):
        w=np.random.default_rng(23).normal(size=(100,2))
        a=aggregate(w,0,.7); b=aggregate(w,0,.7,omit_feedback_for_weights=True)
        np.testing.assert_array_equal(a['alpha'],b['alpha'])

    def test_observable_scores(self):
        rng=np.random.default_rng(82)
        for n in (250,1000):
            data,truth=simulate_smart(n,1.,2.,rng)
            out=fit_smart(data)
            self.assertFalse(out['stabilized'])
            self.assertEqual(out['actions_differ_from_reference'],0)
            self.assertLess(out['response_score_identity_error'],2e-12)
            self.assertAlmostEqual(float(out['alpha'].sum()),1,places=13)
            self.assertTrue(np.isfinite(out['coefficients1']).all())

    def test_logrss_matches_independent_akaike_fit(self):
        from scipy.special import expit
        data, truth=simulate_smart(300,1.,2.,np.random.default_rng(219))
        fit=fit_smart(data,criterion='logrss')
        x,a,z,d,y1,y2=(data[k] for k in ('X1','A1','X2','A2','Y1','Y2'))
        n=len(x)
        n2=np.column_stack((np.ones(n),z,a,d,d*z)); w2=np.column_stack((n2,d*z*z))
        c20=np.linalg.lstsq(n2,y2,rcond=None)[0];c21=np.linalg.lstsq(w2,y2,rcond=None)[0]
        r20=y2-n2@c20;r21=y2-w2@c21
        weight2=expit(n*np.log((r20@r20)/(r21@r21))/2-1)
        c2=np.r_[c20,0]*(1-weight2)+c21*weight2
        value2=c2[0]+c2[1]*z+c2[2]*a+np.abs(c2[3]+c2[4]*z+c2[5]*z*z)
        response=y1+value2
        n1=np.column_stack((np.ones(n),x,a,a*x));w1=np.column_stack((n1,x*x))
        c10=np.linalg.lstsq(n1,response,rcond=None)[0];c11=np.linalg.lstsq(w1,response,rcond=None)[0]
        r10=response-n1@c10;r11=response-w1@c11
        weight1=expit(n*np.log((r10@r10)/(r11@r11))/2-1)
        c1=np.r_[c10,0]*(1-weight1)+c11*weight1
        np.testing.assert_allclose(fit['candidate_coefficients1'][:,4],c1,atol=2e-12)
        np.testing.assert_allclose(fit['candidate_coefficients2'][:,4],c2,atol=2e-12)
        self.assertLess(fit['response_score_identity_error'],2e-12)

    def test_degenerate_design_guard(self):
        data,truth=simulate_smart(40,0,0,np.random.default_rng(71))
        data['A1'][:]=1
        out=fit_smart(data)
        self.assertTrue(out['stabilized'])
        np.testing.assert_array_equal(out['coefficients1'],np.zeros(5))

    def test_input_validation(self):
        with self.assertRaises(ValueError):
            aggregate([0,1],.5,.6,0)
        with self.assertRaises(ValueError):
            aggregate([0,1],.5,.6,prior=[1,0])
        with self.assertRaises(ValueError):
            smooth_map([0,1],float('inf'))

if __name__=='__main__':
    unittest.main(verbosity=2)
