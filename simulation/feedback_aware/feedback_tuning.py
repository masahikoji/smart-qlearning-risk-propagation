"""Feedback-aware soft tuning for smooth, complete two-stage procedures.

This is a reference implementation accompanying the manuscript and Supplementary Material.
The Gaussian risk identities are exact. SMART guarantees are first-order,
under the E3-prime conditions specified in the Supplementary Material. This is not an estimator
of clinical policy value and does not provide universal dominance over AIC.

Only NumPy and SciPy are required. No unknown local parameter is an input to
fit_smart or to aggregate. The temperature and library must be prespecified.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Mapping, Sequence
import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import expit, logsumexp

FloatArray = NDArray[np.float64]

@dataclass(frozen=True)
class Rule:
    upstream: float | None
    downstream: float | None
    name: str


def default_rules() -> tuple[Rule, ...]:
    """Nine smooth penalty pairs and the all-wide reference, fixed in advance."""
    vals = (0.5, 1.0, 2.0)
    return tuple(Rule(a, b, f"p1={a:g},p2={b:g}")
                 for a in vals for b in vals) + (Rule(None, None, "wide"),)


def smooth_map(x: ArrayLike, penalty: float | None) -> tuple[FloatArray, FloatArray, FloatArray]:
    """g(x), g'(x), g''(x); None denotes the identity, not hard selection."""
    x = np.asarray(x, dtype=float)
    if penalty is None:
        return x.copy(), np.ones_like(x), np.zeros_like(x)
    if not np.isfinite(penalty):
        raise ValueError("A penalty must be finite, or None for the identity.")
    w = expit(0.5 * x * x - float(penalty))
    # expit(-u) avoids catastrophic cancellation in 1-expit(u).
    v = expit(float(penalty) - 0.5 * x * x)
    wv = w * v
    return x*w, w + x*x*wv, 3*x*wv + x*x*x*wv*(1-2*w)


def target_matrix(chi: float, kappa: float) -> FloatArray:
    if not np.isfinite(chi) or not np.isfinite(kappa):
        raise ValueError("Geometry must be finite.")
    return np.array([[1., chi], [0., kappa]])


def candidate_maps(w: ArrayLike, chi: float, kappa: float,
                   rules: Sequence[Rule] | None = None) -> dict[str, FloatArray]:
    """Evaluate the local experiment at one point or a batch of points.

    w has final dimension 2. Output has an additional library index before
    output/derivative dimensions. grad_score differentiates with respect to w.
    """
    w = np.asarray(w, dtype=float)
    if w.shape[-1:] != (2,) or not np.isfinite(w).all():
        raise ValueError("w must be finite with final dimension 2.")
    rules = tuple(default_rules() if rules is None else rules)
    if not rules:
        raise ValueError("The library cannot be empty.")
    shape = w.shape[:-1] + (len(rules),)
    f = np.empty(shape+(2,))
    jac = np.zeros(shape+(2,2))
    score = np.empty(shape)
    gs = np.empty(shape+(2,))
    feedback = np.empty(shape)
    gfeedback = np.empty(shape+(2,))
    A = target_matrix(chi, kappa)
    aw = np.einsum('ab,...b->...a', A, w)
    for j, rule in enumerate(rules):
        u, b, bb = smooth_map(w[..., 1], rule.downstream)
        q = w[..., 0] + chi*u
        v, a, aa = smooth_map(q, rule.upstream)
        r1 = v-aw[..., 0]
        h2 = u-w[..., 1]
        f[..., j, 0] = v
        f[..., j, 1] = kappa*u
        jac[..., j, 0, 0] = a
        jac[..., j, 0, 1] = chi*a*b
        jac[..., j, 1, 1] = kappa*b
        score[..., j] = (r1*r1 + kappa*kappa*h2*h2
                          + 2*a*(1+chi*chi*b)+2*kappa*kappa*b
                          - (1+chi*chi+kappa*kappa))
        gs[..., j, 0] = 2*r1*(a-1)+2*aa*(1+chi*chi*b)
        gs[..., j, 1] = (2*r1*chi*(a*b-1)+2*kappa*kappa*h2*(b-1)
                          +2*aa*chi*b*(1+chi*chi*b)
                          +2*(chi*chi*a+kappa*kappa)*bb)
        feedback[..., j] = 2*chi*chi*a*b
        gfeedback[..., j, 0] = 2*chi*chi*aa*b
        gfeedback[..., j, 1] = 2*chi*chi*(chi*aa*b*b+a*bb)
    return dict(f=f, jac=jac, score=score, grad_score=gs,
                feedback=feedback, grad_feedback=gfeedback)


def aggregate(w: ArrayLike, chi: float, kappa: float, temperature: float = 2.,
              rules: Sequence[Rule] | None = None, prior: ArrayLike | None = None,
              omit_feedback_for_weights: bool = False) -> dict[str, FloatArray]:
    """Entropy-regularized soft tuning, including its exact Gaussian risk score.

    omit_feedback_for_weights is an explicitly misspecified ablation. Its
    returned sure is nonetheless the CORRECT risk score for that ablation:
    it differentiates the weight criterion that was actually used. The oracle
    bound proved for the main method is not claimed for this ablation.
    """
    w = np.asarray(w, dtype=float)
    rules = tuple(default_rules() if rules is None else rules)
    J = len(rules)
    if not np.isfinite(temperature) or temperature <= 0:
        raise ValueError("temperature must be positive and finite.")
    if prior is None:
        prior = np.full(J, 1/J)
    prior = np.asarray(prior, dtype=float)
    if prior.shape != (J,) or not np.isfinite(prior).all() or np.any(prior <= 0):
        raise ValueError("prior must have one strictly positive finite entry per rule.")
    prior = prior / prior.sum()
    vals = candidate_maps(w, chi, kappa, rules)
    f, jac, score, gs = (vals[k] for k in ('f','jac','score','grad_score'))
    weight_score = score.copy()
    weight_grad = gs.copy()
    if omit_feedback_for_weights:
        weight_score -= vals['feedback']
        weight_grad -= vals['grad_feedback']
    logw = np.log(prior)-weight_score/temperature
    lz = logsumexp(logw, axis=-1)
    alpha = np.exp(logw-lz[..., None])
    fbar = np.sum(alpha[..., None]*f, axis=-2)
    diff = f-fbar[..., None, :]
    gradbar = np.sum(alpha[..., None]*weight_grad, axis=-2)
    gradalpha = -alpha[..., None]*(weight_grad-gradbar[..., None, :])/temperature
    jacbar = (np.sum(alpha[..., None, None]*jac, axis=-3)
              + np.einsum('...ja,...jb->...ab', f, gradalpha))
    A = target_matrix(chi, kappa)
    Ag = np.einsum('ab,...jb->...ja', A, weight_grad)
    correction = -2/temperature*np.sum(alpha*np.sum(diff*Ag, axis=-1), axis=-1)
    dispersion = np.sum(alpha*np.sum(diff*diff, axis=-1), axis=-1)
    frozen = np.sum(alpha*score, axis=-1)-dispersion
    sure = frozen+correction
    Agbar = np.sum(alpha[..., None]*Ag, axis=-2)
    grad_variance = np.sum(alpha*np.sum((Ag-Agbar[..., None, :])**2, axis=-1), axis=-1)
    log_partition_risk = -temperature*lz
    upper = log_partition_risk+grad_variance/(temperature*temperature)
    return dict(**vals, alpha=alpha, fbar=fbar, jacbar=jacbar, sure=sure,
                frozen_sure=frozen, tuning_correction=correction,
                dispersion=dispersion, gradient_variance=grad_variance,
                risk_upper_score=upper, log_partition_risk=log_partition_risk,
                weight_score=weight_score)


class StabilizationError(ValueError):
    """The prespecified Gram-matrix safeguard was triggered."""


def _ols_map(x: FloatArray, gram_floor: float) -> tuple[FloatArray, FloatArray]:
    n = x.shape[0]
    G = x.T@x/n
    if np.linalg.eigvalsh(G)[0] < gram_floor:
        raise StabilizationError("A required empirical Gram matrix failed the safeguard.")
    inv = np.linalg.inv(G)
    return inv, inv@x.T/n


def fit_smart(data: Mapping[str, ArrayLike], temperature: float = 2.,
              rules: Sequence[Rule] | None = None, prior: ArrayLike | None = None,
              gram_floor: float = 1e-8, scale_floor: float = 1e-6,
              geometry_cap: float = 20.,
              omit_feedback_for_weights: bool = False,
              criterion: str = 'wald') -> dict:
    """Fit the E3-prime procedure using observed data only.

    Required arrays: X1, A1, X2, A2, Y1, Y2, all of the same length.
    Coefficient order: stage 1 (1,X1,A1,A1*X1,X1^2), stage 2
    (1,X2,A1,A2,A2*X2,A2*X2^2). Q1_FA is the mixture of COMPLETE
    recursive fits, not a refit after first mixing the stage-2 Q functions.

    criterion='logrss' uses the actual n*log(RSS0/RSS1) in each recursive
    candidate and agrees with ordinary Akaike weighting at penalty pair (1,1).
    The tuning score remains a first-order Gaussian score in either option.
    A failure of the fixed Gram check returns a zero-coefficient fallback
    with stabilized=True. Positive variance floors and geometry caps keep
    the algorithm defined outside the regular event; asymptotic claims
    require that the true scales/geometry are in the interior of these guards.
    """
    rules = tuple(default_rules() if rules is None else rules)
    if criterion not in ('wald', 'logrss'):
        raise ValueError("criterion must be 'wald' or 'logrss'.")
    if not rules:
        raise ValueError("Empty library.")
    arrays = {key:np.asarray(data[key], dtype=float) for key in
              ('X1','A1','X2','A2','Y1','Y2')}
    n = len(arrays['X1'])
    if n < 12 or any(a.shape != (n,) for a in arrays.values()):
        raise ValueError("All data arrays must be one-dimensional, equal length, n>=12.")
    if any(not np.isfinite(a).all() for a in arrays.values()):
        raise ValueError("Missing/nonfinite data are not silently removed.")
    if any(not np.isin(arrays[a], [-1.,1.]).all() for a in ('A1','A2')):
        raise ValueError("Treatment variables must be coded -1/+1.")
    if gram_floor <= 0 or scale_floor <= 0 or geometry_cap <= 0:
        raise ValueError("Safeguards must be positive.")
    x,a,z,d,y1,y2 = (arrays[k] for k in ('X1','A1','X2','A2','Y1','Y2'))
    X10 = np.column_stack((np.ones(n),x,a,a*x))
    X20 = np.column_stack((np.ones(n),z,a,d,d*z))
    q1 = x*x
    q2 = d*z*z
    try:
        G1inv, L1 = _ols_map(X10, gram_floor)
        G2inv, L2 = _ols_map(X20, gram_floor)
        _ols_map(np.column_stack((X10,q1)), gram_floor)
        _ols_map(np.column_stack((X20,q2)), gram_floor)
        proj1, proj2 = L1@q1, L2@q2
        r1, r2 = q1-X10@proj1, q2-X20@proj2
        v1, v2 = r1@r1, r2@r2
        g1hat, g2hat = v1/n, v2/n
        if min(g1hat,g2hat) < gram_floor:
            raise StabilizationError("A residualised added block failed the safeguard.")
    except StabilizationError as exc:
        return dict(stabilized=True, reason=str(exc), coefficients1=np.zeros(5),
                    coefficients2=np.zeros(6), n=n)
    common2 = L2@y2
    c2 = (r2@y2)/v2
    beta2wide = np.r_[common2-proj2*c2,c2]
    resid2 = y2-X20@common2-r2*c2
    sig2 = max(float(np.sqrt(resid2@resid2/n)), scale_floor)
    w2 = np.sqrt(n*g2hat)*c2/sig2
    # Use actual maximisation for every complete recursive fit.
    def stage2_value(c):
        cc = common2[:,None]-proj2[:,None]*np.asarray(c)[None,:]
        prognostic = cc[0]+z[:,None]*cc[1]+a[:,None]*cc[2]
        contrast = cc[3]+z[:,None]*cc[4]+z[:,None]**2*np.asarray(c)[None,:]
        return prognostic+np.abs(contrast), np.where(contrast>=0,1.,-1.)
    wide_value, ref_action = stage2_value(np.array([c2]))
    ref_action = ref_action[:,0]
    ref_response = y1+wide_value[:,0]
    common1_ref = L1@ref_response
    c1ref = (r1@ref_response)/v1
    resid1 = ref_response-X10@common1_ref-r1*c1ref
    sig1 = max(float(np.sqrt(resid1@resid1/n)), scale_floor)
    eta1wide = np.sqrt(n*g1hat)*c1ref/sig1
    X20ref = np.column_stack((np.ones(n),z,a,ref_action,ref_action*z))
    vref = (ref_action*z*z-X20ref@proj2)/np.sqrt(g2hat)
    chi_raw = (sig2/sig1)*float((r1/np.sqrt(g1hat))@vref/n)
    narrow_transport = X10@(L1@vref)
    kappa_raw = (sig2/sig1)*float(np.sqrt(narrow_transport@narrow_transport/n))
    chi = float(np.clip(chi_raw, -geometry_cap, geometry_cap))
    kappa = float(np.clip(kappa_raw, 0., geometry_cap))
    w1 = eta1wide-chi*w2
    w = np.array([w1,w2])
    local = aggregate(w, chi, kappa, temperature, rules, prior,
                      omit_feedback_for_weights=omit_feedback_for_weights)
    lambda2 = float(n*np.log1p(c2*c2*v2/max(float(resid2@resid2),n*scale_floor**2)))
    terminal_weights = np.array([1. if rule.downstream is None else
        expit((w2*w2 if criterion=='wald' else lambda2)/2-rule.downstream) for rule in rules])
    c2s = terminal_weights*c2
    responses, actions = stage2_value(c2s)
    responses += y1[:,None]
    common1s = L1@responses
    c1wides = (r1@responses)/v1
    eta1s = np.sqrt(n*g1hat)*c1wides/sig1
    stage1_resid = responses-X10@common1s-r1[:,None]*c1wides[None,:]
    rss1 = np.maximum(np.sum(stage1_resid*stage1_resid,axis=0),n*scale_floor**2)
    lambda1s = n*np.log1p(v1*c1wides*c1wides/rss1)
    upstream_weights = np.array([1. if rule.upstream is None else
        expit((eta1s[j]**2 if criterion=='wald' else lambda1s[j])/2-rule.upstream)
        for j,rule in enumerate(rules)])
    c1s = upstream_weights*c1wides
    coefs1 = np.vstack((common1s-proj1[:,None]*c1s[None,:],c1s))
    coefs2 = np.vstack((common2[:,None]-proj2[:,None]*c2s[None,:],c2s))
    alpha = local['alpha']
    B0 = X10.T@X20ref/n
    common_risk = X10.shape[1]*sig1*sig1+sig2*sig2*np.trace(G1inv@B0@G2inv@B0.T)
    theoretical_scores = w1+chi*terminal_weights*w2
    return dict(stabilized=False, n=n, criterion=criterion, coefficients1=coefs1@alpha,
                coefficients2=coefs2@alpha, candidate_coefficients1=coefs1,
                candidate_coefficients2=coefs2, alpha=alpha,
                w=w, chi=chi, kappa=kappa, sigma1=sig1, sigma2=sig2,
                common_risk_estimate=common_risk,
                scaled_risk_estimate=common_risk+sig1*sig1*float(local['sure']),
                local=local, g1=g1hat, g2=g2hat, terminal_weights=terminal_weights,
                upstream_weights=upstream_weights, lambda2=lambda2, lambda1s=lambda1s,
                response_score_identity_error=float(np.max(np.abs(eta1s-theoretical_scores))),
                nonpositive_reference_actions=int(np.sum(ref_action<0)),
                actions_differ_from_reference=int(np.sum(actions!=ref_action[:,None])),
                geometry_clipped=(chi!=chi_raw or kappa!=kappa_raw))


def population_geometry(rho: float, omega: float, sigma2: float=1.,
                        barsigma1: float=1.) -> dict[str, float]:
    su = 1-rho
    m2=(rho*rho+su*su)/3+omega*omega
    m4=(rho**4+su**4)/5+omega**4+2*rho*rho*su*su/3+2*(rho*rho+su*su)*omega*omega
    g2=m4-m2*m2
    if g2 <= 0:
        raise ValueError("Degenerate geometry.")
    chi=rho*rho*sigma2*np.sqrt(4/45)/(barsigma1*np.sqrt(g2))
    kappa=2*rho*omega*sigma2/(np.sqrt(3)*barsigma1*np.sqrt(g2))
    G1=np.diag([1.,1/3,1.,1/3])
    G2=np.diag([1.,m2,1.,1.,m2]); G2[1,2]=G2[2,1]=omega
    B=np.array([[1,0,0,1,0],[0,rho/3,0,0,rho/3],
                [0,omega,1,0,omega],[0,0,0,0,0]],dtype=float)
    rcommon=4*barsigma1**2+sigma2**2*np.trace(np.linalg.solve(G1,B)@np.linalg.solve(G2,B.T))
    return dict(chi=float(chi),kappa=float(kappa),g1=4/45,g2=float(g2),
                common_risk=float(rcommon),m2=float(m2))


def simulate_smart(n: int, delta1: float, delta2: float, rng: np.random.Generator,
                   rho: float=.7, omega: float=.3) -> tuple[dict, FloatArray]:
    """E3-prime DGP. Local parameters enter this simulator ONLY, never the fitter."""
    geom=population_geometry(rho,omega)
    b1=delta1/np.sqrt(geom['g1']); b2=delta2/np.sqrt(geom['g2'])
    c1,c2=b1/np.sqrt(n),b2/np.sqrt(n)
    beta1=np.array([.1,.2,.3,-.1]); beta20,beta21,beta22=1.,.4,.3
    tau20,tau21=5.,.5
    su=1-rho
    sigma1sq=1-(beta21+tau21)**2*su*su/3
    if sigma1sq <= 0:
        raise ValueError("The chosen parameters do not allow barsigma1=1.")
    # Guarantee the exact no-switch population target used by the evaluator.
    zmax=abs(rho)+abs(su)+abs(omega)
    if tau20-abs(tau21)*zmax-abs(c2)*zmax*zmax <= 0:
        raise ValueError("True treatment separation is not guaranteed for this DGP.")
    x=rng.uniform(-1,1,n); a=rng.choice([-1.,1.],n)
    z=rho*x+su*rng.uniform(-1,1,n)+omega*a
    d=rng.choice([-1.,1.],n)
    N1=np.column_stack((np.ones(n),x,a,a*x))
    y1=N1@beta1+c1*x*x+rng.normal(0,np.sqrt(sigma1sq),n)
    y2=beta20+beta21*z+beta22*a+d*(tau20+tau21*z+c2*z*z)+rng.normal(size=n)
    truth=np.array([beta1[0]+beta20+tau20+c2*(omega*omega+su*su/3),
                    beta1[1]+(beta21+tau21)*rho,
                    beta1[2]+beta22+(beta21+tau21)*omega,
                    beta1[3]+2*rho*omega*c2,c1+rho*rho*c2])
    return dict(X1=x,A1=a,X2=z,A2=d,Y1=y1,Y2=y2), truth


def stage1_risk(coef: ArrayLike, truth: ArrayLike) -> FloatArray:
    """Exact population MSE for the five-dimensional E3-prime stage-1 span."""
    G=np.diag([1.,1/3,1.,1/3,1/5]); G[0,4]=G[4,0]=1/3
    e=np.asarray(coef,dtype=float)-np.asarray(truth,dtype=float)
    return np.einsum('...i,ij,...j->...',e,G,e)
