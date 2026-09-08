import numpy as np
from numpy.polynomial.legendre import leggauss

# Independent check of the finite-horizon performance-difference estimator used
# in v2.4.  This is deliberately written without using any R implementation.

rng = np.random.default_rng(260907)

rho = 0.30
xi = 0.00
h = 1.0 - rho
n = 250
# Representative two-stage smooth-boundary construction.
beta0 = np.array([0.25, 0.25])
beta1 = np.array([0.50, 0.50])
beta_past2 = 0.20
tau0 = np.array([1.35, 0.0])
tau1 = np.array([0.20, 1.0])
# Approximate calibrated local coefficient from the frozen v2 scenarios.
c = np.array([0.0, 4.89924 / np.sqrt(n)])

nodes, weights = leggauss(192)
pw = weights / 2.0

def r1(x, a):
    return beta0[0] + beta1[0]*x + a*(tau0[0] + tau1[0]*x + c[0]*x*x)

def r2(x, a1, a2):
    return beta0[1] + beta1[1]*x + beta_past2*a1 + a2*(tau0[1] + tau1[1]*x + c[1]*x*x)

def q2(x, a1, a2):
    return r2(x, a1, a2)

def v2(x, a1):
    return np.maximum(q2(x,a1,1), q2(x,a1,-1))

def q1(x1, a1):
    cen = rho*x1 + xi*a1
    x2 = cen[:,None] + h*nodes[None,:]
    a1m = np.broadcast_to(a1[:,None], x2.shape)
    ev2 = np.sum(v2(x2, a1m)*pw[None,:], axis=1)
    return r1(x1,a1) + ev2

def v1(x1):
    return np.maximum(q1(x1,np.ones_like(x1)), q1(x1,-np.ones_like(x1)))

# A deliberately imperfect deterministic policy.
def d1(x):
    return np.where(x > -0.15, 1.0, -1.0)

def d2(x2, a1):
    return np.where(x2 > 0.12, 1.0, -1.0)

N = 500_000
x1 = rng.uniform(-1,1,N)
u2 = rng.uniform(-1,1,N)
a1 = d1(x1)
x2 = rho*x1 + h*u2 + xi*a1
a2 = d2(x2,a1)

# Performance-difference form.
q1p = q1(x1,np.ones_like(x1))
q1m = q1(x1,-np.ones_like(x1))
gap1 = np.maximum(q1p,q1m) - np.where(a1>0,q1p,q1m)
q2p = q2(x2,a1,1)
q2m = q2(x2,a1,-1)
gap2 = np.maximum(q2p,q2m) - np.where(a2>0,q2p,q2m)
regret_gap = np.mean(gap1 + gap2)

# Direct value difference using the same initial states/transitions for the
# fitted policy and an independent high-accuracy truth evaluation for V1*.
value_star = np.mean(np.maximum(q1p,q1m))
value_policy = np.mean(r1(x1,a1) + r2(x2,a1,a2))
regret_direct = value_star - value_policy

# The two Monte Carlo estimators use the same finite sample and are not exactly
# identical sample-by-sample because Q1 integrates over U2, while the rollout
# uses realized U2.  Their difference should be well within the rollout MCSE.
diff = regret_gap - regret_direct
rollout_term = np.maximum(q1p,q1m) - (r1(x1,a1) + r2(x2,a1,a2))
se_direct = np.std(rollout_term, ddof=1)/np.sqrt(N)

# Oracle policy accumulated gap is algebraically zero at every stage.
a1_star = np.where(q1p >= q1m, 1.0, -1.0)
x2_star = rho*x1 + h*u2 + xi*a1_star
q2ps = q2(x2_star,a1_star,1)
q2ms = q2(x2_star,a1_star,-1)
a2_star = np.where(q2ps >= q2ms, 1.0, -1.0)
oracle_gap = (np.maximum(q1p,q1m) - np.where(a1_star>0,q1p,q1m)) + \
             (np.maximum(q2ps,q2ms) - np.where(a2_star>0,q2ps,q2ms))

print(f"gap_regret={regret_gap:.10f}")
print(f"direct_regret={regret_direct:.10f}")
print(f"difference={diff:.10f}")
print(f"direct_mcse={se_direct:.10f}")
print(f"max_oracle_gap={np.max(np.abs(oracle_gap)):.3e}")

assert regret_gap >= 0
assert np.max(np.abs(oracle_gap)) < 1e-12
assert abs(diff) < 5*se_direct
print("all reference v2.4 regret checks passed")
