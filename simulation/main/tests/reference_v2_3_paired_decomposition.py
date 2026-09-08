import numpy as np

# Independent check of the v2.3 decomposition idea in a two-stage setting.
# Local zeta_1 is Rao-Blackwellized exactly; zeta_2 and rho_2 are evaluated on
# two conditionally independent descendants and paired cross-products estimate
# products of their transported conditional means.
rng = np.random.default_rng(20260907)
rho, xi = 0.5, 0.2
h = 1-rho
N = 500_000
x1 = rng.uniform(-1,1,N)
a1 = np.where(rng.random(N)<0.5,-1.0,1.0)

# Arbitrary fitted Q-functions. True stage-2 optimal action is designated +1.
def q1hat(x,a):
    return 0.2 + 0.1*x + 0.05*x*x + a*(0.8+0.2*x+0.1*x*x)
def p2(x,a1):
    return -0.1 + 0.15*x + 0.03*x*x + 0.07*a1
def c2hat(x):
    return -0.05 + 0.30*x + 0.08*x*x
def q2hat(x,a1,a2):
    return p2(x,a1) + a2*c2hat(x)
def mu1(x,a):
    return 0.1 + 0.2*x + a*(0.7+0.1*x)
def mu2(x,a1,a2):
    return -0.05 + 0.1*x + 0.05*a1 + a2*(0.6+0.08*x)

def roots(a0,a1,a2):
    d=a1*a1-4*a2*a0
    if abs(a2)<1e-14: return np.array([-a0/a1]) if abs(a1)>1e-14 else np.array([])
    if d<0:return np.array([])
    s=np.sqrt(max(d,0));return np.sort(np.array([(-a1-s)/(2*a2),(-a1+s)/(2*a2)]))
def anti(x,a0,a1,a2):return a0*x+a1*x*x/2+a2*x**3/3
def mean_abs_quad(center,h,a0,a1,a2):
    center=np.asarray(center); lo=center-h; hi=center+h
    rr=roots(a0,a1,a2); cuts=np.r_[-np.inf,rr,np.inf]
    out=np.zeros_like(center)
    for L,U in zip(cuts[:-1],cuts[1:]):
        l=np.maximum(lo,L); u=np.minimum(hi,U); ok=u>l
        if not np.any(ok):continue
        mid=(l[ok]+u[ok])/2
        sg=np.sign(a0+a1*mid+a2*mid*mid);sg[sg==0]=1
        out[ok]+=sg*(anti(u[ok],a0,a1,a2)-anti(l[ok],a0,a1,a2))
    return out/(2*h)

def mean_quad(center,h,a0,a1,a2):
    return a0+a1*center+a2*(center*center+h*h/3)

center=rho*x1+xi*a1
# p2 coefficients depend on a1 through the intercept only.
evhat2 = (-0.1 + 0.07*a1) + 0.15*center + 0.03*(center*center+h*h/3) + \
          mean_abs_quad(center,h,-0.05,0.30,0.08)
zeta1 = q1hat(x1,a1)-mu1(x1,a1)-evhat2

# Exact transported stage-2 source and nonlinear remainder.
# zeta2(x,+1)=q2hat-mu2 is quadratic in x.
# Coefficients: intercept (-.1+.07a1)-(-.05+.05a1)-.6, linear .15+.30-(.1+.08), quadratic .03+.08
z20 = (-0.1+0.07*a1-0.05)-(-0.05+0.05*a1)-0.6
z21 = 0.15+0.30-0.1-0.08
z22 = 0.03+0.08
L2 = z20 + z21*center + z22*(center*center+h*h/3)
mean_chat = mean_quad(center,h,-0.05,0.30,0.08)
N2 = mean_abs_quad(center,h,-0.05,0.30,0.08)-mean_chat

# Two independent descendants.
def desc(u):
    x2=center+h*u
    z2=q2hat(x2,a1,1.0)-mu2(x2,a1,1.0)
    rr=np.maximum(np.abs(c2hat(x2))-c2hat(x2),0)
    return z2,rr
z2a,r2a=desc(rng.uniform(-1,1,N)); z2b,r2b=desc(rng.uniform(-1,1,N))

exact={
    'diag_s1':np.mean(zeta1*zeta1),
    'diag_s2':np.mean(L2*L2),
    'cov_s1_s2':np.mean(zeta1*L2),
    'nonlinear_sq':np.mean(N2*N2),
    'linear_nonlinear_cross2':2*np.mean((zeta1+L2)*N2),
    'full':np.mean((zeta1+L2+N2)**2),
}
paired={
    'diag_s1':np.mean(zeta1*zeta1),
    'diag_s2':np.mean(z2a*z2b),
    'cov_s1_s2':0.5*np.mean(zeta1*z2b+z2a*zeta1),
    'nonlinear_sq':np.mean(r2a*r2b),
    'linear_nonlinear_cross2':np.mean((zeta1+z2a)*r2b+r2a*(zeta1+z2b)),
    'full':np.mean((zeta1+z2a+r2a)*(zeta1+z2b+r2b)),
}
for k in exact:
    print(f"{k}: exact={exact[k]:.8f} paired={paired[k]:.8f} diff={paired[k]-exact[k]:.3e}")
print('max_abs_diff=',max(abs(paired[k]-exact[k]) for k in exact))
