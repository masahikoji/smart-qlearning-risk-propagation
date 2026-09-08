import math
import numpy as np
from scipy.integrate import quad

TOL=1e-12

def qv(x,a0,a1,a2): return a0+a1*x+a2*x*x

def qanti(x,a0,a1,a2): return a0*x+0.5*a1*x*x+(a2/3.0)*x**3

def roots(a0,a1,a2,tol=TOL):
    if abs(a2)<tol:
        if abs(a1)<tol: return []
        return [-a0/a1]
    disc=a1*a1-4*a2*a0
    if disc < -tol: return []
    if abs(disc)<=tol: return [-a1/(2*a2)]
    d=math.sqrt(max(0,disc))
    return sorted([(-a1-d)/(2*a2),(-a1+d)/(2*a2)])

def mean_quad(c,h,a0,a1,a2):
    return a0+a1*c+a2*(c*c+h*h/3)

def mean_abs(c,h,a0,a1,a2):
    lo,hi=c-h,c+h
    rs=roots(a0,a1,a2)
    cuts=[-math.inf]+rs+[math.inf]
    val=0.0
    for x0,x1 in zip(cuts[:-1],cuts[1:]):
        l=max(lo,x0); u=min(hi,x1)
        if not u>l: continue
        mid=(l+u)/2
        sg=np.sign(qv(mid,a0,a1,a2)) or 1
        val += max(0.0,sg*(qanti(u,a0,a1,a2)-qanti(l,a0,a1,a2)))
    return val/(2*h)

def mean_signed(c,h,chat,cstar,tie=1):
    if max(abs(x) for x in cstar)<TOL:
        return tie*mean_quad(c,h,*chat)
    lo,hi=c-h,c+h
    rs=roots(*cstar)
    cuts=[-math.inf]+rs+[math.inf]
    val=0.0
    for x0,x1 in zip(cuts[:-1],cuts[1:]):
        l=max(lo,x0); u=min(hi,x1)
        if not u>l: continue
        if math.isfinite(x0) and math.isfinite(x1): mid=(x0+x1)/2
        elif not math.isfinite(x0) and math.isfinite(x1): mid=x1-1
        elif math.isfinite(x0) and not math.isfinite(x1): mid=x0+1
        else: mid=0
        sg=np.sign(qv(mid,*cstar)) or tie
        val += sg*(qanti(u,*chat)-qanti(l,*chat))
    return val/(2*h)

def numeric_signed(c,h,chat,cstar,tie=1):
    f=lambda x: (tie if abs(qv(x,*cstar))<1e-14 else np.sign(qv(x,*cstar)))*qv(x,*chat)
    pts=[r for r in roots(*cstar) if c-h<r<c+h]
    val=quad(f,c-h,c+h,points=pts,epsabs=1e-12,epsrel=1e-12,limit=200)[0]/(2*h)
    return val

rng=np.random.default_rng(260907)
errs=[]
for _ in range(1000):
    c=rng.uniform(-1,1); h=rng.uniform(.05,1)
    chat=tuple(rng.normal(size=3))
    if rng.uniform()<.1:
        cstar=(0.,0.,0.)
    else:
        cstar=tuple(rng.normal(size=3))
    a=mean_signed(c,h,chat,cstar)
    b=numeric_signed(c,h,chat,cstar)
    errs.append(abs(a-b))
print('max signed-quadratic error',max(errs))

# Validate terminal zeta identity for random coefficients/histories using a K=3-like terminal model.
# p_true = beta0 + beta1*x + bp1*a1 + bp2*a2; C_true quadratic.
def mean_poly3(c,h,coef):
    c0,c1,c2,c3=coef
    return c0+c1*c+c2*(c*c+h*h/3)+c3*(c**3+c*h*h)

def analytic_zeta(c,h,a1,a2,bhat,btrue,cstar):
    # bhat: int,X,X2,X3,past1,past2,A,AX,AX2
    p0=(bhat[0]-btrue[0])+(bhat[4]-btrue[2])*a1+(bhat[5]-btrue[3])*a2
    ep=p0+mean_poly3(c,h,(0,bhat[1]-btrue[1],bhat[2],bhat[3]))
    chat=(bhat[6],bhat[7],bhat[8])
    return ep+mean_signed(c,h,chat,cstar)-mean_abs(c,h,*cstar)

def numeric_zeta(c,h,a1,a2,bhat,btrue,cstar):
    def f(x):
        ct=qv(x,*cstar)
        s=1 if ct>=0 else -1
        ph=bhat[0]+bhat[1]*x+bhat[2]*x*x+bhat[3]*x**3+bhat[4]*a1+bhat[5]*a2
        ch=qv(x,bhat[6],bhat[7],bhat[8])
        pt=btrue[0]+btrue[1]*x+btrue[2]*a1+btrue[3]*a2
        qt=pt+s*ct
        qh=ph+s*ch
        return qh-qt
    pts=[r for r in roots(*cstar) if c-h<r<c+h]
    return quad(f,c-h,c+h,points=pts,epsabs=1e-12,epsrel=1e-12,limit=200)[0]/(2*h)

errs=[]
for _ in range(1000):
    c=rng.uniform(-1,1); h=rng.uniform(.1,1)
    a1=rng.choice([-1,1]); a2=rng.choice([-1,1])
    bhat=rng.normal(size=9)
    btrue=rng.normal(size=4)
    if rng.uniform()<.1: cstar=(0.,0.,0.)
    else: cstar=tuple(rng.normal(size=3))
    a=analytic_zeta(c,h,a1,a2,bhat,btrue,cstar)
    b=numeric_zeta(c,h,a1,a2,bhat,btrue,cstar)
    errs.append(abs(a-b))
print('max terminal-zeta error',max(errs))
assert max(errs)<1e-9
