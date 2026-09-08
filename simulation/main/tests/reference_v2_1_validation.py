import math
import numpy as np
from scipy.integrate import quad
from scipy.stats import norm


def moments(rho, xi, K):
    m2=[1/3]; m4=[1/5]; h=1-rho
    for _ in range(1,K):
        a2,a4=m2[-1],m4[-1]
        m2.append(rho**2*a2+h**2/3+xi**2)
        m4.append(rho**4*a4+h**4/5+xi**4+
                  6*rho**2*h**2*a2/3+6*rho**2*xi**2*a2+6*h**2*xi**2/3)
    return np.array(m2),np.array(m4)


def bar_delta(delta):
    def err(mode):
        def f(z):
            xi=delta+z
            w=(1.0 if xi*xi>2 else 0.0) if mode=='sel' else 1/(1+math.exp(1-xi*xi/2))
            return (w*xi-delta)**2*norm.pdf(z)
        return quad(f,-10,10,epsabs=1e-12,epsrel=1e-12,limit=2000)[0]
    return err('sel')-err('ma')


def roots_quad(a0,a1,a2):
    if abs(a2)<1e-14:
        return [] if abs(a1)<1e-14 else [-a0/a1]
    d=a1*a1-4*a2*a0
    if d<0:return []
    if abs(d)<1e-14:return [-a1/(2*a2)]
    s=math.sqrt(d)
    return sorted([(-a1-s)/(2*a2),(-a1+s)/(2*a2)])


def int_quad(x,a0,a1,a2):return a0*x+a1*x*x/2+a2*x**3/3

def exact_rho(center,h,ch,cs):
    roots=sorted(set(roots_quad(*ch)+roots_quad(*cs)))
    cuts=[-math.inf]+roots+[math.inf]
    total=0
    lo_all=center-h;hi_all=center+h
    for lo0,hi0 in zip(cuts[:-1],cuts[1:]):
        lo=max(lo_all,lo0);hi=min(hi_all,hi0)
        if hi<=lo:continue
        if math.isfinite(lo0) and math.isfinite(hi0):mid=(lo0+hi0)/2
        elif math.isfinite(hi0):mid=hi0-1
        elif math.isfinite(lo0):mid=lo0+1
        else:mid=0
        vh=ch[0]+ch[1]*mid+ch[2]*mid*mid
        vs=cs[0]+cs[1]*mid+cs[2]*mid*mid
        sh=1 if vh>=0 else -1; ss=1 if vs>=0 else -1
        total+=(sh-ss)*(int_quad(hi,*ch)-int_quad(lo,*ch))
    return max(total/(2*h),0)


def numeric_rho(center,h,ch,cs):
    def f(x):
        vh=ch[0]+ch[1]*x+ch[2]*x*x
        vs=cs[0]+cs[1]*x+cs[2]*x*x
        ss=1 if vs>=0 else -1
        return max(abs(vh)-ss*vh,0)
    return quad(f,center-h,center+h,epsabs=1e-12,epsrel=1e-12,points=[r for r in roots_quad(*ch)+roots_quad(*cs) if center-h<r<center+h])[0]/(2*h)

for rho,xi,K in [(0.2,0.1,2),(0.5,0.3,2),(0.3,0,2),(0.5,0.3,3),(0.3,0,3)]:
    m2,m4=moments(rho,xi,K); g=m4[-1]-m2[-1]**2
    assert g>0
    print('moments',rho,xi,K,'m2=',m2[-1],'g=',g,'b(delta=1)=',1/math.sqrt(g))

refs={0.5:0.2654727577,1.0:0.3829918740,2.4:-0.0167525400,3.0:-0.2118179436}
for d,ref in refs.items():
    got=bar_delta(d)
    assert abs(got-ref)<2e-8,(d,got,ref)
    print('barDelta',d,got)

ch=(0.03,0.9,0.15);cs=(0,1.0,0.1);h=0.7
for c in [-0.2,0,0.25]:
    a=exact_rho(c,h,ch,cs);b=numeric_rho(c,h,ch,cs)
    assert abs(a-b)<1e-10,(c,a,b)
    print('rho integral',c,a,b)

# E3 constants used as an independent regression check.
rho=.5;xi=.3;m2,m4=moments(rho,xi,2);g=m4[-1]-m2[-1]**2
assert abs(g-89/900)<1e-12
lam=(4/3*rho*rho*xi*xi)/g
assert abs(lam-27/89)<1e-12
print('E3 g2',g,'lambda',lam)
print('all reference v2.1 checks passed')
