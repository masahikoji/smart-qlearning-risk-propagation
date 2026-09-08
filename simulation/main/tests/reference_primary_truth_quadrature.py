import numpy as np
from numpy.polynomial.legendre import leggauss

def real_roots(a0,a1,a2,tol=1e-12):
    if abs(a2)<tol:
        return np.array([]) if abs(a1)<tol else np.array([-a0/a1])
    d=a1*a1-4*a2*a0
    if d<-tol:return np.array([])
    if abs(d)<=tol:return np.array([-a1/(2*a2)])
    s=np.sqrt(max(d,0)); return np.sort(np.array([(-a1-s)/(2*a2),(-a1+s)/(2*a2)]))

def anti(x,a0,a1,a2):return a0*x+.5*a1*x*x+a2*x**3/3

def mean_abs(center,h,a0,a1,a2):
    center=np.asarray(center); l=center-h; u=center+h
    roots=real_roots(a0,a1,a2); cuts=np.r_[-np.inf,roots,np.inf]
    integ=np.zeros(center.shape)
    for lo0,hi0 in zip(cuts[:-1],cuts[1:]):
        lo=np.maximum(l,lo0); hi=np.minimum(u,hi0); ok=hi>lo
        if not np.any(ok):continue
        mid=(lo[ok]+hi[ok])/2
        sg=np.sign(a0+a1*mid+a2*mid*mid); sg[sg==0]=1
        integ[ok]+=sg*(anti(hi[ok],a0,a1,a2)-anti(lo[ok],a0,a1,a2))
    return integ/(2*h)

def moments(rho,xi,K=3):
    m2=np.zeros(K);m4=np.zeros(K);m2[0]=1/3;m4[0]=1/5;h=1-rho
    for s in range(1,K):
        m2[s]=rho**2*m2[s-1]+h*h/3+xi*xi
        m4[s]=rho**4*m4[s-1]+h**4/5+xi**4+6*rho**2*h**2*m2[s-1]/3+6*rho**2*xi**2*m2[s-1]+6*h*h*xi*xi/3
    return m2,m4

def scenario(name,rho,xi,boundary='separated',delta=1,tau0=None,tau1=None,bmanual=None,fixed=None):
    K=3; m2,m4=moments(rho,xi,K); g=m4[-1]-m2[-1]**2
    b=np.zeros(3) if bmanual is None else np.array(bmanual,float)
    if delta is not None: b[-1]=delta/np.sqrt(g)
    if tau0 is None:tau0=np.array([1.35]*3)
    else:tau0=np.array(tau0,float)
    if tau1 is None:tau1=np.array([.2]*3)
    else:tau1=np.array(tau1,float)
    return dict(name=name,rho=rho,xi=xi,boundary=boundary,b=b,tau0=tau0,tau1=tau1,fixed=fixed)

sc=[]
sc.append(scenario('weak',.2,.1))
sc.append(scenario('strong',.5,.3))
sc.append(scenario('multi',.5,.3,bmanual=[0,2,0]))
sc.append(scenario('smooth',.3,0,'smooth',delta=1,tau0=[1.5,1.5,0],tau1=[.15,.15,1]))
sc.append(scenario('tie',.3,0,'exact_tie',delta=None,bmanual=[0,0,0],tau0=[1.5,1.5,0],tau1=[.15,.15,0]))
sc.append(scenario('delta05',.5,.3,delta=.5))
sc.append(scenario('delta24',.5,.3,delta=2.4))
sc.append(scenario('delta30',.5,.3,delta=3.0))
# beta0=.25 beta1=.5; past stage2=.2; stage3=.15,.15

def cvec(S,n):
    c=S['b']/np.sqrt(n)
    if S['boundary']=='exact_tie':c[-1]=0
    return c

def q2(S,n,x2,a1,a2):
    c=cvec(S,n); rho=S['rho'];xi=S['xi'];h=1-rho
    y2=.25+.5*x2+.2*a1+a2*(S['tau0'][1]+S['tau1'][1]*x2+c[1]*x2*x2)
    cen=rho*x2+xi*a2
    ev3=.25+.5*cen+.15*a1+.15*a2+mean_abs(cen,h,S['tau0'][2],S['tau1'][2],c[2])
    return y2+ev3

def q1(S,n,x1,a1,order):
    c=cvec(S,n);rho=S['rho'];xi=S['xi'];h=1-rho
    y1=.25+.5*x1+a1*(S['tau0'][0]+S['tau1'][0]*x1+c[0]*x1*x1)
    nodes,w=leggauss(order); cen=rho*x1+xi*a1
    X=cen[:,None]+h*nodes[None,:]
    A1=np.broadcast_to(a1[:,None],X.shape)
    vp=q2(S,n,X.ravel(),A1.ravel(),np.ones(X.size)).reshape(X.shape)
    vm=q2(S,n,X.ravel(),A1.ravel(),-np.ones(X.size)).reshape(X.shape)
    return y1+(np.maximum(vp,vm)*(w/2)).sum(1)

xs=np.linspace(-1,1,401)
maxerr=0
for S in sc:
  for n in [250,500,1000]:
    for a in [-1,1]:
      aa=np.full(xs.shape,a,float)
      q48=q1(S,n,xs,aa,48)
      q192=q1(S,n,xs,aa,192)
      e=np.max(np.abs(q48-q192));maxerr=max(maxerr,e)
      print(S['name'],n,a,e)
print('MAX',maxerr)
