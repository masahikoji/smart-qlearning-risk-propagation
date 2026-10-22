import numpy as np
from dataclasses import dataclass

@dataclass
class Scenario:
    K:int
    rho:float
    xi:float
    delta:float

BETA0=0.25
BETA1=0.50
SIGMA_Y=1.0
TAU0_2=1.35
TAU1_2=0.20

def beta_past(K,s):
    if s==1: return []
    if s==2: return [0.20]
    if s==3: return [0.15,0.15]
    raise ValueError

def tau0(K,s):
    return 1.35

def tau1(K,s):
    return 0.20

def feature_wide(x,past,a):
    # [1,x,x2,past...,a,a*x,a*x2]
    cols=[np.ones_like(x),x,x*x]
    cols.extend([p for p in past])
    cols.extend([a,a*x,a*x*x])
    return np.column_stack(cols)

def feature_narrow(x,past,a):
    return feature_wide(x,past,a)[:,:-1]

def feature_value_wide(x,past):
    a=np.ones_like(x)
    return feature_wide(x,past,a)

def baseline_q_coefs(K,rho,xi):
    # list 1-indexed of wide coefficient vectors and V history coefficients dict
    q=[None]*(K+1)
    v=[None]*(K+1) # each [int,x,past...]
    for s in range(K,0,-1):
        bp=beta_past(K,s)
        if s==K:
            qint=BETA0; qx=BETA1; qpast=np.array(bp,float)
            qa=tau0(K,s); qax=tau1(K,s)
        else:
            vn=v[s+1]
            v0,vx=vn[0],vn[1]
            vp=np.array(vn[2:],float) # length s, for A1..As
            qint=BETA0+v0
            qx=BETA1+rho*vx
            bpv=np.array(bp,float) if bp else np.zeros(0)
            qpast=bpv+vp[:s-1]
            qa=tau0(K,s)+xi*vx+vp[s-1]
            qax=tau1(K,s)
        vec=np.concatenate([[qint,qx,0.0],qpast,[qa,qax,0.0]])
        q[s]=vec
        v[s]=np.concatenate([[qint+qa,qx+qax],qpast])
    return q,v

def state_moments_terminal(K,rho,xi,N=3000000,seed=123):
    rng=np.random.default_rng(seed)
    # use chunks to compute E x2/x4 terminal
    sm2=0.;sm4=0.;count=0
    chunk=500000
    for start in range(0,N,chunk):
        n=min(chunk,N-start)
        x=rng.uniform(-1,1,n)
        for s in range(1,K):
            a=rng.choice(np.array([-1.,1.]),n)
            u=rng.uniform(-1,1,n)
            x=rho*x+(1-rho)*u+xi*a
        sm2 += np.sum(x*x); sm4 += np.sum(x**4); count += n
    return sm2/count, sm4/count

def calibrate_b(K,rho,xi,delta):
    # analytic even moments from R code
    m2=1/3; m4=1/5
    h=1-rho
    for s in range(2,K+1):
        m2_new=rho**2*m2+h**2/3+xi**2
        m4_new=(rho**4*m4+h**4/5+xi**4+6*rho**2*h**2*m2/3+
                6*rho**2*xi**2*m2+6*h**2*xi**2/3)
        m2,m4=m2_new,m4_new
    g=m4-m2*m2
    return delta/np.sqrt(g),g

def estimate_matrices(sc:Scenario,N=4000000,seed=20261016):
    K,rho,xi=sc.K,sc.rho,sc.xi
    rng=np.random.default_rng(seed)
    dims={s:s+5 for s in range(1,K+1)}
    G={s:np.zeros((dims[s],dims[s])) for s in range(1,K+1)}
    M={s:np.zeros((dims[s],dims[s+1])) for s in range(1,K)}
    count=0
    chunk=200000
    for st in range(0,N,chunk):
        n=min(chunk,N-st)
        x=rng.uniform(-1,1,n)
        past=[]
        Wstage={}
        # simulate all states/actions, store current wide features
        for s in range(1,K+1):
            a=rng.choice(np.array([-1.,1.]),n)
            W=feature_wide(x,past,a)
            Wstage[s]=W
            G[s]+=W.T@W
            if s<K:
                u=rng.uniform(-1,1,n)
                xnext=rho*x+(1-rho)*u+xi*a
                past_next=past+[a]
                Vnext=feature_value_wide(xnext,past_next)
                M[s]+=W.T@Vnext
                x=xnext; past=past_next
        count+=n
    for s in G: G[s]/=count
    for s in M: M[s]/=count
    return G,M

def prepare(sc,Npop=4000000):
    G,M=estimate_matrices(sc,N=Npop)
    K=sc.K
    q0,v0=baseline_q_coefs(K,sc.rho,sc.xi)
    b,g=calibrate_b(K,sc.rho,sc.xi,sc.delta)
    # residual variances of Bellman responses
    resvar={}
    for s in range(1,K+1):
        if s==K: resvar[s]=1.0
        else:
            vx=v0[s+1][1]
            resvar[s]=1.0+(vx*(1-sc.rho))**2/3
    # narrow/wide matrices and residualization
    info={}
    for s in range(1,K+1):
        Gw=G[s]
        Gn=Gw[:-1,:-1]
        grx=Gw[-1,:-1][None,:]
        a=(grx@np.linalg.inv(Gn)).reshape(-1)
        omega=Gw[-1,-1]-a@Gw[:-1,-1]
        info[s]={'Gw':Gw,'Gn':Gn,'a':a,'omega':omega,'resvar':resvar[s]}
        if s<K:
            Mw=M[s]
            info[s]['Mw']=Mw
            info[s]['Mn']=Mw[:-1,:]
            info[s]['mR']=Mw[-1,:]-a@Mw[:-1,:]
    # local score drifts, qstar coefs
    muw={s:np.zeros(G[s].shape[0]) for s in G}
    # terminal local perturbation r = A X^2 is last feature
    muw[K]=b*G[K][:,-1]
    mustar={}
    mustar[K]=np.zeros(G[K].shape[0]); mustar[K][-1]=b
    for s in range(K-1,0,-1):
        mustar[s]=np.linalg.solve(G[s],M[s]@mustar[s+1])
    return info,muw,mustar,b,g,q0,v0

def simulate_limit(sc,ndraw=2000000,Npop=4000000,seed=777,chunk=100000):
    info,muw,qstar,b,g,q0,v0=prepare(sc,Npop=Npop)
    K=sc.K
    rng=np.random.default_rng(seed)
    sums={('SEL',s):0.0 for s in range(1,K+1)}
    sums.update({('MA',s):0.0 for s in range(1,K+1)})
    sums2={k:0.0 for k in sums}
    gap_sums={s:0.0 for s in range(1,K+1)}
    gap_sums2={s:0.0 for s in range(1,K+1)}
    cnt=0
    # Cholesky stage score covariance = resvar * Gw
    chol={s:np.linalg.cholesky(info[s]['resvar']*info[s]['Gw']+1e-14*np.eye(info[s]['Gw'].shape[0])) for s in range(1,K+1)}
    for st in range(0,ndraw,chunk):
        n=min(chunk,ndraw-st)
        U={s:rng.standard_normal((n,info[s]['Gw'].shape[0]))@chol[s].T for s in range(1,K+1)}
        C={'SEL':None,'MA':None}
        Cstage={'SEL':{},'MA':{}}
        # backwards
        for s in range(K,0,-1):
            for proc in ['SEL','MA']:
                down=np.zeros((n, info[s]['Gw'].shape[0])) # placeholder no
                if s==K:
                    dn=np.zeros((n,info[s]['Gn'].shape[0]))
                    dw=np.zeros((n,info[s]['Gw'].shape[0]))
                    dj=np.zeros(n)
                else:
                    Cnext=Cstage[proc][s+1]
                    dn=Cnext@info[s]['Mn'].T
                    dw=Cnext@info[s]['Mw'].T
                    dj=Cnext@info[s]['mR']
                mu=muw[s]
                un=U[s][:,:-1]
                uw=U[s]
                bn=np.linalg.solve(info[s]['Gn'],(un+mu[:-1]+dn).T).T
                bw=np.linalg.solve(info[s]['Gw'],(uw+mu+dw).T).T
                j= (U[s][:,-1]-U[s][:,:-1]@info[s]['a']) + (mu[-1]-info[s]['a']@mu[:-1]) + dj
                lam=j*j/(info[s]['resvar']*info[s]['omega'])
                if proc=='SEL':
                    w=(lam>2.0).astype(float)
                else:
                    w=1/(1+np.exp(1-lam/2))
                bne=np.zeros_like(bw); bne[:,:-1]=bn
                c=bne+w[:,None]*(bw-bne)
                Cstage[proc][s]=c
        risks={'SEL':{},'MA':{}}
        for proc in ['SEL','MA']:
            for s in range(1,K+1):
                d=Cstage[proc][s]-qstar[s]
                risk=np.einsum('bi,ij,bj->b',d,info[s]['Gw'],d)
                risks[proc][s]=risk
                sums[(proc,s)]+=risk.sum(); sums2[(proc,s)]+=(risk*risk).sum()
        for s in range(1,K+1):
            gap=risks['SEL'][s]-risks['MA'][s]
            gap_sums[s]+=gap.sum(); gap_sums2[s]+=(gap*gap).sum()
        cnt+=n
    out=[]
    for s in range(1,K+1):
        msel=sums[('SEL',s)]/cnt; mma=sums[('MA',s)]/cnt
        gap=gap_sums[s]/cnt
        gap_var=max(gap_sums2[s]/cnt-gap*gap,0.0)
        gap_se=(gap_var/cnt)**0.5
        out.append((s,msel,mma,gap,gap_se))
    return out, {'b':b,'g':g,'resvar':{s:info[s]['resvar'] for s in info},'omega':{s:info[s]['omega'] for s in info}}

if __name__=='__main__':
    for K in [2,3]:
        for delta in [1.0,3.0]:
            sc=Scenario(K=K,rho=0.5,xi=0.3,delta=delta)
            out,meta=simulate_limit(sc,ndraw=1000000,Npop=3000000,seed=777+K*10+int(delta*10),chunk=50000)
            print('K',K,'delta',delta,'b',meta['b'],'g',meta['g'])
            for row in out: print(' stage %d SEL %.6f MA %.6f gap %.6f SE %.6f'%row)
