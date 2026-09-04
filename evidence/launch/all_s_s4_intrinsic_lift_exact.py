"""Exact Fraction/interval certificate for the s=4 intrinsic 4x4 lift.

It reproduces the intrinsic first-phase epsilon-series recursion on the
complete certified Hopf q interval and critical-alpha envelope, encloses
every 9x9 solve by a midpoint-residual
Neumann argument, and certifies the separate 4x4 lift matrix.
"""

if not __debug__:
    raise RuntimeError(
        "Verification assertions are disabled. "
        "Run without -O/-OO and unset PYTHONOPTIMIZE."
    )

from contextlib import redirect_stdout
from fractions import Fraction as F
from io import StringIO
from pathlib import Path

env={"__name__":"__intrinsic_lift_exact__"}
with redirect_stdout(StringIO()):
    exec(compile(
        Path("evidence/launch/all_s_s4_transverse_hopf_exact.py").read_text(),
        "transverse_hopf_exact","exec"),env)

I,ia,isu,im=env["I"],env["ia"],env["isu"],env["im"]
idv,iz=env["idv"],env["iz"]
ipa,ipm,ipdivmonic=env["ipa"],env["ipm"],env["ipdivmonic"]
P=ipa(env["P"])
Q=env["family_intervals"](env["hopf_q"])[0]
Delta=[z[0] for z in env["hd"]["Delta"]]
Hpoly=[z[0] for z in env["hd"]["Hpoly"]]
C=I(env["C"])
alpha=env["idec"]("0.2952194200019594","0.2952194200028030")
Z=I(0);O=I(1)
s=4;r=1;N=6;m=N+s-1

def mid(x): return (x[0]+x[1])/2
def au(x): return max(abs(x[0]),abs(x[1]))
def neg(x): return (-x[1],-x[0])
def addp(p,q):
    n=max(len(p),len(q))
    return [ia(p[i] if i<len(p) else Z,q[i] if i<len(q) else Z)
            for i in range(n)]
def subp(p,q):
    n=max(len(p),len(q))
    return [isu(p[i] if i<len(p) else Z,q[i] if i<len(q) else Z)
            for i in range(n)]
def pad(p,n): return list(p)+[Z]*(n-len(p))

def mat_inverse(A):
    n=len(A)
    aug=[list(A[i])+[F(i==j) for j in range(n)] for i in range(n)]
    for j in range(n):
        pivot=next(i for i in range(j,n) if aug[i][j])
        aug[j],aug[pivot]=aug[pivot],aug[j]
        z=aug[j][j]
        aug[j]=[x/z for x in aug[j]]
        for i in range(n):
            if i==j: continue
            z=aug[i][j]
            if z:
                aug[i]=[aug[i][k]-z*aug[j][k] for k in range(2*n)]
    return [row[n:] for row in aug]

def matvec_exact_interval(A,x):
    return [sum_intervals([im(I(A[i][j]),x[j]) for j in range(len(x))])
            for i in range(len(A))]

def matvec_interval(A,x):
    return [sum_intervals([im(A[i][j],x[j]) for j in range(len(x))])
            for i in range(len(A))]

def sum_intervals(xs):
    z=Z
    for x in xs:z=ia(z,x)
    return z

def inf_exact(A): return max(sum(abs(x) for x in row) for row in A)
def inf_interval(A): return max(sum(au(x) for x in row) for row in A)

def midpoint_inverse_data(A):
    n=len(A)
    A0=[[mid(A[i][j]) for j in range(n)] for i in range(n)]
    R=mat_inverse(A0)
    RA=[[sum_intervals([im(I(R[i][k]),A[k][j]) for k in range(n)])
         for j in range(n)] for i in range(n)]
    E=[[isu(I(F(i==j)),RA[i][j]) for j in range(n)] for i in range(n)]
    eta=inf_interval(E)
    assert eta<1, float(eta)
    return R,E,eta

def interval_solve(A,b):
    n=len(A);R,E,eta=midpoint_inverse_data(A)
    A0=[[mid(A[i][j]) for j in range(n)] for i in range(n)]
    b0=[mid(x) for x in b]
    x0=[sum(R[i][j]*b0[j] for j in range(n)) for i in range(n)]
    Ax0=[sum_intervals([im(A[i][j],I(x0[j])) for j in range(n)])
         for i in range(n)]
    res=[isu(b[i],Ax0[i]) for i in range(n)]
    Rres=matvec_exact_interval(R,res)
    err=max(au(x) for x in Rres)/(1-eta)
    box=[(-err,err) for _ in range(n)]
    for _ in range(1):
        image=[ia(Rres[i],sum_intervals([im(E[i][j],box[j])
                                         for j in range(n)]))
               for i in range(n)]
        narrowed=[]
        for old,new in zip(box,image):
            z=(max(old[0],new[0]),min(old[1],new[1]))
            assert z[0]<=z[1]
            narrowed.append(z)
        box=narrowed
    return [(x0[i]+box[i][0],x0[i]+box[i][1])
            for i in range(n)],R,eta,max(au(x) for x in box)

def shifted(poly,j,total=9):
    return [Z]*j+list(poly)+[Z]*(total-j-len(poly))

def one_series(Ps,Ss,next0,solve_log):
    G=[]
    for o in range(3):
        z=[Z]*m
        for aa in range(o+1):
            cc=ipm(Ps[aa],Ss[o-aa])
            for j,v in enumerate(cc):
                if j<m:z[j]=ia(z[j],v)
        G.append(z)
    def matrix(Gv,blocks=True):
        cols=[]
        for j in range(s):cols.append(shifted(Gv,j))
        for j in range(s-1):
            cols.append(shifted([neg(x) for x in Delta],j) if blocks
                        else [Z]*m)
        for j in range(r+1):
            col=[Z]*m
            if blocks:col[j]=I(-1)
            cols.append(col)
        return [[cols[j][i] for j in range(m)] for i in range(m)]
    M0=matrix(G[0]);M1=matrix(G[1],False);M2=matrix(G[2],False)
    def rhs(Gv):
        b=[Z]*m
        for k,v in enumerate(Gv):
            if k+s<m:b[k+s]=isu(b[k+s],v)
        return b
    b0=rhs(G[0])
    for k,v in enumerate(Delta):
        if k+s-1<m:b0[k+s-1]=ia(b0[k+s-1],v)
    # On the exact fixed two-cycle, S A A_+ = Delta(SH)+C S.
    # Insert this identity instead of interval-solving the order-zero
    # equation; this prevents dependency inflation from contaminating the
    # next phase while retaining an exact residual check.
    Sigma0=ipm(Ss[0],Hpoly)
    hS0=[im(C,x) for x in Ss[0]]
    x0=list(next0[:s])+list(Sigma0[:s-1])+list(hS0[:r+1])
    residual0=subp(b0,matvec_interval(M0,x0))
    assert all(iz(x) for x in residual0)
    _,_,eta0=midpoint_inverse_data(M0);err0=F(0)
    b1=subp(rhs(G[1]),matvec_interval(M1,x0))
    x1,_,eta1,err1=interval_solve(M0,b1)
    b2=subp(subp(rhs(G[2]),matvec_interval(M1,x1)),
            matvec_interval(M2,x0))
    x2,_,eta2,err2=interval_solve(M0,b2)
    solve_log.append((eta0,max(err0,err1,err2)))
    xs=[x0,x1,x2]
    Qs=[xs[o][:s]+([O] if o==0 else [Z]) for o in range(3)]
    Avec=[[xs[o][s+s-1+j] for j in range(r+1)] for o in range(3)]
    hh=[Avec[o][r] for o in range(3)]
    assert not iz(hh[0])
    Sout=[[Z]*(r+1) for _ in range(3)]
    for j in range(r+1):
        Sout[0][j]=Ss[0][j]
        Sout[1][j]=idv(isu(Avec[1][j],im(hh[1],Sout[0][j])),hh[0])
        Sout[2][j]=idv(isu(isu(Avec[2][j],im(hh[1],Sout[1][j])),
                              im(hh[2],Sout[0][j])),hh[0])
    Sig=[xs[o][s:s+s-1]+([O] if o==0 else [Z]) for o in range(3)]
    return Qs,Sout,Sig,hh

def state(Ss,Sig):
    Hs=[];Us=[]
    for o in range(3):
        if o==0:
            Hs.append(list(Hpoly));Us.append([Z]);continue
        target=list(Sig[o])
        for aa in range(1,o+1):
            target=subp(target,ipm(Ss[aa],Hs[o-aa]))
        hq,rem=ipdivmonic(target,Ss[0])
        Hs.append(pad(hq,3));Us.append(rem)
    return Hs,Us

zeroP=[Z]*5;zeroS=[Z]*2
def orbit(p):
    log=[]
    S=[neg(alpha),O];Ps=[P,list(p)+[Z],zeroP];Ss=[S,zeroS,zeroS]
    Qs,S1,Sig0,h1=one_series(Ps,Ss,Q,log);H0,U0=state(Ss,Sig0)
    Rs,S2,Sig1,h2=one_series(Qs,S1,P,log);H1,U1=state(S1,Sig1)
    Q3,S3,Sig2,h3=one_series(Rs,S2,Q,log);H2,U2=state(S2,Sig2)
    return (Ss,H0,U0,h1),(S2,H2,U2,h3),log

def first_state(p):
    log=[]
    S=[neg(alpha),O];Ps=[P,list(p)+[Z],zeroP];Ss=[S,zeroS,zeroS]
    Qs,S1,Sig0,h1=one_series(Ps,Ss,Q,log);H0,U0=state(Ss,Sig0)
    return (Ss,H0,U0,h1),log

basis=[]
for j in range(4):
    p=[Z]*4;p[j]=O;basis.append(p)
cols=[];all_logs=[]
for p in basis:
    st0,logs=first_state(p);all_logs.extend(logs)
    cols.append([st0[3][1],st0[1][1][0],st0[1][1][1],st0[2][1][0]])
J=[[cols[j][i] for j in range(4)] for i in range(4)]
published_J_boxes = [
    [(12, 13), (-222, -221), (452, 453), (-2582, -2581)],
    [(-729, -728), (1693, 1694), (-9071, -9070), (28934, 28935)],
    [(485, 486), (-1124, -1123), (6038, 6039), (-19226, -19225)],
    [(-4, -3), (0, 1), (-26, -25), (45, 46)],
]
for row, boxes in zip(J, published_J_boxes):
    for value, (lower, upper) in zip(row, boxes):
        assert F(lower) <= value[0] <= value[1] <= F(upper)
R4,E4,eta4=midpoint_inverse_data(J)
inv4=inf_exact(R4)/(1-eta4)
assert eta4<F(1,4)
assert inv4<2**20
lift,_,_,lift_err=interval_solve(J,[Z,Z,Z,O])
assert all(not iz(x) for x in lift)
assert -3<lift[0][0] and lift[0][1]<-2
assert -2<lift[1][0] and lift[1][1]<-F(1,2)
assert F(1,2)<lift[2][0] and lift[2][1]<1
assert F(1,10)<lift[3][0] and lift[3][1]<F(1,4)
assert max(x[0] for x in all_logs)<F(1,100000)

# Execute the complete P,Q,P recursion used in (C.3.14e), rather than merely
# defining it.  The first-order return is the identity; the second-order
# differences are the certified interval representations of beta and q.
initial_state, returned_state, return_logs = orbit(lift)
S0, H0, U0, C0 = initial_state
S2, H2, U2, C2 = returned_state
return_linear = [
    isu(U2[1][0], U0[1][0]),
    isu(C2[1], C0[1]),
    isu(H2[1][0], H0[1][0]),
    isu(H2[1][1], H0[1][1]),
]
assert all(iz(x) for x in return_linear)
return_beta = isu(U2[2][0], U0[2][0])
return_q = [
    isu(C2[2], C0[2]),
    isu(H2[2][0], H0[2][0]),
    isu(H2[2][1], H0[2][1]),
]
assert return_beta[0] <= return_beta[1]
assert all(x[0] <= x[1] for x in return_q)
assert len(return_logs) == 3
assert max(x[0] for x in return_logs) < F(1,100000)

from decimal import Decimal,getcontext
getcontext().prec=16
def show(x):
    return (Decimal(x[0].numerator)/Decimal(x[0].denominator),
            Decimal(x[1].numerator)/Decimal(x[1].denominator))

print("PASS: exact s=4 intrinsic lift on the Hopf critical envelope")
print("4x4 lift inverse infinity norm < 2^20")
print("4x4 Neumann eta =",show((eta4,eta4)))
print("lift coefficients =",[show(x) for x in lift])
print("maximum 9x9 Neumann eta =",show((max(x[0] for x in all_logs),)*2))
print("complete P,Q,P quadratic beta/q extraction executed")
