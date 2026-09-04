"""Exact rational-interval certificate for the declared transverse s=4 Hopf point."""

if not __debug__:
    raise RuntimeError(
        "Verification assertions are disabled. "
        "Run without -O/-OO and unset PYTHONOPTIMIZE."
    )

from decimal import Decimal, getcontext
from fractions import Fraction as F

Z=F(0); O=F(1)

def dec(s):
    s=str(s); neg=s.startswith("-")
    if neg: s=s[1:]
    if "." in s: a,b=s.split(".")
    else: a,b=s,""
    x=F(int(a+b),10**len(b))
    return -x if neg else x

def trim(p):
    p=list(p)
    while len(p)>1 and p[-1]==0: p.pop()
    return p

def add(p,q):
    n=max(len(p),len(q))
    return trim([(p[i] if i<len(p) else Z)+(q[i] if i<len(q) else Z)
                 for i in range(n)])

def sub(p,q):
    n=max(len(p),len(q))
    return trim([(p[i] if i<len(p) else Z)-(q[i] if i<len(q) else Z)
                 for i in range(n)])

def scale(p,a): return trim([a*x for x in p])

def mul(p,q):
    out=[Z]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q): out[i+j]+=a*b
    return trim(out)

def divmodp(p,q):
    p=trim(p);q=trim(q)
    if q==[Z]: raise ZeroDivisionError
    if len(p)<len(q): return [Z],p
    out=[Z]*(len(p)-len(q)+1)
    while len(p)>=len(q) and p!=[Z]:
        k=len(p)-len(q);a=p[-1]/q[-1];out[k]=a
        for j,b in enumerate(q):p[k+j]-=a*b
        p=trim(p)
    return trim(out),p

def deriv(p):
    return trim([F(i)*p[i] for i in range(1,len(p))]) if len(p)>1 else [Z]

def peval(p,x):
    y=Z
    for a in reversed(p):y=y*x+a
    return y

def proots(rs):
    p=[O]
    for r in rs:p=mul(p,[-r,O])
    return p

def sturm(p):
    seq=[trim(p),deriv(p)]
    while True:
        _,r=divmodp(seq[-2],seq[-1])
        if r==[Z]:return seq
        seq.append(scale(r,-O))

def variations_at(seq,x):
    s=[]
    for p in seq:
        y=peval(p,x)
        if y:s.append(1 if y>0 else -1)
    return sum(s[i]!=s[i-1] for i in range(1,len(s)))

def root_count(p,a,b):
    s=sturm(p)
    return variations_at(s,a)-variations_at(s,b)

# Closed rational intervals.
def I(a,b=None):
    a=F(a);return (a,a) if b is None else (a,F(b))
def idec(a,b):return dec(a),dec(b)
def ia(x,y):return x[0]+y[0],x[1]+y[1]
def ine(x):return -x[1],-x[0]
def isu(x,y):return ia(x,ine(y))
def im(x,y):
    z=[a*b for a in x for b in y];return min(z),max(z)
def ii(x):
    assert not x[0]<=0<=x[1]
    return O/x[1],O/x[0]
def idv(x,y):return im(x,ii(y))
def isc(x,a):return im(I(a),x)
def iz(x):return x[0]<=0<=x[1]

def ipa(p):return [I(x) for x in p]
def ipe(p,x):
    y=I(0)
    for a in reversed(p):y=ia(im(y,x),a)
    return y
def ipd(p):return [isc(p[i],i) for i in range(1,len(p))]
def ipm(p,q):
    out=[I(0)]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):out[i+j]=ia(out[i+j],im(a,b))
    return out
def ipadd(p,q):
    n=max(len(p),len(q));out=[]
    for i in range(n):
        out.append(ia(p[i] if i<len(p) else I(0),
                      q[i] if i<len(q) else I(0)))
    return out
def iproot(p,r):
    n=len(p)-1;q=[I(0)]*n;q[n-1]=p[n]
    for k in range(n-1,0,-1):q[k-1]=ia(p[k],im(r,q[k]))
    return q,ia(p[0],im(r,q[0]))
def ipdivmonic(p,q):
    p=list(p);out=[I(0)]*(len(p)-len(q)+1)
    while len(p)>=len(q):
        k=len(p)-len(q);a=p[-1];out[k]=a
        for j,b in enumerate(q):p[k+j]=isu(p[k+j],im(a,b))
        p.pop()
    return out,p

# Five-component jet (v, d_a, d_aa, d_q, d_aq).
def K(v,a=None,aa=None,q=None,aq=None):
    return (v,I(0) if a is None else a,I(0) if aa is None else aa,
            I(0) if q is None else q,I(0) if aq is None else aq)
def ka(x,y):return tuple(ia(x[i],y[i]) for i in range(5))
def kn(x):return tuple(ine(z) for z in x)
def ks(x,y):return ka(x,kn(y))
def km(x,y):
    v=im(x[0],y[0])
    a=ia(im(x[1],y[0]),im(x[0],y[1]))
    aa=ia(ia(im(x[2],y[0]),isc(im(x[1],y[1]),2)),im(x[0],y[2]))
    q=ia(im(x[3],y[0]),im(x[0],y[3]))
    aq=ia(ia(ia(im(x[4],y[0]),im(x[1],y[3])),
                 im(x[3],y[1])),im(x[0],y[4]))
    return v,a,aa,q,aq
def ki(x):
    v=ii(x[0]);v2=im(v,v);v3=im(v2,v)
    a=ine(im(x[1],v2))
    aa=isu(isc(im(im(x[1],x[1]),v3),2),im(x[2],v2))
    q=ine(im(x[3],v2))
    aq=isu(isc(im(im(x[1],x[3]),v3),2),im(x[4],v2))
    return v,a,aa,q,aq
def kd(x,y):return km(x,ki(y))
def ksc(x,a):return tuple(isc(z,a) for z in x)
def kpe(p,x):
    y=K(I(0))
    for a in reversed(p):y=ka(km(y,x),a)
    return y
def kpd(p):return [ksc(p[i],i) for i in range(1,len(p))]
def kpm(p,q):
    out=[K(I(0))]*(len(p)+len(q)-1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):out[i+j]=ka(out[i+j],km(a,b))
    return out
def kpdivmonic(p,q):
    p=list(p);out=[K(I(0))]*(len(p)-len(q)+1)
    while len(p)>=len(q):
        k=len(p)-len(q);a=p[-1];out[k]=a
        for j,b in enumerate(q):p[k+j]=ks(p[k+j],km(a,b))
        p.pop()
    return out,p
def kproot(p,r):
    n=len(p)-1;q=[K(I(0))]*n;q[n-1]=p[n]
    for j in range(n-1,0,-1):q[j-1]=ka(p[j],km(r,q[j]))
    return q,ka(p[0],km(r,q[0]))

# Fixed rational family.
P=proots([dec("-3.812196375334383"),dec("-0.3408184613601878"),
          dec("1.5054785540322202"),dec("2.190808973956946")])
Qbar=proots([dec("-2.3304189906225"),dec("0.6733921795528344"),
             dec("4.323444838404325")])
C=dec("0.31355365301369736")
q_global=idec("-1.858","-1.851")

def Qexact(q):return mul([-q,O],Qbar)
def rational_family(q):
    Q=Qexact(q);pq=mul(P,Q)
    mi=list(pq);mi[0]-=C
    pl=list(pq);pl[0]+=C
    return Q,mi,pl

# Global root collars. Endpoint Sturm counts label the branches; uniform
# boundary signs and nonzero t-derivatives continue each root for all q.
E_global=[
 idec("-3.812224","-3.812222"),idec("-2.329779","-2.329767"),
 idec("-1.859053","-1.852044"),idec("-0.339465","-0.339456"),
 idec("0.671419","0.671428"),idec("1.507016","1.507022"),
 idec("2.190298","2.190302"),idec("4.323453","4.323455")]
N_global=[
 idec("-3.812171","-3.812168"),idec("-2.331069","-2.331056"),
 idec("-1.856952","-1.849955"),idec("-0.342178","-0.342169"),
 idec("0.675358","0.675367"),idec("1.503932","1.503938"),
 idec("2.191315","2.191319"),idec("4.323434","4.323437")]

for q in q_global:
    _,mi,pl=rational_family(q)
    assert all(root_count(mi,*b)==1 for b in E_global)
    assert all(root_count(pl,*b)==1 for b in N_global)

def family_intervals(qbox):
    Qi=ipm([ine(qbox),I(1)],ipa(Qbar))
    pqi=ipm(ipa(P),Qi)
    mi=list(pqi);mi[0]=isu(mi[0],I(C))
    pl=list(pqi);pl[0]=ia(pl[0],I(C))
    return Qi,mi,pl

def sign_interval(x):
    if x[0]>0:return 1
    if x[1]<0:return -1
    return 0

def family_value(kind,t,qbox):
    base=im(im(ipe(ipa(P),t),isu(t,qbox)),ipe(ipa(Qbar),t))
    return ia(base,I(C)) if kind=="plus" else isu(base,I(C))

def family_t_derivative(t,qbox):
    p=ipe(ipa(P),t);pd=ipe(ipa(deriv(P)),t)
    qb=ipe(ipa(Qbar),t);qbd=ipe(ipa(deriv(Qbar)),t)
    return ia(im(im(pd,isu(t,qbox)),qb),
              im(p,ia(qb,im(isu(t,qbox),qbd))))

def certify_global(kind,boxes):
    for j,b in enumerate(boxes):
        sl=sign_interval(family_value(kind,I(b[0]),q_global))
        sr=sign_interval(family_value(kind,I(b[1]),q_global))
        assert sl*sr==-1,(j,b)
        assert not iz(family_t_derivative(b,q_global)),(j,b)

Qi,mi,pl=family_intervals(q_global)
certify_global("minus",E_global);certify_global("plus",N_global)

def shrink(kind,qbox,box):
    lo,hi=box
    slo=sign_interval(family_value(kind,I(lo),qbox))
    shi=sign_interval(family_value(kind,I(hi),qbox))
    assert slo*shi==-1
    for _ in range(90):
        mid=(lo+hi)/2;sm=sign_interval(family_value(kind,I(mid),qbox))
        if sm==slo:lo=mid
        elif sm==shi:hi=mid
        else:break
    out=(lo,hi)
    assert not iz(family_t_derivative(out,qbox))
    return out

def root_boxes(qbox):
    return ([shrink("minus",qbox,b) for b in E_global],
            [shrink("plus",qbox,b) for b in N_global])

def root_jet(z,Qk):
    # F_q=-P Qbar and F_t=(P Q)' for both signs.
    Pk=[K(I(a)) for a in P];Qbark=[K(I(a)) for a in Qbar]
    pq=kpm(Pk,Qk)
    ft=kpe(kpd(pq),K(z))[0]
    fq=ine(ipe(ipa(mul(P,Qbar)),z))
    zq=ine(idv(fq,ft))
    return K(z,q=zq)

def kernel(R,Rd,B,Bd,t,g,same,alpha):
    beta=kd(kpe(B,alpha),kpe(R,alpha))
    def U(z):
        rz=kpe(R,z);bz=kpe(B,z);num=ks(km(beta,rz),bz)
        za=ks(z,alpha);value=kd(num,ksc(za,C))
        nt=ks(km(beta,kpe(Rd,z)),kpe(Bd,z))
        ut=kd(ks(km(nt,za),num),ksc(km(za,za),C))
        return value,ut
    Ut,Utt=U(t);Ug,_=U(g);Rt=kpe(R,t);Rg=kpe(R,g)
    if same:return ks(km(kpe(Rd,t),Ut),km(Utt,Rt))
    return kd(ks(km(Rt,Ug),km(Ut,Rg)),ks(t,g))

def build(qbox,E,N,abox,full=True):
    qk=K(qbox,q=I(1));Qk=kpm([kn(qk),K(I(1))],[K(I(a)) for a in Qbar])
    Pk=[K(I(a)) for a in P]
    Ej=[root_jet(z,Qk) for z in E];Nj=[root_jet(z,Qk) for z in N]
    hL,hR=Ej[2],Ej[5]
    pq=kpm(Pk,Qk);minus=list(pq);minus[0]=ks(minus[0],K(I(C)))
    q7,r1=kproot(minus,hL);Delta,r2=kproot(q7,hR)
    assert iz(r1[0]) and iz(r2[0])
    H=kpm([kn(hL),K(I(1))],[kn(hR),K(I(1))])
    AP,BP=kpdivmonic(Delta,Pk);AQ,BQ=kpdivmonic(Delta,Qk)
    Pd=kpd(Pk);Qd=kpd(Qk);BPd=kpd(BP);BQd=kpd(BQ)
    alpha=K(abox,a=I(1));labels=[Nj[1],Nj[6]]
    normals=[]
    for g in labels:
        Pg=kpe(Pk,g);Qg=kpe(Qk,g)
        normals.append(ksc(km(kd(Pg,ks(g,alpha)),
                              ks(kd(Pg,kpe(Pk,alpha)),
                                 kd(Qg,kpe(Qk,alpha)))),4))
    T=[[None,None],[None,None]]
    for ih,h in enumerate(labels):
        for ig,g in enumerate(labels):
            Pg=kpe(Pk,g);Ph=kpe(Pk,h);Qh=kpe(Qk,h)
            KP=kernel(Pk,Pd,BP,BPd,h,g,ih==ig,alpha)
            KQ=kernel(Qk,Qd,BQ,BQd,h,g,ih==ig,alpha)
            pg=ksc(km(Pg,KP),-4);qg=ksc(km(Pg,KQ),4)
            VC=ksc(ka(km(Pg,Pg),K(I(C))),4)
            rq=kd(ks(kpe(Delta,h),kpe(Delta,alpha)),ks(h,alpha))
            br=ks(ka(ka(km(pg,Qh),km(Ph,qg)),VC),km(normals[ig],rq))
            T[ih][ig]=ksc(br,-F(2)/C)
    vn=[normals[1],kn(normals[0])]
    rates=[ka(km(T[i][0],vn[0]),km(T[i][1],vn[1])) for i in range(2)]
    G=ks(rates[0],rates[1])
    if not full:
        return {"G":G,"normals":normals,"rates":rates}
    c=[kd(ksc(vn[0],2),rates[0]),kd(ksc(vn[1],2),rates[0])]
    kappa=ksc(kpe(Delta,alpha),F(2)/C)
    ell=[ksc(kd(ks(kpe(Delta,h),kpe(Delta,alpha)),ks(h,alpha)),
             -F(2)/C) for h in labels]
    # p and g=T_alpha c, with correct first alpha/q derivatives.
    def alpha_derivative(f):
        return K(f[1],a=f[2],q=f[4])
    p=ka(km(alpha_derivative(normals[0]),c[0]),
         km(alpha_derivative(normals[1]),c[1]))
    gv=[ka(km(alpha_derivative(T[i][0]),c[0]),
           km(alpha_derivative(T[i][1]),c[1])) for i in range(2)]
    r=km(c[0],normals[0])
    sigma=ka(km(c[1],T[0][1]),km(c[0],T[1][0]))
    dt=ks(K(I(2)),sigma)
    de=ks(ell[0],ell[1]);dg=ks(gv[0],gv[1])
    trace=ka(ka(km(kappa,p),km(r,de)),dt)
    q0=km(kappa,ks(km(dt,p),km(r,dg)))
    disc=ks(km(trace,trace),ksc(q0,4))
    Hopf=ks(ksc(trace,F(1,2)),ksc(disc,F(1,4)))
    nodes=[Ej[i] for i in (0,1,3,4,6,7)]
    weights=[]
    minusd=kpd(minus)
    for z in nodes:
        den=kd(kpe(minusd,z),kpe(H,z))
        weights.append(kd(km(ks(z,alpha),kpe(Qk,z)),den))
        weights.append(kd(km(ks(z,alpha),kpe(Pk,z)),den))
    return {"G":G,"H":Hopf,"disc":disc,"q0":q0,"trace":trace,
            "Delta":Delta,"Hpoly":H,"Tmatrix":T,"ell":ell,
            "kappa":kappa,
            "c":c,"weights":weights,"normals":normals,"rates":rates}

# Localized rational IFT rectangle.  Boundary signs and G_alpha<0 prove one
# continuous simple critical branch alpha(q) across the whole q interval.
qHL=dec("-1.851888497035");qHR=dec("-1.851888497033")
hopf_q=(qHL,qHR)
EH,NH=root_boxes(hopf_q)
hopf_a=idec("0.2952193","0.29521955")
lower_face=build(hopf_q,EH,NH,I(hopf_a[0]),False)
upper_face=build(hopf_q,EH,NH,I(hopf_a[1]),False)
assert lower_face["G"][0][0]>0
assert upper_face["G"][0][1]<0
hd=build(hopf_q,EH,NH,hopf_a)
assert hd["G"][1][1]<0
assert dec("-161517.159") < hd["G"][1][0]
assert hd["G"][1][1] < dec("-161516.321")
assert hd["disc"][0][1]<dec("-55.40")
assert min(w[0][0] for w in hd["weights"])>dec("0.0000027")
assert hd["normals"][0][0][1]<0<hd["normals"][1][0][0]
assert hd["rates"][0][0][0]>0
assert min(c[0][0] for c in hd["c"])>dec("0.00003")
assert hd["c"][0][0][0] > F(1, 2**15)
assert hd["c"][1][0][0] > F(1, 4)

# Implicit transversality on the complete localized branch.
alpha_q=ine(idv(hd["G"][3],hd["G"][1]))
hopf_slope=ia(hd["H"][3],im(hd["H"][1],alpha_q))
assert hopf_slope[0]>dec("109.89")
assert hopf_slope[1]<dec("110.82")
crossing_lower = hopf_slope[0] / (1 + 4 * F(931, 250)**2)
assert crossing_lower > F(194, 100)

# Tight exact alpha collars at the two rational q endpoints.  Their G face
# signs locate the endpoint critical roots, and H has opposite signs there.
def endpoint_data(q,abox):
    qb=I(q);E,N=root_boxes(qb)
    lo=build(qb,E,N,I(abox[0]),False)
    hi=build(qb,E,N,I(abox[1]),False)
    assert lo["G"][0][0]>0 and hi["G"][0][1]<0
    data=build(qb,E,N,abox)
    assert data["G"][1][1]<0 and data["disc"][0][1]<dec("-55")
    return data

hL=endpoint_data(qHL,idec("0.2952194200028029",
                           "0.2952194200028030"))
hR=endpoint_data(qHR,idec("0.2952194200019594",
                           "0.2952194200019595"))
assert hL["H"][0][1]<0<hR["H"][0][0]

getcontext().prec=18
def show(x):
    return (Decimal(x[0].numerator)/Decimal(x[0].denominator),
            Decimal(x[1].numerator)/Decimal(x[1].denominator))

print("PASS: exact physical s=4 transverse Hopf point")
print("Hopf q bracket =",qHL,qHR)
print("G lower alpha-face interval =",show(lower_face["G"][0]))
print("G upper alpha-face interval =",show(upper_face["G"][0]))
print("localized G_alpha interval =",show(hd["G"][1]))
print("Hopf bracket H intervals =",show(hL["H"][0]),show(hR["H"][0]))
print("localized implicit alpha_q interval =",show(alpha_q))
print("localized implicit dH/dq interval =",show(hopf_slope))
print("localized discriminant interval =",show(hd["disc"][0]))
print("localized minimum phase-weight lower bound =",
      min(show(w[0])[0] for w in hd["weights"]))
print("localized critical amplitudes =",
      [show(c[0]) for c in hd["c"]])
print("localized alpha collar =",show(hopf_a))
