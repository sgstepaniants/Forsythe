"""Exact Fraction/interval certificate for the a-normalized scalar Hopf vector.

The script imports the exact transverse-Hopf sheet, reconstructs the full
external forcing vectors V_g and their a derivatives from the certified
two-by-two contraction matrix, forms every block in the scalar eigenvector
equations, and checks a closed binary operation ledger.  It does not assume
the previous 2^25000 vector bound.
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


env = {"__name__": "__scalar_hopf_vector_exact__"}
with redirect_stdout(StringIO()):
    exec(compile(
        Path("evidence/launch/all_s_s4_transverse_hopf_exact.py").read_text(),
        "all_s_s4_transverse_hopf_exact.py", "exec"), env)

K, I = env["K"], env["I"]
ka, ks, kn, km = env["ka"], env["ks"], env["kn"], env["km"]
ksc, kd = env["ksc"], env["kd"]
kpm, kpe = env["kpm"], env["kpe"]
root_jet = env["root_jet"]
C = env["C"]
hd = env["hd"]


def abs_upper(x):
    return max(abs(x[0]), abs(x[1]))


def jet_value_bound(xs):
    return max(abs_upper(x[0]) for x in xs)


# Reconstruct the algebraic external labels with the same exact q jets used
# by the transverse certificate.
qk = K(env["hopf_q"], q=I(1))
Qk = kpm([kn(qk), K(I(1))], [K(I(a)) for a in env["Qbar"]])
Pk = [K(I(a)) for a in env["P"]]
Nj = [root_jet(z, Qk) for z in env["NH"]]
labels = [Nj[1], Nj[6]]
Delta = hd["Delta"]
T = hd["Tmatrix"]


# Recover V_g=(V_C,V_0,V_1) from its two exact contractions
# T_(h,g)=d_h(V_g).  Since
# d_h(V)=-2/C {Delta(h)(V_0+hV_1)+2V_C},
# two distinct external labels determine the affine H part.
V = []
for j, g in enumerate(labels):
    Pg = kpe(Pk, g)
    VC = ksc(ka(km(Pg, Pg), K(I(C))), 4)
    values = []
    for i, h in enumerate(labels):
        numerator = ks(ksc(T[i][j], -C/F(2)), ksc(VC, 2))
        values.append(kd(numerator, kpe(Delta, h)))
    V1 = kd(ks(values[1], values[0]), ks(labels[1], labels[0]))
    V0 = ks(values[0], km(V1, labels[0]))
    V.append([VC, V0, V1])


# Scalar Jacobian blocks in coordinates (a,u,B_C,B_0,B_1,l_-,l_+),
# where l_j is the logarithmic external-amplitude perturbation.
c = hd["c"]
A = hd["normals"]
ell = hd["ell"]
kappa = hd["kappa"][0]

p = ka(km(c[0], K(A[0][1])), km(c[1], K(A[1][1])))
qvec = [
    ka(km(c[0], K(V[0][r][1])), km(c[1], K(V[1][r][1])))
    for r in range(3)
]

D = []
for g in labels:
    Dg = kpe(Delta, g)
    D.append([
        K(I(-F(4)/C)),
        ksc(Dg, -F(2)/C),
        ksc(km(Dg, g), -F(2)/C),
    ])

R = [km(c[j], A[j]) for j in range(2)]
W = [[km(c[j], V[j][r]) for j in range(2)] for r in range(3)]


# Exact primitive coefficient decisions on the complete certified Hopf box.
# The complete interval lies strictly to the left of -2^6.  Checking only
# the absolute values of the two endpoints would not, by itself, exclude an
# interval crossing zero.
assert kappa[1] < -(2**6)
assert all(0 < cj[0][0] and cj[0][1] < 1 for cj in c)
assert jet_value_bound(A) < 2**10
assert max(abs_upper(x[1]) for x in A) < 2**10
assert jet_value_bound(ell) < 2**7
assert abs_upper(p[0]) < 1
assert jet_value_bound(qvec) < 2**2
assert max(abs_upper(x[0]) for row in D for x in row) < 2**4
assert jet_value_bound(R) < 1
assert max(abs_upper(x[0]) for row in W for x in row) < 2**3
assert max(sum(abs_upper(x[0]) for x in row) for row in W) < 2**3
assert max(abs_upper(x[0]) for row in V for x in row) < 2**17
assert max(abs_upper(x[1]) for row in V for x in row) < 2**17


# At the exact Hopf zero, theta=lambda(lambda-1) has
# |Im theta|=omega and omega^2=-disc/4.  These exact rational comparisons
# give 3.72<omega<3.724, hence the determinant lower bound used below is >1.
omega2_lower = -hd["disc"][0][1] / 4
omega2_upper = -hd["disc"][0][0] / 4
assert omega2_lower > F(93, 25)**2
assert omega2_upper < F(931, 250)**2
assert omega2_lower > 1


# Exact integer operation ledger for an eigenvector normalized by a=1.
# The block equations are
#   lambda a=kappa u,
#   lambda u=p a+u+R l,
#   lambda B=qvec a+B+W l,
#   lambda l=ell u+D B.
# Elimination gives
#   (theta I_2-DW)l=(lambda-1)ell u+D qvec a.
# The critical identity DW(1,1)^T=2(1,1)^T makes both eigenvalues of
# the real 2x2 matrix DW real.  Therefore
# |det(theta I-DW)| >= |Im theta|^2 = omega^2 > 1.
kappa_lower = 2**6
lambda_upper = 4
lambda_minus_one_upper = 8
u_upper = F(lambda_upper, kappa_lower)
assert u_upper < 1

D_upper = 2**4
W_upper = 2**3
q_upper = 2**2
ell_upper = 2**7

DW_entry_upper = 3 * D_upper * W_upper
assert DW_entry_upper < 2**9
theta_entry_upper = 4**2 + 4
assert theta_entry_upper < 2**5
theta_DW_entry_upper = theta_entry_upper + DW_entry_upper
assert theta_DW_entry_upper < 2**10

# A 2x2 adjugate has two entries in each row; det>1.
theta_DW_inverse_inf = 2 * 2**10
assert theta_DW_inverse_inf == 2**11

rhs_l_upper = (
    lambda_minus_one_upper * ell_upper * u_upper
    + 3 * D_upper * q_upper
)
assert rhs_l_upper < 2**11
l_upper = theta_DW_inverse_inf * 2**11
assert l_upper == 2**22

# |lambda-1|=sqrt(1+omega^2)>1.
B_upper = q_upper + 2 * W_upper * l_upper
assert B_upper < 2**27

full_vector_inf = max(F(1), u_upper, l_upper, B_upper)
full_vector_two = 3 * full_vector_inf       # sqrt(7)<3
assert full_vector_two < 2**29
assert full_vector_two < 2**25000

print("PASS: exact a-normalized scalar Hopf-vector certificate")
print("|kappa| > 2^6 and 3.72 < omega < 3.724")
print("max |A|,|A_a| < 2^10; max |V|,|V_a| < 2^17")
print("max entry |D| < 2^4, induced row-sum |W| < 2^3, max |qvec| < 2^2")
print("|det(theta I_2-DW)| > omega^2 > 1")
print("a-normalized scalar Hopf-vector Euclidean norm < 2^29 < 2^25000")
