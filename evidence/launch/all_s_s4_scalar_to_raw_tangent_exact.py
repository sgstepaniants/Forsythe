"""Exact coarse certificate for the scalar-to-raw weighted tangent map.

This script does not sample a map.  It imports the exact rational interval
certificates for the two 9-by-9 phase solves and the intrinsic 4-by-4 lift,
checks the primitive margins used in the reconstruction, and checks the
integer exponent arithmetic in the forward tangent majorant.
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


def run_quiet(path, name):
    env = {"__name__": name}
    with redirect_stdout(StringIO()):
        exec(compile(Path(path).read_text(), path, "exec"), env)
    return env


local = run_quiet(
    "evidence/launch/all_s_s4_local_majorants_exact.py",
    "__scalar_to_raw_local__",
)
lift = run_quiet(
    "evidence/launch/all_s_s4_intrinsic_lift_exact.py",
    "__scalar_to_raw_lift__",
)


def abs_upper(x):
    return max(abs(x[0]), abs(x[1]))


# Certified primitive inverse and geometry bounds.
assert max(local["inverse_bounds"]) < 2**40
assert lift["inv4"] < 2**20
assert lift["eta4"] < F(1, 4)
assert local["core_gap"] > 1
assert local["gram_inverse_bound"] < 2**124
assert local["hd"]["weights"][0][1][0] > F(1, 16)
assert all(0 < w[0][0] and w[0][1] < 1
           for w in local["hd"]["weights"])
assert all(0 < c[0][0] and c[0][1] < 1
           for c in local["hd"]["c"])

# J has output order (delta C,delta H0,delta H1,delta U).  Therefore the
# scalar input order (u,B_C,B_0,B_1) uses the permutation
# (B_C,B_0,B_1,u); its u column is exactly the certified p_U.
p_u = lift["lift"]
assert -3 < p_u[0][0] and p_u[0][1] < -2
assert -2 < p_u[1][0] and p_u[1][1] < -F(1, 2)
assert F(1, 2) < p_u[2][0] and p_u[2][1] < 1
assert F(1, 10) < p_u[3][0] and p_u[3][1] < F(1, 4)

P = lift["P"]
Q = lift["Q"]
core = local["core"]
alpha = lift["alpha"]
assert max(abs_upper(x) for x in P) < 2**10
assert max(abs_upper(x) for x in Q) < 2**10
assert max(max(abs(x[0]), abs(x[1])) for x in core) < 5
assert max(max(abs(x[0]), abs(x[1])) for x in local["external"]) < 5
assert max(abs(alpha[0]), abs(alpha[1])) < 1
assert alpha[0] - local["external"][0][1] > 1
assert local["external"][1][0] - alpha[1] > 1

# Exact operation-count majorants from the displayed reconstruction.
# Norms on coefficient vectors and finite matrices are infinity norms.
M_inv = 2**40
J_inv = 2**20
poly_coeff = 2**10

# b=-S p Q has fewer than 2*4 convolution summands per coefficient.
b_per_p = 2**14
x_per_p = M_inv * b_per_p
assert x_per_p == 2**54

# M_a has only the four shifted P columns.  The deliberately weakened
# bound 2^14 covers its row sum.  Differentiating Mx=b at fixed p gives
# x_a=M^{-1}(b_a-M_a x).
Ma = 2**14
ba_per_p = 2**12
xa_per_p = M_inv * (ba_per_p + Ma * x_per_p)
assert xa_per_p < 2**109

# Extraction of (delta C,delta H0,delta H1,delta U) from x costs less
# than 2^3, including its a derivative. Much weaker powers suffice below.
J_per_p = 2**57
Ja_per_p = 2**113
assert 8 * x_per_p <= J_per_p
assert 8 * (xa_per_p + x_per_p) <= Ja_per_p

p_per_c = J_inv
pa_per_c = J_inv * Ja_per_p * p_per_c
assert pa_per_c < 2**154

# Inverse spectral formula:
#   delta w_i=w_i(delta C/C-p(x_i)/P(x_i)).
# A lower-degree quartic perturbation has degree at most three, so its
# evaluation factor on |x|<5 is 1+5+25+125=156<2^8.
eval3 = 1 + 5 + 25 + 125
assert eval3 < 2**8

# On the sheet P(x_i)Q(x_i)=C.  Thus |Q(x_i)|<2^19 and C>1/4
# give 1/|P(x_i)|<2^21.  All physical core weights are below one.
eval_Q = 5**4 + poly_coeff * eval3
assert eval_Q < 2**19
assert local["env"]["C"] > F(1, 4)
inv_P_at_core = 2**21
rho_per_c = 4 + 2**8 * p_per_c * inv_P_at_core
rhoa_per_c = (
    2**19 * rho_per_c
    + 2**8 * pa_per_c * inv_P_at_core
)
assert rho_per_c < 2**50
assert rhoa_per_c < 2**184

# |Delta'(x_i)|>1 and s_1>2^-4 give |s_i/s_1|<2^23.
s_ratio = 2**23
K_per_c = (1 + s_ratio) * rho_per_c
Ka_per_c = (1 + s_ratio) * rhoa_per_c
assert K_per_c < 2**74
assert Ka_per_c < 2**208

# Direct kernel count for V_h.  In the certified monomial basis, four
# Gram-inverse contractions cost <2^137.  Multiplication by the external
# P value, then by a quartic (with its coefficient bound deliberately
# weakened to 2^24), costs the following bounds.  The
# external-to-alpha gap is >1, so the A_h divided-difference term is
# smaller.  Two-step monic division by Delta costs one factor <2^24.
kernel_coeff = 4 * 2**124 * 2**10
assert kernel_coeff < 2**137
kernel_forcing = 4 * 2**19 * kernel_coeff
assert kernel_forcing < 2**158
product_forcing = 5 * kernel_forcing * 2**24
assert product_forcing < 2**185
A_external = (4 * 2**19) * (2**19)
assert A_external == 2**40
delta_divdiff = 2**40
forcing_numerator = 2 * product_forcing + 2**41 + A_external * delta_divdiff
assert forcing_numerator < 2**187
V_bound = 2**24 * forcing_numerator
assert V_bound < 2**212

# At the critical equilibrium u_*=0 and B_*=-sum_h c_h V_h, with two
# amplitudes 0<c_h<1.
critical_c = 2**213
v_star_upper = K_per_c * critical_c + 2
assert v_star_upper < 2**500001
raw_inf_forward = 1 + K_per_c + Ka_per_c * critical_c
assert raw_inf_forward < 2**422

# Passing from the seven-dimensional infinity norm to Euclidean operator
# norm costs less than 2^3.  This intentionally leaves a very large margin.
assert 2**3 * raw_inf_forward < 2**5000
hopf_vector = run_quiet(
    "evidence/launch/all_s_s4_scalar_hopf_vector_exact.py",
    "__scalar_to_raw_hopf_vector__",
)
assert hopf_vector["full_vector_two"] < 2**29
assert 2**5000 * hopf_vector["full_vector_two"] < 2**5029
assert 2**5029 < 2**30000

print("PASS: exact coarse scalar-to-raw weighted tangent certificate")
print("first-phase inverse infinity norm < 2^40")
print("intrinsic lift inverse infinity norm < 2^20")
print("scalar order (u,B_C,B_0,B_1) maps to lift order (B_C,B_0,B_1,u)")
print("first P-phase weight slope s_1 > 1/16")
print("scalar-to-raw Euclidean forward tangent norm < 2^5000")
print("proved scalar Hopf-vector norm < 2^29")
print("a-normalized scalar Hopf vector transfers with raw norm < 2^5029 < 2^30000")
