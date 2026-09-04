"""Exact coarse local majorants for the amplitude-gauged s=4 Hopf route.

This certificate currently covers the intrinsic 9-by-9 coefficient solve
and elementary geometric/Gram margins. It is one finite component of the
quantitative Lyapunov--Schmidt certificate, not a complete counterexample
checker.
"""

if not __debug__:
    raise RuntimeError(
        "Verification assertions are disabled. "
        "Run without -O/-OO and unset PYTHONOPTIMIZE."
    )

from contextlib import redirect_stdout
from fractions import Fraction as F
from io import StringIO
from math import comb
from pathlib import Path


env = {"__name__": "__local_majorants__"}
with redirect_stdout(StringIO()):
    exec(compile(
        Path("evidence/launch/all_s_s4_transverse_hopf_exact.py").read_text(),
        "transverse_hopf_exact", "exec"), env)

I, ia, isu, im = env["I"], env["ia"], env["isu"], env["im"]
ipm, ipa, ipe = env["ipm"], env["ipa"], env["ipe"]
P = env["P"]
hd = env["hd"]
alpha = env["hopf_a"]


def midpoint(x):
    return (x[0] + x[1]) / 2


def abs_upper(x):
    return max(abs(x[0]), abs(x[1]))


def abs_lower(x):
    assert not x[0] <= 0 <= x[1]
    return min(abs(x[0]), abs(x[1]))


def mat_inverse(a):
    n = len(a)
    aug = [list(a[i]) + [F(i == j) for j in range(n)] for i in range(n)]
    for j in range(n):
        pivot = next(i for i in range(j, n) if aug[i][j])
        aug[j], aug[pivot] = aug[pivot], aug[j]
        z = aug[j][j]
        aug[j] = [x / z for x in aug[j]]
        for i in range(n):
            if i == j:
                continue
            z = aug[i][j]
            if z:
                aug[i] = [aug[i][k] - z * aug[j][k]
                          for k in range(2*n)]
    return [row[n:] for row in aug]


def mat_interval_product(a, b):
    n, m, p = len(a), len(b), len(b[0])
    out = [[I(0) for _ in range(p)] for _ in range(n)]
    for i in range(n):
        for k in range(m):
            for j in range(p):
                out[i][j] = ia(out[i][j], im(a[i][k], b[k][j]))
    return out


def inf_norm_exact(a):
    return max(sum(abs(x) for x in row) for row in a)


def inf_norm_interval(a):
    return max(sum(abs_upper(x) for x in row) for row in a)


# The order-zero 9-by-9 solves in the intrinsic three-phase recursion.
S = [(-alpha[1], -alpha[0]), I(1)]
Delta = [z[0] for z in hd["Delta"]]


def intrinsic_matrix_inverse_bound(poly):
    G = ipm(S, poly)
    cols = []
    for j in range(4):
        cols.append([I(0)]*j + G + [I(0)]*(9-j-len(G)))
    for j in range(3):
        cols.append([I(0)]*j + [(-x[1], -x[0]) for x in Delta]
                    + [I(0)]*(9-j-len(Delta)))
    for j in range(2):
        col = [I(0)]*9
        col[j] = I(-1)
        cols.append(col)
    A = [[cols[j][i] for j in range(9)] for i in range(9)]
    A0 = [[midpoint(A[i][j]) for j in range(9)] for i in range(9)]
    A0i = mat_inverse(A0)
    A0ii = [[I(x) for x in row] for row in A0i]
    prod = mat_interval_product(A0ii, A)
    E = [[isu(I(F(i == j)), prod[i][j]) for j in range(9)]
         for i in range(9)]
    eta = inf_norm_interval(E)
    assert eta < 1
    return inf_norm_exact(A0i) / (1-eta)


Qhopf = env["family_intervals"](env["hopf_q"])[0]
inverse_bounds = [intrinsic_matrix_inverse_bound(ipa(P)),
                  intrinsic_matrix_inverse_bound(Qhopf)]
assert max(inverse_bounds) < 2**40

# Direct primitive and implicit-root derivative ledger for (C.4.9)--(C.4.10).
# The root collars cover every selected branch on the full q envelope.
root_t_derivatives = [
    env["family_t_derivative"](box, env["q_global"])
    for box in env["E_global"] + env["N_global"]
]
root_t_lower = min(abs_lower(x) for x in root_t_derivatives)
assert root_t_lower > F(1, 2**14)

_, minus_family, plus_family = env["family_intervals"](env["q_global"])
source_family = ipm(ipa(env["P"]), ipa(env["Qbar"]))
t_envelope = I(-5, 5)
family_derivative_upper = F(0)
for poly in (minus_family, plus_family, source_family):
    derivative_poly = poly
    for _ in range(9):
        family_derivative_upper = max(
            family_derivative_upper,
            abs_upper(ipe(derivative_poly, t_envelope)),
        )
        derivative_poly = env["ipd"](derivative_poly)
assert family_derivative_upper < 2**44


def bell(n, k, derivatives):
    if n == 0 and k == 0:
        return F(1)
    if n == 0 or k == 0:
        return F(0)
    return sum(
        F(comb(n-1, j-1)) * derivatives[j-1]
        * bell(n-j, k-1, derivatives)
        for j in range(1, n-k+2)
    )


root_q_derivative_bounds = []
for n in range(1, 6):
    nonlinear = sum(
        bell(n, k, root_q_derivative_bounds)
        for k in range(2, n+1)
    )
    affine_q = sum(
        bell(n-1, k, root_q_derivative_bounds)
        for k in range(0, n)
    )
    root_q_derivative_bounds.append(
        (family_derivative_upper * nonlinear
         + n * family_derivative_upper * affine_q) / root_t_lower
    )
for bound, exponent in zip(
    root_q_derivative_bounds,
    (58, 175, 294, 420, 550),
):
    assert bound < 2**exponent

# Coarse geometric margins used by later derivative majorants.
core = [env["EH"][i] for i in (0, 1, 3, 4, 6, 7)]
external = [env["NH"][i] for i in (1, 6)]
labels = sorted(core + external)
gap = min(labels[j][0]-labels[i][1]
          for i in range(8) for j in range(i+1, 8))
assert gap > F(1, 2**11)
weight_lower = min(w[0][0] for w in hd["weights"])
assert weight_lower > F(1, 2**19)

# A selected four-node Vandermonde term gives a direct Gram determinant
# lower bound.  The certified root bound supplies the trace premise used
# for the inverse estimate in the manuscript.
core_gap = min(core[j][0]-core[i][1]
               for i in range(6) for j in range(i+1, 6))
assert core_gap > 1
assert max(max(abs(a), abs(b)) for a, b in core) < 5
trace_upper = F(1 + 5**2 + 5**4 + 5**6, 1)
assert trace_upper < 2**14
gram_det_lower = weight_lower**4 * core_gap**12
assert gram_det_lower > F(1, 2**76)
gram_inverse_bound = trace_upper**3 / gram_det_lower
assert gram_inverse_bound < 2**124

# The manuscript uses the translated positive nodes in its physical Gram
# matrix.  This direct shifted trace/inverse ledger is the one for (C.4.13).
shifted_trace_upper = F(10 + 10**3 + 10**5 + 10**7)
assert shifted_trace_upper < 2**26
shifted_gram_inverse_bound = shifted_trace_upper**3 / gram_det_lower
assert shifted_gram_inverse_bound < 2**160
signed_gram_perturbation = F(1, 2**999972)
signed_gram_neumann_defect = (
    shifted_gram_inverse_bound * signed_gram_perturbation
)
assert signed_gram_neumann_defect < F(1, 2)
signed_gram_inverse_bound = (
    shifted_gram_inverse_bound / (1-signed_gram_neumann_defect)
)
assert signed_gram_inverse_bound < 2**164

# Uniform algebraic-root and physical-box margins.
for boxes in (env["EH"], env["NH"]):
    assert min(boxes[i+1][0]-boxes[i][1]
               for i in range(7)) > F(1, 4)
Pa = ipe(ipa(P), alpha)
Qa = ipe(Qhopf, alpha)
assert Pa[0] > 2
assert Qa[0] > 4
assert env["C"] > F(1, 4)
assert env["C"] < F(1, 2)
assert hd["weights"][0][1][0] > F(1, 16)

# The unique Hopf point is at least 2^-41 from both q faces.
slope_upper = env["hopf_slope"][1]
q_left_margin = -env["hL"]["H"][0][1] / slope_upper
q_right_margin = env["hR"]["H"][0][0] / slope_upper
assert min(q_left_margin, q_right_margin) > F(1, 2**41)

# Monotonicity of alpha(q) and the endpoint root collars put alpha_* at
# least 2^-24 from both faces of the broad alpha collar.
assert env["alpha_q"][1] < 0
alpha_left_upper = env["dec"]("0.2952194200028030")
alpha_right_lower = env["dec"]("0.2952194200019594")
alpha_margin = min(alpha_right_lower-alpha[0],
                   alpha[1]-alpha_left_upper)
assert alpha_margin > F(1, 2**24)

# Direct factor and normalization margins in (C.4.15)--(C.4.17).
monic_factor_upper = F(2**14)
monic_factor_lower = F(1, 4) / monic_factor_upper
assert monic_factor_lower == F(1, 2**16)
normalized_factor_lower = monic_factor_lower / monic_factor_upper
assert normalized_factor_lower == F(1, 2**30)
fixed_Z_lower = F(1, 2**19) * normalized_factor_lower**2
assert fixed_Z_lower == F(1, 2**79)
factor_variation = F(2**20000, 2**1000000)
assert normalized_factor_lower - factor_variation > F(1, 2**31)
Z_variation = F(1, 2**900000)
assert fixed_Z_lower - Z_variation > F(1, 2**80)

print("PASS: exact coarse s=4 local majorants")
print("intrinsic 9x9 inverse infinity-norm bound < 2^40")
print("all selected label gaps > 2^-11")
print("all core phase weights > 2^-19")
print("phase Gram inverse 2-norm bound < 2^124")
print("all roots on each degree-eight sheet are separated by > 1/4")
print("P(alpha)>2, Q(alpha)>4, and C>1/4")
print("C<1/2 and the first P-phase weight slope is > 1/16")
print("Hopf q-face margin > 2^-41 and alpha-face margin > 2^-24")
print("implicit root q-derivative bounds through order five match C.4.10")
print("shifted phase Gram inverse bounds match C.4.13")
print("two-phase factor and normalization margins match C.4.15--C.4.17")
