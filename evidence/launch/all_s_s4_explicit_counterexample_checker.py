"""Fail-closed exact checker for the finite s=4 counterexample inputs.

The infinite-dimensional Banach contractions are proved in the accompanying
human proof.  This checker certifies their finite algebraic/interval premises
and every dyadic exponent comparison used to select the two contractions.
It uses no floating-point arithmetic.
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


def run_quiet(filename, module_name):
    env = {"__name__": module_name}
    with redirect_stdout(StringIO()):
        exec(compile(Path(filename).read_text(), filename, "exec"), env)
    return env


# This imports, executes, and checks the transverse-Hopf, local-majorant,
# intrinsic-lift, and scalar-to-raw tangent Fraction certificates.
tangent = run_quiet(
    "evidence/launch/all_s_s4_scalar_to_raw_tangent_exact.py",
    "__explicit_s4_counterexample__",
)
local = tangent["local"]
hopf = local["env"]

assert hopf["hL"]["H"][0][1] < 0 < hopf["hR"]["H"][0][0]
assert hopf["hopf_slope"][0] > hopf["dec"]("109.89")
assert hopf["hopf_slope"][1] < hopf["dec"]("110.82")
assert hopf["crossing_lower"] > F(194, 100)
assert hopf["hd"]["disc"][0][1] < hopf["dec"]("-55.40")
assert min(w[0][0] for w in hopf["hd"]["weights"]) > F(1, 2**19)
assert hopf["hd"]["c"][0][0][0] > F(1, 2**15)
assert hopf["hd"]["c"][1][0][0] > F(1, 4)
hopf_vector = tangent["hopf_vector"]
assert hopf_vector["omega2_lower"] > F(93, 25)**2
assert hopf_vector["omega2_upper"] < F(931, 250)**2
assert max(
    hopf_vector["abs_upper"](x[0])
    for row in hopf_vector["D"] for x in row
) < 2**4
assert max(
    sum(hopf_vector["abs_upper"](x[0]) for x in row)
    for row in hopf_vector["W"]
) < 2**3
assert tangent["v_star_upper"] < 2**500001
assert tangent["raw_inf_forward"] * 8 < 2**5000


# Binary exponents in the explicit table.
D = 3_000_000
V = 200_000
P0 = 18_100_000
Z1 = 21_200_000
Z2 = 64_000_000
MH = 70_000_000
BH = LH = 70_000_000
r0 = 500_000_000
eps = 300_000_000
KP = 80_000_000
K = 100_000_000
delta_x_inv = 81_000_000
L = 250_000_000
R = 410_000_000
M = 600_000_000
N = 1_200_000_000

# Direct operation ledger for (C.4.10)--(C.4.18).
for bound, exponent in zip(
    local["root_q_derivative_bounds"],
    (58, 175, 294, 420, 550),
):
    assert bound < 2**exponent
assert local["root_t_lower"] > F(1, 2**14)
assert local["family_derivative_upper"] < 2**44

coefficient_derivative_exp = 1024
weight_quotient_exp = 6 * coefficient_derivative_exp + 14 + 10
assert weight_quotient_exp < 10_000
for order in range(6):
    assert 10_000 * (order + 1) <= 60_000

moment_solve_exp = max(
    (words + 1) * 164 + words * 2_000 + 15
    for words in range(6)
)
assert moment_solve_exp < 20_000
assert local["shifted_trace_upper"] < 2**26
assert local["gram_det_lower"] > F(1, 2**76)
assert local["shifted_gram_inverse_bound"] < 2**160
assert local["signed_gram_inverse_bound"] < 2**164
assert local["normalized_factor_lower"] == F(1, 2**30)
assert local["fixed_Z_lower"] == F(1, 2**79)
assert (
    local["normalized_factor_lower"] - local["factor_variation"]
    > F(1, 2**31)
)
assert local["fixed_Z_lower"] - local["Z_variation"] > F(1, 2**80)

for order in range(6):
    reciprocal_exp = 40_020 * order + 500
    one_block_exp = 40_020 + reciprocal_exp + 10
    assert one_block_exp < 50_000 * (order + 1)
composition_exp = max(50_000 * (5 + 2*r + 1) for r in range(6))
assert composition_exp == 800_000
assert composition_exp + 100_000 + 60_000 + 20_000 < 1_000_000
assert 1_000_000 + 3 * 500_001 + 100 < 3_000_000

# Projectors, moving equilibrium, augmented Hopf inverse, and C2 equation
# Lipschitz bound (which uses the certified C3 field bound).
projector_inf_exp = 6 * (D + 1)
# For a 7 x 7 matrix, ||A||_2 <= sqrt(7) ||A||_inf < 2^2 ||A||_inf.
projector_two_exp = projector_inf_exp + 2
assert projector_two_exp < P0
assert P0 + 3 < MH
assert P0 + 3 + D < Z1
assert P0 + 3 + D + 2 + 2 * Z1 < Z2
assert 2 * P0 + D + Z1 + 2 * V + 1 < MH
assert 20 + D + max(Z1, 2 * Z1, Z2) + 3 * V < BH

# Lyapunov--Schmidt radius and contraction.
R0_signed_exp = 1 + MH + BH - r0
assert R0_signed_exp == -359_999_999
ls_contraction_exp = MH + LH + (MH + BH + 2) - r0
assert ls_contraction_exp == -219_999_998
assert ls_contraction_exp < -2
assert R0_signed_exp + Z1 < -338_000_000
assert D - 330_000_000 < -eps

# Floquet roughness, tube transport, exact-map Taylor bounds.
assert 4 + (P0 + 3) - eps < -281_000_000
assert 2 + 2 * (P0 + 9) - eps < -263_000_000
assert KP - delta_x_inv < -100
assert 20 + 3 * KP + D < L
assert 20 + 5 * KP + D < R

# Quantitative Lyapunov--Perron inequalities with beta=3/4, sigma=1/4.
assert M > delta_x_inv
assert 4 + K + R + 1 < M
assert 9 + K + L < M
assert 8 + K + R < N
assert N == 2 * M
assert N > 2_000_000
assert N > 2 * 18_200_000 + 1

print("PASS: exact finite inputs for the explicit s=4 counterexample")
print("unique transverse Hopf branch and positive physical margins")
print("intrinsic and scalar-to-raw tangent certificates passed")
print("raw C3/C1, Lyapunov--Schmidt, Floquet, and Lyapunov--Perron exponents passed")
print("selected amplitude exponent = 500000000")
print("selected shadowing index exponent = 1200000000")
