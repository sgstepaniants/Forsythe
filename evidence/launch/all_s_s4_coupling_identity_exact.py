"""Independent exact-interval check of the Part C coupling signs.

The transverse-Hopf certificate evaluates the expanded formula (C.13.2)
directly.  This checker takes a different route: it constructs the source
polynomials p_g and q_g, divides the left side of (C.3.9) by Delta to obtain
(V_g)^H, evaluates (C.3.10), and compares the result with both (C.13.2) and
the certified coupling matrix.  It shares the certified algebraic root boxes
and rational-interval primitives with the transverse script, but it does not
use that script's T-matrix to construct (V_g)^H.
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


env = {"__name__": "__coupling_identity_inputs__"}
source = Path("evidence/launch/all_s_s4_transverse_hopf_exact.py")
source_text = source.read_text()
definitions, marker, _ = source_text.partition("# Localized rational IFT rectangle.")
assert marker
with redirect_stdout(StringIO()):
    exec(compile(definitions, str(source), "exec"), env)

I, K = env["I"], env["K"]
ia, iz = env["ia"], env["iz"]
ka, ks, kn, km = env["ka"], env["ks"], env["kn"], env["km"]
kd, ksc, kpe = env["kd"], env["ksc"], env["kpe"]
kpm, kpdivmonic, kproot = env["kpm"], env["kpdivmonic"], env["kproot"]
root_jet = env["root_jet"]
C = env["C"]
ZERO, ONE = K(I(0)), K(I(1))


def poly_add(p, q):
    n = max(len(p), len(q))
    return [
        ka(p[i] if i < len(p) else ZERO, q[i] if i < len(q) else ZERO)
        for i in range(n)
    ]


def poly_sub(p, q):
    n = max(len(p), len(q))
    return [
        ks(p[i] if i < len(p) else ZERO, q[i] if i < len(q) else ZERO)
        for i in range(n)
    ]


def poly_scalar(p, scalar):
    return [km(scalar, coefficient) for coefficient in p]


def all_components_contain_zero(polynomial):
    return all(iz(component) for coefficient in polynomial for component in coefficient)


def intervals_overlap(x, y):
    return x[0] <= y[1] and y[0] <= x[1]


def jets_overlap(x, y):
    return all(intervals_overlap(x[i], y[i]) for i in range(5))


def quotient_by_linear(polynomial, root):
    quotient, remainder = kpdivmonic(polynomial, [kn(root), ONE])
    assert all_components_contain_zero(remainder)
    return quotient


def u_polynomial(R, B, alpha):
    beta = kd(kpe(B, alpha), kpe(R, alpha))
    numerator = poly_sub(poly_scalar(R, beta), B)
    numerator = [ksc(coefficient, F(1, 1) / C) for coefficient in numerator]
    return quotient_by_linear(numerator, alpha)


def kernel_polynomial(R, B, alpha, g):
    U = u_polynomial(R, B, alpha)
    numerator = poly_sub(poly_scalar(R, kpe(U, g)), poly_scalar(U, kpe(R, g)))
    return quotient_by_linear(numerator, g)


# Reconstruct the same algebraic sheet from the certified input boxes. Loading
# only the definition prefix above deliberately avoids executing the original
# localized Hopf assertions before this independent comparison can run.
q_left = env["dec"]("-1.851888497035")
q_right = env["dec"]("-1.851888497033")
hopf_q = (q_left, q_right)
EH, NH = env["root_boxes"](hopf_q)
hopf_a = env["idec"]("0.2952193", "0.29521955")
implemented_matrix = env["build"](hopf_q, EH, NH, hopf_a)["Tmatrix"]

# Do not use the implemented matrix, Delta, or H polynomial to construct V_g.
q = K(hopf_q, q=I(1))
Pk = [K(I(coefficient)) for coefficient in env["P"]]
Qk = kpm([kn(q), ONE], [K(I(coefficient)) for coefficient in env["Qbar"]])
Ej = [root_jet(root, Qk) for root in EH]
Nj = [root_jet(root, Qk) for root in NH]

product = kpm(Pk, Qk)
minus = list(product)
minus[0] = ks(minus[0], K(I(C)))
degree_seven, first_remainder = kproot(minus, Ej[2])
Delta, second_remainder = kproot(degree_seven, Ej[5])
assert iz(first_remainder[0]) and iz(second_remainder[0])

_, BP = kpdivmonic(Delta, Pk)
_, BQ = kpdivmonic(Delta, Qk)
alpha = K(hopf_a, a=I(1))
labels = [Nj[1], Nj[6]]
rq_numerator = list(Delta)
rq_numerator[0] = ks(rq_numerator[0], kpe(Delta, alpha))
rq = quotient_by_linear(rq_numerator, alpha)

direct_matrix = [[None, None], [None, None]]
expanded_matrix = [[None, None], [None, None]]
vh_polynomials = []

for source_index, g in enumerate(labels):
    Pg = kpe(Pk, g)
    Qg = kpe(Qk, g)
    KP = kernel_polynomial(Pk, BP, alpha, g)
    KQ = kernel_polynomial(Qk, BQ, alpha, g)
    pg = poly_scalar(KP, ksc(Pg, -4))
    qg = poly_scalar(KQ, ksc(Pg, 4))
    VC = ksc(ka(km(Pg, Pg), K(I(C))), 4)
    Ag = ksc(
        km(
            kd(Pg, ks(g, alpha)),
            ks(kd(Pg, kpe(Pk, alpha)), kd(Qg, kpe(Qk, alpha))),
        ),
        4,
    )

    # Literal left side of (C.3.9): p_g Q + P q_g - VC - A_g rq.
    numerator = poly_add(kpm(pg, Qk), kpm(Pk, qg))
    numerator = poly_sub(numerator, [VC])
    numerator = poly_sub(numerator, poly_scalar(rq, Ag))
    VH, remainder = kpdivmonic(numerator, Delta)
    assert len(VH) == 2
    assert all_components_contain_zero(remainder)
    vh_polynomials.append(VH)

    for target_index, h in enumerate(labels):
        # (C.3.10), evaluated after obtaining V_g from (C.3.9).
        direct = ksc(
            ka(km(kpe(Delta, h), kpe(VH, h)), ksc(VC, 2)),
            -F(2, 1) / C,
        )
        direct_matrix[target_index][source_index] = direct

        # Independently expanded formula printed as (C.13.2).
        expanded_bracket = ka(
            ka(km(kpe(pg, h), kpe(Qk, h)), km(kpe(Pk, h), kpe(qg, h))),
            VC,
        )
        expanded_bracket = ks(expanded_bracket, km(Ag, kpe(rq, h)))
        expanded = ksc(expanded_bracket, -F(2, 1) / C)
        expanded_matrix[target_index][source_index] = expanded

        certified = implemented_matrix[target_index][source_index]
        assert jets_overlap(direct, expanded)
        assert jets_overlap(direct, certified)
        assert jets_overlap(expanded, certified)

# The four value intervals are separated from zero with the sign pattern
# (target rows h_-,h_+) by (source columns h_-,h_+): [-,+; +,-].
expected_signs = ((-1, 1), (1, -1))
for i in range(2):
    for j in range(2):
        for matrix in (direct_matrix, expanded_matrix, implemented_matrix):
            value = matrix[i][j][0]
            assert not iz(value)
            assert (1 if value[0] > 0 else -1) == expected_signs[i][j]

print("PASS: independent exact-interval coupling identity")
print("(C.3.9) polynomial division, (C.3.10), and (C.13.2) agree")
print("coupling order: target rows (h_-,h_+), source columns (h_-,h_+)")
print("coupling sign pattern: [-,+; +,-]")
