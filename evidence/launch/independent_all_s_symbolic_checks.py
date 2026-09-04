#!/usr/bin/env python3
"""Exact symbolic identities and seeded regressions for all-s>=4 elevation.

The manuscript supplies the implicit-function, Lyapunov--Schmidt, and
Lyapunov--Perron arguments. This module supplies their finite algebraic layer.
"""
from __future__ import annotations

if not __debug__:
    raise RuntimeError(
        "Verification assertions are disabled; Python optimization flags "
        "are incompatible with this certificate."
    )

import random
import sympy as sp


def check_orientation_identity() -> None:
    k, c1, c2, s = sp.symbols("k c1 c2 s")
    A1, A2, g1, g2 = sp.symbols("A1 A2 g1 g2")
    t11, t12, t21, t22 = sp.symbols("t11 t12 t21 t22")
    l1, l2 = sp.symbols("l1 l2")

    Jc = sp.Matrix(
        [
            [s, c1 * A1, c2 * A2],
            [g1, t11 * c1, t12 * c2],
            [g2, t21 * c1, t22 * c2],
        ]
    )
    R = sp.Matrix([[c1 * A1, c2 * A2]])
    B = sp.Matrix([[t11 * c1, t12 * c2], [t21 * c1, t22 * c2]])
    ell = sp.Matrix([l1, l2])
    g = sp.Matrix([g1, g2])
    M = sp.zeros(3)
    M[0, 0] = k * s + (R * ell)[0]
    M[0, 1:3] = R
    M[1:3, 0] = k * g + B * ell
    M[1:3, 1:3] = B

    # Jc is D_(a,c)Psi with its two amplitude columns multiplied by c1,c2.
    assert sp.expand(M.det() - k * Jc.det()) == 0


def check_characteristic_reduction() -> None:
    """Identity det(lambda I-J)=(lambda-1)^(d-2)det(theta I_3-M).

    Three seeded exact integer instances cover each listed dimension. The
    manuscript supplies the general polynomial elimination argument.
    """
    lam = sp.symbols("lambda")
    rng = random.Random(20260825)
    for d in (3, 4, 5, 6):
        for _ in range(3):
            k = sp.Integer(rng.randint(-4, 4) or 1)
            p = sp.Integer(rng.randint(-4, 4))
            R = sp.Matrix([[rng.randint(-3, 3), rng.randint(-3, 3)]])
            q = sp.Matrix([rng.randint(-3, 3) for _ in range(d)])
            W = sp.Matrix(d, 2, [rng.randint(-3, 3) for _ in range(2 * d)])
            ell = sp.Matrix([rng.randint(-3, 3), rng.randint(-3, 3)])
            D = sp.Matrix(2, d, [rng.randint(-3, 3) for _ in range(2 * d)])

            size = d + 4
            J = sp.zeros(size)
            ia, iu = 0, 1
            iB = 2
            il = 2 + d
            J[ia, iu] = k
            J[iu, ia] = p
            J[iu, iu] = 1
            J[iu, il : il + 2] = R
            J[iB : iB + d, ia] = q
            J[iB : iB + d, iB : iB + d] = sp.eye(d)
            J[iB : iB + d, il : il + 2] = W
            J[il : il + 2, iu] = ell
            J[il : il + 2, iB : iB + d] = D

            DW = D * W
            M = sp.zeros(3)
            M[0, 0] = k * p + (R * ell)[0]
            M[0, 1:3] = R
            M[1:3, 0] = DW * ell + k * D * q
            M[1:3, 1:3] = DW
            theta = lam * (lam - 1)

            lhs = sp.Poly((lam * sp.eye(size) - J).det(), lam)
            rhs = sp.Poly((lam - 1) ** (d - 2) * (theta * sp.eye(3) - M).det(), lam)
            assert lhs == rhs


def check_similarity_scaling() -> None:
    """Block scaling used in the far-root Hopf persistence."""
    R = sp.symbols("R", nonzero=True)
    a = sp.symbols("a")
    x1, x2, y1, y2 = sp.symbols("x1 x2 y1 y2")
    b11, b12, b21, b22 = sp.symbols("b11 b12 b21 b22")
    old = sp.Matrix([[a, x1, x2], [y1, b11, b12], [y2, b21, b22]])
    raw = sp.Matrix(
        [[a, -R * x1, -R * x2], [-y1 / R, b11, b12], [-y2 / R, b21, b22]]
    )
    D = sp.diag(-R, 1, 1)
    assert sp.simplify(D.inv() * raw * D - old) == sp.zeros(3)


def check_elevation_factor_identity() -> None:
    t, R, d, C = sp.symbols("t R d C", nonzero=True)
    P, Q = sp.symbols("P Q")
    lhs = ((t - R) * P) * ((t - R - d) * Q) / (R * (R + d)) - C
    rhs = P * Q * (1 - t / R) * (1 - t / (R + d)) - C
    assert sp.factor(lhs - rhs) == 0


def check_three_phase_differential() -> None:
    """Identity (6.0f) for j=0,1 using SR0=Delta*u-rho."""
    C, rho, Delta_u, S_R0 = sp.symbols("C rho Delta_u S_R0")
    # The only algebraic relation used is S_R0=Delta_u-rho.
    for j in (0, 1):
        sigma_j = -sp.Rational(j, 1) * rho / C
        sigma_next = -sp.Rational(j + 1, 1) * rho / C
        middle = C * sigma_j + S_R0
        assert sp.expand(middle.subs(S_R0, Delta_u - rho) - (Delta_u + C * sigma_next)) == 0


def check_chebyshev_newton_degree(max_d: int = 8) -> None:
    """Newton-sum degree bound for the listed finite range of d."""
    x, y = sp.symbols("x y")
    for d in range(1, max_d + 1):
        Phi = sp.expand(16 * (1 - (-1) ** d * sp.chebyshevt(d, x - 1)))
        poly = sp.Poly(Phi - y, x)
        # Make monic and use Newton sums through 4d-1.
        monic = sp.Poly(poly.as_expr() / poly.LC(), x)
        roots_power_sums = []
        coeffs = monic.all_coeffs()  # [1,a1,...,ad]
        for m in range(1, 4 * d):
            total = 0
            # p_m + a1 p_(m-1)+...+a_(m-1)p_1 + m a_m =0 (m<=d)
            upper = min(m - 1, d)
            for j in range(1, upper + 1):
                total += coeffs[j] * roots_power_sums[m - j - 1]
            if m <= d:
                total += m * coeffs[m]
            else:
                for j in range(upper + 1, d + 1):
                    # For m>d, include a_j p_(m-j), with p_0=d only when m=j,
                    # but here m>d so m-j>=1 for j<=d.
                    total += coeffs[j] * roots_power_sums[m - j - 1]
            pm = sp.cancel(-total)
            roots_power_sums.append(pm)
            assert sp.Poly(pm, y).degree() <= m // d



def check_binary_exponent_ledger() -> None:
    """Integer exponent comparisons in the aggregate s=4 certificate."""
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

    projector_inf_exp = 6 * (D + 1)
    # For a 7 x 7 matrix, ||A||_2 <= sqrt(7) ||A||_inf < 2^2 ||A||_inf.
    projector_two_exp = projector_inf_exp + 2
    assert projector_two_exp < P0
    assert P0 + 3 + D < Z1
    assert P0 + 3 + D + 2 + 2 * Z1 < Z2
    assert 2 * P0 + D + Z1 + 2 * V + 1 < MH
    assert 20 + D + max(Z1, 2 * Z1, Z2) + 3 * V < BH

    R0_signed_exp = 1 + MH + BH - r0
    assert R0_signed_exp == -359_999_999
    ls_contraction_exp = MH + LH + (MH + BH + 2) - r0
    assert ls_contraction_exp == -219_999_998
    assert ls_contraction_exp < -2
    assert R0_signed_exp + Z1 < -338_000_000
    assert D - 330_000_000 < -eps

    assert 4 + (P0 + 3) - eps < -281_000_000
    assert 2 + 2 * (P0 + 9) - eps < -263_000_000
    assert KP - delta_x_inv < -100
    assert 20 + 3 * KP + D < L
    assert 20 + 5 * KP + D < R

    assert M > delta_x_inv
    assert 4 + K + R + 1 < M
    assert 9 + K + L < M
    assert 8 + K + R < N
    assert N == 2 * M
    assert N > 2_000_000
    assert N > 2 * 18_200_000 + 1

def main() -> None:
    check_elevation_factor_identity()
    check_orientation_identity()
    check_characteristic_reduction()
    check_similarity_scaling()
    check_three_phase_differential()
    check_chebyshev_newton_degree()
    check_binary_exponent_ledger()
    print("PASS: exact all-s elevation identities")
    print("- normalized far-root elevation identity")
    print("- critical Jacobian orientation determinant")
    print("- characteristic reduction on seeded exact rational test instances")
    print("- Hopf-matrix similarity scaling")
    print("- three-phase return differential")
    print("- Chebyshev Newton-sum degree bound for d=1,...,8")
    print("- complete binary-exponent ledger from the aggregate s=4 certificate")


if __name__ == "__main__":
    main()
