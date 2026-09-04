# Certificate-to-manuscript ledger

This ledger maps the finite numerical premises used by the restart-four
certificate to exact assertions. The named programs are in `evidence/launch/`
relative to this directory. `verify.sh` executes them with Python assertions
enabled.

## Data and sign conventions

- Polynomial coefficient vectors are in ascending powers of `t`.
- The five jet components are `(value, d_a, d_aa, d_q, d_aq)`.
- External labels are ordered `(h_-, h_+)` throughout.
- The coupling matrix `Tmatrix` has target labels in its rows and source
  labels in its columns: `Tmatrix[ih][ig] = d_h[V_g]`.
- Its certified value-level sign pattern is `[-,+; +,-]` in that row/column
  order.
- Scalar tangent coordinates use `(u,B_C,B_0,B_1)`, while the intrinsic lift
  returns `(B_C,B_0,B_1,u)`. The tangent certificate uses this explicit
  permutation.
- `C`, all terminating decimal seed data, interval endpoints, and all derived
  certificate arithmetic are exact `fractions.Fraction` values. SymPy is used
  by the supplementary `all_s` stage.

For the coupling signs in particular, (C.3.9) gives
`Delta*(V_g)^H = p_g*Q + P*q_g - (V_g)_C - A_g*rq`. Hence (C.3.10) expands as
`d_h[V_g] = (-2/C)*(p_g(h)*Q(h) + P(h)*q_g(h) + (V_g)_C - A_g*rq(h))`, which is
the expression in (C.13.2) and in the transverse-Hopf program.

## Program inputs, outputs, and dependencies

| Program | Direct inputs | Principal outputs/assertions | Dependencies |
|---|---|---|---|
| `all_s_s4_transverse_hopf_exact.py` | Hard-coded rational `P,Qbar,C,q` boxes | Root isolation, critical branch, `Tmatrix`, positive amplitudes, Hopf sign change and transversality | Self-contained; Python standard library |
| `all_s_s4_coupling_identity_exact.py` | Rational seed and interval primitives from the definition prefix of the transverse program | Constructs `p_g,q_g`, obtains `(V_g)^H` by literal division in (C.3.9), and compares (C.3.10), (C.13.2), and `Tmatrix`, including all five jet components | Loads the definition prefix and constructs `(V_g)^H` by polynomial division; the audited matrix is used at the comparison step |
| `all_s_s4_local_majorants_exact.py` | Transverse-Hopf environment | Root/weight/factor margins, phase systems, Gram bounds, derivative bounds | Transverse program prerequisite |
| `all_s_s4_intrinsic_lift_exact.py` | Transverse-Hopf environment | Three-phase recursion, quadratic coefficients, `4x4` lift and inverse bound | Transverse program prerequisite; phase solves rebuilt within the lift |
| `all_s_s4_scalar_hopf_vector_exact.py` | `Tmatrix`, roots, and scalar data from the transverse program | Scalar Jacobian blocks, reconstructed `V_g`, resolvent and Hopf-vector bounds | Transverse program prerequisite; `V_g` reconstruction from `Tmatrix` |
| `all_s_s4_scalar_to_raw_tangent_exact.py` | Local-majorant, intrinsic-lift, and scalar-Hopf-vector environments | Coordinate-order permutation, raw tangent and visible-vector bounds | Transitively depends on the transverse program through its three inputs |
| `all_s_s4_explicit_counterexample_checker.py` | Scalar-to-raw environment and integer exponent ledger | Aggregate finite restart-four assertions and contraction exponent comparisons | Upstream certificate chain |
| `independent_all_s_symbolic_checks.py` | Fixed seed `20260825`; SymPy 1.14.0 | General symbolic identities where stated, characteristic determinant on seeded exact integer instances for `d=3,4,5,6`, Newton-sum checks for `d=1,...,8`, and exponent arithmetic | Self-contained supplementary regression stage |

The coupling stage reconstructs the polynomial identity and all five parameter
jets from the rational seed and root-isolation primitives. Stored matrix
entries enter at the final comparison step. The result independently confirms
the coupling signs in target-row/source-column order at the polynomial level.

| Manuscript item | Certified premise | Program and direct assertion |
|---|---|---|
| C.3.3--C.3.4 | One simple root in each stated collar, uniform continuation, and the selected core/external labels | all_s_s4_transverse_hopf_exact.py: root counts, opposite endpoint signs, and nonzero derivative intervals |
| C.3.9--C.3.10 and C.13.2 | Polynomial-division definition of `(V_g)^H`, the `+2(V_g)_C` radial contribution, and equality with the expanded coupling formula in target-row/source-column order | all_s_s4_coupling_identity_exact.py: independent polynomial construction, zero-containing division remainders, pairwise overlap of all five jet components, and the certified matrix sign pattern |
| C.3.11 | Unique critical branch in the displayed alpha strip | all_s_s4_transverse_hopf_exact.py: lower/upper G face signs and the strict interval for the alpha derivative of G |
| C.3.15 | Negative theta discriminant and 3.72 < omega < 3.724 | all_s_s4_transverse_hopf_exact.py: discriminant < -55.40; all_s_s4_scalar_hopf_vector_exact.py: exact comparisons of the omega-square interval with (93/25)^2 and (931/250)^2 |
| C.3.16 | Hopf crossing speed greater than 1.94 | all_s_s4_transverse_hopf_exact.py: 109.89 < hopf_slope < 110.82 and crossing_lower > 194/100 |
| C.3.17 | Core-weight and external-amplitude lower bounds | all_s_s4_transverse_hopf_exact.py: direct weight, c-minus, and c-plus assertions |
| C.3.17 and C.4.3 | Hopf parameter and alpha face margins | all_s_s4_local_majorants_exact.py: q-star margin > 2^-41 and alpha margin > 2^-24 |
| C.4.1 | First P-phase weight slope s_1 > 1/16 | all_s_s4_local_majorants_exact.py: direct first-weight slope assertion |
| C.4.2aa and C.5.4 | Every displayed lift-matrix entry, Neumann defect < 1/4, and lift inverse < 2^20 | all_s_s4_intrinsic_lift_exact.py: published_J_boxes, eta4, and inv4 assertions |
| C.3.14 and C.4.21--C.4.22 | Both phase 9 x 9 systems and their inverse bound | all_s_s4_local_majorants_exact.py: midpoint-residual Neumann bounds for the P and Q matrices; all_s_s4_intrinsic_lift_exact.py: the same matrices in the coefficient recursion |
| C.3.14e | Rigorous enclosures for the three-phase quadratic coefficients beta and q | all_s_s4_intrinsic_lift_exact.py: executed P,Q,P recursion and explicit return_beta/return_q interval extraction; the coefficient identities themselves are derived algebraically in the manuscript |
| C.4.10 | Root derivatives through order five | all_s_s4_local_majorants_exact.py: direct lower bound for Phi_t, direct derivative bound for the polynomial family, Bell-polynomial implicit-differentiation recurrence, and five asserted output bounds |
| C.4.11--C.4.12 | Weight-chart derivatives through order five | all_s_s4_explicit_counterexample_checker.py: exact integer operation ledger from coefficient differentiation, quotient differentiation, and the affine chart formulas |
| C.4.13 | Determinant, trace, and inverse bounds for the translated phase Gram matrix, including the signed perturbation | all_s_s4_local_majorants_exact.py: gram_det_lower, shifted_trace_upper, shifted_gram_inverse_bound, and the Neumann-corrected signed_gram_inverse_bound |
| C.4.15--C.4.17 | Fixed-sheet factor lower bound, normalisation lower bound, and persistence in the physical box | all_s_s4_local_majorants_exact.py: normalized_factor_lower, fixed_Z_lower, factor_variation, and Z_variation |
| C.4.4, C.4.8, and C.4.18 | Fifth derivative of the exact return map and the C3/C1 rescaled-map bounds | all_s_s4_explicit_counterexample_checker.py: direct integer exponents for reciprocal differentiation, one-block derivatives, two-block composition, chart derivatives, and rescaling |
| C.4.2e | Critical uncentred state bound | all_s_s4_scalar_to_raw_tangent_exact.py: v_star_upper < 2^500001 |
| C.5.1--C.5.2 | Scalar-to-weight tangent and visible Hopf-vector bounds | all_s_s4_scalar_to_raw_tangent_exact.py: raw_inf_forward, Euclidean norm conversion, and transferred vector assertions |
| C.5.13 | Scalar block bounds | all_s_s4_scalar_hopf_vector_exact.py: entrywise D, induced row-sum W, A, V, ell, qvec, and kappa assertions |
| C.5.14 | Resolvent and full scalar Hopf-vector bounds | all_s_s4_scalar_hopf_vector_exact.py: determinant lower bound, explicit dimension factors, and final full-vector assertion |
| C.6.4--C.6.14 | Projector, moving-equilibrium, augmented inverse, Sobolev, radius, and Lyapunov--Schmidt contraction exponents | all_s_s4_explicit_counterexample_checker.py: projector bound including the infinity-to-Euclidean norm conversion, derivative bounds, and radius/contraction exponent assertions |
| C.7.4--C.7.15 | Floquet roughness, nonlinear/remainder bounds, tube size, and all four Lyapunov--Perron inequalities | all_s_s4_explicit_counterexample_checker.py: Floquet, K, L, R, M, and N exponent assertions |

The `all_s` SymPy stage provides symbolic elevation identities and regression
coverage for signs and block placement. The characteristic-reduction
regression uses three seeded exact integer instances at each of `d=3,4,5,6`;
the finite Newton-sum range is `d=1,...,8`. The manuscript supplies the general
characteristic argument and the analytic contractions in C.6, C.7, C.15, and
C.16.
