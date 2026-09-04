# Forsythe conjecture: exact certificate for restart length s = 4

This repository contains the exact computational certificate for an explicit
counterexample to the Forsythe conjecture at conjugate-gradient restart length
s = 4.

The manuscript and this certificate together form a computer-assisted proof.
The manuscript develops the periodic-orbit, contraction, shadowing, and
infinite-time arguments. The certificate supplies the finite algebraic and
interval premises.

## Certified statements

Forsythe's conjecture predicts that the normalized residual directions of
fixed-length restarted conjugate gradients approach an asymptotic two-cycle:
the even and odd subsequences should each converge. The construction supported
by this certificate violates that conclusion for s = 4.

The certificate includes:

- the physical transverse Hopf point and its unique critical branch;
- the coupling identity and its sign convention;
- local majorants, phase solves, the intrinsic lift, and tangent-coordinate
  bounds;
- scalar Hopf-vector and resolvent bounds; and
- every finite exponent comparison used by the Lyapunov--Schmidt, Floquet,
  and Lyapunov--Perron estimates.

The proof-bearing computations use integers, `fractions.Fraction`, closed
rational intervals, polynomial division, and Sturm sequences. Decimal
certificate data enter the computation directly as exact base-ten rational
numbers.

The `all_s` stage complements the s = 4 certificate with general
degree-elevation identities and seeded exact regression examples. The
manuscript contains the corresponding arbitrary-degree arguments.

[The certificate ledger](CERTIFICATE_LEDGER.md) gives the precise
correspondence between code and manuscript equations.

## Reproduction

### Isolated environment

`reproduce.sh` creates a temporary virtual environment, installs the pinned
dependencies, and executes every certificate stage:

```bash
bash reproduce.sh
```

Requirements: Python 3, Bash, and network access for the package installation.
Successful output ends with:

```text
PASS: computer-validated certificate replay completed
PASS: isolated certificate replay completed
```

### Existing environment

```bash
python3 -m pip install -r requirements-verification.txt
bash verify.sh
```

`FORSYTHE_PYTHON` selects a specific interpreter:

```bash
FORSYTHE_PYTHON=/path/to/python bash verify.sh
```

The replay runs with Python assertions enabled. Optimized Python modes (`-O`,
`-OO`, or an equivalent `PYTHONOPTIMIZE` setting) cause `verify.sh` to exit.

## Certificate stages

`verify.sh <stage>` selects one stage. With no stage argument, all four execute
in the order shown below.

| Stage | Contents | Dependency |
|---|---|---|
| `coupling` | Independent reconstruction of (C.3.9) and comparison with (C.3.10), (C.13.2), and the certified coupling matrix | Python standard library |
| `scalar` | Focused s = 4 scalar Hopf-vector certificate | Python standard library |
| `aggregate` | Finite inputs and exponent ledger for the explicit s = 4 counterexample | Python standard library |
| `all_s` | Supplementary degree-elevation identities and seeded exact regressions | SymPy 1.14.0 and mpmath 1.3.0 |

Single-stage form:

```bash
bash verify.sh aggregate
```

The `coupling` stage uses the shared rational seed and root-isolation
primitives. Its outputs include the reconstructed polynomial identity behind
(C.3.9)--(C.3.10) and comparisons of all five jet components. The ledger gives
the full construction and dependency boundary.

## Files

- [`CERTIFICATE_LEDGER.md`](CERTIFICATE_LEDGER.md) maps every certified premise
  to a manuscript equation and a direct assertion.
- [`evidence/launch/`](evidence/launch/) contains the exact certificate
  programs.
- [`verification-lock.json`](verification-lock.json) records interpreter,
  package, and stage requirements in machine-readable form.
- [`VERIFICATION_ENVIRONMENT.md`](VERIFICATION_ENVIRONMENT.md) describes the
  reference environments and exact-arithmetic design.
- [`verification-transcript.txt`](verification-transcript.txt) records a
  successful end-to-end replay.
- [`SHA256SUMS`](SHA256SUMS) authenticates the seven proof-bearing s = 4
  programs.

## Integrity and deterministic archive

Manifest validation:

```bash
# Linux
sha256sum -c SHA256SUMS

# macOS
shasum -a 256 -c SHA256SUMS
```

Manifest and deterministic ZIP regeneration:

```bash
python3 build_archive.py
```

The archive builder fixes member timestamps, ordering, and POSIX modes. The
checksum manifest covers the seven proof-bearing s = 4 programs.
