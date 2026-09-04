# Forsythe conjecture: exact certificate for restart length s = 4

This repository contains the machine-checkable certificate for an explicit
counterexample to the Forsythe conjecture at conjugate-gradient restart length
s = 4.

> [!IMPORTANT]
> This is the computational part of a computer-assisted proof, not a
> standalone proof. The programs rigorously verify the proof's finite
> algebraic and interval premises. The accompanying manuscript supplies the
> analytic periodic-orbit, contraction, shadowing, and infinite-time
> arguments that turn those premises into a counterexample.

## What is being checked?

Forsythe's conjecture predicts that the normalized residual directions of
fixed-length restarted conjugate gradients approach an asymptotic two-cycle:
the even and odd subsequences should each converge. The construction supported
by this certificate violates that conclusion for s = 4.

The exact checks in this repository cover:

- the physical transverse Hopf point and its unique critical branch;
- the coupling identity and its sign convention;
- local majorants, phase solves, the intrinsic lift, and tangent-coordinate
  bounds;
- scalar Hopf-vector and resolvent bounds; and
- every finite exponent comparison used by the Lyapunov--Schmidt, Floquet,
  and Lyapunov--Perron estimates.

The proof-bearing computations use integers, `fractions.Fraction`, closed
rational intervals, polynomial division, and Sturm sequences. Decimal
certificate data are converted directly to rational numbers; they are not
passed through binary floating-point arithmetic.

This repository does **not** by itself verify:

- the manuscript's analytic periodic-orbit, shadowing, or infinite-time
  arguments;
- a counterexample for every restart length; or
- an independent differentiation of the full eight-node weight-simplex
  return map.

The stage named `all_s` checks some general degree-elevation identities and
seeded exact regression examples. It is supplementary and is not a proof for
arbitrary restart length or arbitrary degree.

For the precise correspondence between code and manuscript equations, see
[the certificate ledger](CERTIFICATE_LEDGER.md).

## Reproduce the certificate

The simplest end-to-end replay creates a temporary virtual environment,
installs the pinned dependencies, and runs every verification stage:

```bash
bash reproduce.sh
```

Requirements are Python 3, Bash, and network access for the package
installation. A successful run ends with:

```text
PASS: computer-validated certificate replay completed
PASS: isolated certificate replay completed
```

To use an existing Python environment instead:

```bash
python3 -m pip install -r requirements-verification.txt
bash verify.sh
```

Set `FORSYTHE_PYTHON` when the desired interpreter has another name or path:

```bash
FORSYTHE_PYTHON=/path/to/python bash verify.sh
```

Python assertions must be enabled. `verify.sh` rejects `-O`, `-OO`, and any
`PYTHONOPTIMIZE` setting that disables assertions.

## Verification stages

Run one stage by passing its name to `verify.sh`; with no argument, all four
stages run in the order shown below.

| Stage | What it checks | Dependency |
|---|---|---|
| `coupling` | Independently reconstructs (C.3.9), then compares (C.3.10), (C.13.2), and the certified coupling matrix | Python standard library |
| `scalar` | Focused s = 4 scalar Hopf-vector certificate | Python standard library |
| `aggregate` | Finite inputs and exponent ledger for the explicit s = 4 counterexample | Python standard library |
| `all_s` | Supplementary degree-elevation identities and seeded exact regression checks | SymPy 1.14.0 and mpmath 1.3.0 |

For example:

```bash
bash verify.sh aggregate
```

The `coupling` stage independently reconstructs the polynomial identity behind
(C.3.9)--(C.3.10), after loading only the shared rational seed and
root-isolation primitives. It does not re-differentiate the full return map;
the ledger explains this independence boundary in detail.

## Repository guide

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

## Integrity and deterministic archive creation

Check the proof-program manifest with either:

```bash
# Linux
sha256sum -c SHA256SUMS

# macOS
shasum -a 256 -c SHA256SUMS
```

After changing the bundle, regenerate the manifest and the deterministic ZIP
archive with:

```bash
python3 build_archive.py
```

The archive builder fixes member timestamps, ordering, and POSIX modes. The
checksum manifest deliberately covers only the seven proof-bearing s = 4
programs; documentation, replay tooling, environment metadata, transcripts,
and the supplementary `all_s` program remain outside it.
