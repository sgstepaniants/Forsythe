# Verification environment

The original reference replay for the numerical certificate was run on
2026-09-02 with:

- macOS 14.6.1 on arm64;
- CPython 3.10.13 in an isolated virtual environment;
- SymPy 1.14.0 and mpmath 1.3.0, installed from
  `computer_validated_certificate/requirements-verification.txt`.

The verification environment was created and exercised with:

```bash
python3 -m venv /tmp/forsythe-referee-verify-venv
/tmp/forsythe-referee-verify-venv/bin/python -m pip install \
  -r computer_validated_certificate/requirements-verification.txt
FORSYTHE_PYTHON=/tmp/forsythe-referee-verify-venv/bin/python \
  bash ./computer_validated_certificate/verify.sh
```

The same isolated sequence is automated by
`bash ./computer_validated_certificate/reproduce.sh`, which creates a fresh
temporary virtual environment, installs the exact versions in
`computer_validated_certificate/requirements-verification.txt`, and runs
every active verification stage. The ancillary package contains no manuscript
build driver and makes no claim to rebuild the manuscript.

`verify.sh` first checks that Python assertions are enabled and checks the
exact SymPy version before any stage that imports SymPy. The proof-bearing
restart-four stages use Python integers, `fractions.Fraction`, rational
closed intervals, polynomial division, and Sturm sequences. Every decimal
literal used as certificate data is parsed digit-by-digit into a base-ten
`Fraction`; it is never routed through a binary floating-point value.

The output from an end-to-end certificate-driver replay after the package
repair is included as `verification-transcript.txt`. It was recorded
on 2026-09-03 with Windows, CPython 3.11.9, SymPy 1.14.0, and mpmath 1.3.0;
all four active stages passed. The `coupling` stage was added specifically to
reconstruct (C.3.9) independently of the stored coupling matrix and compare it
with (C.3.10) and (C.13.2).

`verification-lock.json` records these package and stage requirements in
machine-readable form. `requirements-verification.txt` remains the direct pip
input for the isolated replay.

No permanent external archive identifier has been assigned to this working
package; one should be added by the publication repository when the ancillary
archive is deposited.
