# Computer-validated certificate

This directory is the active certificate bundle for the groups of finite
algebraic and interval premises listed in `CERTIFICATE_LEDGER.md`. It does not
attempt to certify the analytic periodic-orbit, shadowing, or infinite-time
arguments in Part C of the manuscript.

## Replay

From the repository root, install the pinned symbolic dependency and run the
driver through a POSIX shell:

```bash
python3 -m pip install -r computer_validated_certificate/requirements-verification.txt
bash computer_validated_certificate/verify.sh
```

Set `FORSYTHE_PYTHON` if the desired Python executable is not named `python3`
or `python`. For example:

```bash
FORSYTHE_PYTHON=/path/to/python bash computer_validated_certificate/verify.sh
```

The available stages are:

| Stage | Purpose | Extra dependency |
|---|---|---|
| `coupling` | Independently reconstruct (C.3.9), then compare (C.3.10), (C.13.2), and the certified coupling matrix | Python standard library |
| `scalar` | Restart-four scalar Hopf-vector certificate | Python standard library |
| `aggregate` | Restart-four finite inputs and exponent ledger | Python standard library |
| `all_s` | Supplementary degree-elevation identities and seeded exact regression tests | SymPy 1.14.0 and mpmath 1.3.0 |

Pass a stage name as the first argument, for example:

```bash
bash computer_validated_certificate/verify.sh coupling
```

With no argument, the driver runs all four stages. Python assertions must be
enabled; the driver rejects `-O`, `-OO`, and `PYTHONOPTIMIZE`.

The `coupling` stage is an independent polynomial-division reconstruction of
(C.3.9)--(C.3.10) after loading the shared rational seed and root-isolation
primitives. It is not a second differentiation of the full eight-node
weight-simplex return map.

For an isolated certificate replay, including creation of a temporary virtual
environment and installation of the pinned packages, run:

```bash
bash computer_validated_certificate/reproduce.sh
```

That command verifies this certificate only. A manuscript build driver and
the manuscript build itself are intentionally outside this ancillary bundle.

## Integrity and archive creation

`verification-lock.json` records the interpreter/package expectations and
stage dependencies in machine-readable form. `SHA256SUMS` covers exactly the
seven restart-four programs whose exact calculations supply finite premises
to the proof. Documentation, replay and archive tooling, environment metadata,
transcripts, and the supplementary `all_s` regression program are deliberately
outside the checksum manifest.

After changing the extracted directory, regenerate both the manifest and the
repository-root archive with:

```bash
python3 computer_validated_certificate/build_archive.py
```

The builder uses a fixed member timestamp, deterministic directory and file
lists, and POSIX modes. In particular, `verify.sh` and `reproduce.sh` are
stored as executable in the ZIP.
The README nevertheless uses `bash ...` invocations so replay does not depend
on an extraction tool preserving mode bits.

See `CERTIFICATE_LEDGER.md` for manuscript mappings, coordinate/sign
conventions, program inputs and outputs, and the import-dependency graph.
