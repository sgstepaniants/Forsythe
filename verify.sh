#!/usr/bin/env bash
set -euo pipefail
script_path=${BASH_SOURCE[0]}
if [[ "$script_path" == */* ]]; then
  script_directory=${script_path%/*}
else
  script_directory=.
fi
cd -- "$script_directory"
export PYTHONHASHSEED=0
export PYTHONUNBUFFERED=1
export PYTHONDONTWRITEBYTECODE=1

if [[ -n "${FORSYTHE_PYTHON:-}" ]]; then
  :
elif command -v python3 >/dev/null 2>&1; then
  FORSYTHE_PYTHON=$(command -v python3)
elif command -v python >/dev/null 2>&1; then
  FORSYTHE_PYTHON=$(command -v python)
else
  echo 'missing verification dependency: Python 3' >&2
  exit 2
fi

if ! "$FORSYTHE_PYTHON" -c \
  'import sys; sys.exit(0 if __debug__ else 1)' >/dev/null 2>&1; then
  echo 'Verification requires Python assertions; do not use -O/-OO.' >&2
  echo 'Unset PYTHONOPTIMIZE and retry.' >&2
  exit 4
fi

require_sympy() {
  if ! "$FORSYTHE_PYTHON" -c \
    'import mpmath, sympy, sys; sys.exit(0 if sympy.__version__ == "1.14.0" and mpmath.__version__ == "1.3.0" else 1)' \
      >/dev/null 2>&1; then
    echo 'SymPy 1.14.0 and mpmath 1.3.0 are required for this stage.' >&2
    echo 'Install the pinned dependency with:' >&2
    echo '  python3 -m pip install -r requirements-verification.txt' >&2
    echo 'or set FORSYTHE_PYTHON to a Python executable that provides it.' >&2
    exit 3
  fi
}

run_stage() {
  case "$1" in
    coupling)
      echo '[coupling] independent C.3.9/C.3.10/C.13.2 identity check'
      "$FORSYTHE_PYTHON" evidence/launch/all_s_s4_coupling_identity_exact.py ;;
    scalar)
      echo '[scalar] focused restart-four Hopf-vector certificate'
      "$FORSYTHE_PYTHON" evidence/launch/all_s_s4_scalar_hopf_vector_exact.py ;;
    aggregate)
      echo '[aggregate] restart-four finite inputs and exponent ledger'
      "$FORSYTHE_PYTHON" evidence/launch/all_s_s4_explicit_counterexample_checker.py ;;
    all_s)
      echo '[all_s] degree-elevation identities and seeded exact regression tests'
      require_sympy
      "$FORSYTHE_PYTHON" evidence/launch/independent_all_s_symbolic_checks.py ;;
    *) echo "unknown verification stage: $1" >&2; exit 2 ;;
  esac
}

stage="${1:-all}"
if [[ "$stage" == all ]]; then
  for s in coupling scalar aggregate all_s; do
    run_stage "$s"
    echo
  done
  echo 'PASS: computer-validated certificate replay completed'
else
  run_stage "$stage"
fi
