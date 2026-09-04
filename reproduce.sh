#!/usr/bin/env bash
set -euo pipefail
script_path=${BASH_SOURCE[0]}
if [[ "$script_path" == */* ]]; then
  script_directory=${script_path%/*}
else
  script_directory=.
fi
cd -- "$script_directory"
export PYTHONDONTWRITEBYTECODE=1

bootstrap_python="${FORSYTHE_BOOTSTRAP_PYTHON:-}"
if [[ -z "$bootstrap_python" ]]; then
  if command -v python3 >/dev/null 2>&1; then
    bootstrap_python=$(command -v python3)
  elif command -v python >/dev/null 2>&1; then
    bootstrap_python=$(command -v python)
  else
    echo 'missing replay dependency: Python 3' >&2
    exit 2
  fi
fi

temporary_root="${TMPDIR:-/tmp}"
if [[ ! -d "$temporary_root" ]]; then
  echo "temporary directory does not exist: $temporary_root" >&2
  exit 2
fi
work=$(mktemp -d "$temporary_root/forsythe-certificate-replay.XXXXXX")
trap 'rm -rf "$work"' EXIT

"$bootstrap_python" -m venv "$work/venv"
if [[ -x "$work/venv/bin/python" ]]; then
  venv_python="$work/venv/bin/python"
elif [[ -x "$work/venv/Scripts/python.exe" ]]; then
  venv_python="$work/venv/Scripts/python.exe"
else
  echo 'virtual environment did not provide a Python executable' >&2
  exit 2
fi

"$venv_python" -m pip install --disable-pip-version-check \
  -r requirements-verification.txt
FORSYTHE_PYTHON="$venv_python" bash ./verify.sh

printf 'PASS: isolated certificate replay completed\n'
