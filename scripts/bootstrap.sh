#!/usr/bin/env bash
# Rebuild the full working environment for this project from nothing.
#
#   git clone https://github.com/jparkhill/2rdm_noneq_rdm.git
#   cd 2rdm_noneq_rdm && ./scripts/bootstrap.sh
#
# Idempotent: safe to re-run.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT="$PWD"

echo "==> Python dependencies"
python3 -m pip install -q --upgrade pip
python3 -m pip install -q numpy scipy sympy pyscf h5py matplotlib pytest
python3 - <<'PY'
import numpy, scipy, sympy, pyscf
print(f"    numpy {numpy.__version__}  scipy {scipy.__version__}  "
      f"sympy {sympy.__version__}  pyscf {pyscf.__version__}")
PY

echo "==> Upstream reference repos (read-only sources, never committed here)"
mkdir -p refs && cd refs
# The five that matter, in order of relevance to the rebuild.
for r in TightBindingTCL2 RealtimePySCF MukherjeeKutzelniggNO \
         Symbolic-Many-Fermion-Codes CorrelatedPolaron \
         SupportingInformation_PerovskiteDynamics \
         Transient_Absorption_Paper_SupportingInformation ; do
  if [ -d "$r/.git" ]; then
    echo "    $r (present)"
  else
    echo "    cloning $r"
    git clone -q "https://github.com/jparkhill/$r.git" "$r"
  fi
done
cd "$ROOT"

echo "==> Regenerating derived artifacts from MKNO.nb"
if [ -f refs/MukherjeeKutzelniggNO/MKNO.nb ]; then
  ( ulimit -v 4000000 || true
    python3 tools/nbdump.py refs/MukherjeeKutzelniggNO/MKNO.nb docs/MKNO_dump.md ) || \
      echo "    (nbdump skipped)"
else
  echo "    MKNO.nb missing; skipping"
fi

echo "==> Self-check"
python3 verify_1rdm.py || echo "    verify_1rdm.py FAILED - investigate before continuing"

echo
echo "Bootstrap complete. Read docs/CONTINUITY.md first, then docs/02_project_plan.md."
