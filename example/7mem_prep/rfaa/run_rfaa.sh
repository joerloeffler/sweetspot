#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")"

# Run from your RoseTTAFold-All-Atom environment.
python -m rf2aa.run_inference -cd . --config-name rfaa
