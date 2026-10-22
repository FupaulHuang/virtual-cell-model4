#!/usr/bin/env bash
# Use an explicitly selected VCC executable, or the installed CLI on PATH.
set -euo pipefail
exec "${VCC_BIN:-vcc}" "$@"
