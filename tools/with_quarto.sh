#!/usr/bin/env bash
set -euo pipefail

if [[ -z "${SVCFIT_QUARTO_PREFIX:-}" ]] && command -v quarto >/dev/null 2>&1; then
  exec "$@"
fi

quarto_prefix="${SVCFIT_QUARTO_PREFIX:-$HOME/.local/share/svcfit-quarto}"
quarto_bin="$quarto_prefix/bin"
quarto_share="$quarto_prefix/share/quarto"

if [[ ! -x "$quarto_bin/quarto" ]]; then
  echo "Quarto is not installed at $quarto_prefix" >&2
  echo "Set SVCFIT_QUARTO_PREFIX to the user-local Quarto environment." >&2
  exit 1
fi

export PATH="$quarto_bin:$PATH"
export QUARTO_SHARE_PATH="$quarto_share"
export QUARTO_DENO="$quarto_bin/deno"
export QUARTO_DENO_DOM="$quarto_prefix/lib/deno_dom.dylib"
export QUARTO_PANDOC="$quarto_bin/pandoc"
export QUARTO_ESBUILD="$quarto_bin/esbuild"
export QUARTO_TYPST="$quarto_bin/typst"
export QUARTO_DART_SASS="$quarto_bin/sass"

exec "$@"
