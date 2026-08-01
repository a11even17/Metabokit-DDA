#!/usr/bin/env bash

set -Eeuo pipefail

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
readonly TAURI_CONFIG="crates/metabokit-app/tauri.conf.json"

universal=false
run_tests=true

usage() {
  cat <<'EOF'
Build MetaboKit DDA for macOS.

Usage:
  ./build-macos.sh [--universal] [--skip-tests]

Options:
  --universal   Build one app for both Apple Silicon and Intel Macs.
  --skip-tests  Skip the workspace test suite before packaging.
  -h, --help    Show this help text.

Build products are written under target/release/bundle/ (or the universal
target's equivalent when --universal is used).
EOF
}

fail() {
  printf 'error: %s\n' "$*" >&2
  exit 1
}

while (($# > 0)); do
  case "$1" in
    --universal)
      universal=true
      ;;
    --skip-tests)
      run_tests=false
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      usage >&2
      fail "unknown option: $1"
      ;;
  esac
  shift
done

[[ "$(uname -s)" == "Darwin" ]] || fail "this script must be run on macOS"
command -v cargo >/dev/null 2>&1 || fail "Rust is missing. Install it from https://rustup.rs and run this script again."
command -v rustc >/dev/null 2>&1 || fail "rustc is missing. Reinstall or repair the Rust toolchain."
command -v xcode-select >/dev/null 2>&1 || fail "Xcode Command Line Tools are missing. Run: xcode-select --install"
xcode-select -p >/dev/null 2>&1 || fail "Xcode Command Line Tools are missing. Run: xcode-select --install"

cd "$SCRIPT_DIR"

printf '\nMetaboKit DDA macOS build\n'
printf '  Host: %s (%s)\n' "$(sw_vers -productVersion)" "$(uname -m)"
printf '  Rust: %s\n\n' "$(rustc --version)"

if ! cargo tauri --version >/dev/null 2>&1; then
  printf 'Tauri CLI was not found; installing the current Tauri 2 CLI...\n'
  cargo install tauri-cli --version '^2' --locked
fi

if [[ "$run_tests" == true ]]; then
  printf '\nRunning tests...\n'
  cargo test --workspace
fi

build_args=(tauri build --config "$TAURI_CONFIG")

if [[ "$universal" == true ]]; then
  command -v rustup >/dev/null 2>&1 || fail "a universal build requires rustup: https://rustup.rs"
  printf '\nInstalling the Apple Silicon and Intel Rust targets...\n'
  rustup target add aarch64-apple-darwin x86_64-apple-darwin
  build_args+=(--target universal-apple-darwin)
  bundle_dir="$SCRIPT_DIR/target/universal-apple-darwin/release/bundle"
else
  bundle_dir="$SCRIPT_DIR/target/release/bundle"
fi

printf '\nBuilding the optimized macOS app and disk image...\n'
cargo "${build_args[@]}"

printf '\nBuild complete. Products:\n'

found=false
while IFS= read -r -d '' product; do
  printf '  %s\n' "$product"
  found=true
done < <(find "$bundle_dir" -maxdepth 3 \( -name '*.app' -o -name '*.dmg' \) -print0 2>/dev/null)

if [[ "$found" != true ]]; then
  printf '  %s\n' "$bundle_dir"
fi
