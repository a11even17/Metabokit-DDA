#!/usr/bin/env bash

set -Eeuo pipefail

readonly SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
readonly TAURI_CONFIG="crates/metabokit-app/tauri.conf.json"

universal=false
run_tests=true
release=false
sign_only=false
notary_profile="${METABOKIT_NOTARY_PROFILE:-}"

usage() {
  cat <<'EOF'
Build MetaboKit DDA for macOS.

Usage:
  ./build-macos.sh [--universal] [--skip-tests]
  ./build-macos.sh --sign-only [--universal] [--skip-tests]
  ./build-macos.sh --release --notary-profile PROFILE [--universal] [--skip-tests]

Modes:
  (default)     Build an unsigned DMG for testing on this Mac.
  --sign-only   Sign with the installed Developer ID Application certificate,
                but do not notarize. This is not ready for public distribution.
  --release     Sign, submit the DMG to Apple's notary service, staple its
                ticket, and verify it with Gatekeeper.

Options:
  --notary-profile PROFILE
                Keychain profile created by `xcrun notarytool store-credentials`.
                METABOKIT_NOTARY_PROFILE may be used instead.
  --universal   Build one app for both Apple Silicon and Intel Macs.
  --skip-tests  Skip the workspace test suite before packaging.
  -h, --help    Show this help text.

Products are written under target/release/bundle/ (or the universal target's
equivalent when --universal is used).
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
    --release)
      release=true
      ;;
    --sign-only)
      sign_only=true
      ;;
    --notary-profile)
      (($# >= 2)) || fail "--notary-profile requires a profile name"
      notary_profile="$2"
      shift
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

[[ "$release" != true || "$sign_only" != true ]] || fail "choose either --release or --sign-only, not both"
[[ "$(uname -s)" == "Darwin" ]] || fail "this script must be run on macOS"

# Prefer the full Xcode installation without changing the user's global
# xcode-select setting. The standalone Command Line Tools do not contain every
# distribution and notarization component.
if ! xcode-select -p 2>/dev/null | grep -q '/Xcode.app/'; then
  if [[ -d /Applications/Xcode.app/Contents/Developer ]]; then
    export DEVELOPER_DIR=/Applications/Xcode.app/Contents/Developer
  fi
fi

command -v cargo >/dev/null 2>&1 || fail "Rust is missing. Install it from https://rustup.rs and run this script again."
command -v rustc >/dev/null 2>&1 || fail "rustc is missing. Reinstall or repair the Rust toolchain."
command -v xcrun >/dev/null 2>&1 || fail "Xcode Command Line Tools are missing."
xcodebuild -version >/dev/null 2>&1 || fail "full Xcode is unavailable or its license/setup has not been completed"

# `C.UTF-8` is common in development shells but is not a locale shipped by
# macOS. Tauri's DMG helper invokes Perl, which aborts under that locale.
if [[ "${LANG:-}" == "C.UTF-8" || "${LC_ALL:-}" == "C.UTF-8" ]]; then
  export LANG=en_US.UTF-8
  export LC_ALL=en_US.UTF-8
fi

cd "$SCRIPT_DIR"

printf '\nMetaboKit DDA macOS build\n'
printf '  Host: %s (%s)\n' "$(sw_vers -productVersion)" "$(uname -m)"
printf '  Xcode: %s\n' "$(xcodebuild -version | paste -sd ' ' -)"
printf '  Rust: %s\n' "$(rustc --version)"

if ! cargo tauri --version >/dev/null 2>&1; then
  fail "Tauri CLI is missing. Install it with: cargo install tauri-cli --version '^2' --locked"
fi

if [[ "$release" == true || "$sign_only" == true ]]; then
  if [[ -z "${APPLE_SIGNING_IDENTITY:-}" ]]; then
    signing_identities=()
    while IFS= read -r line; do
      if [[ "$line" =~ \"(Developer\ ID\ Application:[^\"]+)\" ]]; then
        signing_identities+=("${BASH_REMATCH[1]}")
      fi
    done < <(security find-identity -v -p codesigning)

    ((${#signing_identities[@]} > 0)) || fail "no valid Developer ID Application identity was found in the keychain"
    ((${#signing_identities[@]} == 1)) || fail "multiple Developer ID Application identities were found; set APPLE_SIGNING_IDENTITY explicitly"
    export APPLE_SIGNING_IDENTITY="${signing_identities[0]}"
  fi
  printf '  Signing: %s\n' "$APPLE_SIGNING_IDENTITY"
fi

if [[ "$release" == true ]]; then
  [[ -n "$notary_profile" ]] || fail "--release requires --notary-profile PROFILE (see MACOS_DISTRIBUTION.md)"
  printf '  Notary profile: %s\n' "$notary_profile"
  if ! xcrun notarytool history --keychain-profile "$notary_profile" >/dev/null 2>&1; then
    fail "notary profile '$notary_profile' is missing or invalid; create it as described in MACOS_DISTRIBUTION.md"
  fi
fi

if [[ "$run_tests" == true ]]; then
  printf '\nRunning tests...\n'
  cargo test --workspace
fi

build_args=(tauri build --config "$TAURI_CONFIG" --bundles dmg)

if [[ "$release" != true && "$sign_only" != true ]]; then
  build_args+=(--no-sign)
fi

if [[ "$universal" == true ]]; then
  command -v rustup >/dev/null 2>&1 || fail "a universal build requires rustup: https://rustup.rs"
  printf '\nEnsuring Apple Silicon and Intel Rust targets are installed...\n'
  rustup target add aarch64-apple-darwin x86_64-apple-darwin
  build_args+=(--target universal-apple-darwin)
  bundle_dir="$SCRIPT_DIR/target/universal-apple-darwin/release/bundle"
else
  bundle_dir="$SCRIPT_DIR/target/release/bundle"
fi

printf '\nBuilding the optimized macOS app and disk image...\n'

# A release uses notarytool's Keychain profile below. Keep Tauri from starting
# a second notarization submission if unrelated Apple variables are present in
# the caller's shell.
env -u APPLE_ID \
  -u APPLE_PASSWORD \
  -u APPLE_TEAM_ID \
  -u APPLE_API_KEY \
  -u APPLE_API_ISSUER \
  -u APPLE_API_KEY_PATH \
  cargo "${build_args[@]}"

app_path="$(find "$bundle_dir" -type d -name 'MetaboKit DDA.app' -print -quit 2>/dev/null || true)"
dmg_path="$(find "$bundle_dir" -type f -name '*.dmg' -print -quit 2>/dev/null || true)"

[[ -n "$dmg_path" ]] || fail "the DMG was not found under $bundle_dir"

if [[ "$release" == true || "$sign_only" == true ]]; then
  printf '\nVerifying the Developer ID signature...\n'
  codesign --verify --strict --verbose=2 "$dmg_path"
fi

if [[ "$release" == true ]]; then
  printf '\nSubmitting the signed DMG to Apple for notarization...\n'
  xcrun notarytool submit "$dmg_path" --keychain-profile "$notary_profile" --wait

  printf '\nStapling and validating the notarization ticket...\n'
  xcrun stapler staple "$dmg_path"
  xcrun stapler validate "$dmg_path"

  printf '\nChecking the finished DMG with Gatekeeper...\n'
  spctl --assess --type open --context context:primary-signature --verbose=2 "$dmg_path"
fi

printf '\nBuild complete. Products:\n'
if [[ -n "$app_path" ]]; then
  printf '  App: %s\n' "$app_path"
fi
printf '  DMG: %s\n' "$dmg_path"

if [[ "$sign_only" == true ]]; then
  printf '\nThe products are signed but not notarized; do not distribute them publicly yet.\n'
elif [[ "$release" == true ]]; then
  printf '\nThe DMG is signed, notarized, stapled, and ready to distribute.\n'
else
  printf '\nThis unsigned build is for local testing only.\n'
fi
