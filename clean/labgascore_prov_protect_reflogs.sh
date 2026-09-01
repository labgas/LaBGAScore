#!/usr/bin/env bash
#
# labgascore_prov_protect_reflogs.sh
#
# Disables git reflog expiry on every repository cloned under a given root.
#
# The reflog records which commit each clone actually had checked out at
# any past moment, and is the only way to reconstruct the dependency
# versions behind an analysis that was run months ago. Git's defaults
# (gc.reflogExpire = 90 days, gc.reflogExpireUnreachable = 30 days) prune
# it silently during garbage collection. This makes both settings "never".
#
# Companion to LaBGAScore_prov_protect_reflogs.m, which archives the
# reflogs and reports which repositories are unprotected, but deliberately
# does not shell out to change anything.
#
# Idempotent, non-destructive, and touches nothing but those two config
# keys. Undo for one repo with:
#   git -C <repo> config --unset gc.reflogExpire
#   git -C <repo> config --unset gc.reflogExpireUnreachable
#
# Usage: bash labgascore_prov_protect_reflogs.sh [githubrootdir]
#
# ---------------------------------------------------------------------------
# Lukas Van Oudenhove & Claude Opus 5, KU Leuven, August 2026
# ---------------------------------------------------------------------------

set -euo pipefail

ROOT="${1:-/data/master_github_repos}"

if [ ! -d "$ROOT" ]; then
    echo "error: $ROOT does not exist" >&2
    exit 1
fi

total=0
changed=0
failed=0

for d in "$ROOT"/*/; do
    repo="$(basename "${d%/}")"
    [ -e "$d/.git" ] || continue
    total=$((total + 1))

    before="$(git -C "$d" config --get gc.reflogExpire 2>/dev/null || true)"

    if git -C "$d" config gc.reflogExpire never 2>/dev/null \
    && git -C "$d" config gc.reflogExpireUnreachable never 2>/dev/null; then
        if [ "$before" != "never" ]; then
            changed=$((changed + 1))
            printf '  %-32s protected\n' "$repo"
        fi
    else
        failed=$((failed + 1))
        printf '  %-32s FAILED\n' "$repo" >&2
    fi
done

echo
echo "$total git repositories under $ROOT"
echo "$changed newly protected, $((total - changed - failed)) already protected, $failed failed"

if [ "$failed" -gt 0 ]; then
    exit 1
fi
