#!/usr/bin/env bash
#
# labgascore_rotate_github_pat.sh
#
# Updates every place the GitHub personal access token is stored, from one file,
# after validating the token against the GitHub API.
#
#   usage: labgascore_rotate_github_pat.sh <path-to-file-containing-the-new-PAT>
#   e.g.   labgascore_rotate_github_pat.sh ~/tokens/Github_PAT_07_10_26
#
# WHY THIS EXISTS. On the LaBGAS server the PAT lives in three places, and all
# three must agree or git and gh fail differently, with a cause that is not
# obvious from either error message:
#
#   1. ~/tokens/<file>         the copy you keep              (source of truth)
#   2. ~/.git-credentials      credential.helper=store     -> git fetch/push
#   3. ~/.config/gh/hosts.yml  gh                          -> gh commands
#
# WHAT AN EXPIRED TOKEN LOOKS LIKE, because the error names the wrong thing:
#
#   remote: Invalid username or token. Password authentication is not supported
#   fatal: Authentication failed for 'https://github.com/labgas/<repo>.git/'
#
# That message appears only AFTER git has already given up on the stored
# credential and prompted interactively, because credential.helper=store ERASES
# a credential the server rejects - so ~/.git-credentials is left EMPTY (0
# bytes) and what failed is the password typed at the prompt, which can never
# work: GitHub removed password authentication for git operations in 2021. The
# real cause is one layer up, and `gh auth status` names it directly:
#
#   The token in ~/.config/gh/hosts.yml is invalid.
#
# Note also that with http.<github>.proactiveAuth set (see README.md, "GitHub
# authentication"), an expired token breaks READS as well as pushes: git demands
# a credential for public repos instead of falling back to anonymous, so
# headless fetches hang or fail rather than quietly working.
#
# SCOPES. Pushing needs only 'repo'. gh additionally requires 'read:org' and
# refuses the token without it, which is why step 3 below is non-fatal: a token
# that is perfectly good for git should not be reported as a failure. Neither
# 'admin:org' nor 'delete_repo' is needed by anything in this workflow.
#
# SAFETY PROPERTIES, all deliberate:
#   - The token is VALIDATED FIRST. Nothing is written if GitHub rejects it, so
#     a typo or a half-copied paste cannot replace working credentials with
#     broken ones.
#   - The token is never echoed (only its length and 4-character prefix).
#   - It is never passed as a command-line argument: argv is world-readable via
#     ps on a shared server. The API call takes it through a curl config on
#     stdin; gh takes it on stdin.
#   - It is never written anywhere except ~/.git-credentials itself, created
#     mode 600 under umask 077. No temporary file ever holds it.
#   - Non-github.com lines in ~/.git-credentials are preserved.
#
# Server-specific by design: it assumes the three locations above and gh
# installed user-local on PATH. It is idempotent - re-running it on the same
# token is safe and simply re-verifies.
#
# Exit status: 0 all good; 1 token rejected or git verification failed;
# 2 usage/unreadable file; 3 git fixed but gh not (missing scope).
#
# -----------------------------------------------------------------------------
#
# author: Lukas Van Oudenhove
# date:   October, 2026
#
# LaBGAS, KU Leuven
#
# version: 1.0

set -euo pipefail

VERIFY_REPO="${LABGAS_VERIFY_REPO:-https://github.com/labgas/LaBGAScore.git}"

TOKFILE="${1:-}"
[ -n "$TOKFILE" ] || { echo "usage: $(basename "$0") <path-to-file-containing-the-new-PAT>" >&2; exit 2; }
[ -r "$TOKFILE" ] || { echo "ERROR: cannot read $TOKFILE" >&2; exit 2; }

umask 077

TOK=$(tr -d ' \t\r\n' < "$TOKFILE")
[ -n "$TOK" ] || { echo "ERROR: $TOKFILE is empty" >&2; exit 2; }
printf 'token file: %s\n  %d chars, prefix %s...\n' \
       "$TOKFILE" "${#TOK}" "$(printf '%s' "$TOK" | cut -c1-4)"

if [ "$(stat -c %a "$TOKFILE")" != "600" ]; then
    echo "  NOTE: $TOKFILE is mode $(stat -c %a "$TOKFILE"); 600 is expected."
    echo "        (Harmless if its parent directory is 700, which blocks traversal"
    echo "        regardless of file mode - check with: stat -c %a \"\$(dirname $TOKFILE)\")"
fi

# ---- 1. validate BEFORE changing anything ----------------------------------
# Header and body go to files; the TOKEN does not - it reaches curl through a
# config read from stdin, so it stays out of argv.
HDR=$(mktemp); BODY=$(mktemp)
trap 'rm -f "$HDR" "$BODY"' EXIT

CODE=$(printf 'header = "Authorization: Bearer %s"\nsilent\nwrite-out = "%%{http_code}"\ndump-header = "%s"\noutput = "%s"\nurl = "https://api.github.com/user"\n' \
        "$TOK" "$HDR" "$BODY" | curl --config -)

if [ "$CODE" != "200" ]; then
    echo "ERROR: GitHub rejected this token (HTTP $CODE). NOTHING WAS CHANGED." >&2
    grep -o '"message": *"[^"]*"' "$BODY" | head -1 | sed 's/^/  /' >&2 || true
    exit 1
fi

LOGIN=$(grep -o  '"login": *"[^"]*"' "$BODY" | head -1 | cut -d'"' -f4)
SCOPES=$(grep -i '^x-oauth-scopes:' "$HDR" | cut -d: -f2- | tr -d '\r' | sed 's/^ *//')
EXPIRY=$(grep -i '^github-authentication-token-expiration:' "$HDR" | cut -d: -f2- | tr -d '\r' | sed 's/^ *//')

echo "  valid, account: $LOGIN"
echo "  scopes        : ${SCOPES:-<fine-grained token; scopes not reported by the API>}"
echo "  expires       : ${EXPIRY:-<no expiry set>}"

case ",$(printf '%s' "$SCOPES" | tr -d ' ')," in
  *,admin:org,*|*,delete_repo,*)
    echo "  NOTE: this token carries admin:org and/or delete_repo. Nothing in this"
    echo "        workflow needs either; 'repo' + 'read:org' is sufficient." ;;
esac
case ",$(printf '%s' "$SCOPES" | tr -d ' ')," in
  *,repo,*) ;;
  *) echo "  WARNING: no 'repo' scope - git push to a private repo will fail." ;;
esac

# ---- 2. ~/.git-credentials, preserving any other host ----------------------
# Read the lines to keep into a variable, then write the file ONCE. The token
# therefore never lands in a temp file - only in $CRED, created mode 600.
CRED="$HOME/.git-credentials"
KEEP=""
if [ -f "$CRED" ]; then
    KEEP=$(grep -v '@github\.com$' "$CRED" || true)
fi
: > "$CRED"
chmod 600 "$CRED"
[ -n "$KEEP" ] && printf '%s\n' "$KEEP" >> "$CRED"
printf 'https://%s:%s@github.com\n' "$LOGIN" "$TOK" >> "$CRED"
echo "  wrote $CRED ($(grep -c . "$CRED") line(s), mode $(stat -c %a "$CRED"))"

# ---- 3. gh -- NON-FATAL ON PURPOSE -----------------------------------------
# gh refuses a classic token lacking 'read:org' even though git push needs only
# 'repo'. Aborting here would hide a working git setup behind a gh-only problem
# and skip the verification below. Measured 2026-10-07: a repo+user token
# rotated git correctly and failed here with
#   error validating token: missing required scope 'read:org'
export PATH="$HOME/.local/bin:$PATH"
GH_OK=skipped
if command -v gh >/dev/null 2>&1; then
    GHERR=$(mktemp)
    trap 'rm -f "$HDR" "$BODY" "$GHERR"' EXIT
    if printf '%s\n' "$TOK" | gh auth login --hostname github.com --with-token 2>"$GHERR"; then
        GH_OK=yes
        echo "  updated gh ($(gh --version | head -1))"
    else
        GH_OK=no
        echo "  gh NOT updated:"
        sed 's/^/    /' "$GHERR" | head -2
        if grep -q 'read:org' "$GHERR"; then
            echo "    -> add 'read:org' to this token on GitHub. A classic token can be"
            echo "       edited in place, no need to regenerate: Settings > Developer"
            echo "       settings > Personal access tokens > the token > tick read:org >"
            echo "       Update token. Then re-run this script. git is unaffected."
        fi
    fi
else
    echo "  gh not on PATH - hosts.yml not touched"
fi

# ---- 4. prove it, end to end -----------------------------------------------
echo
echo "verifying:"
if command -v gh >/dev/null 2>&1; then gh auth status 2>&1 | sed 's/^/  /' || true; fi
if GIT_TERMINAL_PROMPT=0 git ls-remote --heads "$VERIFY_REPO" >/dev/null 2>&1; then
    echo "  authenticated ls-remote OK - GIT FETCH AND PUSH WILL WORK"
else
    echo "  authenticated ls-remote FAILED against $VERIFY_REPO" >&2
    echo "  check the scopes reported above" >&2
    exit 1
fi

if [ "$GH_OK" = "no" ]; then
    echo
    echo "SUMMARY: git is FIXED; gh is NOT (missing scope, see above)."
    exit 3
fi
echo
echo "SUMMARY: git and gh both updated."
