#!/bin/bash
#
# labgascore_run_headless.sh -- publish LaBGAS analysis reports without a display
#
# USAGE
#
#   labgascore_run_headless.sh -d <projdir> -s <setup_script> <script> [<script> ...]
#
#   ALL OPTIONS MUST COME BEFORE THE SCRIPT NAMES. Anything after the first
#   script name is treated as another script.
#
#   -d <projdir>       superdataset root to cd into before anything else.
#                      REQUIRED: prep_s0_define_directories derives every path
#                      from pwd, so running from the wrong directory silently
#                      points the whole analysis at the wrong tree.
#   -s <setup_script>  script that defines htmlsavedir (and DAT). Normally your
#                      <proj>_secondlevel_m<M>_s0_a_set_up_paths_always_run_first.
#                      Repeat -s to run several setup scripts in order, e.g. s0,
#                      then prep_1, then prep_1b.
#   -a <artefact>      optional. Expected output of the correspondingly numbered
#                      script, checked after publishing. Repeat once per script,
#                      using "-" for scripts with nothing to check.
#   -p <dir>           optional. addpath(genpath(<dir>)) before running, appended
#                      in the order given. Repeat as needed. Use this if your
#                      startup.m does not already put CanlabCore, LaBGAScore, the
#                      CANlab_help_examples fork and your study's code subdataset
#                      on the path. SPM12 should be first on the path, so pass it
#                      with -P instead.
#   -P <dir>           optional. addpath(<dir>,'-begin'), i.e. PREPEND without
#                      genpath. This is what SPM12 wants: /opt/KUL_apps/spm12.
#   -l <logfile>       optional. Default: headless_<timestamp>.log in the cwd.
#   -h                 this help.
#
# EXAMPLES
#
#   # one script
#   labgascore_run_headless.sh -d /data/proj_xxx \
#       -s proj_secondlevel_m1_s0_a_set_up_paths_always_run_first \
#       proj_secondlevel_m1_s4_prep_2_load_image_data_and_save
#
#   # a whole chain, with artefact checks
#   labgascore_run_headless.sh -d /data/proj_xxx \
#       -s proj_secondlevel_m1_s0_a_set_up_paths_always_run_first \
#       -s proj_secondlevel_m1_s2_prep_1_set_conditions_contrasts \
#       -s proj_secondlevel_m1_s3_prep_1b_behavioral_data \
#       -a data_objects.mat -a contrast_data_objects.mat -a - \
#       proj_secondlevel_m1_s4_prep_2_load_image_data_and_save \
#       proj_secondlevel_m1_s5_prep_3_calc_univariate_contrasts \
#       proj_secondlevel_m1_s6_prep_3a_run_second_level_regression
#
#   # long chain in the background, survives logout
#   setsid nohup labgascore_run_headless.sh -d /data/proj_xxx -s ... > run.log 2>&1 < /dev/null &
#
# WHY THE MATLAB INVOCATION LOOKS LIKE THIS
#
#   * "matlab -batch" CANNOT be used. publish() is unsupported there and fails
#     with "Unable to run the 'publish' function, because it is not supported
#     for this ...". Use -nodisplay with -r, as below.
#
#   * "< /dev/null" is required. With -r, MATLAB keeps reading commands from
#     stdin after running them; if stdin is an open terminal the session hangs
#     around, and if it is a closed pipe MATLAB can exit before running anything.
#     Redirecting from /dev/null makes the behaviour deterministic.
#
#   * The trailing "exit" is inside the -r string, so MATLAB always terminates.
#
#   * Figures are NOT capped by the virtual screen. Headless, publish() prints
#     figures rather than screen-capturing them, so the 1024x768 / 72 dpi
#     headless "screen" does not limit figure size the way an X2go session does.
#
# SEE ALSO
#
#   LaBGAScore_run_reports.m   -- does the publishing and the failure detection
#   LaBGAScore_check_display.m -- for interactive X2go sessions instead
#   LaBGAS_fMRI_analysis_workflow.md, section "Running scripts and publishing reports"
#
# -----------------------------------------------------------------------------
# Lukas Van Oudenhove, KU Leuven, September 2026
# -----------------------------------------------------------------------------

set -u

PROJDIR=""
SETUP=()
ARTEFACTS=()
ADDPATHS=()
PREPATHS=()
LOGFILE=""

usage () { sed -n '2,/^# ----/p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

while getopts ":d:s:a:p:P:l:h" opt; do
    case $opt in
        d) PROJDIR="$OPTARG" ;;
        s) SETUP+=("$OPTARG") ;;
        a) ARTEFACTS+=("$OPTARG") ;;
        p) ADDPATHS+=("$OPTARG") ;;
        P) PREPATHS+=("$OPTARG") ;;
        l) LOGFILE="$OPTARG" ;;
        h) usage 0 ;;
        \?) echo "unknown option -$OPTARG" >&2; usage 1 ;;
        :)  echo "option -$OPTARG needs an argument" >&2; usage 1 ;;
    esac
done
shift $((OPTIND - 1))

SCRIPTS=("$@")

# All options must precede the script list: getopts stops at the first
# non-option, so a flag written after the scripts is silently swallowed as a
# script name. Catch that here instead of failing later with a confusing count.
for s in ${SCRIPTS[@]+"${SCRIPTS[@]}"}; do
    case "$s" in
        -*) echo "ERROR: option '$s' appears after the script list." >&2
            echo "       All options (-d -s -a -p -P -l) must come BEFORE the script names." >&2
            exit 1 ;;
    esac
done

[ -z "$PROJDIR" ]        && { echo "ERROR: -d <projdir> is required" >&2; usage 1; }
[ ${#SETUP[@]} -eq 0 ]   && { echo "ERROR: at least one -s <setup_script> is required" >&2; usage 1; }
[ ${#SCRIPTS[@]} -eq 0 ] && { echo "ERROR: no scripts to publish" >&2; usage 1; }
[ -d "$PROJDIR" ]        || { echo "ERROR: no such directory: $PROJDIR" >&2; exit 1; }

if [ ${#ARTEFACTS[@]} -gt 0 ] && [ ${#ARTEFACTS[@]} -ne ${#SCRIPTS[@]} ]; then
    echo "ERROR: -a given ${#ARTEFACTS[@]} time(s) for ${#SCRIPTS[@]} script(s);" >&2
    echo "       give one -a per script, using '-' where there is nothing to check" >&2
    exit 1
fi

[ -z "$LOGFILE" ] && LOGFILE="headless_$(date +%Y%m%d_%H%M%S).log"

# ---- build the MATLAB command ------------------------------------------------

mat_cellstr () {   # quote each argument as a MATLAB char in a cellstr
    local out="{" first=1
    for a in "$@"; do
        [ $first -eq 0 ] && out="$out, "
        out="$out'${a//\'/\'\'}'"
        first=0
    done
    echo "$out}"
}

SETUP_CELL=$(mat_cellstr "${SETUP[@]}")
SCRIPT_CELL=$(mat_cellstr "${SCRIPTS[@]}")

if [ ${#ARTEFACTS[@]} -gt 0 ]; then
    ART_CELL="{"
    for i in "${!ARTEFACTS[@]}"; do
        [ "$i" -gt 0 ] && ART_CELL="$ART_CELL, "
        if [ "${ARTEFACTS[$i]}" = "-" ]; then
            ART_CELL="$ART_CELL[]"
        else
            ART_CELL="$ART_CELL'${ARTEFACTS[$i]}'"
        fi
    done
    ART_CELL="$ART_CELL}"
    ART_ARG=", 'artefacts', $ART_CELL"
else
    ART_ARG=""
fi

PATHCODE=""
for d in ${PREPATHS[@]+"${PREPATHS[@]}"}; do
    PATHCODE="$PATHCODE addpath('${d//\'/\'\'}','-begin');"
done
for d in ${ADDPATHS[@]+"${ADDPATHS[@]}"}; do
    PATHCODE="$PATHCODE addpath(genpath('${d//\'/\'\'}'),'-end');"
done

# The MATLAB code is written to a temporary .m file and run(), rather than
# squeezed into -r as one line. Flattening newlines into a single -r string is
# what breaks such commands: a MATLAB "%" comment then comments out everything
# after it, and "for"/"if" headers run into their bodies.
MFILE=$(mktemp "${TMPDIR:-/tmp}/labgascore_headless_XXXXXX.m")
if [ "${LABGASCORE_KEEP_MFILE:-0}" = "1" ]; then
    echo "    mfile   : $MFILE (kept: LABGASCORE_KEEP_MFILE=1)"
else
    trap 'rm -f "$MFILE"' EXIT
fi

PROJDIR_ESC="${PROJDIR//\'/\'\'}"

# The static MATLAB code goes through a QUOTED heredoc so the shell performs no
# expansion or escape processing on it. Building these lines with printf in
# double quotes does not work: "\\n" is collapsed to a real newline, which then
# splits a MATLAB char literal across two lines and yields
# "Character vector is not terminated properly".
{
    printf '%s\n' "$PATHCODE"
    printf "cd('%s');\n" "$PROJDIR_ESC"
    printf 'setup = %s;\n' "$SETUP_CELL"
    cat <<'MEOF'
for k = 1:numel(setup)
    fprintf('--- setup: %s ---\n', setup{k});
    eval([setup{k} ';']);
end
if ~exist('htmlsavedir','var') || isempty(htmlsavedir)
    error('htmlsavedir was not defined by the setup script(s)');
end
MEOF
    # prep_s0 derives every path from pwd, and publish() can change directory
    printf "cd('%s');\n" "$PROJDIR_ESC"
    printf 'results = LaBGAScore_run_reports(%s, htmlsavedir%s);\n' "$SCRIPT_CELL" "$ART_ARG"
    printf 'if ~all(results.ok), exit(1); end\n'
} > "$MFILE"

echo "=== headless run $(date '+%F %T') ==="
echo "    project : $PROJDIR"
echo "    paths   : ${PREPATHS[*]+${PREPATHS[*]}} ${ADDPATHS[*]+${ADDPATHS[*]}}"
echo "    setup   : ${SETUP[*]}"
echo "    scripts : ${SCRIPTS[*]}"
echo "    log     : $LOGFILE"
echo

# -batch cannot publish(); -r with stdin from /dev/null is the supported route.
matlab -nodisplay -nosplash \
    -r "try, run('$MFILE'); catch e, disp(getReport(e)); exit(1); end, exit(0)" \
    < /dev/null 2>&1 | tee "$LOGFILE"

STATUS=${PIPESTATUS[0]}

echo
if [ "$STATUS" -eq 0 ]; then
    echo "=== all reports published $(date '+%F %T') ==="
else
    echo "=== FAILURES: see $LOGFILE, and grep for '>>>' ===" >&2
fi
exit "$STATUS"
