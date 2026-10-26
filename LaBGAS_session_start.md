# LaBGAS: read this at the start of every Claude session

Server-wide instructions for every lab member's Claude Code session on the LaBGAS
server. `/data/CLAUDE.md` is a symlink to this file, so any session started under
`/data` loads it automatically. Edit it here, in LaBGAScore, and commit.

## 1. Read the ground-truth docs first

Before doing anything else in a session, read:

- **LaBGAScore**: `/data/master_github_repos/LaBGAScore/CLAUDE.md` and `README.md`.
- **CANlab_help_examples (LaBGAS fork)**: only
  `/data/master_github_repos/CANlab_help_examples/Second_level_analysis_template_scripts/`
  (`CLAUDE.md`, `README.md`). That is the only folder with LaBGAS-edited scripts;
  ignore the rest of the fork.
- **The subfolder you work in**: its `README*.md` and `CLAUDE.md` (e.g.
  `LaBGAScore/secondlevel/`, `firstlevel/`, `clean/`).
- **A study project**: the project root `CLAUDE.md`, then every `CLAUDE.md`,
  `README*.md` and `BUILD_DECISIONS.md` on the path down to the model directory you
  touch.

| project | root |
|---|---|
| DISCOVERIE | `/data/proj_discoverie` |
| CFS | `/data/proj_cfs` |
| MOODBUGS WP2 | `/data/proj_moodbugs/proj_moodbugs_wp2` (nested) |
| bitter-reward | `/data/proj_bitter-reward` |

Numbers in these docs are dated measurements. Re-derive them before relying on them.

## 2. How to work

1. Ask for clarification when the request is underspecified; do not guess.
2. Once the task is clear, show a plan and get approval before executing.
3. Never delete files, drop annexed content or push to GIN without explicit approval.
   Read the project `CLAUDE.md` on its GIN policy first: content missing from GIN is
   often missing by decision.
4. Never annex files a CONN project reads or writes (unzipped anat/probseg `.nii`,
   `esub-*`, `dssub-*`, `*.QC_*.mat` in `derivatives`) until the CONN analysis is
   signed off: read-only symlinks can break a CONN re-run. Check the project's
   `CONN_x.Setup` file lists first. The `ssub-*` rest images are the deliberate
   exception.
5. Be concise.
