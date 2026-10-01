# clean — infrastructure: provenance, dependencies, checkers, housekeeping

The only folder here that is not part of an analysis. Nothing in it produces a
result; it records how results were produced, catches the mistakes that waste a
day, and does the chores.

Almost everything is a **function called from elsewhere** rather than a script
you run by hand — 16 of the 17 `.m` files are functions, the exception being
`LaBGAScore_smart_parallel_pool_setup.m`. The two `.py` checkers and the two
`.sh` wrappers are run directly.

**[`README_provenance.md`](README_provenance.md) is authoritative for the
`prov_*` and `dep_*` tooling** — the design, the six rules that stop the call
graph inventing dependencies, the negative controls, and the per-dataset state.
This file is an index, and does not repeat it.

---

## What is here

### Provenance — which version of everything produced a result

| file | |
|---|---|
| `LaBGAScore_prov_publish.m` | drop-in replacement for the bare `publish()` in every script header; adds a provenance table to the HTML report plus a `.tsv`/`.mat` sidecar |
| `LaBGAScore_prov_snapshot.m` | records which version of every dependency was in place at the time of a run |
| `LaBGAScore_prov_resolve_retrospective.m` | reconstructs the same information for analyses **already** run, from each artifact's embedded date and each clone's reflog. Writes sidecars; never modifies existing files |
| `LaBGAScore_prov_gitinfo.m` | reads commit, branch, remote and dirty state straight from the plain-text files under `.git` — nothing shells out |
| `LaBGAScore_prov_gitstatus.m` | reimplements git's own dirty check; verified to reproduce `git status --porcelain` exactly |
| `LaBGAScore_prov_protect_reflogs.m`, `labgascore_prov_protect_reflogs.sh` | disables `gc.reflogExpire`, which otherwise prunes the perishable evidence the retrospective depends on after 90 days, silently |

### Dependencies — the generated documentation

| file | |
|---|---|
| `LaBGAScore_dep_build_index.m` | indexes every callable file under a root |
| `LaBGAScore_dep_map.m` | maps which dependency files a script actually reaches |
| `LaBGAScore_dep_report.m` | **generates** `DEPENDENCIES.md`, `dependencies.tsv` and `dependencies.yml` — never hand-edit those |

### Checkers — run all three before launching a long chain

| file | catches | |
|---|---|---|
| `LaBGAScore_check_all_scripts.m` | parse errors, plus style/performance suggestions listed separately | wraps MATLAB's `checkcode` over a whole tree |
| `use_before_def.py` | an option **read above the line that defines it** | a gate — dies at runtime with *"Unrecognized function or variable"* |
| `set_after_use.py` | an option **set below the line that already consumed it** | **advisory, not a gate** — legitimate reuse (`savefilename`, `figtitle`) trips it; expect a handful per model and read them |

`checker_positive_controls/` holds one control file per Python checker plus
`run_controls.sh`, which asserts that each checker flags its own control and
ignores the other's. **`checkcode` reports zero messages on either control
file** — measured, not assumed, which is exactly why these exist alongside it.
If a change to a checker makes its control pass, the checker is broken.

Between them these three cover parse errors and the two ordering failures. They
do **not** catch undefined variables used at runtime, calls to functions that do
not exist, or logic bugs — those still need reading.

### Running and publishing

| file | |
|---|---|
| `LaBGAScore_run_reports.m` | publishes a chain of scripts and **fails loudly when one errors**. `publish()` catches a script error into the HTML and returns normally, so a crashed run otherwise raises nothing and exits 0 — the single biggest source of silently wrong results in this pipeline |
| `labgascore_run_headless.sh` | wraps the headless invocation, including the traps: `matlab -batch` cannot `publish()`, stdin must be redirected from `/dev/null`, and the MATLAB code goes into a temporary `.m` file rather than a one-line `-r` string |
| `LaBGAScore_check_display.m` | reports whether the current X2go session is big enough for a full-size report figure, and what to change if not. Applies to the **interactive** route only |

### One-off corrections

| file | |
|---|---|
| `LaBGAScore_stats_rederive_storey_q.m` | re-derives Storey q-values in **already saved** second-level results and reports what changes, without re-running any analysis. Possible because `q = pi0 * q_BH` floored at `p` is a scalar transform of p-values that the results tables already store next to the q-values — everything expensive produced those p-values. Report-only by default; `'write', true` with `'output','sidecar'` (the default) saves corrected tables as `*_storeyfix.mat` and leaves the originals alone, `'output','inplace'` overwrites them, and `'csv',true` writes a flat `.csv` beside each |

**Reading its report: `n_sig_BH` first, not `n_sig_old`.** The stored `q_Storey`
column is not always a correction — where the old pi0 collapsed towards 0 the
`q >= p` floor handed the raw p-values back under an FDR name, which was the case
in 47 of 73 families across the eight final models. Scoring a corrected q against
that baseline makes every genuine improvement look like a lost result. `q_BH` is
the reference that is always valid, and an upper bound: since pi0 ≤ 1 the new q
can never be *stricter* than BH, so `n_sig_new > n_sig_BH` is power gained and
`n_sig_old > n_sig_BH` is significance the old column claimed without support.

**`'output','inplace'` refuses to write through a git-annex symlink**, and says
so rather than failing quietly. git-annex deduplicates by content hash, so one
object can back the identical table in several models; saving over the link
would edit that shared object and corrupt every model pointing at it. Run
`git annex unlock <file>` first — which leaves the previous version in annex
history, so the overwrite stays reversible.

> **The eight final models were corrected in place on 2026-10-01** and each
> carries `results/README_storey_correction_2026-10-01.md` recording that the
> tables now **supersede the `q_Storey` columns in the published
> `results/html/` reports**, which were not regenerated. Everything else in
> those reports — p-values, `q_BH`, effect sizes, all voxelwise results — is
> unaffected.

### Housekeeping

| file | |
|---|---|
| `LaBGAScore_clean_gzip_all_nii.m` | gzips `.nii` recursively under a folder |
| `LaBGAScore_clean_sourcedata.m` | cleans the `sourcedata` subdataset, branching on whether DICOMs are gitignored or annexed |
| `LaBGAScore_move_repos_matlabpath.m` | moves one repo below another on the MATLAB path and calls `savepath` — for when path order decides which of two same-named functions wins |
| `LaBGAScore_smart_parallel_pool_setup.m` | sizes a `parpool` to a fraction of available cores (~60%). The one **script** here; called from the second-level scripts before bootstrapping or permutation |

## Conventions

Two small things specific to this folder:

- **These functions use the `*USAGE*` header style** of the repo's *scripts*,
  not the MATLAB H1-line convention that `*/functions/` files use. Harmless, but
  it means `help <name>` gives you a section heading rather than a one-line
  summary, and the repo root's description of the two conventions does not quite
  cover this folder.
- **`.dep_index_*.mat` is a machine-specific cache** and is gitignored
  (`.gitignore:11`), so the one you will find here belongs to whichever machine
  last built an index.

## Dependencies

No CanlabCore or SPM dependency — that is the point of keeping this folder
separate. The provenance functions read `.git` as plain text rather than shelling
out to `git`, so they work on a clone with no `git` binary on `PATH`.
`matlab.codetools.requiredFilesAndProducts` is **not** usable in this tree and
should not be re-attempted: it aborts entirely on a syntax error anywhere in the
transitive closure, and there are 33 such files across the cloned repos.
`LaBGAScore_dep_map.m` records those as `unparseable` and carries on.

## Status

This folder is recent work and has been exercised in anger: the retrospective has
been run over two studies at both levels, `prov_protect_reflogs` over all 29
clones on the server, and the checkers against their own positive controls. Code
and headers were reviewed as they were written, so no audit findings are recorded
here — unlike `prep/` and `decoding_toolbox/`, which needed one.

The gap worth knowing about: `LaBGAScore_dep_map.m` walks `.m` files only, so
`stats_tools/sas_macros/` is covered by nothing.
