# clean — infrastructure: provenance, dependencies, checkers, housekeeping

The only folder here that is not part of an analysis. Nothing in it produces a
result; it records how results were produced, catches the mistakes that waste a
day, and does the chores.

Almost everything is a **function called from elsewhere** rather than a script
you run by hand — 16 of the 17 `.m` files are functions, the exception being
`LaBGAScore_smart_parallel_pool_setup.m`. The two `.py` checkers and the three
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

`checker_positive_controls/` holds `run_controls.sh` plus **three** checker
control files — two for `use_before_def.py`, one per guard spelling, and one for
`set_after_use.py`. The runner asserts that each checker flags its own
control(s) and ignores the other's. (`control_mvpa_reg_cov_permutation.m` also
lives there but is a different animal: a runnable statistical control for
`prep_3a`'s permutation block, not a checker control. Both checkers scan it and
neither should flag it, which is why the runner reports four scripts.) **`checkcode` reports zero messages on any
control file** — measured, not assumed, which is exactly why these exist
alongside it. If a change to a checker makes its control pass, the checker is
broken.

**Why `use_before_def.py` needs two controls.** Until 2026-10-05 it built its
set of option-ish names from bare `name = value` assignments only. The templates
overwhelmingly use the one-line guard instead —

```matlab
if ~exist('x','var'), x = []; end
```

— whose assignment is not at the start of a line, so those options were never
recognised as options, recorded **zero uses**, and passed unconditionally. In
`prep_3a` that was **all 11 of 11** guarded options: the check could only ever
return PASS, and did. It was hiding a live case — `cv_seed_mvpa_reg_cov`, read
by the `tuned_seed` default four lines above its own guard, which would have
died on *"Unrecognized function or variable"* in any model script that did not
set `cv_seed_mvpa_reg_cov` itself (`a2_set_default_options` always does, which is
what masked it). The fix adds guard-form names to that set;
`control_use_before_def_guard.m` is the regression test, and the original
control covers the three-line spelling, whose assignment the old pattern did
see. Re-run across every template and study model directory after the fix: no
other live case.

**The lesson generalises past this checker.** A checker that silently examines
nothing reports PASS, which is indistinguishable from a clean result and worse
than no checker at all. A positive control per *idiom the checker must handle*,
not per *failure it was written for*, is what catches that.

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

### GitHub authentication, and rotating the PAT

| file | |
|---|---|
| `labgascore_rotate_github_pat.sh` | updates every place the GitHub PAT is stored, from one file, **after validating it against the API** so a bad paste cannot replace working credentials with broken ones |

All clones under `/data/master_github_repos` use HTTPS remotes and share one
PAT, stored in **three** places that must agree — `~/tokens/<file>` (the copy
you keep), `~/.git-credentials` (`credential.helper=store`, used by git), and
`~/.config/gh/hosts.yml` (used by `gh`). The script writes the last two from the
first and verifies with an authenticated `ls-remote`.

**Two things make an expired token hard to recognise.** First, the error names
the wrong thing:

```
remote: Invalid username or token. Password authentication is not supported
fatal: Authentication failed for 'https://github.com/labgas/<repo>.git/'
```

That appears only *after* git has given up on the stored credential and
prompted, because `credential.helper=store` **erases** a credential the server
rejects — so `~/.git-credentials` is left **empty (0 bytes)** and what actually
failed is the password typed at the prompt, which can never work (GitHub
removed password auth for git operations in 2021). `gh auth status` names the
real cause directly: *"The token in ~/.config/gh/hosts.yml is invalid."*

Second, **an expired token breaks reads too**, not just pushes. `~/.gitconfig`
carries

```
git config --global 'http.https://github.com/.proactiveAuth' basic
```

which exists because git sends no credentials when fetching a *public* repo —
GitHub answers the anonymous request with 200, so git never sees the 401 that
would make it consult the credential helper, and every fetch went against the
server's IP-wide anonymous quota, which GitHub began throttling on 2026-09-03.
(Neither `gh auth login` nor `gh auth setup-git` fixes that: they register a
helper, and a helper that is never consulted changes nothing. Requires git
≥ 2.46.) The side effect is that with no valid credential git now *demands* one
for public repos instead of falling back to anonymous, so headless fetches hang
or fail rather than quietly working.

**Scopes:** pushing needs only `repo`. `gh` additionally requires `read:org`
and refuses a classic token without it — which is why the script's `gh` step is
non-fatal and exits 3 rather than 1: a token that is fine for git should not be
reported as a total failure. Nothing here needs `admin:org` or `delete_repo`.

### Reclaiming space with `git annex drop` — nine ways it goes wrong

None of this is in the DataLad docs in one place, and all five have cost real
time on this server. Measured examples are from `proj_bitter-reward`, 2026-10-02.

**1. `dropunused` works from a single shared list, not one list per remote.**
`git annex unused` writes `.git/annex/unused`, and **every invocation overwrites
it**. `git annex dropunused --from <remote> all` then operates on whatever that
file currently holds, regardless of the `--from` you pass to the *drop*. So this
sequence silently reclaims nothing:

```bash
git annex unused --from gin          # 319 keys, 163.89 GB  -> written to the list
git annex unused                     # 4 local keys         -> OVERWRITES the list
git annex dropunused --force --from gin all   # drops those 4, not the 319
```

It exits 0 and reports success. Always regenerate the list for the remote
*immediately* before the drop, in one chain, with nothing in between:

```bash
git annex unused --from gin && git annex dropunused --force --from gin all
```

**2. Push the deletion before dropping.** Until the commit that removed the files
is on the remote, that remote's `master` still references their content, so the
content is not "unused" there and the drop reclaims nothing — again silently.

**3. Diff the key sets before any forced drop, because annex deduplicates by
content hash.** One object can back identical files in several models. In
`proj_bitter-reward/firstlevel`, all **84** `mask.nii` files across four models
shared a single key — SPM's analysis mask does not depend on the motion
regressors, so the defective and corrected fits produced byte-identical masks. A
path-scoped `git annex drop --force model_1_food_images model_2_FID` would have
destroyed the mask for the two models being *kept*.

**Prefer `dropunused` precisely because its key list is computed from
reachability rather than from paths** — a key still referenced by anything cannot
appear in it. Verify anyway, since it is two commands:

```bash
git annex unused --from gin | grep -oE '(MD5E|SHA256E)-[^ ]+' | sort -u > /tmp/unused.txt
git annex find --include '*' --format='${key}\n' | sort -u > /tmp/used.txt
comm -12 /tmp/unused.txt /tmp/used.txt     # MUST be empty
```

**4. `numcopies` will refuse, and `--force` means the content is gone.** Where
the remote holds the only copy — normal here, since local content is routinely
dropped — the drop fails with *"Could not verify the existence of the 1 necessary
copy"* on every key. Overriding with `--force` destroys the last copy: git keeps
the filenames in history, so checking out the pre-removal commit afterwards gives
broken annex symlinks with nothing behind them. For a superseded model that is
the intent; know that refitting is the only way back.

Related: **drop from the remote first, while you still hold a local copy.** Then
`numcopies` is satisfied by the local copy and that step needs no `--force` at
all — only the final copy does. Doing it the other way round forces both steps
and gives up the check earlier than necessary.

**5. After dropping, push again — the location log changed even though nothing
else did.** A drop does not touch the working tree or `master`; `git status` stays
clean and `master` stays at the same commit, so it looks like there is nothing to
save or push. But the drop records "this remote no longer has these keys" in the
location log, which lives on the **`git-annex` branch**. Measured after the two
drops above: `master` 0 ahead, `git-annex` **3 ahead** in `derivatives` and **2
ahead** in `firstlevel`.

Until that is pushed, the remote's own copy of the log still claims it holds
content it has deleted, so a later `datalad get` from a fresh clone is told the
content is available, tries to fetch it, and fails confusingly instead of
reporting cleanly that it is gone. `datalad save` is **not** what you want here —
nothing in the tree changed — just:

```bash
datalad push --to gin      # transfers no content, only the log commits
```

A worked sequence that gets the first five right:

```bash
cd <subdataset>
git rm -r <superseded_model> && datalad save -m "remove ..."   # 1. remove
datalad push --to gin                                          # 2. push FIRST
git annex unused --from gin                                    # 3. list, then verify
#    ... run the comm check above ...
git annex dropunused --force --from gin all                    # 4. drop, nothing in between
datalad push --to gin                                          # 5. push the location log
```

Verify afterwards that nothing still needed was taken — this should print 0:

```bash
git annex find --not --in gin | wc -l
```

The five above all concern **unused** keys — content no longer referenced by the
tree. The next four were found on 2026-10-08 while reclaiming ~1.1 TB across
`proj_cfs` and `proj_discoverie`, and three of them only bite once you start
dropping content that is **still tracked**.

**6. Dropping tracked content is a different operation, and there `datalad push`
DOES undo it.** Step 5 above is safe for unused keys precisely because nothing in
the tree references them, so the push moves no content. But a *selective* drop —
keeping a model's files in the tree while removing their content from the remote
— is the opposite case: `datalad push` exists to make the sibling hold everything
the tree references, so it re-uploads exactly what you just dropped. Measured on
`proj_cfs/secondlevel`: 6.24 GB dropped from five tracked models, then a
`datalad push --to gin` intended only to sync the location log began restoring it
and was cut off by a 2-minute timeout, leaving models 6/10/11 fully back and
model_7 half back — a partial-transfer fingerprint, not a bookkeeping glitch.
After a selective drop of **tracked** content, publish refs only:

```bash
git push git@gin.g-node.org:/labgas/<repo>.git master:master git-annex:git-annex
```

That carries the location log and no content. Verify against the remote itself
rather than the log: `git annex checkpresentkey <key> gin` queries it directly and
signals through its **exit status** (0 present, 1 absent) while printing nothing,
so capture `$?`.

**7. `git annex unused` counts the branch's whole history, so it understates
badly.** A version superseded in place is still reachable from the commit that
held it, and `unused` therefore calls it used. To count only the branch tip:

```bash
git annex unused --used-refspec='+refs/heads/master'
```

On `proj_discoverie/derivatives` the default reported **18801 keys / 83.90 GB**;
the tips-only form reported **22504 keys / 500.71 GB**. The 416.80 GB difference
was almost entirely 431 unzipped `s6-*task-MIST*.nii` images superseded in place
by their `.nii.gz` twins. The refspec needs a `+`/`-` prefix — `--used-refspec='HEAD'`
errors with *"bad refspec item"*, which reads like nothing to do.

Dropping that extra set means old commits can no longer be checked out **with
content**; the current tree is untouched. Decide that deliberately.

**Every `unused` run rewrites `.git/annex/unused`, which is the numbering
`dropunused` consumes.** Run the default scan between a tips-only scan and the
drop and `dropunused 1-22504` silently addresses a different, smaller set. Re-run
the exact scan you intend to drop from, immediately before dropping, and do not
interleave another.

**8. A local unused list cannot see keys that were never local.** `git annex
unused` examines the local object store, so content that only ever existed on the
remote is invisible to it. Scan the remote instead:

```bash
git annex unused --from gin
```

This found **138.18 GB** in `proj_discoverie/firstlevel` (the whole
`model_1_basic` first-level output, long since removed from the tree) and
**55.23 GB** in `derivatives`, neither of which appeared in any local list. Such
keys are single-copy by definition, so dropping them needs `--force` and is
permanent — diff against the current tree first:

```bash
git annex unused --from gin | awk '/^ +[0-9]+ +[A-Z]/{print $2}' | sort -u > /tmp/u.txt
git annex find --format='${key}\n' | sort -u > /tmp/cur.txt
comm -12 /tmp/cur.txt /tmp/u.txt | wc -l        # must be 0
```

**9. The `gin` remote's fetch URL is https while its push URL is ssh.** In 14 of
18 subdatasets across `proj_cfs` and `proj_discoverie`. Push works; anything that
**reads** from GIN blocks forever on a username prompt that never arrives, with
no error. This stalled a `git annex drop` in `proj_discoverie/firstlevel` for 68
minutes at zero I/O, on paths that did not even exist. Two consequences: address
GIN by its ssh URL explicitly, and treat every `gin/master` remote-tracking ref as
stale — `git rev-list --count gin/master..HEAD` is then meaningless. Get the real
tip with:

```bash
git ls-remote "$(git remote get-url --push gin)" refs/heads/master
```

#### Telling a re-gzip from a real result

Two places in this workflow produce a file that differs from its recorded key
while containing the same data, and both look like new output:

- **A typechange whose size matches the committed key exactly but whose md5
  differs** is a recompression. Only the gzip header's embedded mtime changed. 38
  such `s6-*rest*.nii.gz` files (53.65 GB) turned up in
  `proj_discoverie/derivatives`; saving them would have added 53.65 GB of
  duplicate content locally and the same again to the push backlog. `git checkout`
  restored the symlinks instead. Confirm by comparing **decompressed** md5s before
  deciding.
- **A `.nii` orphan alongside a tracked `.nii.gz`**: because annex keys are
  `MD5E-s<size>--<md5>` — size plus MD5 of the content — you can derive the key
  the unzipped twin *would* have and look it up directly:

```bash
sz=$(gunzip -c file.nii.gz | wc -c); md5=$(gunzip -c file.nii.gz | md5sum | cut -d' ' -f1)
grep -F "MD5E-s${sz}--${md5}.nii" /tmp/unused_list.txt
```

A hit proves `gunzip` reproduces the orphan bit-for-bit. This replaced
`git annex whereused`, which managed 88 of 19238 keys in several minutes and had
to be abandoned.

#### One more on #3, because it recurred

The key-set diff in #3 is needed on **every** drop, including ones that look safe
because the dropped paths all have local copies — local copies protect the
*dropped* model, not the *kept* ones. Skipped on `proj_discoverie/secondlevel`, a
13.75 GB drop of models 20–30 removed **9 keys shared with the kept models**: the
per-ROI mask files (`Amyg_L.nii`, `aINS_R.nii`, `Tha_L.nii` …) are byte-identical
across model directories and so share one key. 49 files in `model_2h`–`model_2l`
became unavailable on GIN; found by re-running `--not --in=gin` on the kept models
afterwards, repaired by copying the 9 keys back. Note that after such a repair the
restored keys show up under the dropped models' paths too — one key serves every
path that references it, so that is not a failed drop.


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
