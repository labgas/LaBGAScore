# atlas_mask_tools — masks, ROIs, and the bundled mask library

Two scripts that build masks and ROI sets out of CANlab atlases, plus the
ready-made gray-matter and brain masks that the rest of the pipeline defaults to.

**The naming rule, since it is the first thing anyone needs:** for both mask
families the trailing number is a **threshold**, so a **higher number means a
smaller, sparser, more conservative mask** — never a more liberal one. All masks
within a family are strictly nested. The measured voxel counts are in
[The bundled masks](#the-bundled-masks) below.

---

## Contents

- [What this is](#what-this-is)
- [The scripts](#the-scripts)
  - [`LaBGAScore_atlas_binary_mask_from_atlas.m`](#labgascore_atlas_binary_mask_from_atlasm)
  - [`LaBGAScore_atlas_rois_from_atlas.m`](#labgascore_atlas_rois_from_atlasm)
  - [Which one do I want](#which-one-do-i-want)
- [Who consumes the output](#who-consumes-the-output)
- [The bundled masks](#the-bundled-masks)
  - [Gray-matter masks](#gray-matter-masks)
  - [Brain masks](#brain-masks)
  - [Brain templates](#brain-templates)
- [Choosing a mask](#choosing-a-mask)
- [Dependencies](#dependencies)
- [Notes](#notes)

---

## What this is

`atlas_mask_tools/` holds two **scripts** — no `functions/` subfolder — and three
directories of ready-made image files:

```
LaBGAScore_atlas_binary_mask_from_atlas.m    one binary mask from atlas parcels
LaBGAScore_atlas_rois_from_atlas.m           a SET of ROIs from atlas parcels
gray_matter_masks/                           8 GM masks + the source probseg
brain_masks/                                 3 whole-brain masks
brain_templates/                             1 MNI152 template for display
```

Both scripts are **worked examples rather than generic templates**, and both say
so in their own headers: the directory setup has been genericized, but the atlas
names, label strings and ROI definitions in the body are those of one specific
study (`bit_rew_m1m`, reward regions) and must be rewritten for yours. What
generalizes is the structure — `select_atlas_subset` on one or more atlases, merge,
save in the shapes the second-level scripts expect.

## The scripts

### `LaBGAScore_atlas_binary_mask_from_atlas.m`

Builds **one binary mask** by combining regions from one or more CANlab atlases,
and writes it into the model's `maskdir` as an `fmri_mask_image` plus a `.nii`.

Options control how much of the provenance is kept alongside it:

- `save_original_atlas_obj` — the atlas object *before* the selected parcels are
  merged, i.e. one index per parcel. Keep this if you want to label voxel-wise
  results with the parcel names later.
- `save_merged_atlas_obj` — the atlas object *after* merging, one index for the
  whole mask. Needed to extract ROI averages over the mask as a unit.
- `single_roi` — set true when the mask is one contiguous region.

### `LaBGAScore_atlas_rois_from_atlas.m`

Builds **a set of ROIs**, each one itself a combination of atlas parcels, and
saves them as atlas objects in `maskdir` — optionally also as `.nii` files, and
optionally as a single "flat" object carrying one index per ROI.

The distinction between *original* and *flat* matters and is worth getting right
first time:

| object | what it holds | what it is for |
|---|---|---|
| original | one index per **parcel** that went into each ROI | labelling voxel-wise results with parcel names |
| flat | one index per **ROI** | extracting ROI averages — always saved, because the second-level ROI machinery needs it |

### Which one do I want

- Restricting an analysis to a region — an SVM inside a reward mask, a
  whole-brain GLM masked to gray matter: **binary mask**.
- Reporting one number per region across several regions — an ROI analysis with
  eight a priori targets: **ROI set**.

## Who consumes the output

The ROI script's output is read by `prep_3a_run_second_level_regression_and_save.m`
in the [`CANlab_help_examples`](https://github.com/labgas/CANlab_help_examples)
LaBGAS fork when `doroi_analysis = true`. It loads

```
<maskdir>/<roi_modelname>_rois_<roi_set_name>.mat
```

and **does not create it**, so this script has to run first, with the same
`roi_modelname` and `roi_set_name` as that model's `a2_set_default_options`.

Two traps recorded in the script header, both easy to hit:

- `prep_3a` loads from the `maskdir` of **the model it is running**, whatever
  `roi_modelname` says — `roi_modelname` is only a filename prefix. To reuse an
  ROI set across models, **copy the `.mat` into the other model's `maskdir`**.
- The masks and ROI objects are reached from the second-level options through
  option *strings* (`atlasname_glm`, `roi_names`, `roi_modelname`,
  `roi_set_name`) rather than direct function calls, which is why these two
  scripts do not appear in the generated `dependencies.tsv` even though the
  second level cannot run its ROI analyses without them.

## The bundled masks

Every figure below was **measured from the files in this directory**, not read off
the names. `ml` is in-mask volume in millilitres (nonzero voxels × voxel volume).

### Gray-matter masks

The trailing number is the **gray-matter probability threshold × 100**, applied to
the bundled `tpl-MNI152NLin2009cAsym_res-01_label-GM_probseg.nii`. Verified: the
voxel count of each mask matches thresholding that probseg at the stated value,
exactly, at every threshold.

**1 mm grid** (193×229×193, 1 mm³ voxels) — thresholds of the fMRIPrep template GM probability map:

| mask | GM prob > | voxels | ml | vs the 0.10 mask |
|---|---|---|---|---|
| `gm_mask_fmriprep_20_0_10.nii` | 0.10 | 1,557,143 | 1557.1 | — |
| `gm_mask_fmriprep_20_0_15.nii` | 0.15 | 1,469,809 | 1469.8 | −5.6% |
| `gm_mask_fmriprep_20_0_20.nii` | 0.20 | 1,396,443 | 1396.4 | −10.3% |
| `gm_mask_fmriprep_20_0_25.nii` | 0.25 | 1,330,521 | 1330.5 | −14.6% |
| `gm_mask_fmriprep_20_0_30.nii` | 0.30 | 1,267,998 | 1268.0 | −18.6% |
| `gm_mask_fmriprep_20_0_35.nii` | 0.35 | 1,205,208 | 1205.2 | −22.6% |

**2 mm grid** (97×115×97, 8 mm³ voxels) — the analysis-space variants, named for the
`canlab2023_coarse` atlas:

| mask | GM prob > | voxels | ml |
|---|---|---|---|
| `gm_mask_canlab2023_coarse_fmriprep20_0_20.nii` | 0.20 | 150,630 | 1205.0 |
| `gm_mask_canlab2023_coarse_fmriprep20_0_25.nii` | 0.25 | 138,833 | 1110.7 |

`0_25` is a strict subset of `0_20`. These are the masks second-level analyses
actually run in, because the data are resampled to 2 mm — 150,630 is the voxel
count you will see reported by a whole-brain analysis masked with the `0_20` one.

> **Do not compare `ml` across grids.** `gm_mask_canlab2023_coarse_…_0_20`
> (1205.0 ml) and `gm_mask_fmriprep_20_0_35` (1205.2 ml) have nearly identical
> volumes by coincidence. They are different masks at different thresholds on
> different grids, and the 2 mm one has 8× fewer voxels — which is what matters
> for multiple comparisons, not the volume.

The source probability map, `tpl-MNI152NLin2009cAsym_res-01_label-GM_probseg.nii`,
is also bundled: continuous, 0–1, 2,107,582 nonzero voxels. Threshold it yourself if
you need a value that is not shipped.

### Brain masks

Also thresholds, but of a different quantity, and with a far smaller effect. The
trailing number is a **T1w intensity threshold in the template's own units**
(the fMRIPrep `MNI152NLin2009cAsym` `res-01` T1w runs 0–10000), applied **within
the template brain mask**:

| mask | construction | voxels | ml | removed vs `_0` |
|---|---|---|---|---|
| `brain_mask_fmriprep20_template_0.nii` | the template brain mask itself, no intensity threshold | 1,886,574 | 1886.6 | — |
| `brain_mask_fmriprep20_template_1000.nii` | `_0` ∩ (T1w ≥ 1000) | 1,886,414 | 1886.4 | 160 voxels (0.008%) |
| `brain_mask_fmriprep20_template_1800.nii` | `_0` ∩ (T1w ≥ 1800) | 1,882,234 | 1882.2 | 4,340 voxels (0.23%) |

How this was established, since the names do not say it:

- `_0` is **byte-for-byte identical** to templateflow's
  `tpl-MNI152NLin2009cAsym_res-01_desc-brain_mask` — same 1,886,574 voxels, zero
  differing. So `_0` means "no intensity threshold applied".
- `_1000` is **exactly** `_0 ∩ (T1w ≥ 1000)`; the minimum T1w value inside it is
  precisely 1000.
- `_1800` matches `_0 ∩ (T1w ≥ 1800)` to within **2 voxels out of 1.88 million**
  (its minimum interior T1w value is 1799), so 1800 is plainly the intended
  threshold and the residual is rounding.
- All three are strictly nested: `_1800` ⊂ `_1000` ⊂ `_0`.

**So higher is sparser here too — but the difference is negligible.** Going from
`_0` to `_1800` removes under a quarter of one percent of voxels, at the rim
where the template is darkest. Do not expect a mask choice within this family to
change a result; choose one and be consistent. `_1000` is the LaBGAS default and
is what `pet/scripts/LaBGAScore_pet_a2_set_default_options.m` sets as
`maskname_brain`.

### Brain templates

`mni152_MRIcroGL.nii` — an MNI152 T1 at 0.4 mm³ (207×256×215), continuous, range
0–91.5. A **display underlay**, not a mask: it is not binary and not on any
analysis grid.

## Choosing a mask

- **Gray matter, 2 mm, second level** — `gm_mask_canlab2023_coarse_fmriprep20_0_20`
  is the usual choice, and the one whose 150,630-voxel count appears in this
  repo's own searchlight and TFCE results.
- **Gray matter, 1 mm** — pick by how much partial-volume tissue you are willing
  to carry. 0.20 is the common middle; 0.10 is permissive enough to include a
  good deal of white matter and CSF boundary, 0.35 tight enough to start losing
  cortex in thin gyri.
- **Whole brain** — `brain_mask_fmriprep20_template_1000`, the LaBGAS default.
- **Anything region-specific** — build it with the two scripts above rather than
  hand-editing an existing mask, so the parcel provenance is saved with it.

**Resample deliberately, not incidentally.** Masks here are at 1 mm, 2 mm and
0.4 mm. CANlab's `resample_space` will happily match a mask to your data, but a
mask resampled from 1 mm to 2 mm is not the same mask as the 2 mm one shipped
here, and the voxel count that ends up in your multiple-comparison correction is
whichever one you actually used. The first-level and decoding scripts in this
repo resample the mask to the data and write the resampled copy out, so what was
used is on disk.

## Dependencies

- **CanlabCore** and **Neuroimaging_Pattern_Masks** on the MATLAB path, for
  `atlas`, `select_atlas_subset`, `fmri_mask_image`, `resample_space`. Clone from
  <https://github.com/canlab> if absent. Neither is vendored, and
  Neuroimaging_Pattern_Masks is ~18 GB.
- **SPM12** for image reading/writing.
- Both scripts expect to be run **from the superdataset root** and reach the
  study's directories through the standard `prep_s0_define_directories` setup.

Useful reading: `help atlas.select_atlas_subset` and CANlab's
[using CANlab atlases](https://canlab.github.io/_pages/using_canlab_atlases/using_canlab_atlases.html).

## Notes

- All 13 bundled `.nii` files are ordinary git-tracked files — LaBGAScore is a
  plain git repo, not a DataLad/git-annex dataset — so they are present in every
  clone with no `datalad get` and no broken symlinks. This is the one place in the
  LaBGAS workflow where image files behave that way; everything under a
  `proj_*` superdataset is annexed.
- Neither script has an automated test. The cheap check after editing one is to
  load the mask you produced and confirm its nonzero voxel count and its
  dimensions are what you expected — the same measurements tabulated above,
  which is how the naming conventions in this README were established rather
  than assumed.
