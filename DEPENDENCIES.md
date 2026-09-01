# LaBGAScore — dependency overview

Core scripts and templates for the LaBGAS standard neuroimaging analysis workflow.

> **Generated file — do not edit.** Regenerate with
> `LaBGAScore_dep_report('/data/master_github_repos/LaBGAScore')`
> (see `clean/LaBGAScore_dep_report.m` in LaBGAScore).
> Generated 2026-09-01 by MATLAB 2021a.

This document records which **external** functions each file calls and which
repository those live in. Calls that resolve back into this repository, and
MathWorks' own functions, are omitted.

Direct calls only (depth 1). The transitive closure of a second-level
script runs to thousands of files and reduces to "CanlabCore and SPM";
it is what the provenance tooling uses, for a different purpose.

## How to read this

Resolution is static, and MATLAB makes some of it genuinely undecidable.
Every edge carries a confidence, and the uncertain ones are reported rather
than guessed:

| Confidence | Meaning |
|---|---|
| `resolved` | exactly one definition exists for this name |
| `ambiguous` | several classes define it — e.g. `threshold` exists in `@atlas`, `@glm_map`, `@image_vector` and `@statistic_image`. Deciding which one runs needs type inference these workspace-chained scripts do not support. All candidates are listed. |
| `dotcall` | called as `obj.name(...)`, matched to a class method by name. Could also be a struct field. |
| `dynamic` | the file uses `feval`/`eval`/`str2func`, so its real call set cannot be recovered statically |
| `unparseable` | the file has a syntax error and could not be walked |

## Summary

| Script | Domain | Depends on | Direct calls | Caveats |
|---|---|---|---:|---:|
| `LaBGAScore_atlas_binary_mask_from_atlas` | atlas_mask_tools | CANlab_help_examples, CanlabCore | 6 | 4 |
| `LaBGAScore_atlas_rois_from_atlas` | atlas_mask_tools | CANlab_help_examples, CanlabCore | 10 | 4 |
| `LaBGAScore_check_all_scripts` | clean | — | 0 | 0 |
| `LaBGAScore_check_display` | clean | — | 0 | 0 |
| `LaBGAScore_clean_gzip_all_nii` | clean | — | 0 | 0 |
| `LaBGAScore_clean_sourcedata` | clean | — | 0 | 0 |
| `LaBGAScore_dep_build_index` | clean | — | 0 | 0 |
| `LaBGAScore_dep_map` | clean | — | 0 | 0 |
| `LaBGAScore_dep_report` | clean | — | 0 | 0 |
| `LaBGAScore_move_repos_matlabpath` | clean | — | 0 | 0 |
| `LaBGAScore_prov_gitinfo` | clean | — | 0 | 0 |
| `LaBGAScore_prov_gitstatus` | clean | — | 0 | 0 |
| `LaBGAScore_prov_protect_reflogs` | clean | — | 0 | 0 |
| `LaBGAScore_prov_publish` | clean | — | 0 | 1 |
| `LaBGAScore_prov_resolve_retrospective` | clean | — | 0 | 0 |
| `LaBGAScore_prov_snapshot` | clean | — | 0 | 1 |
| `LaBGAScore_smart_parallel_pool_setup` | clean | — | 0 | 0 |
| `LaBGAScore_cosmomvpa_searchlight_rsa` | cosmomvpa | CANlab_help_examples, CanlabCore, CoSMoMVPA | 25 | 1 |
| `LaBGAScore_decoding_SVM_between_subjects` | decoding_toolbox | CANlab_help_examples, CanlabCore, CoSMoMVPA, spm12 | 10 | 1 |
| `LaBGAScore_decoding_template_xclass_acc` | decoding_toolbox | — | 0 | 1 |
| `canlabCmap` | figures | — | 0 | 0 |
| `cluster_surface_plots` | figures | spm12 | 2 | 0 |
| `save_all_open_figures_smart` | figures | — | 0 | 0 |
| `LaBGAScore_firstlevel_s1_options_dsgn_struct` | firstlevel | — | 0 | 0 |
| `LaBGAScore_firstlevel_s1a_options_dsgn_multisess_multitask` | firstlevel | — | 0 | 0 |
| `LaBGAScore_firstlevel_s1b_fit_phMRI_model` | firstlevel | spm12 | 4 | 1 |
| `LaBGAScore_firstlevel_s2_fit_model` | firstlevel | CanlabCore, CanlabPrivate, spm12 | 13 | 42 |
| `LaBGAScore_firstlevel_s2a_fit_model_multisess_multitask` | firstlevel | CanlabCore, CanlabPrivate, spm12 | 13 | 42 |
| `LaBGAScore_firstlevel_s2b_phMRI_covar_contrasts` | firstlevel | spm12 | 2 | 2 |
| `LaBGAScore_firstlevel_s3_diagnose_model` | firstlevel | CanlabCore | 11 | 11 |
| `LaBGAScore_firstlevel_s3b_phMRI_signature_responses` | firstlevel | CanlabCore | 3 | 1 |
| `LaBGAScore_firstlevel_s3c_phMRI_neurotransmitter_maps` | firstlevel | CanlabCore, canlab_single_trials | 5 | 3 |
| `canlab_glm_subject_levels_old` | firstlevel | CanlabCore, spm12 | 5 | 2 |
| `canlab_glm_subject_levels_run1subject_old` | firstlevel | CanlabCore, spm12 | 10 | 11 |
| `LaBGAScore_prep_graphvar_input_from_conn` | graphvar | CanlabCore | 3 | 2 |
| `LaBGAScore_juspace_corr_behav` | juspace | CANlab_help_examples | 1 | 0 |
| `LaBGAScore_prep_juspace_input` | juspace | CANlab_help_examples, CanlabCore | 8 | 4 |
| `LaBGAScore_mrs_osprey_jobfile_GE` | mrs | — | 0 | 0 |
| `LaBGAScore_mrs_osprey_jobfile_Philips` | mrs | — | 0 | 0 |
| `LaBGAScore_mrs_osprey_single_sess_jobfile_GE` | mrs | — | 0 | 0 |
| `LaBGAScore_mrs_osprey_single_sess_jobfile_Philips` | mrs | — | 0 | 0 |
| `LaBGAScore_mrs_run_osprey_GE` | mrs | osprey | 8 | 1 |
| `LaBGAScore_mrs_run_osprey_Philips` | mrs | osprey | 8 | 1 |
| `LaBGAScore_prep_mrs2bids` | mrs | — | 0 | 0 |
| `LCN12_PET_preprocessing` | pet | spm12 | 5 | 5 |
| `LCN12_analyze_headmovement_PET` | pet | — | 0 | 0 |
| `LCN12_generate_qc_from4D` | pet | — | 0 | 0 |
| `LCN12_read_image` | pet | spm12 | 5 | 2 |
| `LCN12_realign2_dy_PET` | pet | spm12 | 5 | 6 |
| `LCN12_realign_dy_PET` | pet | spm12 | 5 | 6 |
| `LCN12_smooth` | pet | spm12 | 1 | 0 |
| `LCN12_write_image` | pet | spm12 | 2 | 2 |
| `LCN_3Dimage_dilate` | pet | — | 0 | 0 |
| `LCN_LOGAN` | pet | — | 0 | 0 |
| `LCN_calc2_model_2T4k` | pet | — | 0 | 0 |
| `LCN_calc2_model_2T4k_Vb` | pet | — | 0 | 0 |
| `LCN_calc_intact_tracer_hill` | pet | — | 0 | 0 |
| `LCN_calc_model_2T4k_vasc1k` | pet | — | 0 | 0 |
| `LCN_calc_model_2T4k_vasc2k` | pet | — | 0 | 0 |
| `LCN_calc_model_selection` | pet | — | 0 | 0 |
| `LCN_check_filename` | pet | — | 0 | 0 |
| `LCN_cost_intact_tracer_hill` | pet | — | 0 | 0 |
| `LCN_string` | pet | — | 0 | 0 |
| `LaBGAScore_pet_a2_set_default_options` | pet | — | 0 | 0 |
| `LaBGAScore_pet_a_set_up_paths_always_run_first` | pet | — | 0 | 0 |
| `LaBGAScore_pet_dcm2bids` | pet | spm12 | 1 | 0 |
| `LaBGAScore_pet_model_TSPO_DPA714` | pet | CanlabCore, spm12 | 5 | 3 |
| `LaBGAScore_pet_prep_1_dat_behavioral_data` | pet | — | 0 | 0 |
| `LaBGAScore_pet_prep_2_load_image_data_and_save` | pet | CANlab_help_examples, CanlabCore | 8 | 0 |
| `LaBGAScore_pet_prep_4_apply_signatures_and_save` | pet | CanlabCore, MasksPrivate, Neuroimaging_Pattern_Masks | 3 | 0 |
| `LaBGAScore_pet_preprocess_data` | pet | — | 0 | 0 |
| `LaBGAScore_pet_run_plot_PLS_ENet_pipeline` | pet | CanlabCore | 6 | 2 |
| `holmthreshold` | power | — | 0 | 0 |
| `medianIQR_to_meanSD` | power | — | 0 | 0 |
| `print_LEAR_matrix` | power | — | 0 | 0 |
| `sd_from_ci` | power | — | 0 | 0 |
| `LaBGAScore_prep_parrec2bids` | prep | spm12 | 2 | 0 |
| `LaBGAScore_prep_s0_define_directories` | prep | — | 0 | 0 |
| `LaBGAScore_prep_s1_write_events_tsv` | prep | — | 0 | 1 |
| `LaBGAScore_prep_s1_write_events_tsv_multisess_multitask` | prep | — | 0 | 1 |
| `LaBGAScore_prep_s2_smooth` | prep | spm12 | 2 | 1 |
| `LaBGAScore_prep_s2_smooth_multisess` | prep | spm12 | 2 | 1 |
| `ENet_neuroimaging_pipeline` | secondlevel | CanlabCore | 1 | 0 |
| `LaBGAScore_secondlevel_MS_mat_pipeline` | secondlevel | spm12 | 1 | 0 |
| `LaBGAScore_secondlevel_extractparcels_sessions` | secondlevel | CANlab_help_examples, CanlabCore | 4 | 1 |
| `LaBGAScore_secondlevel_mvpa_beta_maps_conn` | secondlevel | CanlabCore, ooFmriDataObjML | 18 | 16 |
| `LaBGAScore_secondlevel_ooFmriDataObjML_example` | secondlevel | CanlabCore, ooFmriDataObjML | 10 | 2 |
| `LaBGAScore_secondlevel_roi_run_plot_PLS_ENet_pipeline` | secondlevel | CANlab_help_examples, CanlabCore | 4 | 1 |
| `PLSDA_neuroimaging_pipeline` | secondlevel | CanlabCore | 1 | 0 |
| `PLSDA_paired_neuroimaging_pipeline` | secondlevel | — | 0 | 0 |
| `PLSR_neuroimaging_pipeline` | secondlevel | CanlabCore | 1 | 0 |
| `ProgressTracker` | secondlevel | — | 0 | 0 |
| `applyScaling` | secondlevel | — | 0 | 0 |
| `bootstrapOOB_ENet` | secondlevel | CanlabCore | 1 | 0 |
| `bootstrapOOB_PLSDA` | secondlevel | CanlabCore | 1 | 0 |
| `bootstrapOOB_PLSR` | secondlevel | CanlabCore | 1 | 0 |
| `capLV` | secondlevel | — | 0 | 0 |
| `dice_statistic_image` | secondlevel | — | 0 | 0 |
| `dice_statistic_image_by_roi` | secondlevel | — | 0 | 0 |
| `enetLambdaGrid` | secondlevel | — | 0 | 0 |
| `foldPreprocess` | secondlevel | — | 0 | 0 |
| `globalBaselineCV` | secondlevel | CanlabCore | 1 | 0 |
| `group_tfce_from_subject_maps` | secondlevel | CanlabCore | 2 | 0 |
| `logitSafe` | secondlevel | — | 0 | 0 |
| `makeGroupedFolds` | secondlevel | — | 0 | 0 |
| `maskToSignificant` | secondlevel | — | 0 | 0 |
| `plot_ENet_diagnostics_neuroimaging` | secondlevel | spm12 | 3 | 0 |
| `plot_PLSDA_diagnostics_neuroimaging` | secondlevel | spm12 | 3 | 0 |
| `plot_PLSR_diagnostics_neuroimaging` | secondlevel | spm12 | 3 | 0 |
| `quickCV_ENet` | secondlevel | CanlabCore | 1 | 0 |
| `quickCV_ENet_PR` | secondlevel | CanlabCore | 1 | 0 |
| `quickCV_PLSDA` | secondlevel | — | 0 | 0 |
| `quickCV_PLSDA_PR` | secondlevel | — | 0 | 0 |
| `quickCV_PLSDA_core` | secondlevel | CanlabCore | 1 | 0 |
| `quickCV_PLSR` | secondlevel | CanlabCore | 1 | 0 |
| `quickGroupedCV` | secondlevel | — | 0 | 0 |
| `residualizeFold` | secondlevel | — | 0 | 0 |
| `residualizeY` | secondlevel | — | 0 | 0 |
| `selectENetHyperparams` | secondlevel | CanlabCore | 1 | 0 |
| `setParforStream` | secondlevel | — | 0 | 0 |
| `swapWithinSubjectLabels` | secondlevel | — | 0 | 0 |
| `tfce_one_fmri_dat` | secondlevel | — | 0 | 0 |
| `tfce_transform_3d` | secondlevel | — | 0 | 0 |
| `tfce_volume` | secondlevel | — | 0 | 0 |
| `thresholded_fmri_data_from_pval_nii` | secondlevel | CanlabCore, canlab_single_trials | 6 | 5 |
| `thresholded_fmri_data_from_statistic_image` | secondlevel | CanlabCore, canlab_single_trials | 2 | 0 |
| `validateAtlasLabels` | secondlevel | — | 0 | 0 |
| `validateCovariates` | secondlevel | — | 0 | 0 |
| `warnUnknownOptions` | secondlevel | — | 0 | 0 |
| `LaBGAScore_Storey_FDR` | stats_tools | — | 0 | 0 |

## Dependencies by repository

| Repository | Call edges | Distinct functions |
|---|---:|---:|
| CanlabCore | 635 | 53 |
| spm12 | 179 | 29 |
| ooFmriDataObjML | 18 | 11 |
| CoSMoMVPA | 17 | 17 |
| osprey | 16 | 8 |
| CANlab_help_examples | 9 | 2 |
| CanlabPrivate | 6 | 2 |
| canlab_single_trials | 4 | 2 |
| MasksPrivate | 1 | 1 |
| Neuroimaging_Pattern_Masks | 1 | 1 |

## Per-script detail

### `LaBGAScore_atlas_binary_mask_from_atlas`

`atlas_mask_tools/LaBGAScore_atlas_binary_mask_from_atlas.m`

**CANlab_help_examples**

- `a_set_up_paths_always_run_first`

**CanlabCore**

- `fmri_mask_image` *(@fmri_mask_image)*
- `load_atlas`
- `merge_atlases` *(@atlas)*
- `select_atlas_subset` *(@atlas)*
- `threshold` *(@atlas)* — `dotcall`, 4 candidates

### `LaBGAScore_atlas_rois_from_atlas`

`atlas_mask_tools/LaBGAScore_atlas_rois_from_atlas.m`

**CANlab_help_examples**

- `a_set_up_paths_always_run_first`

**CanlabCore**

- `addblobs` *(@fmridisplay)*
- `atlas2region` *(@atlas)*
- `canlab_results_fmridisplay`
- `load_atlas`
- `merge_atlases` *(@atlas)*
- `scn_standard_colors`
- `select_atlas_subset` *(@atlas)*
- `threshold` *(@atlas)* — `dotcall`, 4 candidates
- `title_montage` *(@fmridisplay)*

### `LaBGAScore_check_all_scripts`

`clean/LaBGAScore_check_all_scripts.m`

No external dependencies.

### `LaBGAScore_check_display`

`clean/LaBGAScore_check_display.m`

No external dependencies.

### `LaBGAScore_clean_gzip_all_nii`

`clean/LaBGAScore_clean_gzip_all_nii.m`

No external dependencies.

### `LaBGAScore_clean_sourcedata`

`clean/LaBGAScore_clean_sourcedata.m`

No external dependencies.

### `LaBGAScore_dep_build_index`

`clean/LaBGAScore_dep_build_index.m`

No external dependencies.

### `LaBGAScore_dep_map`

`clean/LaBGAScore_dep_map.m`

No external dependencies.

### `LaBGAScore_dep_report`

`clean/LaBGAScore_dep_report.m`

No external dependencies.

### `LaBGAScore_move_repos_matlabpath`

`clean/LaBGAScore_move_repos_matlabpath.m`

No external dependencies.

### `LaBGAScore_prov_gitinfo`

`clean/LaBGAScore_prov_gitinfo.m`

No external dependencies.

### `LaBGAScore_prov_gitstatus`

`clean/LaBGAScore_prov_gitstatus.m`

No external dependencies.

### `LaBGAScore_prov_protect_reflogs`

`clean/LaBGAScore_prov_protect_reflogs.m`

No external dependencies.

### `LaBGAScore_prov_publish`

`clean/LaBGAScore_prov_publish.m`

No external dependencies.

### `LaBGAScore_prov_resolve_retrospective`

`clean/LaBGAScore_prov_resolve_retrospective.m`

No external dependencies.

### `LaBGAScore_prov_snapshot`

`clean/LaBGAScore_prov_snapshot.m`

No external dependencies.

### `LaBGAScore_smart_parallel_pool_setup`

`clean/LaBGAScore_smart_parallel_pool_setup.m`

No external dependencies.

### `LaBGAScore_cosmomvpa_searchlight_rsa`

`cosmomvpa/LaBGAScore_cosmomvpa_searchlight_rsa.m`

**CANlab_help_examples**

- `a_set_up_paths_always_run_first`

**CanlabCore**

- `canlab_results_fmridisplay`
- `fmri_data` *(@fmri_data)*
- `fmri_mask_image` *(@fmri_mask_image)*
- `load_atlas`
- `montage` *(@region)* — `ambiguous`
- `region` *(@region)*
- `resample_space` *(@image_vector)* — `ambiguous_within_repo`, 2 candidates
- `threshold` *(@atlas)* — `ambiguous_within_repo`, 4 candidates
- `title_montage` *(@fmridisplay)*

**CoSMoMVPA**

- `cosmo_cluster_neighborhood`
- `cosmo_disp`
- `cosmo_fmri_dataset`
- `cosmo_fx`
- `cosmo_map2fmri`
- `cosmo_montecarlo_cluster_stat`
- `cosmo_pdist`
- `cosmo_randomize_targets`
- `cosmo_remove_useless_data`
- `cosmo_searchlight`
- `cosmo_slice`
- `cosmo_spherical_neighborhood`
- `cosmo_squareform`
- `cosmo_stack`
- `cosmo_target_dsm_corr_measure`

### `LaBGAScore_decoding_SVM_between_subjects`

`decoding_toolbox/LaBGAScore_decoding_SVM_between_subjects.m`

**CANlab_help_examples**

- `a_set_up_paths_always_run_first`

**CanlabCore**

- `fmri_mask_image` *(@fmri_mask_image)*
- `load_atlas`
- `train` *(@algorithm)* — `dotcall`

**CoSMoMVPA**

- `load_nii`
- `save_nii`

**spm12**

- `decoding`
- `decoding_defaults`
- `make_design_cv`
- `make_design_permutation`

### `LaBGAScore_decoding_template_xclass_acc`

`decoding_toolbox/LaBGAScore_decoding_template_xclass_acc.m`

No external dependencies.

### `canlabCmap`

`figures/canlabCmap.m`

No external dependencies.

### `cluster_surface_plots`

`figures/cluster_surface_plots.m`

**spm12**

- `spm_read_vols` — `ambiguous_within_repo`, 2 candidates
- `spm_vol` — `ambiguous_within_repo`, 2 candidates

### `save_all_open_figures_smart`

`figures/save_all_open_figures_smart.m`

No external dependencies.

### `LaBGAScore_firstlevel_s1_options_dsgn_struct`

`firstlevel/LaBGAScore_firstlevel_s1_options_dsgn_struct.m`

No external dependencies.

### `LaBGAScore_firstlevel_s1a_options_dsgn_multisess_multitask`

`firstlevel/LaBGAScore_firstlevel_s1a_options_dsgn_multisess_multitask.m`

No external dependencies.

### `LaBGAScore_firstlevel_s1b_fit_phMRI_model`

`firstlevel/LaBGAScore_firstlevel_s1b_fit_phMRI_model.m`

**spm12**

- `cfg_dep` *(@cfg_dep)*
- `spm` — `ambiguous_within_repo`, 2 candidates
- `spm_jobman`
- `spm_select` — `ambiguous_within_repo`, 2 candidates

### `LaBGAScore_firstlevel_s2_fit_model`

`firstlevel/LaBGAScore_firstlevel_s2_fit_model.m`

**CanlabCore**

- `display` *(@inplace)* — `dotcall`, 11 candidates
- `onsets2fmridesign`
- `param` *(@param)* — `dotcall`
- `plotDesign`
- `plot_matrix_cols`
- `read_hdr`
- `scn_session_spike_id`

**CanlabPrivate**

- `canlab_spm_contrast_job_single_trials_lukasvo`
- `display` *(@InteractiveViewer)* — `dotcall`, 2 candidates

**spm12**

- `conditions` *(@meeg)* — `dotcall`
- `display` *(@file_array)* — `dotcall`, 19 candidates
- `spm_hrf`
- `time` *(@meeg)* — `dotcall`
- `tr` *(@mardo_2)* — `dotcall`, 2 candidates
- `type` *(@meeg)* — `dotcall`, 5 candidates

### `LaBGAScore_firstlevel_s2a_fit_model_multisess_multitask`

`firstlevel/LaBGAScore_firstlevel_s2a_fit_model_multisess_multitask.m`

**CanlabCore**

- `display` *(@inplace)* — `dotcall`, 11 candidates
- `onsets2fmridesign`
- `param` *(@param)* — `dotcall`
- `plotDesign`
- `plot_matrix_cols`
- `read_hdr`
- `scn_session_spike_id`

**CanlabPrivate**

- `canlab_spm_contrast_job_single_trials_lukasvo`
- `display` *(@InteractiveViewer)* — `dotcall`, 2 candidates

**spm12**

- `conditions` *(@meeg)* — `dotcall`
- `display` *(@file_array)* — `dotcall`, 19 candidates
- `spm_hrf`
- `time` *(@meeg)* — `dotcall`
- `tr` *(@mardo_2)* — `dotcall`, 2 candidates
- `type` *(@meeg)* — `dotcall`, 5 candidates

### `LaBGAScore_firstlevel_s2b_phMRI_covar_contrasts`

`firstlevel/LaBGAScore_firstlevel_s2b_phMRI_covar_contrasts.m`

**spm12**

- `contrasts` *(@mardo)* — `dotcall`
- `spm_jobman`

### `LaBGAScore_firstlevel_s3_diagnose_model`

`firstlevel/LaBGAScore_firstlevel_s3_diagnose_model.m`

**CanlabCore**

- `activate_figures` *(@fmridisplay)*
- `addblobs` *(@fmridisplay)*
- `canlab_results_fmridisplay`
- `create_figure`
- `display` *(@inplace)* — `dotcall`, 11 candidates
- `fmri_mask_image` *(@fmri_mask_image)*
- `resample_space` *(@image_vector)* — `ambiguous_within_repo`, 2 candidates
- `scn_spm_design_check`
- `statistic_image` *(@statistic_image)*
- `threshold` *(@atlas)* — `ambiguous_within_repo`, 4 candidates
- `title_montage` *(@fmridisplay)*

### `LaBGAScore_firstlevel_s3b_phMRI_signature_responses`

`firstlevel/LaBGAScore_firstlevel_s3b_phMRI_signature_responses.m`

**CanlabCore**

- `apply_mask` *(@image_vector)*
- `fmri_data` *(@fmri_data)*
- `load_image_set`

### `LaBGAScore_firstlevel_s3c_phMRI_neurotransmitter_maps`

`firstlevel/LaBGAScore_firstlevel_s3c_phMRI_neurotransmitter_maps.m`

**CanlabCore**

- `fmri_data` *(@fmri_data)*
- `get_wh_image` *(@image_vector)* — `dotcall`
- `image_similarity_plot` *(@image_vector)*
- `load_image_set`

**canlab_single_trials**

- `fmri_data_st` *(@fmri_data_st)*
- `get_wh_image` *(@fmri_data_st)* — `dotcall`

### `canlab_glm_subject_levels_old`

`firstlevel/functions/canlab_glm_subject_levels_old.m`

**CanlabCore**

- `batch_t_histograms`
- `canlab_glm_publish`
- `filenames`

**spm12**

- `conditions` *(@meeg)* — `dotcall`
- `spm_jobman`

### `canlab_glm_subject_levels_run1subject_old`

`firstlevel/functions/canlab_glm_subject_levels_run1subject_old.m`

**CanlabCore**

- `expand_4d_filenames`
- `filenames`

**spm12**

- `conditions` *(@meeg)* — `dotcall`
- `contrasts` *(@mardo)* — `dotcall`
- `nifti` *(@nifti)* — `ambiguous_within_repo`, 2 candidates
- `spm_get_defaults` — `ambiguous_within_repo`, 2 candidates
- `spm_jobman`
- `time` *(@meeg)* — `dotcall`
- `tr` *(@mardo_2)* — `dotcall`, 2 candidates
- `type` *(@meeg)* — `dotcall`, 5 candidates

### `LaBGAScore_prep_graphvar_input_from_conn`

`graphvar/LaBGAScore_prep_graphvar_input_from_conn.m`

**CanlabCore**

- `data` *(@data)* — `ambiguous`
- `names`
- `size` *(@inplace)* — `dotcall`

### `LaBGAScore_juspace_corr_behav`

`juspace/LaBGAScore_juspace_corr_behav.m`

**CANlab_help_examples**

- `a_set_up_paths_always_run_first`

### `LaBGAScore_prep_juspace_input`

`juspace/LaBGAScore_prep_juspace_input.m`

**CANlab_help_examples**

- `a_set_up_paths_always_run_first`

**CanlabCore**

- `apply_mask` *(@image_vector)*
- `fmri_data` *(@fmri_data)*
- `fmri_mask_image` *(@fmri_mask_image)*
- `get_wh_image` *(@image_vector)*
- `load_atlas`
- `resample_space` *(@image_vector)* — `ambiguous_within_repo`, 2 candidates
- `threshold` *(@atlas)* — `dotcall`, 4 candidates

### `LaBGAScore_mrs_osprey_jobfile_GE`

`mrs/LaBGAScore_mrs_osprey_jobfile_GE.m`

No external dependencies.

### `LaBGAScore_mrs_osprey_jobfile_Philips`

`mrs/LaBGAScore_mrs_osprey_jobfile_Philips.m`

No external dependencies.

### `LaBGAScore_mrs_osprey_single_sess_jobfile_GE`

`mrs/LaBGAScore_mrs_osprey_single_sess_jobfile_GE.m`

No external dependencies.

### `LaBGAScore_mrs_osprey_single_sess_jobfile_Philips`

`mrs/LaBGAScore_mrs_osprey_single_sess_jobfile_Philips.m`

No external dependencies.

### `LaBGAScore_mrs_run_osprey_GE`

`mrs/LaBGAScore_mrs_run_osprey_GE.m`

**osprey**

- `OspreyCoreg`
- `OspreyFit`
- `OspreyJob`
- `OspreyLoad`
- `OspreyOverview`
- `OspreyProcess`
- `OspreyQuantify`
- `OspreySeg`

### `LaBGAScore_mrs_run_osprey_Philips`

`mrs/LaBGAScore_mrs_run_osprey_Philips.m`

**osprey**

- `OspreyCoreg`
- `OspreyFit`
- `OspreyJob`
- `OspreyLoad`
- `OspreyOverview`
- `OspreyProcess`
- `OspreyQuantify`
- `OspreySeg`

### `LaBGAScore_prep_mrs2bids`

`mrs/LaBGAScore_prep_mrs2bids.m`

No external dependencies.

### `LCN12_PET_preprocessing`

`pet/functions/LCN12_PET_preprocessing.m`

**spm12**

- `dim` *(@file_array)* — `dotcall`, 2 candidates
- `estimate` *(@mardo_2)* — `dotcall`, 3 candidates
- `spm_jobman`
- `spm_vol` — `ambiguous_within_repo`, 2 candidates
- `spm_write_vol` — `ambiguous_within_repo`, 3 candidates

### `LCN12_analyze_headmovement_PET`

`pet/functions/LCN12_analyze_headmovement_PET.m`

No external dependencies.

### `LCN12_generate_qc_from4D`

`pet/functions/LCN12_generate_qc_from4D.m`

No external dependencies.

### `LCN12_read_image`

`pet/functions/LCN12_read_image.m`

**spm12**

- `dim` *(@file_array)* — `dotcall`, 2 candidates
- `spm_matrix` — `ambiguous_within_repo`, 2 candidates
- `spm_read_vols` — `ambiguous_within_repo`, 2 candidates
- `spm_slice_vol` — `ambiguous_within_repo`, 5 candidates
- `spm_vol` — `ambiguous_within_repo`, 2 candidates

### `LCN12_realign2_dy_PET`

`pet/functions/LCN12_realign2_dy_PET.m`

**spm12**

- `dim` *(@file_array)* — `dotcall`, 2 candidates
- `estimate` *(@mardo_2)* — `dotcall`, 3 candidates
- `spm_imatrix`
- `spm_jobman`
- `spm_vol` — `ambiguous_within_repo`, 2 candidates

### `LCN12_realign_dy_PET`

`pet/functions/LCN12_realign_dy_PET.m`

**spm12**

- `dim` *(@file_array)* — `dotcall`, 2 candidates
- `estimate` *(@mardo_2)* — `dotcall`, 3 candidates
- `spm_imatrix`
- `spm_jobman`
- `spm_vol` — `ambiguous_within_repo`, 2 candidates

### `LCN12_smooth`

`pet/functions/LCN12_smooth.m`

**spm12**

- `spm_smooth`

### `LCN12_write_image`

`pet/functions/LCN12_write_image.m`

**spm12**

- `dim` *(@file_array)* — `dotcall`, 2 candidates
- `spm_write_vol` — `ambiguous_within_repo`, 3 candidates

### `LCN_3Dimage_dilate`

`pet/functions/LCN_3Dimage_dilate.m`

No external dependencies.

### `LCN_LOGAN`

`pet/functions/LCN_LOGAN.m`

No external dependencies.

### `LCN_calc2_model_2T4k`

`pet/functions/LCN_calc2_model_2T4k.m`

No external dependencies.

### `LCN_calc2_model_2T4k_Vb`

`pet/functions/LCN_calc2_model_2T4k_Vb.m`

No external dependencies.

### `LCN_calc_intact_tracer_hill`

`pet/functions/LCN_calc_intact_tracer_hill.m`

No external dependencies.

### `LCN_calc_model_2T4k_vasc1k`

`pet/functions/LCN_calc_model_2T4k_vasc1k.m`

No external dependencies.

### `LCN_calc_model_2T4k_vasc2k`

`pet/functions/LCN_calc_model_2T4k_vasc2k.m`

No external dependencies.

### `LCN_calc_model_selection`

`pet/functions/LCN_calc_model_selection.m`

No external dependencies.

### `LCN_check_filename`

`pet/functions/LCN_check_filename.m`

No external dependencies.

### `LCN_cost_intact_tracer_hill`

`pet/functions/LCN_cost_intact_tracer_hill.m`

No external dependencies.

### `LCN_string`

`pet/functions/LCN_string.m`

No external dependencies.

### `LaBGAScore_pet_a2_set_default_options`

`pet/scripts/LaBGAScore_pet_a2_set_default_options.m`

No external dependencies.

### `LaBGAScore_pet_a_set_up_paths_always_run_first`

`pet/scripts/LaBGAScore_pet_a_set_up_paths_always_run_first.m`

No external dependencies.

### `LaBGAScore_pet_dcm2bids`

`pet/scripts/LaBGAScore_pet_dcm2bids.m`

**spm12**

- `spm_jsonread` — `ambiguous_within_repo`, 3 candidates

### `LaBGAScore_pet_model_TSPO_DPA714`

`pet/scripts/LaBGAScore_pet_model_TSPO_DPA714.m`

**CanlabCore**

- `downsample_parcellation` *(@atlas)*
- `load_atlas`

**spm12**

- `dim` *(@file_array)* — `dotcall`, 2 candidates
- `spm_imatrix`
- `spm_vol` — `ambiguous_within_repo`, 2 candidates

### `LaBGAScore_pet_prep_1_dat_behavioral_data`

`pet/scripts/LaBGAScore_pet_prep_1_dat_behavioral_data.m`

No external dependencies.

### `LaBGAScore_pet_prep_2_load_image_data_and_save`

`pet/scripts/LaBGAScore_pet_prep_2_load_image_data_and_save.m`

**CANlab_help_examples**

- `plugin_get_options_for_analysis_script`

**CanlabCore**

- `check_valid_imagename`
- `create_figure`
- `enforce_variable_types` *(@image_vector)*
- `fmri_data` *(@fmri_data)*
- `fmri_mask_image` *(@fmri_mask_image)*
- `qc_metrics_second_level` *(@image_vector)*
- `remove_empty` *(@image_vector)*

### `LaBGAScore_pet_prep_4_apply_signatures_and_save`

`pet/scripts/LaBGAScore_pet_prep_4_apply_signatures_and_save.m`

**CanlabCore**

- `ste`

**MasksPrivate**

- `apply_nps`

**Neuroimaging_Pattern_Masks**

- `apply_all_signatures`

### `LaBGAScore_pet_preprocess_data`

`pet/scripts/LaBGAScore_pet_preprocess_data.m`

No external dependencies.

### `LaBGAScore_pet_run_plot_PLS_ENet_pipeline`

`pet/scripts/LaBGAScore_pet_run_plot_PLS_ENet_pipeline.m`

**CanlabCore**

- `downsample_parcellation` *(@atlas)*
- `fmri_data` *(@fmri_data)*
- `load_atlas`
- `remove_atlas_region` *(@atlas)* — `dotcall`
- `resample_space` *(@image_vector)* — `ambiguous_within_repo`, 2 candidates
- `write` *(@image_vector)* — `dotcall`

### `holmthreshold`

`power/holmthreshold.m`

No external dependencies.

### `medianIQR_to_meanSD`

`power/medianIQR_to_meanSD.m`

No external dependencies.

### `print_LEAR_matrix`

`power/print_LEAR_matrix.m`

No external dependencies.

### `sd_from_ci`

`power/sd_from_ci.m`

No external dependencies.

### `LaBGAScore_prep_parrec2bids`

`prep/LaBGAScore_prep_parrec2bids.m`

**spm12**

- `spm_jsonread` — `ambiguous_within_repo`, 3 candidates
- `spm_jsonwrite`

### `LaBGAScore_prep_s0_define_directories`

`prep/LaBGAScore_prep_s0_define_directories.m`

No external dependencies.

### `LaBGAScore_prep_s1_write_events_tsv`

`prep/LaBGAScore_prep_s1_write_events_tsv.m`

No external dependencies.

### `LaBGAScore_prep_s1_write_events_tsv_multisess_multitask`

`prep/LaBGAScore_prep_s1_write_events_tsv_multisess_multitask.m`

No external dependencies.

### `LaBGAScore_prep_s2_smooth`

`prep/LaBGAScore_prep_s2_smooth.m`

**spm12**

- `spm_jobman`
- `spm_select` — `ambiguous_within_repo`, 2 candidates

### `LaBGAScore_prep_s2_smooth_multisess`

`prep/LaBGAScore_prep_s2_smooth_multisess.m`

**spm12**

- `spm_jobman`
- `spm_select` — `ambiguous_within_repo`, 2 candidates

### `ENet_neuroimaging_pipeline`

`secondlevel/functions/ENet_neuroimaging_pipeline.m`

**CanlabCore**

- `training` *(@chain)* — `ambiguous_within_repo`, 38 candidates

### `LaBGAScore_secondlevel_MS_mat_pipeline`

`secondlevel/scripts/LaBGAScore_secondlevel_MS_mat_pipeline.m`

**spm12**

- `MA_model_space`

### `LaBGAScore_secondlevel_extractparcels_sessions`

`secondlevel/scripts/LaBGAScore_secondlevel_extractparcels_sessions.m`

**CANlab_help_examples**

- `a_set_up_paths_always_run_first`

**CanlabCore**

- `downsample_parcellation` *(@atlas)*
- `group` *(@group)* — `dotcall`
- `load_atlas`

### `LaBGAScore_secondlevel_mvpa_beta_maps_conn`

`secondlevel/scripts/LaBGAScore_secondlevel_mvpa_beta_maps_conn.m`

**CanlabCore**

- `addblobs` *(@fmridisplay)*
- `canlab_results_fmridisplay`
- `fit` *(@glm_map)* — `dotcall`, 2 candidates
- `fmri_data` *(@fmri_data)*
- `fmri_mask_image` *(@fmri_mask_image)*
- `pipeline` *(@pipeline)* — `ambiguous`
- `plot` *(@fmri_data)* — `dotcall`, 11 candidates
- `predict` *(@fmri_data)*
- `region` *(@region)*
- `test` *(@algorithm)* — `dotcall`
- `title_montage` *(@fmridisplay)*

**ooFmriDataObjML**

- `bayesOptCV`
- `crossValScore`
- `cvpartition2`
- `fmri2VxlFeatTransformer`
- `get_mse`
- `pcrRegressor`
- `pipeline` — `ambiguous`
- `plsRegressor`

### `LaBGAScore_secondlevel_ooFmriDataObjML_example`

`secondlevel/scripts/LaBGAScore_secondlevel_ooFmriDataObjML_example.m`

**CanlabCore**

- `pipeline` *(@pipeline)* — `ambiguous`

**ooFmriDataObjML**

- `bayesOptCV`
- `crossValScore`
- `cvpartition2`
- `fmri2VxlFeatTransformer`
- `functionTransformer`
- `get_mse`
- `gridSearchCV`
- `pipeline` — `ambiguous`
- `plsRegressor`
- `zscoreVxlTransformer`

### `LaBGAScore_secondlevel_roi_run_plot_PLS_ENet_pipeline`

`secondlevel/scripts/LaBGAScore_secondlevel_roi_run_plot_PLS_ENet_pipeline.m`

**CANlab_help_examples**

- `a_set_up_paths_always_run_first`

**CanlabCore**

- `fmri_data` *(@fmri_data)*
- `resample_space` *(@image_vector)* — `ambiguous_within_repo`, 2 candidates
- `write` *(@image_vector)* — `dotcall`

### `PLSDA_neuroimaging_pipeline`

`secondlevel/functions/PLSDA_neuroimaging_pipeline.m`

**CanlabCore**

- `training` *(@chain)* — `ambiguous_within_repo`, 38 candidates

### `PLSDA_paired_neuroimaging_pipeline`

`secondlevel/functions/PLSDA_paired_neuroimaging_pipeline.m`

No external dependencies.

### `PLSR_neuroimaging_pipeline`

`secondlevel/functions/PLSR_neuroimaging_pipeline.m`

**CanlabCore**

- `training` *(@chain)* — `ambiguous_within_repo`, 38 candidates

### `ProgressTracker`

`secondlevel/classes/ProgressTracker.m`

No external dependencies.

### `applyScaling`

`secondlevel/functions/applyScaling.m`

No external dependencies.

### `bootstrapOOB_ENet`

`secondlevel/functions/bootstrapOOB_ENet.m`

**CanlabCore**

- `training` *(@chain)* — `ambiguous_within_repo`, 38 candidates

### `bootstrapOOB_PLSDA`

`secondlevel/functions/bootstrapOOB_PLSDA.m`

**CanlabCore**

- `training` *(@chain)* — `ambiguous_within_repo`, 38 candidates

### `bootstrapOOB_PLSR`

`secondlevel/functions/bootstrapOOB_PLSR.m`

**CanlabCore**

- `training` *(@chain)* — `ambiguous_within_repo`, 38 candidates

### `capLV`

`secondlevel/functions/capLV.m`

No external dependencies.

### `dice_statistic_image`

`secondlevel/functions/dice_statistic_image.m`

No external dependencies.

### `dice_statistic_image_by_roi`

`secondlevel/functions/dice_statistic_image_by_roi.m`

No external dependencies.

### `enetLambdaGrid`

`secondlevel/functions/enetLambdaGrid.m`

No external dependencies.

### `foldPreprocess`

`secondlevel/functions/foldPreprocess.m`

No external dependencies.

### `globalBaselineCV`

`secondlevel/functions/globalBaselineCV.m`

**CanlabCore**

- `training` *(@chain)* — `ambiguous_within_repo`, 38 candidates

### `group_tfce_from_subject_maps`

`secondlevel/functions/group_tfce_from_subject_maps.m`

**CanlabCore**

- `group` *(@group)*
- `statistic_image` *(@statistic_image)*

### `logitSafe`

`secondlevel/functions/logitSafe.m`

No external dependencies.

### `makeGroupedFolds`

`secondlevel/functions/makeGroupedFolds.m`

No external dependencies.

### `maskToSignificant`

`secondlevel/functions/maskToSignificant.m`

No external dependencies.

### `plot_ENet_diagnostics_neuroimaging`

`secondlevel/functions/plot_ENet_diagnostics_neuroimaging.m`

**spm12**

- `spm_read_vols` — `ambiguous_within_repo`, 2 candidates
- `spm_vol` — `ambiguous_within_repo`, 2 candidates
- `spm_write_vol` — `ambiguous_within_repo`, 3 candidates

### `plot_PLSDA_diagnostics_neuroimaging`

`secondlevel/functions/plot_PLSDA_diagnostics_neuroimaging.m`

**spm12**

- `spm_read_vols` — `ambiguous_within_repo`, 2 candidates
- `spm_vol` — `ambiguous_within_repo`, 2 candidates
- `spm_write_vol` — `ambiguous_within_repo`, 3 candidates

### `plot_PLSR_diagnostics_neuroimaging`

`secondlevel/functions/plot_PLSR_diagnostics_neuroimaging.m`

**spm12**

- `spm_read_vols` — `ambiguous_within_repo`, 2 candidates
- `spm_vol` — `ambiguous_within_repo`, 2 candidates
- `spm_write_vol` — `ambiguous_within_repo`, 3 candidates

### `quickCV_ENet`

`secondlevel/functions/quickCV_ENet.m`

**CanlabCore**

- `training` *(@chain)* — `ambiguous_within_repo`, 38 candidates

### `quickCV_ENet_PR`

`secondlevel/functions/quickCV_ENet_PR.m`

**CanlabCore**

- `training` *(@chain)* — `ambiguous_within_repo`, 38 candidates

### `quickCV_PLSDA`

`secondlevel/functions/quickCV_PLSDA.m`

No external dependencies.

### `quickCV_PLSDA_PR`

`secondlevel/functions/quickCV_PLSDA_PR.m`

No external dependencies.

### `quickCV_PLSDA_core`

`secondlevel/functions/quickCV_PLSDA_core.m`

**CanlabCore**

- `training` *(@chain)* — `ambiguous_within_repo`, 38 candidates

### `quickCV_PLSR`

`secondlevel/functions/quickCV_PLSR.m`

**CanlabCore**

- `training` *(@chain)* — `ambiguous_within_repo`, 38 candidates

### `quickGroupedCV`

`secondlevel/functions/quickGroupedCV.m`

No external dependencies.

### `residualizeFold`

`secondlevel/functions/residualizeFold.m`

No external dependencies.

### `residualizeY`

`secondlevel/functions/residualizeY.m`

No external dependencies.

### `selectENetHyperparams`

`secondlevel/functions/selectENetHyperparams.m`

**CanlabCore**

- `training` *(@chain)* — `ambiguous_within_repo`, 38 candidates

### `setParforStream`

`secondlevel/functions/setParforStream.m`

No external dependencies.

### `swapWithinSubjectLabels`

`secondlevel/functions/swapWithinSubjectLabels.m`

No external dependencies.

### `tfce_one_fmri_dat`

`secondlevel/functions/tfce_one_fmri_dat.m`

No external dependencies.

### `tfce_transform_3d`

`secondlevel/functions/tfce_transform_3d.m`

No external dependencies.

### `tfce_volume`

`secondlevel/functions/tfce_volume.m`

No external dependencies.

### `thresholded_fmri_data_from_pval_nii`

`secondlevel/functions/thresholded_fmri_data_from_pval_nii.m`

**CanlabCore**

- `apply_mask` *(@image_vector)*
- `autolabel_regions_using_atlas` *(@region)* — `dotcall`
- `region` *(@region)*
- `statistic_image` *(@statistic_image)*
- `threshold` *(@atlas)* — `dotcall`, 4 candidates

**canlab_single_trials**

- `fmri_data_st` *(@fmri_data_st)*

### `thresholded_fmri_data_from_statistic_image`

`secondlevel/functions/thresholded_fmri_data_from_statistic_image.m`

**CanlabCore**

- `region` *(@region)*

**canlab_single_trials**

- `fmri_data_st` *(@fmri_data_st)*

### `validateAtlasLabels`

`secondlevel/functions/validateAtlasLabels.m`

No external dependencies.

### `validateCovariates`

`secondlevel/functions/validateCovariates.m`

No external dependencies.

### `warnUnknownOptions`

`secondlevel/functions/warnUnknownOptions.m`

No external dependencies.

### `LaBGAScore_Storey_FDR`

`stats_tools/functions/LaBGAScore_Storey_FDR.m`

No external dependencies.

