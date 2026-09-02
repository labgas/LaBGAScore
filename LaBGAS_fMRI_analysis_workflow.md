# LaBGAS fMRI analysis workflow

**NOTE:** this can be considered the default LaBGAS fMRI analysis workflow from scanner to paper, based on the example of the Discoverie and Erythritol_4a projects, but there may obviously be study- or project-specific deviations from this default. Discuss this thoroughly with the project team before you start organizing your data with DataLad and other tools!

## Contents

- [Before you start](#before-you-start)
- [Data storage](#data-storage)
- [DICOM to NIfTI conversion \& BIDS organization](#dicom-to-nifti-conversion--bids-organization)
- [Set up QC and preprocessing pipeline](#set-up-qc-and-preprocessing-pipeline)
- [Quality control](#quality-control)
- [Preprocessing](#preprocessing)
- [Creating events.tsv files (and a phenotype file)](#creating-eventstsv-files-and-a-phenotype-file)
- [Smoothing](#smoothing)
- [Clean up sourcedata to save space (optional)](#clean-up-sourcedata-to-save-space-optional)
- [First-level analysis](#first-level-analysis)
- [Second-level analysis](#second-level-analysis)
- [Publish your dataset on GIN and/or Github](#publish-your-dataset-on-gin-andor-github)

## Before you start

### 1. Connect to the LaBGAS Linux server on your computer

Follow the instructions in *Getting started with the LaBGAS Linux server*, or mail [linuxteam.gbiomed@kuleuven.be](mailto:linuxteam.gbiomed@kuleuven.be) (or, better, create a ticket via the ICTS helpdesk) with Lukas in cc, to request access to the server if you do not have it yet.

### 2. Set up your X2go display for publishing figures

MATLAB's `publish` function, which produces the HTML reports for first- and second-level analyses, **captures figures from the screen**. Your X2go display settings therefore decide how the figures in your reports come out, and mismatched settings are why the same script can produce different-looking reports for two people.

Two settings matter, and they trade against each other:

- **Window size** — how many pixels are available.
- **Display DPI** (*Session preferences > Input/Output > "Set display DPI"*) — how those pixels convert to inches. Figure sizes are specified in inches, because MATLAB font sizes are in points (1 pt = 1/72 inch), so **raising the DPI shrinks the figure you can produce** at a given window size. High DPI makes menus comfortable and figures small.

The analysis scripts target a **12 × 7.5 inch** figure, which needs `12 × DPI` by `7.5 × DPI` pixels of window. Find your screen and use these settings:

| Your screen | X2go window | Set display DPI | Figure in report | Still works up to |
|---|---|---|---|---|
| 1366 × 768 laptop | fullscreen / maximum available | **96** | 1152 × 720 px | 96 (at the limit) |
| 1600 × 900 laptop | fullscreen / maximum available | **96** | 1152 × 720 px | 112 |
| 1920 × 1080 | fullscreen, or custom ≥ 1470 × 960 | **120** | 1440 × 900 px | 135 |
| 2560 × 1440 or larger | fullscreen | **144** | 1728 × 1080 px | 180 |
| 3440 × 1440 ultrawide, half width | custom ≈ 1720 × 1400 | **133** | 1596 × 998 px | 140 |

Set the value in the **"Set display DPI"** column; that is the single number to type. The last column is only there in case your window ends up smaller than assumed — anything at or below it still produces the full-size figure.

Why the recommended DPI is not simply the lowest that works: **as long as the target fits, a higher DPI is better.** The figure stays 12 × 7.5 inches either way, so it stays comparable with everyone else's, but it is captured at more pixels and is therefore sharper in the report — and on-screen text is more comfortable to work with. The only cost is that a higher DPI needs a bigger window.

Rules of thumb:

- **On a laptop, use fullscreen or "maximum available" and leave DPI at 96.** A small window at high DPI is the one combination that cannot produce a usable figure.
- **"Custom" width and height only set the size the session *starts* at.** If you drag the X2go window afterwards, the session resizes with it, and the new size is what `publish` captures. What counts is the window size at the moment you run the script.

To check your own session rather than read from the table:

```matlab
LaBGAScore_check_display        % from LaBGAScore/clean
```

It reports what your session can currently produce, whether it meets the target, and — if not — both fixes (lower the DPI to a stated value, or enlarge the window to a stated size). Run it once after changing your X2go settings.

If your session is too small, nothing breaks: `plugin_set_figure_size` scales the figure down proportionally, so the aspect ratio is always correct and nothing is cropped. The figure is simply smaller than a colleague's, and text looks relatively larger. `LaBGAScore_prov_publish` records your screen size, DPI and the resulting figure dimensions in every report, so such differences are visible afterwards instead of merely puzzling.

### 3. Study the following documents

Use the learning resources in them to familiarize yourself with the tools we will be using:

1. *Getting started with the LaBGAS Linux server*
2. *Getting started with Git, Github, GIN, and DataLad*
3. *Getting started with BIDS, mriqc, and fmriprep*

### 4. Take an fMRI course/tutorials which focuses on SPM

For example:

1. [Principles of fMRI](https://www.youtube.com/channel/UC_BIby85hZmcItMrkAlc8eA) (Tor Wager, Martin Lindquist)
2. Andrew Jahn
   - [Andy's Brain Blog](https://www.andysbrainblog.com/)
   - [Andy's Brain Book](https://andysbrainbook.readthedocs.io/en/latest/#)
3. [Designing and analysing fMRI experiments](https://ucl.podia.com/view/courses/designing-and-analysing-fmri-experiments/1228405-week-0-getting-started/3741966-welcome) (UCL)

### 5. Familiarize yourself with Matlab, SPM, and CANlab tools

1. [Matlab learning resources](https://matlabacademy.mathworks.com/?s_tid=pl_learn)
2. [Tutorials with Matlab live scripts](https://canlab.github.io/tutorials/)
3. [Walkthroughs](https://canlab.github.io/walkthroughs/)

### 6. Configure your Git identity

Follow [these instructions](http://handbook.datalad.org/en/latest/intro/installation.html#initial-configuration) in the Datalad Handbook.

**NOTE:** it is recommended (but not compulsory) to use the same username and e-mail as the one you use for your Github and GIN accounts!

## Data storage

### 1. Export your sourcedata and logfiles

From the scanner patient database/stimulus PC to your external hard disk, or to `E:\Export` on the scanner PC, immediately after scanning.

### 2. Transfer your sourcedata and logfiles

To the root folder for your project/study on the LaBGAS J-drive, from your external hard disk (on your own PC) or directly from `E:\Export` on the stimulus (not scanner) PC using WinSCP (choose the LaBGAS export configuration), sticking to BIDS-compliant subject naming and the folder structure of the project.

**NOTE:** the E-drive on the scanner PC is tiny (10GB) and can be finicky — it complains about being 90% full even when it is still well below 90%, preventing you from exporting to that drive! This can be prevented by:

1. cleaning the E-drive after copying to the J-drive, and/or getting rid of old stuff from other people who have not cleaned up their stuff;
2. exporting all your scans in one go after you completed them all, rather than doing interim exports during scanning — but if one of them has more than 16384 slices, the export of all of them will fall back to classic DICOM rather than the more convenient enhanced DICOM format!

Hence, always take your external hard drive with you as a fallback option!

**NOTE:** root folders for new projects on the J-drive need to be created by Lukas, and access needs to be assigned to team members by IT.

### 3. Create the corresponding root directory for your project or study under `/data` on the LaBGAS server as a DataLad superdataset

We use the [`datalad create`](http://docs.datalad.org/en/stable/generated/man/datalad-create.html) command:

```bash
cd /data
datalad create --description "superdataset for discoverie WP4 on LaBGAS server" -c text2git proj_discoverie
```

#### Adapt the superdataset configuration

This is done to prevent datalad from storing anything in `/code` subdirs, as well as files with a specific extension and small files, in the superdataset and its subdatasets in git-annex rather than git — see [this section](http://handbook.datalad.org/en/latest/basics/101-123-config2.html#gitattributes) of the Datalad handbook.

```bash
cd /data/proj_discoverie
nano .gitattributes
```

Add the following lines:

```
**/code/** annex.largefiles=nothing
** annex.largefiles=(largerthan=100kb)
```

Change this line:

```
* annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

into:

```
** annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

`Ctrl+O` to save the file, press enter, then `Ctrl+X` to exit the nano editor.

Save this change to the superdataset:

```bash
datalad save -m "adapted rules for annexing in .gitattributes"
```

#### Make the git (annex) repo writable by our user group

```bash
git config core.sharedrepository world
cd .git
chmod -R a+rwX ./*
```

#### Create a .gitignore file

This is done to keep certain files/directories completely out of git or git-annex, such as temp, work, or log directories, to prevent them from clogging the repos and keep the dataset lightweight!

```bash
cd <superdataset>   # do not stay in the .git subfolder
nano .gitignore
```

```
KUL_LOG/**
work/**
```

`Ctrl+O` to save the file, press enter, then `Ctrl+X` to exit the nano editor.

Save this change to the superdataset:

```bash
datalad save -m "created .gitignore file"
```

### 4. Create the sourcedata subdirectory as a subdataset

```bash
cd /data/proj_discoverie
datalad create --description "sourcedata subdataset for discoverie WP4 on LaBGAS server" -c text2git -d . sourcedata
```

#### Adapt the subdataset configuration

Same rationale as for the superdataset above — see [this section](http://handbook.datalad.org/en/latest/basics/101-123-config2.html#gitattributes) of the Datalad handbook.

```bash
cd /data/proj_discoverie/sourcedata
nano .gitattributes
```

Add the following line:

```
** annex.largefiles=(largerthan=100kb)
```

Change this line:

```
* annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

into:

```
** annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

`Ctrl+O` to save the file, press enter, then `Ctrl+X` to exit the nano editor.

Save this change to the subdataset, and one level higher to the superdataset:

```bash
datalad save -m "adapted rules for annexing in .gitattributes"
```

#### Make the git (annex) repo writable by our user group

```bash
git config core.sharedrepository world
cd .git
chmod -R a+rwX ./*
```

#### Create a .gitignore file in sourcedata

This is done to keep certain files/directories completely out of git or git-annex, such as DICOM and DICOMDIR files, to prevent them from clogging the repos and keep the dataset lightweight!

DICOMDIR can always be gitignored; ignoring DICOMs is optional, but a good idea if you have a lot of classic DICOM files (which we try to avoid, but is unavoidable in case of longer runs and/or shorter TRs, as there is a maximal number of volumes per run that can be stored as Philips enhanced DICOM, our default format).

```bash
cd sourcedata   # do not stay in the .git subfolder
nano .gitignore
```

```
**/DICOM
**/DICOMDIR
```

`Ctrl+O` to save the file, press enter, then `Ctrl+X` to exit the nano editor.

Save this change to the subdataset (and again to the superdataset afterwards):

```bash
datalad save -m "created .gitignore file"
```

### 5. Copy your subject directory from the J-drive to the sourcedata directory on the server

**NOTE:** if you know from your notes during scanning (which should be summarized in the first "scanning" sheet of your `trouble_in_paradise.xlsx` file, see below) that you will exclude an entire subject at this stage because too much went wrong during scanning, simply do not copy it over to the server!

#### Option a

1. Make the `/data` directory on the server accessible as a network drive in your Windows system, following the instructions in *Getting started with the LaBGAS Linux server*.
2. Copy the Windows way.
3. Record this change using the [`datalad save`](http://docs.datalad.org/en/stable/generated/man/datalad-save.html) command:
   - in the subdataset:

     ```bash
     cd /data/proj_discoverie/sourcedata
     datalad save -m "copied sourcedata for sub-KUL002 from LaBGAS J-drive on KU Leuven server"
     ```

   - in the superdataset:

     ```bash
     cd /data/proj_discoverie
     datalad save -m "added sub-KUL001 to sourcedata subdataset"
     ```

#### Option b

1. Mount the LaBGAS J-drive in the server's Linux filesystem:

   ```bash
   sudo gbiomed-mount J-drive
   ```

   The J-drive is now accessible in your Linux filesystem under `/media/<your u-number>/J-Drive`.

2. [Copy in your Linux terminal using the rsync command](https://linux.die.net/man/1/rsync):

   ```bash
   rsync -a --chmod=g+w /media/u0115135/J-drive/GBW-0264_TARGID-Brain-Gut-Axis/0018_proj-moodbugs/proj-moodbugs_wp4/sourcedata/sub-INF007 /data/proj_moodbugs/proj_moodbugs_wp4/sourcedata/.
   ```

   *(if you want a specific folder from within a subject, just copy the directory you want, e.g. `J:\GBW-0264_TARGID-Brain-Gut-Axis\0018_proj-moodbugs\proj-moodbugs_wp4\sourcedata\sub-INF001\ses-01\PET\` — make sure all slashes are forward slashes; this is where it ends up on the server: `/data/proj_moodbugs/proj_moodbugs_wp4/sourcedata`)*

   **NOTE:** there's a space between your subject folder and the `/data/...` destination directory.

   **NOTE:** we prefer the `rsync` command here over the more common `cp` command, as it has an option (`--chmod`) that lets you set permissions correctly while copying, and only copies over files not yet present in the destination.

   **TIPS:**
   - if you want to copy an entire directory with subdirectories, add the `-r` flag to the rsync command (`rsync -ar`);
   - if you want to perform a dry run of your copy (i.e. no changes made), add the `-n` flag to the rsync command (`rsync -an`).

   Save to the subdataset:

   ```bash
   cd /data/proj_moodbugs/proj_moodbugs_wp4/sourcedata
   datalad save -m "copied sourcedata for sub-INF007 from LaBGAS J-drive on KU Leuven server"
   ```

   Save the change to the subdataset to the superdataset:

   ```bash
   cd /data/proj_moodbugs/proj_moodbugs_wp4
   datalad save -m "sourcedata subdataset: copied sourcedata for sub-INF003 from LaBGAS J-drive on KU Leuven server"
   ```

**NOTE #1:** in either scenario, do NOT copy your DICOMDIR file over (or delete it after copying it to the server), as it causes trouble in the next step! In case you still do this, `**/DICOMDIR` in your `.gitignore` file (see above) prevents it from clogging the git(-annex) repo.

**NOTE #2:** if you have multiple sessions (for example two study visits), organize your sourcedata in separate folders per session (`ses-01`, or a name for your session).

### 6. Check your DICOMs using fsleyes

File > Add from DICOM > select DICOM directory (but do not open it), choose the scan(s) you want to open from the menu (all of them in one go is easiest — fsleyes allows you to toggle them on and off one by one), and click Open.

**NOTE:** MRIcroGL is a good alternative to fsleyes!

- Log any abnormalities in your `trouble_in_paradise.xlsx` file in your superdataset directory at this and any previous or following QC step (including at the scanner)! [Here](https://www.dropbox.com/s/jwgqn2fter3bb7z/trouble_in_paradise.xlsx?dl=0) is an example of such a file (on a messy dataset that can still be improved — no tab for abnormalities during scanning yet).
- If you are certain you have to **exclude the entire subject** at this early stage of processing (for example because of incomplete brain coverage in all scans, or in the structural scan):
  - create a new subdirectory `sourcedata/off_study`;
  - move your entire subject directory there;
  - save the change to your datalad sub- and superdataset — see [instructions](http://handbook.datalad.org/en/latest/basics/101-136-filesystem.html#moving-files-from-or-into-subdirectories) about moving files/folders.
- If you think you may have to **exclude certain runs**:
  - do not do anything yet (moving, renaming, deleting, etc. raw DICOMs is a bad idea);
  - proceed and check again after DICOM-to-NIfTI conversion.

### 7. Zip your DICOMs on the J-drive

Use 7-zip with these settings, and delete the unzipped files afterwards.

## DICOM to NIfTI conversion & BIDS organization

**NOTE:** we use (a) the non-nested data organization option (i.e. other subdatasets, including sourcedata, derivatives, etc., are NOT nested within the BIDS subdataset — each subdirectory is a datalad subdataset under the root directory/superdataset), and (b) a derivatives subdirectory with one subdirectory per pipeline, as described [here](https://bids-specification.readthedocs.io/en/stable/02-common-principles.html#storage-of-derived-datasets), with each pipeline subdirectory NOT being a datalad subdataset itself, just a subdirectory of the derivatives subdataset.

### 1. Create the code subdirectory as a subdataset

```bash
cd /data/proj_discoverie
datalad create --description "code subdataset for discoverie project on LaBGAS server" -c text2git -d . code
```

Adapt the subdataset configuration:

```bash
cd /data/proj_discoverie/code
nano .gitattributes
```

Change as follows:

```
* annex.backend=MD5E
**/.git* annex.largefiles=nothing
** annex.largefiles=nothing
```

`Ctrl+O` followed by enter to save, `Ctrl+X` to close nano.

Save this change to the subdataset:

```bash
datalad save -m "adapted rules for annexing in .gitattributes"
```

We use a fixed structure for the code directory across all projects — for now, we need the subdirectory `prep`:

```bash
mkdir prep
```

#### Make the git (annex) repo writable by our user group

```bash
git config core.sharedrepository world
cd .git
chmod -R a+rwX ./*
```

### 2. Download KUL_dcm2bids.sh from Github in the code subdataset using the datalad download-url command

**NOTE:** you need the raw file URL — see [this Neurostars thread](https://neurostars.org/t/datalad-download-url-failing-to-download-scripts-from-github-correctly/20353).

```bash
cd /data/proj_discoverie/code
datalad download-url -m "downloading KUL_dcm2bids.sh from Github" -d . --path=prep/KUL_dcm2bids.sh https://raw.githubusercontent.com/treanus/KUL_NIS/master/KUL_dcm2bids.sh
```

Save this modification of the subdataset to your superdataset:

```bash
cd /data/proj_discoverie
datalad save -m "downloaded KUL_dcm2bids.sh from Github into code/prep"
```

### 3. Create the BIDS subdirectory as a subdataset

Output of `KUL_dcm2bids.sh` (.nii.gz images, BIDS-compliant) will be written here.

```bash
datalad create --description "BIDS subdataset for discoverie project" -c text2git -c bids -d . BIDS
```

#### Adapt the subdataset configuration

```bash
cd /data/proj_discoverie/BIDS
nano .gitattributes
```

Modify as follows:

```
* annex.backend=MD5E
**/.git* annex.largefiles=nothing
.bidsignore annex.largefiles=nothing
CHANGES annex.largefiles=nothing
README annex.largefiles=nothing
code/** annex.largefiles=nothing
** annex.largefiles=(largerthan=100kb)
```

Change this line:

```
* annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

into:

```
** annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

`Ctrl+O` to save the file, followed by enter, `Ctrl+X` to exit nano.

Save this change to the subdataset (and again to the superdataset afterwards):

```bash
datalad save -m "adapted rules for annexing in .gitattributes"
```

#### Make the git (annex) repo writable by our user group

From the BIDS folder:

```bash
git config core.sharedrepository world
cd .git
chmod -R a+rwX ./*
```

#### Create a .gitignore file in BIDS

This is done to keep certain files/directories completely out of git or git-annex, such as temp, work, or log directories, to prevent them from clogging the repos and keep the dataset lightweight!

```bash
cd ..
nano .gitignore
```

```
tmp_dcm2bids/**
sourcedata/**
derivatives/**
```

`Ctrl+O` to save the file, followed by enter, `Ctrl+X` to exit nano.

Save this change to the subdataset (and again to the superdataset afterwards):

```bash
datalad save -m "created .gitignore file"
```

### 4. Run KUL_dcm2bids.sh wrapped in a datalad run command

#### Create your sequences.txt file

```bash
cd /data/proj_discoverie
mkdir study_config
```

Follow the instructions in *Getting started with BIDS, mriqc, and fmriprep*.

#### Run the command

`KUL_dcm2bids.sh` is a shell script, so unlike the Matlab scripts later in this workflow it can be wrapped in a `datalad run` command:

```bash
cd /data/proj_moodbugs/proj_moodbugs_wp4
datalad run -m "run KUL_dcm2bids.sh on sub-INF003" --input "sourcedata/sub-INF003/*" --output "BIDS/sub-INF003/*" "KUL_dcm2bids.sh -d sourcedata/sub-INF003 -p INF003 -c study_config/sequences.txt -v"
```

**NOTE #1:** if you have multiple sessions (for example two study visits), organize your sourcedata in separate folders per session (`ses-01`, or a name for your session), and run `KUL_dcm2bids.sh` per session using the following command:

```bash
datalad run -m "run KUL_dcm2bids.sh on sub-INF003 session 1" --input "sourcedata/sub-INF003/ses-01/*" --output "BIDS/sub-INF003/ses-01/*" "KUL_dcm2bids.sh -d sourcedata/sub-INF003/ses-01 -p INF003 -c study_config/sequences.txt -v"
```

**NOTE #2:** if `KUL_dcm2bids.sh` runs correctly (zero exit code in Linux terms), the `datalad run` command will automatically perform a `datalad save` to the subdataset and the superdataset, leaving your working tree clean upon a `datalad status -r`! Cool stuff!

**NOTE #3:** if you have an appropriately named study-specific [`*_events.tsv`](https://bids-specification.readthedocs.io/en/stable/04-modality-specific-files/05-task-events.html) file at this stage already, and place it in your `study_config/` directory under your superdataset, the script will copy it to the right place in the BIDS subdataset if you add the `-e` flag to the above command! This may also work for subject- or run-specific events.tsv files, but that has not been tested yet. See below for more about events.tsv files.

**NOTE #4:** `KUL_dcm2bids.sh` leaves a `tmp_dcm2bids` subdirectory with a bunch of correctly converted but not BIDS-compliant, renamed .nii.gz files corresponding to acquisitions (e.g. revphase, phase1/phase2) that are not in `sequences.txt`. However, thanks to the `.gitignore` file, this directory is not added to git nor git-annex! Hence, these can be deleted without requiring a `datalad save`, but you may want to keep the log:

```bash
cd /data/proj_discoverie/BIDS/tmp_dcm2bids
rm -r sub-INF003
```

### 5. Check your converted .nii.gz files using fsleyes

- Log any abnormalities in your `trouble_in_paradise.xlsx` file in your superdataset directory!
- If you notice **abnormalities that were not present when checking your DICOMs** for the same subject, something went wrong with conversion — you may want to try again, or use another conversion method.
- If you decide to **exclude an entire subject** at this stage (unlikely — typically you will exclude because of fatally bad DICOMs, or after QC and preprocessing):
  - do not proceed with the next steps for this subject;
  - create a new subdirectory `BIDS/off_study`;
  - move (not copy) your entire subject directory there;
  - save the change to your datalad sub- and superdataset — see [instructions](http://handbook.datalad.org/en/latest/basics/101-136-filesystem.html#moving-files-from-or-into-subdirectories) about moving files/directories.
- If you decide to **exclude (a) certain run(s)** at this stage (also unlikely):
  - do proceed with the next steps for this subject;
  - that run will be moved to a subject-specific `off_study` directory after QC, preprocessing, and creation of events.tsv files (see below).

## Set up QC and preprocessing pipeline

**NOTE:** we largely follow the procedures in [this usecase](http://handbook.datalad.org/en/latest/beyond_basics/101-171-enki.html) of the Datalad handbook, which crosslinks to [this](http://handbook.datalad.org/en/latest/basics/101-133-containersrun.html) and [this](http://handbook.datalad.org/en/latest/usecases/reproducible_neuroimaging_analysis.html#usecase-reproduce-neuroimg) relevant section.

### Create the pipeline subdataset

```bash
datalad create --description "pipeline subdataset for mriqc and fmriprep containers for discoverie project" -d . -c text2git pipeline
```

#### Make the git (annex) repo writable by our user group

```bash
cd pipeline
git config core.sharedrepository world
cd .git
chmod -R g+rwX ./*
```

### Add the mriqc container to the pipeline subdataset

**NOTE:** the container with the latest stable mriqc version is available on the [Datalad Repository](http://datasets.datalad.org/?dir=/repronim/containers/images/bids), and we use the [`datalad containers-add`](http://docs.datalad.org/projects/container/en/latest/generated/man/datalad-containers-add.html) command to add it to our subdataset. Hence, change the URL below to the latest stable version!

```bash
cd /data/proj_discoverie/pipeline
datalad containers-add mriqc --url https://datasets.datalad.org/repronim/containers/images/bids/bids-mriqc--24.0.2.sing --call-fmt 'singularity run --cleanenv -B "$PWD,/data/singularity/:/data/singularity/" {img} {cmd}'
```

### Copy the freesurfer license file into the pipeline subdataset

```bash
cp /opt/KUL_apps/freesurfer/license.txt .
datalad save --to-git -m "add freesurfer license-file" license.txt
```

### Add the fmriprep container to the pipeline subdataset

**NOTE:** the container with the latest stable fmriprep version is available on the [Datalad Repository](http://datasets.datalad.org/), and we use the [`datalad containers-add`](http://docs.datalad.org/projects/container/en/latest/generated/man/datalad-containers-add.html) command to add it to our subdataset. However, do not upgrade fmriprep beyond version 24 in the link below, because more recent versions (25 and above) have issues with the fieldmap-less distortion correction we are using!

```bash
datalad containers-add fmriprep --url https://datasets.datalad.org/repronim/containers/images/bids/bids-fmriprep--24.0.1.sing --call-fmt 'singularity run --cleanenv -B "$PWD,/data/singularity/:/data/singularity/" {img} {cmd} --fs-license-file "$PWD/{img_dspath}/license.txt"'
```

## Quality control

### 1. Create the mriqc subdataset

Go back to your superdataset folder:

```bash
datalad create --description "mriqc subdataset for discoverie project" -c text2git -d . mriqc
```

#### Adapt the subdataset configuration

```bash
cd /data/proj_discoverie/mriqc
nano .gitattributes
```

```
* annex.backend=MD5E
**/.git* annex.largefiles=nothing
```

Add this line:

```
** annex.largefiles=(largerthan=100kb)
```

Change this line:

```
* annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

into:

```
** annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

`Ctrl+O` to save, enter, `Ctrl+X` to exit nano.

Save this change to the subdataset (and again to the superdataset afterwards):

```bash
datalad save -m "adapted rules for annexing in .gitattributes"
```

#### Make the git (annex) repo writable by our user group

```bash
git config core.sharedrepository world
cd .git
chmod -R a+rwX ./*
```

### 2. Run mriqc wrapped in a datalad containers-run command

**NOTE:** contrary to the mriqc documentation, we use [singularity containers to run](https://www.nipreps.org/apps/singularity/) mriqc rather than docker containers, since running with docker leads to permission issues.

See the [`datalad containers-run` docs](http://docs.datalad.org/projects/container/en/latest/generated/man/datalad-containers-run.html).

```bash
cd /data/proj_discoverie
datalad containers-run -m "run mriqc on sub-KUL001" -i "BIDS/sub-KUL001/*" -o "mriqc/sub-KUL001/*" -n pipeline/mriqc BIDS mriqc participant --participant-label KUL001 --verbose --verbose-reports --fd_thres 0.9 --n_cpus 24 --mem-gb 128
```

**NOTE:** in addition to correctly putting the output in the mriqc directory, this command currently creates a `work` directory with lots of temporary junk under the root of the superdataset, but thanks to the `.gitignore` file in the superdataset, this directory is not added to git-annex and can simply be deleted without needing a `datalad save` afterwards:

```bash
rm -r work
```

### 3. Check the html output of mriqc

- Log any abnormalities in your `trouble_in_paradise.xlsx` file in the root of your project directory! Some examples of head motion artifacts can be found on the last page of the MR8 scan protocol.
- **Do NOT exclude subjects based on mriqc alone** — always run fmriprep and decide afterwards!

## Preprocessing

**NOTE:** we largely follow the procedures in [this usecase](http://handbook.datalad.org/en/latest/beyond_basics/101-171-enki.html) of the Datalad handbook, which crosslinks to [this](http://handbook.datalad.org/en/latest/basics/101-133-containersrun.html) and [this](http://handbook.datalad.org/en/latest/usecases/reproducible_neuroimaging_analysis.html#usecase-reproduce-neuroimg) relevant section.

### 1. Create the derivatives subdataset

```bash
datalad create --description "derivatives subdataset for moodbugs wp3 project" -c text2git -d . derivatives
```

#### Adapt the subdataset configuration

```bash
cd /data/proj_moodbugs/proj_moodbugs_wp4/derivatives
nano .gitattributes
```

```
* annex.backend=MD5E
**/.git* annex.largefiles=nothing
```

Add this line:

```
** annex.largefiles=(largerthan=100kb)
```

Change this line:

```
* annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

into:

```
** annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

`Ctrl+O` to save, enter, `Ctrl+X` to exit nano.

Save this change to the subdataset (and again to the superdataset afterwards):

```bash
datalad save -m "adapted rules for annexing in .gitattributes"
```

#### Make the git (annex) repo writable by our user group

Make sure you are in the correct `derivatives` folder:

```bash
git config core.sharedrepository world
cd .git
chmod -R a+rwX ./*
```

### 2. Run fmriprep wrapped in a datalad containers-run command

**NOTE:** contrary to the fmriprep documentation, we use [singularity containers to run fmriprep](https://www.nipreps.org/apps/singularity/) rather than docker containers, since running with docker leads to permission issues (see this [Neurostars thread](https://neurostars.org/t/problem-wrapping-bash-commands-in-datalad-run/20232/23))!

```bash
cd /data/proj_discoverie
datalad containers-run -m "run fmriprep on sub-KUL001" -i "BIDS/sub-KUL001/*" -o "derivatives/fmriprep/sub-KUL001/*" -n pipeline/fmriprep BIDS derivatives participant --participant-label KUL001 --fs-no-reconall --fd-spike-threshold 0.9 --dvars-spike-threshold 2 --verbose --output-spaces MNI152NLin2009cAsym:res-2 --ignore fieldmaps --use-syn-sdc --n-cpus 24 --mem-mb 128000
```

**NOTE:** we use the [fieldmap-less susceptibility distortion correction method](https://fmriprep.org/en/0.8.0/api/index.html#sdc-fieldmapless).

**NOTE:** in addition to correctly putting the output in the `derivatives/fmriprep` directory, this command currently creates a `work` directory with lots of temporary junk under the root of the subdataset, but thanks to the `.gitignore` file in the superdataset, this directory is not added to git-annex and can simply be deleted without needing a `datalad save` afterwards:

```bash
rm -r work
```

### 3. Check the html output of fmriprep

- Log any abnormalities in your `trouble_in_paradise.xlsx` file in the root of your project directory! Some examples of head motion artifacts can be found on the last page of the MR8 scan protocol.
- If you decide to **exclude the entire subject** at this stage (based on a priori criteria noted in the respective sheet of your `trouble_in_paradise.xlsx` file):
  - *derivatives subdataset*: create a new subdirectory `derivatives/fmriprep/off_study`, and move (not copy) your entire subject directory there;
  - *BIDS and sourcedata subdatasets*: move your entire subject directory into newly created `BIDS/off_study` and `sourcedata/off_study` directories, to keep the list of subjects in sourcedata, BIDS, and derivatives subdatasets consistent;
  - save the changes to your datalad sub- and superdatasets — see [instructions](http://handbook.datalad.org/en/latest/basics/101-136-filesystem.html#moving-files-from-or-into-subdirectories) about moving files/directories.
- If you decide to **exclude (a) certain run(s)** at this stage (based on a priori criteria noted in the respective sheet of your `trouble_in_paradise.xlsx` file):
  - *derivatives subdataset*: create a new subject-specific subdirectory `derivatives/fmriprep/sub-xxx/func/off_study`, and move (not copy) all files corresponding to the excluded run there;
  - *BIDS subdataset*: create a new subject-specific subdirectory `BIDS/sub-xxx/func/off_study`, and move (not copy) all files corresponding to the excluded run there;
  - *sourcedata subdataset*: create a new subject-specific subdirectory **within the logfiles subdirectory**, `sourcedata/sub-xxx/logfiles/off_study`, and move (not copy) the log files corresponding to the excluded run(s) there — do NOT move nor delete DICOMs!
  - save the changes to your datalad sub- and superdatasets — see [instructions](http://handbook.datalad.org/en/latest/basics/101-136-filesystem.html#moving-files-from-or-into-subdirectories) about moving files/directories.

## Creating events.tsv files (and a phenotype file)

### Case #1: study-specific events.tsv files

For some tasks (like the MIST used in proj_discoverie), the onsets and durations of (blocks of) stimuli are the same for every run and every subject. In this case, it suffices to manually create one `task-MIST_events.tsv` file for the entire study and copy it to the root of your BIDS subdataset (as noted above, `KUL_dcm2bids.sh` will put it there for you if you include it in your `study_config` directory in your superdataset) — all subjects and runs will inherit it due to [BIDS' inheritance principle](https://bids-specification.readthedocs.io/en/stable/02-common-principles.html#the-inheritance-principle).

### Case #2: subject-specific events.tsv files

There may be cases where onsets and durations of (blocks of) stimuli are the same for every run in a subject, but not across subjects, although this case is expected to be rare and has not yet been encountered, so the optimal approach has not been tested. The "manual" approach above may work, as well as a more automated/scripted approach below. In any case, your subject-specific file should be named `sub-KUL001_task-MIST_events.tsv` and placed under `BIDS/sub-KUL001/func` for each subject.

### Case #3: run-specific events.tsv files

For most tasks (like the sweettaste task used in proj_erythritol_4a), the onsets and/or durations are randomized in each run in each subject. In this case, creating the run-specific `sub-KUL001_task-sweettaste_run-x_events.tsv` files manually would be very tedious and error-prone, so you need to write a (Matlab) script that generates these automatically from the logfiles generated by Presentation or other software used during scanning to present stimuli and log responses. This script will hence be highly study-specific, but there is an example in the [LaBGAScore GitHub repo](https://github.com/labgas/LaBGAScore) which you can download, rename, and adapt to your study. This example script also includes an option (turned on by default) to create a single phenotype file (e.g. `ratings.tsv`).

#### 1. Download the example script from the LaBGAScore GitHub repo and rename it

**NOTE:** you need the raw file URL — see [this Neurostars thread](https://neurostars.org/t/datalad-download-url-failing-to-download-scripts-from-github-correctly/20353).

```bash
cd /data/proj_erythritol/proj_erythritol_4a/code
datalad download-url -m "downloading example script to create events.tsv files from LaBGAScore Github repo" -d . --path=prep/ery_4a_prep_s1_write_events_tsv.m https://raw.githubusercontent.com/labgas/LaBGAScore/main/prep/LaBGAScore_prep_s1_write_events_tsv.m
```

Since this script calls the [`LaBGAScore_prep_s0_define_directories`](https://github.com/labgas/LaBGAScore/blob/main/prep/LaBGAScore_prep_s0_define_directories.m) LaBGAScore script under the hood, we want to add that script to our code subdataset in the same way:

```bash
datalad download-url -m "downloading helper script defining directories from LaBGAScore Github repo" -d . --path=prep/ery_4a_prep_s0_define_directories.m https://raw.githubusercontent.com/labgas/LaBGAScore/main/prep/LaBGAScore_prep_s0_define_directories.m
```

Save these modifications of the subdataset to your superdataset:

```bash
cd /data/proj_erythritol/proj_erythritol_4a
datalad save -m "downloaded example script to create events.tsv files, and helper script, from LaBGAScore Github repo"
```

#### 2. Add the code subdataset to your Matlab path

```bash
matlab &
```

Execute the following commands in the Matlab terminal:

```matlab
addpath(genpath('/data/proj_erythritol/proj_erythritol_4a/code'),'-end');
savepath;
```

#### 3. Adapt the script to your study

**NOTE:** to test/debug your script, it is recommended to make a copy of the sourcedata, BIDS, and derivatives directories for one (or a few, if you want to test the option to loop over subjects) subject(s) before running it on your datalad dataset. Make sure this copy is not inside your superdataset. Use the `--dereference` flag in your [`cp` command](http://manpages.ubuntu.com/manpages/trusty/man1/cp.1.html) to avoid copying the symlinks to git-annex outside of the datalad dataset!

#### 4. Run the script

Run it from the Matlab command line, with the root directory of your superdataset as the working directory:

```matlab
ery_4a_prep_s1_write_events_tsv
```

Then save the `events.tsv` files (and `phenotype.tsv`, if you chose that option) it wrote to the BIDS subdataset and the superdataset:

```bash
cd /data/proj_erythritol/proj_erythritol_4a
datalad save -m "ran ery_4a_prep_s1_write_events_tsv.m on all subjects"
```

#### 5. Check the events.tsv files (and phenotype.tsv file if that option is chosen) generated by the script!

Go back to debugging/testing if you encounter errors! There are some sanity checks built into the example script above, which you should keep for your script, but no errors or warnings does not automatically imply that the onsets and durations in your output files are correct!

## Smoothing

Use the [`LaBGAScore_prep_s2_smooth.m`](https://github.com/labgas/LaBGAScore/blob/main/prep/LaBGAScore_prep_s2_smooth.m) script, which creates, saves, and runs SPM smoothing batches.

### 1. Download the script from the LaBGAScore GitHub repo and rename it

**NOTE:** you need the raw file URL — see [this Neurostars thread](https://neurostars.org/t/datalad-download-url-failing-to-download-scripts-from-github-correctly/20353).

```bash
cd /data/proj_erythritol/proj_erythritol_4a/code
datalad download-url -m "downloading smoothing script from LaBGAScore Github repo" -d . --path=prep/ery_4a_prep_s2_smooth.m https://raw.githubusercontent.com/labgas/LaBGAScore/main/prep/LaBGAScore_prep_s2_smooth.m
```

Save this modification of the subdataset to your superdataset:

```bash
cd /data/proj_erythritol/proj_erythritol_4a
datalad save -m "code subdataset: downloaded smoothing script from Github"
```

### 2. Check the script and change options if needed

```bash
cd /data/proj_erythritol/proj_erythritol_4a/code/prep
matlab &
```

```matlab
edit ery_4a_prep_s2_smooth.m
```

**NOTE #1:** the options for smoothing defined in the first section of the script (fwhm for kernel, prefix for smoothed images) can be considered LaBGAS defaults, but there may be study-specific reasons to change them!

**NOTE #2:** the `subjs2smooth` option is turned off by default, which causes the script to loop over all subjects in your `derivatives/fmriprep` directory. We recommend smoothing all of your subjects in one go at the end of data collection (contrary to all the steps above, which you should do ASAP after scanning a subject) — in that case, you can keep the default options in the core script. However, if you want to smooth only one or a few subjects, you can specify them as a cell array in the `subjs2smooth` variable at the beginning of the script. In that case, you may not always want to `datalad save` the minor changes to the script, so you may discard the changes in your workspace (but be careful with this destructive git command — do not use it if you don't understand what it is doing!):

```bash
git reset --hard
```

### 3. Run the script

Run it from the Matlab command line, with the root directory of your superdataset as the working directory:

```matlab
ery_4a_prep_s2_smooth
```

Then save the smoothed images it wrote to the derivatives subdataset and the superdataset:

```bash
cd /data/proj_erythritol/proj_erythritol_4a
datalad save -m "ran ery_4a_prep_s2_smooth.m on all subjects"
```

## Clean up sourcedata to save space (optional)

From this stage onwards, or later, we can save a considerable amount of space by cleaning up the sourcedata subdataset, **provided its structure matches exactly the structure of the sourcedata folder on the J-drive**.

**NOTE:** only do this if you know what you are doing — RTFM and/or talk to Lukas!

There is now a [LaBGAScore Matlab function](https://github.com/labgas/LaBGAScore/blob/main/clean/LaBGAScore_move_repos_matlabpath.m) to do this automatically, but again, know what you are doing!

There are two possible scenarios.

### DICOM files are stored in git-annex

This is the case if `**/DICOM` is not added to your `.gitignore` file in your sourcedata subdataset (see [above](#create-a-gitignore-file-in-sourcedata)).

In this case, we use the very useful [`datalad drop`](https://docs.datalad.org/en/stable/generated/man/datalad-drop.html) command, rather than simply using Linux' [`rm`](https://manpages.ubuntu.com/manpages/bionic/man1/rm.1.html) command. This removes the files irreversibly from the local subdataset, but keeps the symlinks in git (with the purpose of being able to get them back using [`datalad get`](https://docs.datalad.org/en/stable/generated/man/datalad-get.html) if a remote copy is available, for example on GIN — see below), but this is not the case here since we never publish the sourcedata subdataset, as it contains patient info. However, we have a copy of the sourcedata on the J-drive, so dropping the files here is safe.

**NOTE:** if your files are stored in git-annex, this is the only way to save space, since using `rm` will remove the files and their symlinks, but keep them in git-annex's history, from which they can be retrieved by moving back in time — they keep on using space there! This is clearly explained [here](https://handbook.datalad.org/en/latest/basics/101-136-filesystem.html#removing-a-file-but-keeping-content-in-history) in the Datalad handbook.

```bash
cd /data/superdataset/sourcedata
datalad drop --reckless availability
```

Note that the `--reckless` flag is needed to override the default safety setting, which only allows dropping if at least one verified remote copy exists — not the case here, as explained above.

### DICOM files are gitignored

In case `**/DICOM` is added to your `.gitignore` file in your sourcedata subdataset, the files cannot be dropped, only removed using `rm` followed by `datalad save`, which will however save space since the files were never stored in git-annex and hence do not stay there nor take up space.

**NOTE:** do NOT remove the folder structure with one folder per subject (needed for sanity checks in subsequent scripts) — remove content in each subject folder instead (the abovementioned script does this); below is the manual way:

```bash
cd /data/superdataset/sourcedata/sub-001
rm -r *
```

Repeat for every subject folder, and `datalad save` at the end.

## First-level analysis

First-level analysis can be run flexibly using the three scripts in the [`firstlevel` directory of the LaBGAScore Github repo](https://github.com/labgas/LaBGAScore/tree/main/firstlevel), plus `s1a`/`s2a` as the multisession/multitask variant. A separate phMRI (pharmacological challenge) chain (`s1b`/`s2b`/`s3b`/`s3c`) exists for continuous designs in which a substance is administered mid-scan.

The scripts themselves are extensively annotated and documented, and that documentation is consolidated in two reference guides in the same folder: [`firstlevel/README.md`](https://github.com/labgas/LaBGAScore/blob/main/firstlevel/README.md) for the task-fMRI chain (the `LaBGAS_options` reference, the `DSGN` structure and the CANlab documentation for it, the noise model, parametric modulators, single-trial models, what lands on disk, and the contrast-order contract with the second-level scripts) and [`firstlevel/README_phMRI.md`](https://github.com/labgas/LaBGAScore/blob/main/firstlevel/README_phMRI.md) for the phMRI chain. The procedure below — creating the subdataset, downloading and renaming the scripts, and running them — is not repeated there.

Another great resource on first-level analysis using CANlab tools is Michael Sun's Matlab live script walkthrough, which can be found [here](https://drive.google.com/drive/folders/1-M5UvibmsWXVCIrR31-qJNu506pDA_0t).

### 1. Create the first level subdataset

```bash
datalad create --description "firstlevel subdataset for erythritol_4a project" -c text2git -d . firstlevel
```

#### Adapt the subdataset configuration

**NOTE:** this may need to be slightly adapted after thoroughly checking which types of files are annexed and which are not with the default configuration.

```bash
cd /data/proj_erythritol/proj_erythritol_4a/firstlevel
nano .gitattributes
```

```
* annex.backend=MD5E
**/.git* annex.largefiles=nothing
```

Add this line:

```
** annex.largefiles=(largerthan=100kb)
```

Change this line:

```
* annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

into:

```
** annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

`Ctrl+O` to save, enter, `Ctrl+X` to exit nano.

Save this change to the subdataset (and again to the superdataset afterwards):

```bash
datalad save -m "adapted rules for annexing in .gitattributes"
```

#### Make the git (annex) repo writable by our user group

```bash
cd /data/proj_erythritol/proj_erythritol_4a/firstlevel
git config core.sharedrepository world
cd .git
chmod -R a+rwX ./*
```

#### Create a .gitignore file

This is done to keep certain files/directories completely out of git or git-annex, such as temp, work, or log directories, to prevent them from clogging the repos and keep the dataset lightweight!

```bash
cd /data/proj_erythritol/proj_erythritol_4a/firstlevel
nano .gitignore
```

```
**/canlab_glm_logs/*
```

`Ctrl+O` to save the file, followed by enter, `Ctrl+X` to exit nano.

Save this change to the subdataset (and again to the superdataset afterwards):

```bash
datalad save -m "created .gitignore file"
```

### 2. Create a model-specific subdirectory

**IMPORTANT NOTE:** even if you do not anticipate running several first-level models, always create a model-specific directory to keep things consistent and organized, e.g.:

```bash
cd /data/proj_erythritol/proj_erythritol_4a/firstlevel
mkdir model_1_conds_pmods
```

### 3. Download the scripts from the LaBGAScore GitHub repo to your code subdataset and rename them

We want to include all our code in the code subdataset of our datalad superdataset, also the generic scripts, so we consistently download and rename them, like we did with the earlier prep scripts.

**IMPORTANT NOTE:** always create a model-specific directory (with the same name as the corresponding model-specific directory in your firstlevel — and later secondlevel — subdatasets), and include a model index (`m1`, etc.) in your script name, even if you do not anticipate running several first-level models, to keep things consistent and organized.

Example for the first script:

```bash
cd /data/proj_erythritol/proj_erythritol_4a/code
mkdir firstlevel
mkdir model_1_conds_pmods
```

Run the following from the code directory, NOT the subdirectories:

```bash
datalad download-url -m "downloading first level design script from LaBGAScore Github repo" -d . --path=firstlevel/model_1_conds_pmods/ery_4a_firstlevel_m1_s1_options_dsgn_struct.m https://raw.githubusercontent.com/labgas/LaBGAScore/main/firstlevel/LaBGAScore_firstlevel_s1_options_dsgn_struct.m
```

Save this modification of the subdataset to your superdataset:

```bash
cd /data/proj_erythritol/proj_erythritol_4a
datalad save -m "code subdataset: downloaded first level design script from Github"
```

### 4. Check the script(s) and make study-specific adaptations where needed

The scripts themselves are extensively annotated and documented, hence we refer to their documentation — and to [`firstlevel/README.md`](https://github.com/labgas/LaBGAScore/blob/main/firstlevel/README.md), which consolidates it — here. The table below gives a short description of each of the three core scripts, based on the *USAGE* section of their Matlab header documentation.

| # | Script | Description |
|---|---|---|
| 1 | [`LaBGAScore_firstlevel_s1_options_dsgn_struct`](https://github.com/labgas/LaBGAScore/blob/main/firstlevel/LaBGAScore_firstlevel_s1_options_dsgn_struct.m) | Sets the options and creates a CANlab-style DSGN structure variable, which is used by `LaBGAScore_firstlevel_s2_fit_model.m` to fit the first-level models using CANlab and SPM functions. Should be run from the root directory of your superdataset. Highly study-specific, since conditions, contrasts, and other design choices differ across studies — the version in the repo is a worked example from the `proj_erythritol_4a` study, and needs to be downloaded and adapted to your own code subdataset. **NOTE:** Michael Sun's [`promptDSGN.m`](https://github.com/canlab/CanlabCore/blob/master/CanlabCore/GLM_Batch_tools/promptDSGN.m) script can help you build your DSGN structure! |
| 2 | [`LaBGAScore_firstlevel_s2_fit_model`](https://github.com/labgas/LaBGAScore/blob/main/firstlevel/LaBGAScore_firstlevel_s2_fit_model.m) | Fits the first-level model in 8 steps: (1) checks out and adds the necessary CANlab Github repos to your Matlab path; (2) creates the first-level directory structure (**NOTE:** if using datalad, create the firstlevel directory without subdirectories first via datalad commands); (3) reorganizes your derivatives directories and their content per subject into run-specific subdirectories; (4) extracts/creates noise regressors from the fMRIprep output or CANlab functions (by default the global CSF signal, 24 head motion parameters, and dummy spike regressors) and saves them to the run-specific directories; (5) extracts onsets, durations, and parametric modulators from the `events.tsv` files created by `LaBGAScore_prep_s1_write_events_tsv.m` and saves them to the run-specific directories; (6) plots the design (matrix) of conditions and parametric modulators per run and saves the plots to a model-specific subdirectory; (7) defines and estimates the first-level model and saves the SPM batches; (8) runs diagnostics and publishes an html report by calling `LaBGAScore_firstlevel_s3_diagnose_model.m`. `LaBGAScore_firstlevel_s1_options_dsgn_struct.m` should always be run first — this script checks whether that has been done, and runs it if not. Generic script, tested on Ubuntu 20.04.3 and Windows 10; should not need study-specific modification beyond renaming the generic script names to your study-specific ones at the line numbers indicated in the script. |
| 3 | [`LaBGAScore_firstlevel_s3_diagnose_model`](https://github.com/labgas/LaBGAScore/blob/main/firstlevel/LaBGAScore_firstlevel_s3_diagnose_model.m) | Runs diagnostics on the first-level model and publishes an html report using CANlab's `scn_spm_design_check` function, and saves the variance inflation factors as a `.mat` file. Should be called from `LaBGAScore_firstlevel_s2_fit_model.m` — not intended for standalone use. |

### 5. Run the script

Like all scripts in the LaBGAScore repo, it should be run from the Matlab command line, with the root directory of your superdataset as the working directory! You only need to run the second script, since it calls the first and third under the hood!

```matlab
ery_4a_firstlevel_m1_s2_fit_model
```

Then save the output it wrote to the derivatives and firstlevel subdatasets and the superdataset:

```bash
cd /data/proj_erythritol/proj_erythritol_4a
datalad save -m "ran ery_4a_firstlevel_m1_s2_fit_model.m on all subjects"
```

### 6. Check the first-level results

You can obviously check your results using the SPM GUI, but the scripts also produce output which you can use for quality control and scanning the results:

#### Run-specific plots of the design (matrix)

Saved as individual .png image files in `derivdir/fmriprep/sub-xxx/func/run-x/model_x`.

#### Subject-specific reports of design, collinearity, and results per contrast

Saved as an html report and individual .png image files in `firstlevel/model_x/sub-xxx/diagnostics`.

Further quality control will be done at the group level in the early stage of second-level analysis (see below).

## Second-level analysis

Second-level analysis can be run flexibly using the scripts in the [LaBGAS fork of the CANlab_help_examples Github repo](https://github.com/labgas/CANlab_help_examples/tree/master/Second_level_analysis_template_scripts).

**NOTE:** the most important scripts may move/get copied to LaBGAScore in the future, but not for now, as this would imply losing the ability to merge the LaBGAS fork of this repo with the original CANlab fork, which is ultimately the aim.

The scripts themselves are annotated and documented, hence we refer to their documentation here, but also provide an overview of the typical order in which to execute them [below](#4-check-the-scripts-and-make-study-specific-adaptations-where-needed-1). See also the [CANlab documentation](https://canlab.github.io/batch/) about these second-level batch scripts.

### 1. Create the second level subdataset

```bash
datalad create --description "secondlevel subdataset for erythritol_4a project" -c text2git -d . secondlevel
```

#### Adapt the subdataset configuration

**NOTE:** this may need to be slightly adapted after thoroughly checking which types of files are annexed and which are not with the default configuration.

```bash
cd /data/proj_erythritol/proj_erythritol_4a/secondlevel
nano .gitattributes
```

```
* annex.backend=MD5E
**/.git* annex.largefiles=nothing
```

Add this line:

```
** annex.largefiles=(largerthan=100kb)
```

Change this line:

```
* annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

into:

```
** annex.largefiles=((mimeencoding=binary)and(largerthan=0))
```

`Ctrl+O` to save + enter, `Ctrl+X` to exit nano.

Save this change to the subdataset (and again to the superdataset afterwards):

```bash
datalad save -m "adapted rules for annexing in .gitattributes"
```

#### Make the git (annex) repo writable by our user group

```bash
git config core.sharedrepository world
cd .git
chmod -R g+rwX ./*
```

### 2. Create a model-specific subdirectory corresponding to its firstlevel sibling

**IMPORTANT NOTE:** even if you do not anticipate running several first-level models, always create a model-specific directory to keep things organized and consistent with your firstlevel subdataset, e.g.:

```bash
cd /data/proj_erythritol/proj_erythritol_4a/secondlevel
mkdir model_1_conds_pmods
```

### 3. Download the scripts from the LaBGAS fork of the CANlab_help_examples GitHub repo to your code subdataset and rename them

We want to include all our code in the code subdataset of our datalad superdataset, also the generic scripts, so we consistently download and rename them, like we did with the earlier prep and first-level scripts.

Which scripts you will use depends on your study-specific aims and hypotheses, but there are a number of scripts that should always be run first to set up second-level analysis.

**IMPORTANT NOTES:**

1. always create a model-specific directory in your secondlevel directory in your code subdataset, with the same name as the corresponding model-specific directory in your firstlevel subdataset, and include a model index (`m1`, etc.) in your script name;
2. if you want to run several second-level models based on the same firstlevel data (for example with and without second-level covariates, different types of analyses, etc.), you can make different copies of the same script (with letters as index, e.g. `s6a`, `s6b`) and use options built into the scripts to add a suffix to the saved .mat files!

Example for the first script:

```bash
cd /data/proj_erythritol/proj_erythritol_4a/code
mkdir secondlevel
cd secondlevel
mkdir model_1_conds_pmods
cd ..
datalad download-url -m "downloading second level prep script from LaBGAS fork of CANlab_help_examples Github repo" -d . --path=secondlevel/model_1_conds_pmods/ery_4a_secondlevel_m1_s0_a_set_up_paths_always_run_first.m https://raw.githubusercontent.com/labgas/CANlab_help_examples/master/Second_level_analysis_template_scripts/b_copy_to_local_scripts_dir_and_modify/a_set_up_paths_always_run_first.m
```

**NOTE:** we keep the original name of the script in CANlab_help_examples as the last part of the name of our study-specific script, to keep things consistent and easy to trace back. However, sometimes the length of the name may exceed Matlab's 63-character limit, in which case you need to shorten it.

Save this modification of the subdataset to your superdataset:

```bash
cd /data/proj_erythritol/proj_erythritol_4a
datalad save -m "code subdataset: downloaded second level prep script from Github"
```

### 4. Check the script(s) and make study-specific adaptations where needed

When copying the html address of a script from Github, you need the *raw* version (look for the "raw" button next to the file). See the [`README.md`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/README.md) of this script folder for the full script reference and workflow order.

Scripts under `b_copy_to_local_scripts_dir_and_modify/` are always study-specific and need to be copied and adapted; scripts under `core_scripts_to_run_without_modifying/` are largely generic and typically only need small, targeted adaptations (or none at all). The table below gives the typical order of execution and a short description of each script.

| # | Script | Description |
|---|---|---|
| 1 | [`a_set_up_paths_always_run_first`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/b_copy_to_local_scripts_dir_and_modify/a_set_up_paths_always_run_first.m) | Straightforward, mostly generic script defining the path of, and creating subdirs within, your model-specific second-level directory — the second-level equivalent of `LaBGAScore_prep_s0_define_directories.m`. The only study-specific modification needed is the reference to the correct prep and firstlevel scripts for the model, and to the correct secondlevel script to set the options (`<study_name>_prep_s0_define_directories.m`, `<study_name>_firstlevel_<mx>_s1_options_dsgn_struct.m`, `<study_name>_secondlevel_<mx>_s1_a2_set_default_options.m`). Must be run each time before running any of the following scripts; also creates the output dirs of the second-level GLM when you modify the model name. |
| 2 | [`a2_set_default_options`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/b_copy_to_local_scripts_dir_and_modify/a2_set_default_options.m) | Sets default options for all core secondlevel scripts, organized into one section per script; see comments in code for more info, and the [CANlab documentation in repo](https://github.com/labgas/CANlab_help_examples) and [canlab.github.io](https://canlab.github.io/). Automatically called by `a_set_up_paths_always_run_first`. Always copy to your study/model-specific code subdataset before editing — never edit the checked-in repo copy. |
| 3 | [`prep_1_set_conditions_contrasts_colors`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/b_copy_to_local_scripts_dir_and_modify/prep_1_set_conditions_contrasts_colors.m) | Defines `DAT.conditions`, `DAT.contrasts`, and `DAT.colors` (the core design specification used by every later script) and saves them with DSGN to `image_names_and_setup.mat`. A worked example from one LaBGAS study, not a generic template — study-specific modifications will typically be extensive. |
| 4 | [`prep_1b_prep_behavioral_data`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/b_copy_to_local_scripts_dir_and_modify/prep_1b_prep_behavioral_data.m) | Optional script that attaches behavioral/non-imaging data to DAT: reads .tsv phenotype files, computes z-scored condition-/contrast-level ratings into `DAT.BEHAVIOR`, and populates `DAT.BETWEENPERSON` group/condition/contrast tables for between-subject covariates. Like `prep_1_`, a worked example, not a generic template. Prior to running, make sure the `participants.tsv` file is structured correctly and contains all the subjects. |
| 5 | [`prep_2_load_image_data_and_save`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_2_load_image_data_and_save.m) | Loads first-level beta/con images into `fmri_data_st` objects, runs QC (with plots if requested), z-scores them, and saves both raw and scaled versions (`data_objects.mat` / `data_objects_scaled.mat`). Inspect the QC plots in the html report before proceeding to `prep_3_`. Run using Matlab's `publish` function from your terminal (followed by a manual `datalad save`): `publish('prep_2_load_image_data_and_save.m','outputDir',htmlsavedir)`. **NOTE:** to prevent Matlab from snapping the same figure more than once while publishing (which messes up the report), make sure your cursor is in another window than the one running Matlab during execution of the publish command! |
| 6 | [`prep_3_calc_univariate_contrast_maps_and_save`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_3_calc_univariate_contrast_maps_and_save.m) | Calculates within-person contrast images from `prep_2`'s raw and z-scored condition images, plus an l2norm-rescaled variant, runs QC on all three, and saves them to `contrast_data_objects.mat`. Inspect the QC output before proceeding to `prep_3a_` or other downstream scripts. |
| 7 | [`prep_3a_run_second_level_regression_and_save`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_3a_run_second_level_regression_and_save.m) | Runs second-level regression (voxel-wise via `regress()`, or parcel-wise/robust via `robfit_parcelwise()`) for each contrast/condition in DAT, optionally converts t-maps to BayesFactor maps, optionally runs cross-validated MVPA regression on continuous covariates, and saves the results. Publish `c2a_second_level_regression` afterward for thresholded results. Adaptations should only be necessary if you want to run several second-level GLM analyses within the same model (e.g. with/without covariate control) — in that case, make copies with letters as indices (`s6a`, `s6b`, etc.) and add the necessary lines of code to overwrite the options set in `a2_set_default_options.m`. **NOTE:** in that case, run this script followed immediately by the next one for each model, rather than this script first for all models, to avoid confusion of the `regression_stats_results` variable name in the Matlab workspace, or load the correct `regressions_stats_and_maps.mat` file corresponding to the model you want to display using the next script. Within a model, you can use multiple `prep_3a` scripts to run analyses with different covariates, provided they are included in the DAT structure defined in `prep_1` and `prep_1b`; to include other/new covariates, use a new/different model (start again from `a_set_up_paths_always_run_first`). Run the script using Matlab's `publish` function from your Matlab terminal (followed by a manual `datalad save`). |
| 8 | [`c2a_second_level_regression`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/c2a_second_level_regression.m) | Displays the thresholded second-level regression results generated by `prep_3a_run_second_level_regression_and_save.m`. See that script's documentation for available options. |
| 9 | [`prep_3c_run_SVMs_on_contrasts_masked`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_3c_run_SVMs_on_contrasts_masked.m) | Runs second-level support vector machine classification for each contrast in `DAT.contrasts`, plots uncorrected result montages, and saves the results. Publish `c2_SVM_contrasts_masked` afterward for thresholded results (after bootstrapping/searchlight analysis). Adaptations should only be necessary if you want to run several second-level SVM analyses within the same model (e.g. with/without the built-in scaling options, or with different cross-validation settings) — in that case, make copies with letters as indices (`s6a`, `s6b`, etc.) and add the necessary lines of code to overwrite the options set in `a2_set_default_options.m`. |
| 10 | [`c2_SVM_contrasts_masked`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/c2_SVM_contrasts_masked.m) | Displays the thresholded SVM results generated by `prep_3c_run_SVMs_on_contrasts_masked.m`. See that script's documentation for available options. |
| 11 | [`prep_3g_create_fmri_data_runwise_contrast_object`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_3g_create_fmri_data_runwise_contrast_object.m) | Creates and saves an `fmri_data_st` object of runwise contrast images (difference between two specified condition beta images per run) from firstlevel results, matched to QC-passing runs with non-missing behavioral outcomes from the phenotype file, after sanity-checking that runs/ratings/image names line up. |
| 12 | [`prep_3f_create_fmri_data_single_trial_object`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_3f_create_fmri_data_single_trial_object.m) | Creates and saves an `fmri_data_st` object of single-trial con images and per-trial variance inflation factors (VIFs) from firstlevel results, attaches single-trial ratings, excludes high-VIF trials, and runs sanity checks that trial/subject identifiers line up. Inspect the VIF plots in the html report — high-VIF trials indicate multicollinearity with noise regressors and are auto-excluded. |
| 13 | [`c2f_run_MVPA_regression_single_trial`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/c2f_run_MVPA_regression_single_trial.m) | Runs MVPA regression (default cross-validated PCR, adaptable to PLS or other algorithms) on a continuous outcome using the `prep_3f_` single-trial object, via either CANlab's classic `predict` function or Bogdan's object-oriented `ooFmriDataObjML` toolbox; supports optional bootstrapping, permutation testing, and source reconstruction ("structure coefficients"). |
| 14 | [`c2g_run_multivariate_mediation_single_trial`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/c2g_run_multivariate_mediation_single_trial.m) | Runs multivariate (PDM) mediation analysis on a continuous outcome using the `prep_3f_` single-trial object, via CANlab's `multivariateMediation` function; supports optional bootstrapping and source reconstruction. |
| 15 | [`prep_4_apply_signatures_and_save`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/prep_4_apply_signatures_and_save.m) | Calculates selected CANlab signature responses (via `apply_all_signatures()`) for all conditions/contrasts in DAT and saves them to `DAT.SIG_conditions`/`DAT.SIG_contrasts`; additionally computes NPS subregion responses via `apply_nps()` if `'nps'` is among `keyword_sigs`. Choose which signatures to analyze by looking at [this site](https://sites.google.com/dartmouth.edu/canlab-brainpatterns/multivariate-brain-signatures). |
| 16 | [`d_signature_responses_generic`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/d_signature_responses_generic.m) | Plots and tests significance of selected signature responses (from `prep_4_apply_signatures_and_save.m`) for conditions/contrasts in DAT — for individual signatures or, by default, all signatures in `keyword_sigs` — via `plugin_signature_condition_contrast_plot`. |
| 17 | [`d10_signature_riverplots`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/d10_signature_riverplots.m) | Generates cosine-similarity riverplots of signature responses (loaded via `load_image_set`) against condition and contrast images in DAT, showing statistically significant associations only. Unlike `d_signature_responses_generic.m`, works only on signature groups, not individual signatures. |
| 18 | [`h_signature_responses_group_diff`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/h_signature_responses_group_diff.m) | Runs a two-sample t-test for selected signature responses per contrast, plotting group differences and printing between-group test statistics for each (group membership from `DAT.BETWEENPERSON.group`, binarized via median split if continuous). Note: unlike other Group 2 scripts, this script does not call `a_set_up_paths_always_run_first` or reload DAT/DATA_OBJ from saved .mat files itself; it assumes these are already in the workspace from a previous script run earlier in the same Matlab session (e.g. `prep_4_apply_signatures_and_save.m`). Run this with a `publish` command. |
| 19 | [`e1_corr_patterns`](https://github.com/labgas/CANlab_help_examples/blob/master/Second_level_analysis_template_scripts/core_scripts_to_run_without_modifying/e1_corr_patterns.m) | Calculates, thresholds, and plots pairwise searchlight correlation maps between all condition or contrast images in DAT via CANlab's `searchlight_correlation()` function, optionally masked/restricted to an atlas. Independent of `prep_4`/signatures. |

### 5. Run script(s)

**NOTE:** like every other script in this workflow, second-level scripts are run from the Matlab terminal rather than from the Linux command line. Matlab's `publish` function, used by several of them, does not behave the same way when Matlab is driven from a shell. If the script produces any output, perform a manual `datalad save` to the relevant subdataset and superdataset after running it.

**Record which version of the dependencies you ran against.** The scripts in your `code` subdataset are frozen once you copy them, but CanlabCore, the LaBGAS fork of CANlab_help_examples and the other repos under `/data/master_github_repos` keep changing. Instead of the bare `publish` call given in each script header, use:

```matlab
LaBGAScore_prov_publish('<projname>_secondlevel_m<M>_s<N>_<scriptname>', htmlsavedir)
```

This publishes the html report at exactly the same figure resolution as before (it does not set `maxWidth`/`maxHeight`, which would permanently shrink the saved `.png` files), adds a small stylesheet so the report reads correctly on a laptop as well as a large desktop screen, and adds a Provenance section recording the commit of every dependency the script reaches, any uncommitted local changes to files it uses, and the MATLAB and SPM versions. A machine-readable copy is written to `<model>/results/notes/`, small enough to stay in git rather than git-annex.

For analyses you have **already** run, `LaBGAScore_prov_resolve_retrospective('/data/proj_<yourstudy>')` reconstructs the same record after the fact and writes it alongside the existing reports without modifying them. See [`clean/README_provenance.md`](https://github.com/labgas/LaBGAScore/blob/main/clean/README_provenance.md) in LaBGAScore.

## Publish your dataset on GIN and/or Github

**NOTE:** this is work in progress — see this [issue](https://neurostars.org/t/configuring-remotes-on-and-pushing-to-github-gin/20244) opened on Neurostars, and [this rapidly evolving section of the Datalad handbook](http://handbook.datalad.org/en/latest/basics/basics-thirdparty.html), particularly [this walkthrough](http://handbook.datalad.org/en/latest/basics/101-139-gin.html).

**NOTE:** it is highly recommended to configure your GIN siblings and repos early on and push fairly frequently, rather than waiting until your superdataset is complete, after which pushing will take a long time and be more sensitive to errors.

### [GIN](https://gin.g-node.org/G-Node/Info/wiki)

#### Automated approach

##### 1. Create an SSH key pair and a personal access token (PAT) for authentication

Follow [these instructions](https://docs.github.com/en/authentication/connecting-to-github-with-ssh) for SSH, and [these](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) for PAT; you can use the same pair for GIN and Github, but need a separate one for each computer you work on (LaBGAS server, your own laptop, etc.).

##### 2. Create siblings and GIN repos for superdataset and subdatasets

Run the following [`datalad create-sibling-gin`](http://docs.datalad.org/en/stable/generated/man/datalad-create-sibling-gin.html) command from your superdataset:

```bash
datalad create-sibling-gin labgas/proj_erythritol_4a -s gin -r --private
```

This automatically creates GIN siblings locally, as well as remote repos on GIN!

##### 3. Delete the sibling and GIN repo for sourcedata subdataset

The sourcedata subdataset contains non-pseudonymized data and should hence NEVER BE PUBLISHED!

**Delete the repo on GIN:** go to `labgas/proj_erythritol_4a-sourcedata`, click "Settings" (top right), scroll down and click "Delete this repository", type the name of the repo to confirm.

**Remove the sibling from the sourcedata subdataset:**

```bash
datalad siblings -s gin remove
```

##### 4. Add the url of the subrepos for each of the corresponding subdatasets in your superdataset

Run the following command from your superdataset — this is an example for the code subdataset only:

```bash
datalad subdatasets --contains code --set-property url https://gin.g-node.org/labgas/proj_discoverie_code
```

##### 5. Push recursively from your superdataset to GIN

```bash
datalad push --to gin -r
```

**NOTE:** no need to worry about the error for the sourcedata subdataset — it was deleted in the previous step on purpose!

**NOTE:** you may need additional flags in case of push errors, such as `--data anything -f all`. See also the ongoing [Neurostars thread](https://neurostars.org/t/datalad-push-to-gin-errors/24051) about pushing errors to GIN.

##### 6. Use datalad drop (wisely) in large subdatasets to save space

**NOTE:** only do this if you know what you are doing — RTFM and/or talk to Lukas!

Once your BIDS and derivatives (and sometimes also firstlevel) subdatasets are complete, they typically take up a lot of space. However, once pushed to GIN, a remote copy is available, which allows using [`datalad drop`](https://docs.datalad.org/en/stable/generated/man/datalad-drop.html) to remove files entirely on the server. They no longer take up space, but the symlinks are kept, and [`datalad get`](https://docs.datalad.org/en/stable/generated/man/datalad-get.html) allows retrieving them from the GIN remote if needed.

This process is nicely described [here](https://handbook.datalad.org/en/latest/basics/101-136-filesystem.html#removing-annexed-content-entirely) in the Datalad handbook, and contrary to the sourcedata subdataset situation described above, the `--reckless availability` option is not needed since a verified remote copy exists on GIN.

As is clear from the [Datalad handbook](https://handbook.datalad.org/en/latest/basics/101-136-filesystem.html#removing-annexed-content-entirely) and the [`datalad drop` man page](https://docs.datalad.org/en/stable/generated/man/datalad-drop.html), you can choose to drop the entire subdataset, or selected files only.

#### Manual approach — DEPRECATED

##### 1. Create an SSH key pair and a personal access token for authentication

Follow [these instructions](https://docs.github.com/en/authentication/connecting-to-github-with-ssh) for SSH, and [these](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) for PAT; you can use the same pair for GIN and Github, but need a separate one for each computer you work on (LaBGAS server, your own laptop, etc.).

##### 2. Add a GIN "superrepo" as a sibling (and common data source) to your superdataset

**NOTE:** this describes the [manual workflow](http://handbook.datalad.org/en/latest/basics/101-139-hostingservices.html#how-to-add-a-sibling-on-a-git-repository-hosting-site-the-manual-way); a new automatic workflow is described above (using `datalad create-sibling-xxx` commands, which can be combined with [`datalad siblings configure`](http://docs.datalad.org/en/latest/generated/man/datalad-siblings.html) to adapt configuration flexibly).

After creating your empty superrepo on GIN, run the following command from your superdataset:

```bash
datalad siblings add -d . --name gin --pushurl git@gin.g-node.org:/labgas/proj_discoverie.git --url https://gin.g-node.org/labgas/proj_discoverie --as-common-datasrc gin-common
```

Then make sure that annex is supported for this sibling by running:

```bash
git config --unset-all remote.gin.annex-ignore
```

##### 3. Add a GIN "subrepo" as a sibling (and common data source) for each subdataset

**NOTE:** we do **NOT** do this for the sourcedata subdataset, since we do not want it to be "datalad gettable", even for people with access to the private GIN repo!

Essentially, repeat step 2 above in a slightly simplified way for each subdataset. After creating your empty subrepo on GIN, run from your subdataset:

```bash
datalad siblings add -d . --name gin --pushurl git@gin.g-node.org:/labgas/proj_discoverie_code.git --url https://gin.g-node.org/labgas/proj_discoverie_code --as-common-datasrc gin-common-code
```

Then make sure that annex is supported for this sibling:

```bash
git config --unset-all remote.gin.annex-ignore
```

##### 4. Add the url of the subrepos for each of the corresponding subdatasets in your superdataset

Run the following command from your superdataset:

```bash
datalad subdatasets --contains code --set-property url https://gin.g-node.org/labgas/proj_discoverie_code
```

##### 5. Push recursively from your superdataset to GIN

```bash
datalad push --to gin -r
```

**NOTE:** no need to worry about the error for the sourcedata subdataset — no GIN sibling was created for it, on purpose.

##### 6. Clone the entire superdataset or subdatasets wherever you like

**NOTE:** this may also be easier with the new approach described above, but that has not been tested.

**Using datalad commands:**

```bash
datalad clone https://gin.g-node.org/labgas/proj_discoverie
```

If you want a subdataset with annexed files downloaded to your computer:

```bash
datalad get -d . BIDS
```

However, this generates errors due to the fact that the configuration for the siblings of the cloned dataset differs from the original dataset that was pushed to GIN. This can be solved by manually editing `.git/config` of the cloned (sub)dataset to match the info on the remotes from the original dataset on the server — specifically, copy the `annexurl`, `pushurl`, and `annex-uuid` lines of the `gin` remote of the original local dataset to both the `origin` and `gin-common-*` remotes of the cloned dataset, and, importantly, set `annex-ignore` to `false`.

Here is an example for the BIDS subdataset:

```bash
cd BIDS
nano .git/config
```

```ini
[core]
        repositoryformatversion = 0
        filemode = true
        bare = false
        logallrefupdates = true
[remote "origin"]
        url = https://gin.g-node.org/labgas/proj_discoverie_BIDS
        fetch = +refs/heads/*:refs/remotes/origin/*
        annexurl = git@gin.g-node.org:/labgas/proj_discoverie_BIDS.git
        pushurl = git@gin.g-node.org:/labgas/proj_discoverie_BIDS.git
        annex-uuid = de666db2-5da4-475f-a8ce-c34adf5e09e0
        annex-ignore = false
[branch "master"]
        remote = origin
        merge = refs/heads/master
[annex]
        uuid = f6c28cb9-43ae-4b71-b21b-6bda8c91ada8
        version = 8
[filter "annex"]
        smudge = git-annex smudge -- %f
        clean = git-annex smudge --clean -- %f
[remote "gin-common-bids"]
        url = https://gin.g-node.org/labgas/proj_discoverie_BIDS
        fetch = +refs/heads/*:refs/remotes/gin-common-bids/*
        annexurl = git@gin.g-node.org:/labgas/proj_discoverie_BIDS.git
        pushurl = git@gin.g-node.org:/labgas/proj_discoverie_BIDS.git
        annex-uuid = de666db2-5da4-475f-a8ce-c34adf5e09e0
        annex-ignore = false
```

This enables `datalad get`:

```bash
datalad get .
```

However, there should be a more elegant solution, especially as this needs to be done for the superdataset and all subdatasets (at least the ones with annexed content) if you want to `datalad get` the latter into the superdataset.

**Using gin commands:**

An alternative to `datalad clone` is to use [gin commands](https://gin.g-node.org/G-Node/Info/wiki/GIN+CLI+Help), specifically `gin get` (the gin equivalent of `datalad`/`git clone`). `gin get-content` gets all the (annexed) files in your local repo and is hence the equivalent of `datalad get`.

This works perfectly for the subdatasets, and contrary to `datalad clone`, the annex config is correctly preserved! However, when `gin get`ting the superdataset, `gin get-content` for the subdatasets does not immediately work. There is a fairly easy workaround using gin commands and minor edits to `.git/config` of the superdataset:

```bash
cd ..   # superdataset root
rm -r BIDS
gin get labgas/proj_discoverie_BIDS
mv proj_discoverie_BIDS BIDS
nano .git/config
```

Add info on the BIDS submodule:

```ini
[submodule "BIDS"]
        active = true
        url = https://gin.g-node.org/labgas/proj_discoverie_BIDS
        path = ./BIDS
```

`Ctrl+O`, `Ctrl+X`.

**NOTE:** it is not clear whether this step is really needed, since the information on submodules/subdatasets is already correctly stored in the `.gitmodules` file under the root of the superdataset, but it definitely does not harm to make these two files consistent.

```bash
cd BIDS
gin get-content .
```

**NOTE:** this last step is only needed for subdatasets with annexed content, hence not for `code` or `mriqc`, for example.

### Github

**NOTE:** Github does not support large files nor annexed content, so it is less convenient than GIN, but it is more widely known, so we want our dataset — and particularly the code subdataset — available on Github as well, preferably in a clonable way (through a link to the common data source on GIN).

#### 1. Add a Github "superrepo" as a sibling to your superdataset

Like for GIN, currently only the manual approach works reliably, so create an empty (i.e. no license, no README, no .gitignore) data repo on Github first, then run the following command from your superdataset:

```bash
datalad siblings add -d . --name github --url https://github.com/labgas/proj_discoverie.git
```

#### 2. Add a Github "subrepo" as a sibling to your code subdataset

Create an empty repo on Github first, then run the following command from your code subdataset:

```bash
cd code
datalad siblings add -d . --name github --url https://github.com/labgas/proj_discoverie_code.git
```

#### 3. Push recursively from your superdataset to Github

```bash
cd ..
datalad push --to github -r
```

**NOTE:** no need to worry about the errors for most subdatasets — no Github sibling was created for them on purpose, since they are all on GIN anyway, and we do not want them public prior to publication (private repos are not free on Github, contrary to GIN).

**NOTE:** cloning works from Github as from GIN, including the workarounds for the remote config explained above, at least while we are awaiting a more elegant solution — this is being worked on with the Datalad developers.
