# EELS-ZEST GUI Guide

Launch the app from the toolbox root:

```matlab
addpath('EELS-ZEST');
app = run_eels_zest;
```

## Data

Use the Data page to:

- load `.mat` or `.dm3` spectra and datacubes,
- inspect energy-filtered maps and integrated spectra,
- align the ZLP to `0 meV`,
- apply rank-limited SVD/PCA denoising,
- select and export a spatial ROI,
- truncate fitted sparse phonon source intensity and reconvolve it with the fitted ZLP,
- multiply the positive-loss data by the STA factor,
- stack point and rectangular-ROI spectra in a nonblocking window,
- run robust ROI, x-line, or datacube multi-peak fitting,
- load a saved fit result,
- run datacube Richardson-Lucy conversion after fitting.

For full datacubes, crop the energy axis before fitting. `[-250, 250] meV` is
a practical default for the bundled example.

**Send signal to Data** returns the background-subtracted phonon reconstruction,
not an additive ZLP peak. Its per-pixel fitted-ZLP kernels, Stage 1/Stage 2
temperature maps, validity masks, and result provenance are transferred in the
Data metadata. RL uses the displayed cutoff spectrum and the transferred kernel;
STA uses the transferred temperature maps. These fields survive `EELS` export
and reload.

The same operations can be called without dialogs:

```matlab
options = struct('mode', 'Manual cutoff', 'cutoff_meV', 100);
app.sendSignalToData('stage2', options);
app.multiplySTA('stage2');
```

## Stage 1

Stage 1 fits the ZLP and sparse vibrational signal. Run a selected pixel first
to inspect:

- fitted ZLP and ZLP components,
- reconstructed Stokes/anti-Stokes signal,
- sparse mode-energy markers and source intensities,
- residual structure,
- R2 and NRMSE.

After the single-pixel result is satisfactory, run the complete datacube. The
default v2.1 datacube path uses ZLP-only clustering, class templates, and
an edge-aware Stage 1 spatial Huber loss. The additional map views show ZLP
cluster, membership, outliers, template deviation, and ZLP spatial loss.
Serial mode supports checkpoints. Parallel mode can use multiple workers.

## Stage 2

Stage 2 reuses Stage 1 and performs attention-guided thermometry. Inspect:

- hard and soft source-intensity partitions,
- paired attention evidence,
- BE gain/loss linearity,
- effective hard-pair count,
- temperature cost valley,
- uncertainty and reliability.

Only hard attention pairs constrain temperature. Soft source intensities still
participate in the full-spectrum reconstruction.

For datacubes, the default map is **T data valid**. Compare it with:

- **T raw**, the unchanged pixel-local fit;
- **T data valid**, which excludes boundary and unconstrained valleys;
- **T joint**, a ZLP-guided view that keeps every data-valid value fixed;
- **T joint uncertainty**;
- **Spatial inferred**, which marks data-invalid pixels estimated from valid
  same-side neighbors.

There is no temperature-difference smoothing term in v2.1. Sharp interfaces
are scored from ZLP features; Stage 1 coupling and temperature inference are
disabled across a marked interface edge. Inferred pixels do not feed back into
data-valid anchors. Inspect **ZLP interface score** and **ZLP sharp interface**
with the temperature maps.

Full-datacube GUI runs save numerical MAT and summary outputs before any
optional graphics. Automatic large diagnostic-figure export is disabled by
default. Template polish, each spatial color pass, and final Stage 2 report
live progress; Cancel is checked between bounded parallel batches and writes a
checkpoint.

## Phonon-Intensity Analysis

The Phonon Intensity page interpolates fitted Stokes and anti-Stokes
phonon-mode source intensities onto a signed energy axis. It supports:

- source-intensity integral maps,
- selected energy-window export,
- source-intensity NMF,
- source-intensity clustering.

Negative energy denotes gain/anti-Stokes source intensity; positive energy
denotes loss/Stokes source intensity.

## Recommended Sequence

1. Load data.
2. Align ZLP when needed.
3. Apply SVD/PCA denoising.
4. Crop the energy range.
5. Run and inspect one Stage 1 pixel.
6. Run the complete datacube with **Hierarchical ZLP + ZLP interface guidance** enabled.
7. Inspect ZLP classes, soft membership, template deviation, Stage 1 spatial
   loss, and the sharp-interface mask.
8. Inspect hard-pair thermometry and compare Raw, Data-valid, and ZLP-guided maps.
9. Export the fit result and sparse source-intensity products.

All export controls save to files. The optional base-workspace checkbox is off
by default. `Export data` writes the Data page exactly as displayed, including
its processing history, as an `EELS` object; it does not trigger hidden
background subtraction or reconstruction.

## Author, Citation, And License

EELS-ZEST is developed by Ruilin Mao, International Center for Quantum
Materials (ICQM), Peking University. Citation metadata is provided in
`../CITATION.cff`. The software is licensed under `GPL-3.0-only`; see
`../LICENSE`.

## Terminology And MATLAB Field Names

This guide uses **sparse phonon-mode source intensity** for the nonnegative
coefficients fitted before convolution with the ZLP kernel. For compatibility
with existing MATLAB result files, the code abbreviates **phonon-mode
intensity** as **PM**, including fields such as `PM_stokes` and
`PM_antistokes`.
