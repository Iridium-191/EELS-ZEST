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
- load a saved fit result,
- run datacube Richardson-Lucy conversion after fitting.

For full datacubes, crop the energy axis before fitting. `[-250, 250] meV` is
a practical default for the bundled example.

## Stage 1

Stage 1 fits the ZLP and sparse vibrational signal. Run a selected pixel first
to inspect:

- fitted ZLP and ZLP components,
- reconstructed Stokes/anti-Stokes signal,
- sparse mode-energy markers and source intensities,
- residual structure,
- R2 and NRMSE.

After the single-pixel result is satisfactory, run the complete datacube.
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
6. Run the Stage 1 datacube.
7. Run Stage 2 using the saved Stage 1 result.
8. Inspect temperature diagnostics and maps.
9. Export the fit result and sparse source-intensity products.

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
