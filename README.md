# EELS-ZEST

**Zero-loss Estimation and Sparse Thermometry**

EELS-ZEST is a MATLAB toolbox for STEM-EELS spectra and datacubes. It combines
zero-loss peak estimation, sparse-convolution signal reconstruction, sparse
phonon-mode source-intensity analysis, and attention-guided Bose-Einstein
thermometry.

## What Version 2.1 Does

The workflow separates spectral reconstruction from temperature inference:

1. **Stage 1 reconstructs spectral shape.**
   It fits the complete ZLP and uses the fitted total ZLP, centered at `0 meV`
   on an aligned energy axis, as the convolution kernel for sparse Stokes and
   anti-Stokes phonon-mode source intensities.
2. **Datacube Stage 1 learns spatially coherent ZLP families.**
   A ZLP-only probe is clustered with deterministic robust k-means. Each class
   receives a shared ZLP template. Soft membership continuously blends template
   targets and admissible parameter ranges, avoiding hard changes at class
   boundaries. Stage 1 also applies an
   edge-aware eight-neighbor Huber loss on peak center, width, shape, and
   asymmetry. ZLP amplitude and phonon source intensities remain pixel-local.
3. **Stage 2 estimates temperature from selected evidence.**
   Paired modes with strong attention constrain the Bose-Einstein temperature.
   Remaining modes stay in the full-spectrum fit as independent soft
   source-intensity coefficients, so mixed phases, tails, and non-equilibrium
   signal can be fitted without controlling the reported temperature.
4. **ZLP differences define sharp interfaces and safe one-way inference.**
   A robust ZLP-only edge score marks sharp interfaces and disables Stage 1
   spatial coupling across them. Data-valid temperatures are never smoothed or
   moved. A data-invalid pixel may receive a marked, ZLP-guided estimate only
   from at least two valid same-side neighbors; inferred pixels never feed back
   into valid anchors.
5. **Diagnostics explain the result.**
   Outputs include fit residuals, sparse source-intensity distributions,
   hard/soft mode partitions, BE linearity, the temperature cost valley,
   uncertainty, reliability, ZLP clusters, ZLP interface score/mask,
   raw/data-valid/ZLP-guided temperature, and inference masks for datacubes.

The default hard-pair score is:

```text
pair_score(EM) = sqrt(att(+EM) * att(-EM))
```

Pairs with `pair_score >= 0.35` form the main hard-BE set. A `0.10` floor can
supplement a small set. Temperature diagnostics use the selected hard pairs
only.

## Requirements

- MATLAB
- Optimization Toolbox
- Parallel Computing Toolbox for multi-worker datacube fitting

## Quick Start

Launch the GUI:

```matlab
addpath('EELS-ZEST');
app = run_eels_zest;
```

Fit one spectrum:

```matlab
settings = struct();
settings.run_stage2 = true;
settings.crop_range = [-250, 250];
settings.stage2_temperature_algorithm = 'grid_scan';
settings.stage1_hierarchical_enabled = true;
settings.stage2_spatial_map_enabled = true;

result = eels_fit_spectrum(energy_meV, spectrum, settings);
```

Fit the bundled example datacube:

```matlab
settings = struct();
settings.run_stage2 = true;
settings.crop_range = [-250, 250];
settings.datacube_sample_strategy = 'all';
settings.datacube_parallel_enabled = true;
settings.datacube_parallel_workers = 16;
settings.test_output_dir = fullfile(pwd, 'eels_zest_output');

cuberes = eels_fit_datacube( ...
    fullfile('examples', 'data', 'example_data.mat'), settings);
```

Ready-to-run scripts:

```matlab
run(fullfile('examples', 'zest_datacube_full_stage1.m'));
run(fullfile('examples', 'zest_datacube_full_stage2_parallel.m'));
```

## Recommended Datacube Workflow

1. Load a `.mat` or `.dm3` datacube.
2. Align the ZLP to `0 meV` when needed.
3. Apply rank-limited SVD/PCA denoising.
4. Crop to the useful energy interval, normally `[-250, 250] meV`.
5. Run the full datacube with the hierarchical model enabled, then inspect ZLP
   clusters, membership, outliers, template deviation, residual, R2, and NRMSE.
6. Inspect hard-pair BE linearity and the raw temperature cost valley.
7. Use Data-valid as the primary temperature map. Compare Raw and the optional
   ZLP-guided result, and inspect the sharp-interface and inferred masks before
   interpreting spatial trends.

The v2.1 GUI reports template polish, each spatial color pass, and final
Stage 2 separately. Cancel is checked between bounded parallel batches.
Full-cube numerical results are saved before optional plots, and automatic
large diagnostic-figure export is disabled by default.

For a saved v2.1 Stage 1 datacube result:

```matlab
stage2_settings = struct();
stage2_settings.run_stage2 = true;
stage2_settings.stage2_temperature_algorithm = 'grid_scan';

stage2_cuberes = ...
    EELS_PearsonSparse_FullFit.promote_datacube_stage1_to_stage2( ...
        stage1_cuberes, stage2_settings);
```

Reusing a pre-v2.1 Stage 1 result preserves its original ZLP fit; hierarchical
Stage 1 constraints cannot be added retroactively. Run `eels_fit_datacube`
again when the new ZLP clustering and spatial loss are required.

## Main Entry Points

| File | Purpose |
| --- | --- |
| `run_eels_zest.m` | Launch the graphical interface |
| `eels_fit_spectrum.m` | Fit one spectrum |
| `eels_fit_datacube.m` | Fit a spectrum image from a file |
| `EELS_PearsonSparse_FullFit.m` | Core Stage 1/Stage 2 solver and diagnostics |
| `EELS_ZEST_App.m` | GUI implementation |
| `EELS_ZEST_DataOps.m` | sparse-source cutoff, reconvolution, STA, and EELS export |
| `EELS_ZEST_PeakFit.m` | robust Gaussian, Lorentzian, and Pseudo-Voigt multi-peak fitting |
| `examples/` | Runnable scripts |
| `docs/example_data_figures/` | Figure data, CSV panels, SVG/PNG/TIFF figures, and reproduction scripts |
| `index.html` | Project documentation webpage |

## Input Data

Supported inputs include:

- `.mat` files containing an object or struct with `.data` and `.ene`,
- `.mat` files with recognizable datacube and energy-axis variables,
- `.dm3` files loaded through `DM3Import`.

Preferred datacube orientation:

```text
Ny x Nx x Ne
```

where `Ne` is the number of energy channels.

## Data Tools In Version 2.1

- **Send signal to Data** truncates Stage 1 or Stage 2 sparse phonon source
  intensities using a manual energy, a Materials Project phonon JSON cutoff,
  or the Stage 2 hard-attention edge plus a tail. It then reconvolves the
  retained source with the same fitted-ZLP kernel, so the reconstructed
  background-subtracted phonon spectrum retains its natural convolution tail.
  The fitted-ZLP kernel cube and Stage 1/Stage 2 temperature maps are transferred
  with the Data state and retained in `EELS.info`; RL therefore acts on the
  displayed cutoff spectrum, while STA can reuse the corresponding temperature.
- **Multiply STA** keeps the positive-loss side and applies
  `I(E)*E/2*[1-exp(-E/(kB*T))]` using a manual temperature, Stage 1 soft
  temperature, or Stage 2 temperature with a Stage 1 fallback.
- **Stack Spectra** uses a nonblocking map-selection window for points and
  rectangles. **Multi-peak fit** supports ROI, x-line, and full-cube robust
  Gaussian, Lorentzian, and Pseudo-Voigt fits.
- **Export data** writes the currently displayed processing state as an
  `EELS` object. Other export controls save MAT/CSV and, where applicable,
  FIG/SVG/PNG/TIFF files; base-workspace assignment is optional.

## Core Outputs

Single-spectrum results include:

- `raw_data`, `fitted_spectrum`, `residual`
- `zlp_only`, `zlp_components`, `vib_only`
- sparse Stokes and anti-Stokes source-intensity arrays
- `metrics`, `stage2`, `T_fit`, `T_uncertainty`, `reliability_score`

Datacube results include:

- `R2_map`, `NRMSE_map`, `fit_status_map`
- `T_map_raw`, `T_map_data_valid`, `T_map_joint`, `T_map_spatial_refined`
- `T_uncertainty_map`, `T_uncertainty_joint_map`, `T_spatial_inferred_mask`
- `cluster_label_map`, `cluster_membership_map`, `cluster_outlier_mask`
- `template_deviation_map`, `zlp_spatial_loss_map`, `T_joint_shift_map`
- `zlp_interface_score_map`, `zlp_interface_mask`, and directional edge maps
- `reliability_score_map`, `attention_score_map`
- `warning_map`, `bad_pixel_mask`, `interface_mask`

## Example Figure Reproduction

The bundled figure case uses `examples/data/example_data.mat`.

```matlab
cd(fullfile('docs', 'example_data_figures'));
export_example_figure_data;
plot_example_figures;
```

The directory also contains the exported panel CSV files and figure-ready MAT
data, so plotting can be repeated without rerunning the full fit.

## Author And Citation

EELS-ZEST is developed by **Ruilin Mao**, International Center for Quantum
Materials (ICQM), Peking University.

If EELS-ZEST contributes to research, please cite the software using
`CITATION.cff`. GitHub can render this metadata directly through
**Cite this repository**.

## License

Copyright (C) 2026 Ruilin Mao.

EELS-ZEST is free software licensed under the
**GNU General Public License version 3 only** (`GPL-3.0-only`). See `LICENSE`
for the complete terms. Research users are encouraged to notify the author
and cite the software in resulting publications.

## Documentation

- `README_zh-CN.md`: Chinese quick guide
- `docs/GUI_TOOLBOX_GUIDE.md`: GUI workflow
- `docs/LOSS_AND_ATTENTION_THEORY_v1.3.md`: loss and attention model
- `docs/theory.html`: rendered equations and method details
- `index.html`: browsable project guide

## Terminology And MATLAB Field Names

This documentation calls the fitted coefficients **sparse phonon-mode source
intensities**. They are nonnegative source-space weights before convolution
with the fitted ZLP kernel and should not be interpreted as a directly
calibrated phonon DOS.

For compatibility with existing MATLAB result files, the code abbreviates
**phonon-mode intensity** as **PM**. The main fields are `PM_stokes`,
`PM_antistokes`, and `EM_final`.
