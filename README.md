# EELS-ZEST 2.2.0

Zero-loss Estimation and Sparse Thermometry. MATLAB interfaces with C++ Stage1/Stage2 numerical kernels. Released 2026-11-08.

## Start

Open a fresh MATLAB session, change into this directory, and run:

```matlab
app = run_eels_zest;
```

Use `run_reference_tuner` for a standalone tuner or `run_eels_twostage_toolbox` for the alternate toolbox layout. Both GUIs, including dialogs, menus, parameter descriptions and tuner status text, are in English. MATLAB/OS-generated dialogs and errors follow the installed system language.

## Release defaults

| Setting | Default |
| --- | --- |
| Tuner objective | Temperature mean absolute error, in K |
| Maximum residual / baseline | 2.0 |
| Maximum candidate evaluations | 500, including baseline |
| Parallel workers | 8, limited by local capacity |
| Save checkpoint / resume | Off; explicitly enable when needed |
| Search algorithm | Grouped + joint simulated annealing |
| SA temperature | 0.1, dimensionless |
| Initial search radius | 0.12 in normalized parameter space |
| Map sampling | 32 training + 8 heldout points |

500 means candidate evaluations, not 500 full particle-swarm generations. Each candidate fits the training spectra. Small groups use no more workers than spectra. Stage1/2 and the tuner share the process pool.

## Reference tuning

1. Import data in Data, then use Reference tuning to import reference spectra or sample a map. Training points must have the supplied known temperature.
2. Choose parameters and run. Up to 87 numeric controls are available, with dependencies shown in the full catalog. Model choices stay fixed within a run.
3. Inspect training and heldout errors, residuals, zero-peak shapes and temperature valleys. Heldout points never select the strategy.
4. Click Apply to Stage1/2 to activate the best eligible strategy. Import/run alone does not apply it. Without an applied strategy, fitting uses current/default settings.
5. Import `fit_strategy.mat` to continue from a saved strategy. New spectra are always independently fitted; a reference shape or known temperature is not imposed on them.

The outer loss is `mean(abs(raw_valley_temperature_K - known_temperature_K))`. There is no valley-height term, residual soft penalty, logarithmic normalization or worst-point weighting. Residuals remain a separate eligibility gate: any training point above 2 times its initial NRMSE prevents that candidate from becoming the applied best. Failed fits are marked and receive a solver sentinel score. Inner spectral fitting weights and regularizers remain active. Known temperatures are outer tuning labels only; Stage1/Stage2 infer temperature from each spectrum. Configuration import keeps supported parameters and runtime context.


Grouped PSO, all-selected PSO, surrogate/PSO/local refinement and local GP remain available for comparison. Use identical starting settings, spectra and evaluation budgets when comparing solvers. No superiority or cross-material calibration claim is made for this release.

## Outputs

Tuning runs automatically create unique folders under `tuning_runs/`, including strategy MAT/JSON, complete report, candidate/point CSVs, parameter inventory, search stages, plots and benchmark datasets. Tuning exports are independent of the datacube Save checkpoint checkbox. Existing beta tuning runs are not bundled or overwritten.

The dataset records objective definitions, spectral features, actions, rewards, eligibility flags, per-point outcomes, optimizer metadata and code hashes. Keep all trials from one physical acquisition together when splitting benchmark datasets. Historical composite losses cannot be directly compared with the new loss in K.

## Requirements

- MATLAB R2024b on Windows was used for release checks.
- Optimization Toolbox for the inner fitting solvers; Global Optimization Toolbox for SA/PSO/surrogate/pattern search.
- Parallel Computing Toolbox for process workers. Local GP is provided as a toolbox-independent outer search, but inner fitting requirements still apply.
- Some optional image-processing/ROI and analysis tools require their corresponding MATLAB toolboxes or optional external tools.

`eels_native_v22_mex.mexw64` is the included Windows binary. Other systems must build their own binary; see [native setup](README_NATIVE.md). This is a MATLAB source release, not a standalone executable.

## Checks

Run `addpath('tests'); test_temperature_only_loss; test_grouped_sa` for inexpensive numerical-policy checks. `run_release_smoke(outputDirectory)` checks GUI defaults and language. To include a short native tuning run, supply a second argument pointing to a MAT file with `energy_meV`, `spectra` (at least three columns) and `known_T_K`.

See [release notes](RELEASE_NOTES.md), [package contents](PACKAGE_CONTENTS.md), `CITATION.cff` and `LICENSE`.
