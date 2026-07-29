# Examples

The bundled datacube is `data/example_data.mat`.

## Datacube Workflows

- `zest_datacube_full_stage1.m`: align/SVD-ready Stage 1 reconstruction.
- `zest_datacube_full_stage2_parallel.m`: Stage 1 + Stage 2 thermometry with
  parallel execution.
- `zest_reuse_stage1_for_stage2.m`: reuse a saved Stage 1 result.

The full-cube examples use a `[-250, 250] meV` crop. The parallel script
requests 16 workers; reduce this value when memory is limited.

## Single-Spectrum And GUI Workflows

- `zest_gui_launch.m`
- `zest_single_stage1_background_removal.m`
- `zest_single_stage2_thermometry.m`
- `zest_sparse_intensity_analysis_from_result.m`
- `zest_synthetic_smoke.m`

The single-spectrum scripts expect `energy_meV` and `spectrum` in the MATLAB
workspace. `zest_synthetic_smoke.m` is self-contained.

EELS-ZEST is developed by Ruilin Mao, ICQM, Peking University, and licensed
under `GPL-3.0-only`. Citation metadata is provided in `../CITATION.cff`.
