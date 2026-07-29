# EELS-ZEST v2.0.0 Notes

## Attention Hard-BE + Soft-Intensity Stage 2

- Stage 2 uses the fitted total ZLP as its convolution kernel and keeps the
  kernel center at `0 meV` on an already aligned energy axis.
- Hard-BE pairs are selected from paired attention evidence using
  `sqrt(att(+EM) * att(-EM))`. The default hard threshold is `0.35`; a `0.10`
  floor is used only to supplement too-small hard sets.
- Remaining EM centers stay in the spectral fit as independent Stokes and
  anti-Stokes source-intensity coefficients. They represent mixed phases,
  tails, and non-equilibrium signal without controlling the reported
  temperature.
- Temperature linearity, effective-pair count, reliability, and the cost valley
  use only hard attention pairs.
- `grid_scan` is the single public temperature path, with adaptive
  coarse-to-fine scanning and optional local refinement.
- Full datacubes may be preprocessed by ZLP alignment and rank-limited SVD,
  cropped to `[-250, 250] meV`, and processed with parallel workers.
- The cost-profile smoother uses temperature-distance-aware local regression,
  preventing coarse/fine grid gaps from creating artificial valley corners.
- Publication ROI selection requires the fitted-ZLP kernel center to remain
  within the configured energy-step tolerance.

## Bundled Example

- `examples/data/example_data.mat` can be processed directly with the Stage 1
  and Stage 2 scripts under `examples/`.
- Figure-ready data, panel CSV files, and reproduction scripts are under
  `docs/example_data_figures/`.

## Author, Citation, And License

- Author: Ruilin Mao, International Center for Quantum Materials (ICQM),
  Peking University.
- License: GNU General Public License version 3 only (`GPL-3.0-only`).
- Machine-readable citation metadata is provided in `CITATION.cff`.
