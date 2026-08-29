# EELS-ZEST v2.1.0 Notes

## Formal Release

Version 2.1.0 promotes the beta.4 hierarchical Stage 1 and ZLP-guided
temperature policy without changing the Stage 1/Stage 2 fitting model. It adds
file-oriented Data tools, robust multi-peak fitting, and GUI lifecycle and
memory fixes.

## Data Processing And Export

- Stage 1 and Stage 2 **Send signal to Data** actions now use the explicitly
  selected stage. Sparse phonon source intensity can be truncated by a manual
  cutoff, a Materials Project phonon JSON cutoff, or the Stage 2 hard-attention
  edge plus a configurable tail.
- Loss and gain source coefficients above the cutoff are zeroed together and
  reconvolved with that pixel's fitted-ZLP kernel. The final spectrum is not
  sharply clipped at the cutoff.
- The returned Data state remains the background-subtracted phonon
  reconstruction. Per-pixel fitted-ZLP kernels, Stage 1/Stage 2 temperature
  maps, validity masks, and result provenance now travel with the Data and are
  retained by `EELS` export/reload. RL deconvolves the current cutoff spectrum
  rather than silently returning to the pre-cutoff fit result.
- Data-map dimensions and canonical `sample_pixels` now control reconstruction
  placement and temperature transfer. This fixes failed Data refresh and
  `No matching Stage 2 fitted-temperature map is available` after binned fits.
- Programmatic GUI equivalents are available as
  `app.sendSignalToData(stage, options)` and `app.multiplySTA(source, T)`.
- **Multiply STA** keeps the positive-loss side and supports a manual
  temperature, Stage 1 soft temperature, or Stage 2 fitted temperature with a
  per-pixel Stage 1 fallback.
- **Export data** writes the currently displayed processing state as an `EELS`
  object with its energy axis, spatial metadata, masks, and processing history.
  Other export actions write MAT/CSV and figure formats where applicable;
  optional base-workspace assignment is disabled by default.

## Interactive Tools

- **Stack Spectra** is a nonblocking point/rectangle selection window with
  add, rename, delete, clear, reorder, and file export actions.
- **Robust Multi-Peak Fit** supports Gaussian, Lorentzian, and Pseudo-Voigt
  profiles, up to 20 peaks, six multi-starts, nonnegative heights, one robust
  refit, x-line warm starts, and parallel batched datacube fitting.
- Multi-peak exports include amplitudes, integrated areas, positions, FWHM,
  Pseudo-Voigt eta, NRMSE, success masks, long-form CSV, and parameter-map
  SVG/PNG/TIFF files. Fitted datacubes are not retained in GUI memory.

## GUI Reliability

- ROI listeners and owned figures use cell registries and explicit cleanup.
  ROI callbacks resume the waiting operation before deletion, avoiding invalid
  `event.listener` arrays.
- MAT input is inspected with `whos -file`; only candidate EELS, cube, energy,
  or result variables are deserialized. Unrelated saved GUI/listener objects
  are not loaded.
- Status logs are bounded and throttled. Map images and colorbars are reused.
  RL datacube conversion stores only the latent cube and compact quality maps;
  detailed single-pixel diagnostics use a three-entry LRU cache.
- STA and full-map multi-peak processing are block based, and parallel work is
  cancelable between bounded batches.

## Hierarchical Stage 1 And Temperature Policy

- ZLP-only deterministic clustering, soft class-template switching, and the
  eight-neighbor edge-aware Huber Stage 1 loss are unchanged from beta.4.
- ZLP feature differences mark sharp interfaces and disable coupling across
  them.
- `T_map_data_valid` contains only spectra with usable temperature evidence.
  `T_map_joint` copies every data-valid value unchanged and may fill a marked
  invalid pixel from at least two same-side ZLP-compatible valid neighbors.
  There is no temperature-difference smoothing term.

## Compatibility And Packaging

- `run_eels_zest` remains the main entry point. The legacy
  `EELS_TwoStage_Toolbox` remains compatible and has no new controls.
- The public archive contains only `examples/data/example_data.mat`; internal
  validation data and generated test outputs are excluded.
- Author: Ruilin Mao, International Center for Quantum Materials (ICQM),
  Peking University.
- License: GNU General Public License version 3 only (`GPL-3.0-only`).
