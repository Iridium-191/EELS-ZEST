# Native kernels and parallel execution

The MATLAB API and GUI call C++ numerical kernels for Stage1/Stage2. MATLAB retains the outer nonlinear search, model policy, temperature scans, diagnostics and user interface.

The Windows binary `eels_native_v22_mex.mexw64` is included. It was built for MATLAB R2024b with MSVC 2022 and uses MATLAB BLAS/LAPACK. Source: `native/eels_native_mex.cpp`.

Build date: 2026-11-09. `NATIVE_BUILD.json` records the source and binary SHA256 hashes.

## Backend selection

- `native_backend=auto`: use available native kernels with supported MATLAB fallbacks.
- `native_backend=cpp`: require the native binary; report unavailable/failed native paths.
- `native_backend=matlab`: use the reference MATLAB implementation.

Some specialized constraints still use MATLAB. Stage2 can fall back to MATLAB for a temperature point failing native optimality checks. Numerical trajectories can differ between native and MATLAB paths.

## Build on another host

Configure a MATLAB-supported C++ compiler, then run:

```matlab
mex -setup C++
build_eels_native
```

C++17 and MATLAB BLAS/LAPACK are required. Other platforms are not verified by this Windows release test. No standalone executable is provided.

## Parallel policy

Default: 8 process workers, limited by physical cores and local cluster capacity. A tuning group uses no more workers than spectra. Spectra are parallelized; native kernels use one BLAS thread per worker. An existing pool can be reused. The outer tuning solver does not start nested worker parallelism.

Save checkpoint and Resume checkpoint are off by default. Serial and parallel GUI runs both honor the checkbox. Automatic tuning report export remains enabled independently.
