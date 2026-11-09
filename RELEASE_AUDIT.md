# Release audit: 2.2.0

Checked on 2026-11-08 with Windows MATLAB R2024b. All checks below passed.

- MATLAB syntax checks for every top-level runtime source file.
- Pure-temperature objective unit checks: MAE arithmetic, no residual/valley-height/worst-point penalties, separate residual eligibility, failure handling and tied-valley diagnostics.
- Primary and alternate GUIs in compact and wide layouts, plus standalone reference tuner startup.
- Defaults asserted: residual ratio 2, budget 500, workers 8, checkpoint saving/resume disabled.
- Serial Stage1 and Stage2 settings follow the checkpoint checkbox both when disabled and enabled.
- GUI labels, menus, tooltips and table headers checked for untranslated Chinese text. English parameter descriptions/advice and dynamic messages were also scanned in source. A DM3 metadata pattern and localized MATLAB error-recognition string are intentionally preserved; neither is a UI label.
- Translation format placeholders preserved, including dynamic error/status messages.
- Native tuning smoke: two real training spectra, one heldout spectrum, two workers and seven candidate evaluations. Pure MAE matched candidate errors and search scores.
- Complete export, strategy application to Stage1/2 with all 87 numeric controls compared, and strategy reimport passed.
- Primary tuning, main wide view and alternate compact view screenshots inspected.

Build verification on 2026-11-09:

- Recompiled the C++ source with MATLAB R2024b and Microsoft Visual C++ 2022. Source/binary hashes are recorded in `NATIVE_BUILD.json`.
- Checked the configuration schema, runtime context, GUI strategy import/application and all 87 tuning controls.
- Ran one real spectrum through the rebuilt native Stage1/Stage2 backend and exported its strategy. The fitted temperature was 366.25 K; the Stage2 cost grid and fitted spectrum matched the numerical baseline exactly.
- Checked MATLAB syntax and the pure-temperature outer objective. No full-map or long optimizer run was performed.
- Lightweight check: `addpath('tests'); test_settings_contract`. Supply a fixture MAT with `energy_meV`, `spectra` and `known_T_K` to include a native single-spectrum fit.
