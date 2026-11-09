# 2.2.0 release notes

Release date: 2026-11-08. Based on the validated `version2.2beta_reference_tuner` work.

- Integrated reference tuning before Stage1/Stage2, multiple training points and separate heldout validation.
- 87 searchable controls covering ZLP shape, PM/EM regularization, attention selection and fit weights; full configuration inventory and dependency checks.
- Grouped/joint PSO and simulated annealing, plus optional surrogate hybrid search. SA is the current default; PSO remains available.
- Pure temperature MAE objective in K and independent residual eligibility gate.
- Adaptive FWHM-based zero-peak protection and retained MATLAB/C++ numerical interfaces.
- Strategy import/application to both stages, exact preservation of unedited numeric parameters, and benchmark dataset exports.
- English primary GUI, alternate toolbox GUI, standalone tuner, menus, dialogs and parameter advice.
- Defaults: residual ratio 2, 500 candidate evaluations, 8 workers, checkpoint saving/resume off.
- Fixed serial execution overriding the checkpoint checkbox: it now follows the requested setting.
