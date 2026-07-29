% Stage 1 + Stage 2 example: sparse reconstruction and thermometry.
%
% Required variables in the workspace before running:
%   energy_meV : Ne x 1 energy axis
%   spectrum   : Ne x 1 EELS spectrum
%
addpath(fileparts(fileparts(mfilename('fullpath'))));

settings = struct();
settings.run_stage2 = true;
settings.crop_range = [-250, 250];
settings.stage2_temperature_algorithm = 'grid_scan';
settings.use_bg = false;

result = eels_fit_spectrum(energy_meV, spectrum, settings);

fig = EELS_PearsonSparse_FullFit.plot_result(result, 'EELS-ZEST Stage 2 Thermometry');
fprintf('T = %.2f K, uncertainty = %.2f K, reliability = %.3f\n', ...
    result.T_fit, result.T_uncertainty, result.reliability_score);
