% Self-contained EELS-ZEST Stage 1 package smoke test.

root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(root_dir);

energy_meV = (-250:1:250).';
zlp = 1200 .* exp(-0.5 .* (energy_meV ./ 8.5) .^ 2);
loss = 42 .* exp(-0.5 .* ((energy_meV - 62) ./ 8) .^ 2) + ...
    28 .* exp(-0.5 .* ((energy_meV - 94) ./ 10) .^ 2);
gain = 15 .* exp(-0.5 .* ((energy_meV + 62) ./ 8) .^ 2) + ...
    8 .* exp(-0.5 .* ((energy_meV + 94) ./ 10) .^ 2);
spectrum = zlp + loss + gain + 0.25;

settings = struct();
settings.run_stage2 = false;
settings.crop_range = [-250, 250];
settings.stage1_speed_profile = 'v20_fast';
settings.use_bg = false;

result = eels_fit_spectrum(energy_meV, spectrum, settings);
assert(isfield(result, 'metrics') && isfinite(result.metrics.R2), ...
    'Synthetic smoke test did not produce a finite R2.');
assert(isfield(result, 'fitted_spectrum') && ...
    numel(result.fitted_spectrum) == numel(energy_meV), ...
    'Synthetic smoke test returned an unexpected fitted spectrum.');

fprintf('EELS-ZEST synthetic smoke passed: R2=%.8f, NRMSE=%.8f\n', ...
    result.metrics.R2, result.metrics.NRMSE);
