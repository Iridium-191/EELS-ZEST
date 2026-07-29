% Reuse an existing Stage 1 datacube result for Stage 2 thermometry.
%
% Required variable:
%   stage1_cuberes : output from eels_fit_datacube with run_stage2 = false
%
addpath(fileparts(fileparts(mfilename('fullpath'))));

settings = struct();
settings.run_stage2 = true;
settings.stage2_temperature_algorithm = 'grid_scan';
settings.spatial_refit_enabled = true;
settings.test_output_dir = fullfile(pwd, 'eels_zest_reused_stage2_output');
settings.use_bg = false;

stage2_cuberes = EELS_PearsonSparse_FullFit.promote_datacube_stage1_to_stage2(stage1_cuberes, settings);
