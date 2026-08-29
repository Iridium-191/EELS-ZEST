% Example data: Stage 1 + Stage 2 with parallel execution.

example_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(example_dir));
input_file = fullfile(example_dir, 'data', 'example_data.mat');
if ~isfile(input_file)
    error('Missing example data: %s', input_file);
end

settings = struct();
settings.run_stage2 = true;
settings.save_plots = false;
settings.save_diagnostic_plots = false;
settings.crop_range = [-250, 250];
settings.datacube_sample_strategy = 'all';
settings.datacube_parallel_enabled = true;
settings.datacube_parallel_workers = 16;
settings.datacube_checkpoint_enabled = true;
settings.datacube_resume_checkpoint = true;
settings.stage1_hierarchical_enabled = true;
settings.stage1_spatial_zlp_enabled = true;
settings.stage2_spatial_map_enabled = true;
settings.spatial_refit_enabled = false;
settings.spatial_bin = [1, 1];
settings.energy_bin = 1;
settings.stage2_temperature_algorithm = 'grid_scan';
settings.test_output_dir = fullfile(example_dir, 'output', 'example_stage2');
settings.use_bg = false;

cuberes = eels_fit_datacube(input_file, settings);
