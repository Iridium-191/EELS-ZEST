% Example data: Stage 1 reconstruction.

example_dir = fileparts(mfilename('fullpath'));
addpath(fileparts(example_dir));
input_file = fullfile(example_dir, 'data', 'example_data.mat');
if ~isfile(input_file)
    error('Missing example data: %s', input_file);
end

settings = struct();
settings.run_stage2 = false;
settings.save_plots = false;
settings.save_diagnostic_plots = false;
settings.crop_range = [-250, 250];
settings.datacube_sample_strategy = 'all';
settings.datacube_parallel_enabled = false;
settings.stage1_hierarchical_enabled = true;
settings.stage1_spatial_zlp_enabled = true;
settings.spatial_bin = [1, 1];
settings.energy_bin = 1;
settings.test_output_dir = fullfile(example_dir, 'output', 'example_stage1');
settings.use_bg = false;

cuberes = eels_fit_datacube(input_file, settings);
