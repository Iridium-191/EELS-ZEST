function cuberes = eels_fit_datacube(filename, settings)
% Fit a sampled or full EELS datacube with spatial diagnostics.
%
% Example:
%   addpath('pearson_sparse_fullfit');
%   settings = struct();
%   settings.datacube_sample_strategy = 'all';
%   settings.crop_range = [-200, 200];
%   cuberes = eels_fit_datacube('your_datacube.mat', settings);
%
% To process all pixels, set:
%   settings.datacube_sample_strategy = 'all';

if nargin < 2 || isempty(settings)
    settings = struct();
end

cuberes = EELS_PearsonSparse_FullFit.solve_datacube(filename, settings);
end
