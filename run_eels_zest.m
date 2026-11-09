function app = run_eels_zest()
%RUN_EELS_ZEST Launch the EELS-ZEST GUI.
%
% EELS-ZEST = Zero-loss Estimation and Sparse Thermometry.
%
% Usage:
%   addpath('pearson_sparse_fullfit');
%   app = run_eels_zest;

root_dir = fileparts(mfilename('fullpath'));
addpath(root_dir);
app = EELS_ZEST_App();
end
