function app = run_eels_twostage_toolbox()
% Launch the EELS two-stage sparse-convolution GUI toolbox.
%
% Usage:
%   addpath('pearson_sparse_fullfit');
%   app = run_eels_twostage_toolbox;

root_dir = fileparts(mfilename('fullpath'));
addpath(root_dir);
app = EELS_TwoStage_Toolbox;
end
