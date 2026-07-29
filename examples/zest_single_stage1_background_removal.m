% Stage 1 example: ZLP + sparse-convolution signal recovery for one spectrum.
%
% Required variables in the workspace before running:
%   energy_meV : Ne x 1 energy axis
%   spectrum   : Ne x 1 EELS spectrum
%
addpath(fileparts(fileparts(mfilename('fullpath'))));

settings = struct();
settings.run_stage2 = false;
settings.crop_range = [-250, 250];
settings.use_bg = false;

result = eels_fit_spectrum(energy_meV, spectrum, settings);

fig = EELS_PearsonSparse_FullFit.plot_result(result, 'EELS-ZEST Stage 1');
