function result = eels_fit_spectrum(ene, spectrum, settings)
% Fit one EELS spectrum with the Pearson sparse-convolution workflow.
%
% Example:
%   addpath('pearson_sparse_fullfit');
%   settings = struct('run_stage2', true);
%   result = eels_fit_spectrum(eneref, spectrum, settings);

if nargin < 3 || isempty(settings)
    settings = struct();
end

result = EELS_PearsonSparse_FullFit.solve_spectrum(ene, spectrum, settings);

if isfield(settings, 'output_dir') && ~isempty(settings.output_dir)
    if ~isfolder(settings.output_dir), mkdir(settings.output_dir); end
    save(fullfile(settings.output_dir, 'single_spectrum_result.mat'), 'result', '-v7.3');
    if ~isfield(settings, 'save_plots') || settings.save_plots
        fig = EELS_PearsonSparse_FullFit.plot_result(result, 'single spectrum');
        saveas(fig, fullfile(settings.output_dir, 'single_spectrum_diagnostic.png'));
        close(fig);
    end
end
end
