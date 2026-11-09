classdef EELS_ZEST_PeakFit
    %EELS_ZEST_PEAKFIT Robust Gaussian/Lorentzian/Pseudo-Voigt fitting.

    methods (Static)
        function config = defaults()
            config = struct('method', 'Gaussian', 'fit_range_meV', [-20 120], ...
                'center_init_meV', [20 50 75 100], ...
                'center_lower_meV', [10 40 65 90], ...
                'center_upper_meV', [30 60 85 110], ...
                'fwhm_lower_meV', 2 .* ones(1, 4), ...
                'fwhm_upper_meV', 50 .* ones(1, 4), ...
                'eta_init', 0.5 .* ones(1, 4), ...
                'eta_lower', zeros(1, 4), 'eta_upper', ones(1, 4), ...
                'num_starts', 6, 'robust_refit', true, 'random_seed', 23);
        end

        function result = fit_spectrum(ene, spectrum, config, prior_centers)
            if nargin < 3 || isempty(config), config = EELS_ZEST_PeakFit.defaults(); end
            if nargin < 4, prior_centers = []; end
            config = EELS_ZEST_PeakFit.sanitize(config);
            xall = double(ene(:)); yall = double(spectrum(:));
            fit_mask = xall >= config.fit_range_meV(1) & xall <= config.fit_range_meV(2) & ...
                isfinite(xall) & isfinite(yall);
            npeak = numel(config.center_init_meV);
            result = EELS_ZEST_PeakFit.empty_result(xall, config, npeak);
            if nnz(fit_mask) < max(8, 3 .* npeak), result.status = 'insufficient_points'; return; end
            x = xall(fit_mask); y = yall(fit_mask);
            scale = max(abs(y), [], 'omitnan');
            if ~isfinite(scale) || scale <= 0, result.status = 'empty_spectrum'; return; end
            ys = y ./ scale;
            centers = config.center_init_meV;
            if numel(prior_centers) == npeak && all(isfinite(prior_centers))
                centers = min(max(double(prior_centers(:).'), config.center_lower_meV), config.center_upper_meV);
            end
            [p0, lb, ub] = EELS_ZEST_PeakFit.parameter_vectors(x, ys, centers, config);
            stream = RandStream('mt19937ar', 'Seed', config.random_seed);
            best_p = p0; best_cost = Inf; best_exit = -1;
            for ss = 1:max(1, round(config.num_starts))
                trial = p0;
                if ss > 1
                    trial = lb + (ub - lb) .* rand(stream, size(p0));
                    trial(1:EELS_ZEST_PeakFit.parameter_stride(config):end) = ...
                        p0(1:EELS_ZEST_PeakFit.parameter_stride(config):end) .* ...
                        (0.55 + 0.9 .* rand(stream, 1, npeak));
                end
                [p, cost, exitflag] = EELS_ZEST_PeakFit.solve_once(x, ys, trial, lb, ub, config, []);
                if isfinite(cost) && cost < best_cost
                    best_p = p; best_cost = cost; best_exit = exitflag;
                end
            end
            if config.robust_refit && isfinite(best_cost)
                resid = ys - EELS_ZEST_PeakFit.model(best_p, x, config.method);
                med = median(resid, 'omitnan');
                sigma = 1.4826 .* median(abs(resid - med), 'omitnan');
                if isfinite(sigma) && sigma > eps
                    u = (resid - med) ./ (4.685 .* sigma);
                    weights = (abs(u) < 1) .* (1 - min(u .^ 2, 1)) .^ 2;
                    [p, cost, exitflag] = EELS_ZEST_PeakFit.solve_once(x, ys, best_p, lb, ub, config, weights);
                    if isfinite(cost) && cost <= best_cost .* 1.2
                        best_p = p; best_cost = cost; best_exit = exitflag;
                    end
                end
            end
            params = EELS_ZEST_PeakFit.unpack(best_p, config.method);
            fitted_scaled = EELS_ZEST_PeakFit.model(best_p, xall, config.method);
            fitted = fitted_scaled .* scale;
            params.amplitude = params.amplitude .* scale;
            params.area = EELS_ZEST_PeakFit.areas(params, config.method);
            residual = y - fitted(fit_mask);
            nrmse = sqrt(mean(residual .^ 2, 'omitnan')) ./ max(max(y) - min(y), eps);
            result.status = 'ok'; result.success = best_exit > 0 || isfinite(best_cost);
            result.amplitude = params.amplitude; result.area = params.area;
            result.peak_position_meV = params.center; result.fwhm_meV = params.fwhm;
            result.eta = params.eta; result.fitted_spectrum = fitted;
            result.components = EELS_ZEST_PeakFit.components(best_p, xall, config.method) .* scale;
            result.NRMSE = nrmse; result.cost = best_cost; result.exitflag = best_exit;
            result.fit_mask = fit_mask; result.scale = scale;
        end

        function batch = fit_line_profile(cube, ene, y_range, config)
            [ny, nx, ~] = size(cube);
            y_range = min(max(round(y_range), 1), ny);
            ys = min(y_range):max(y_range);
            npeak = numel(config.center_init_meV);
            batch = EELS_ZEST_PeakFit.empty_batch([1 nx], npeak, config, 'line_profile');
            prior = [];
            for xx = 1:nx
                spectrum = squeeze(mean(cube(ys, xx, :), 1, 'omitnan'));
                fit = EELS_ZEST_PeakFit.fit_spectrum(ene, spectrum, config, prior);
                batch = EELS_ZEST_PeakFit.store_batch(batch, fit, 1, xx);
                if fit.success, prior = fit.peak_position_meV; end
            end
            batch.y_range = [min(ys), max(ys)];
        end

        function batch = fit_cube(cube, ene, config, use_parallel, workers, progress_fcn, stop_fcn)
            if nargin < 4, use_parallel = false; end
            if nargin < 5 || isempty(workers), workers = 16; end
            if nargin < 6, progress_fcn = []; end
            if nargin < 7, stop_fcn = []; end
            [ny, nx, ~] = size(cube); npix = ny .* nx;
            npeak = numel(config.center_init_meV);
            batch = EELS_ZEST_PeakFit.empty_batch([ny nx], npeak, config, 'datacube');
            batch_size = max(8, 4 .* max(1, round(workers)));
            for first = 1:batch_size:npix
                if ~isempty(stop_fcn) && stop_fcn(), batch.cancelled = true; break; end
                ids = first:min(npix, first + batch_size - 1);
                fits = cell(size(ids));
                if use_parallel && numel(ids) >= 4
                    EELS_ZEST_PeakFit.ensure_pool(workers);
                    parfor jj = 1:numel(ids)
                        [yy, xx] = ind2sub([ny nx], ids(jj));
                        fits{jj} = EELS_ZEST_PeakFit.fit_spectrum(ene, squeeze(cube(yy, xx, :)), config);
                    end
                else
                    for jj = 1:numel(ids)
                        [yy, xx] = ind2sub([ny nx], ids(jj));
                        fits{jj} = EELS_ZEST_PeakFit.fit_spectrum(ene, squeeze(cube(yy, xx, :)), config);
                    end
                end
                for jj = 1:numel(ids)
                    [yy, xx] = ind2sub([ny nx], ids(jj));
                    batch = EELS_ZEST_PeakFit.store_batch(batch, fits{jj}, yy, xx);
                end
                if ~isempty(progress_fcn), progress_fcn(ids(end), npix); end
            end
        end

        function table_out = long_table(batch)
            [ny, nx, npeak] = size(batch.amplitude);
            [Y, X, P] = ndgrid(1:ny, 1:nx, 1:npeak);
            table_out = table(Y(:), X(:), P(:), batch.amplitude(:), batch.area(:), ...
                batch.peak_position_meV(:), batch.fwhm_meV(:), batch.eta(:), ...
                repmat(batch.NRMSE(:), npeak, 1), repmat(batch.success_mask(:), npeak, 1), ...
                'VariableNames', {'y', 'x', 'peak', 'amplitude', 'area', ...
                'peak_position_meV', 'fwhm_meV', 'eta', 'NRMSE', 'success'});
        end

        function files = export_parameter_maps(batch, output_dir, format)
            assert(isstruct(batch) && isfield(batch, 'amplitude') && isfield(batch, 'mode'), ...
                'EELS:PeakMapRequired', 'A full-map or line-profile peak-fit result is required.');
            if ~isfolder(output_dir), mkdir(output_dir); end
            format = lower(char(string(format)));
            if strcmp(format, 'tiff'), ext = 'tiff'; else, ext = format; end
            map_fields = {'amplitude', 'area', 'peak_position_meV', 'fwhm_meV', 'eta'};
            files = {};
            for pp = 1:size(batch.amplitude, 3)
                for ff = 1:numel(map_fields)
                    name = map_fields{ff};
                    all_values = batch.(name);
                    values = all_values(:, :, pp);
                    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 760 620]);
                    cleanup = onCleanup(@()delete(fig));
                    ax = axes(fig); imagesc(ax, values); axis(ax, 'image');
                    colormap(ax, turbo); colorbar(ax); xlabel(ax, 'x'); ylabel(ax, 'y');
                    title(ax, sprintf('Peak %d %s', pp, strrep(name, '_', ' ')));
                    filename = fullfile(output_dir, sprintf('peak_%02d_%s.%s', pp, name, ext));
                    EELS_ZEST_Export.save_figure(filename, fig, format);
                    files{end + 1} = filename; %#ok<AGROW>
                    clear cleanup
                end
            end
            diagnostics = {'NRMSE', 'success_mask'};
            for ff = 1:numel(diagnostics)
                name = diagnostics{ff};
                fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 760 620]);
                cleanup = onCleanup(@()delete(fig));
                ax = axes(fig); imagesc(ax, batch.(name)); axis(ax, 'image');
                colormap(ax, turbo); colorbar(ax); xlabel(ax, 'x'); ylabel(ax, 'y');
                title(ax, strrep(name, '_', ' '));
                filename = fullfile(output_dir, sprintf('%s.%s', lower(name), ext));
                EELS_ZEST_Export.save_figure(filename, fig, format);
                files{end + 1} = filename; %#ok<AGROW>
                clear cleanup
            end
        end
    end

    methods (Static, Access = private)
        function config = sanitize(config)
            defaults = EELS_ZEST_PeakFit.defaults();
            names = fieldnames(defaults);
            for ii = 1:numel(names)
                if ~isfield(config, names{ii}) || isempty(config.(names{ii}))
                    config.(names{ii}) = defaults.(names{ii});
                end
            end
            config.method = char(string(config.method));
            config.fit_range_meV = sort(double(config.fit_range_meV(1:2)));
            config.center_init_meV = double(config.center_init_meV(:).');
            n = numel(config.center_init_meV);
            fields = {'center_lower_meV', 'center_upper_meV', 'fwhm_lower_meV', ...
                'fwhm_upper_meV', 'eta_init', 'eta_lower', 'eta_upper'};
            for ii = 1:numel(fields)
                value = double(config.(fields{ii})(:).');
                if isscalar(value), value = repmat(value, 1, n); end
                assert(numel(value) == n, 'EELS:PeakConfigSize', '%s must match the peak count.', fields{ii});
                config.(fields{ii}) = value;
            end
            assert(n >= 1 && n <= 20, 'EELS:PeakCount', 'Peak count must be between 1 and 20.');
            config.center_init_meV = min(max(config.center_init_meV, config.center_lower_meV), config.center_upper_meV);
            config.fwhm_lower_meV = max(config.fwhm_lower_meV, eps);
            config.fwhm_upper_meV = max(config.fwhm_upper_meV, config.fwhm_lower_meV + eps);
            config.eta_lower = min(max(config.eta_lower, 0), 1);
            config.eta_upper = max(min(config.eta_upper, 1), config.eta_lower + eps);
            config.eta_init = min(max(config.eta_init, config.eta_lower), config.eta_upper);
        end

        function stride = parameter_stride(method_or_config)
            if isstruct(method_or_config), method = method_or_config.method; else, method = method_or_config; end
            if contains(lower(string(method)), 'pseudo'), stride = 4; else, stride = 3; end
        end

        function [p0, lb, ub] = parameter_vectors(x, y, centers, config)
            n = numel(centers); stride = EELS_ZEST_PeakFit.parameter_stride(config);
            p0 = zeros(1, n .* stride); lb = p0; ub = p0;
            amp_max = max(max(y), eps) .* 3;
            for ii = 1:n
                idx = (ii - 1) .* stride + (1:stride);
                [~, near] = min(abs(x - centers(ii)));
                amp = max(y(near), max(y) ./ max(n, 1));
                base = [amp, centers(ii), 0.5 .* (config.fwhm_lower_meV(ii) + config.fwhm_upper_meV(ii))];
                low = [0, config.center_lower_meV(ii), config.fwhm_lower_meV(ii)];
                high = [amp_max, config.center_upper_meV(ii), config.fwhm_upper_meV(ii)];
                if stride == 4
                    base(4) = config.eta_init(ii); low(4) = config.eta_lower(ii); high(4) = config.eta_upper(ii);
                end
                p0(idx) = base; lb(idx) = low; ub(idx) = high;
            end
        end

        function [p, cost, exitflag] = solve_once(x, y, p0, lb, ub, config, weights)
            if isempty(weights), weights = ones(size(y)); end
            weights = sqrt(max(double(weights(:)), 0));
            model = @(p, xx) EELS_ZEST_PeakFit.model(p, xx, config.method);
            if exist('lsqnonlin', 'file') == 2
                objective = @(p) weights .* (model(p, x) - y);
                opts = optimoptions('lsqnonlin', 'Display', 'off', 'MaxIterations', 700, ...
                    'FunctionTolerance', 1e-8, 'StepTolerance', 1e-8);
                [p, cost, ~, exitflag] = lsqnonlin(objective, p0, lb, ub, opts);
            else
                objective = @(p) sum((weights .* (model(min(max(p, lb), ub), x) - y)) .^ 2);
                opts = optimset('Display', 'off', 'MaxIter', 4000, 'TolX', 1e-7, 'TolFun', 1e-8);
                [p, cost, exitflag] = fminsearch(objective, p0, opts);
                p = min(max(p, lb), ub);
            end
        end

        function y = model(p, x, method)
            components = EELS_ZEST_PeakFit.components(p, x, method);
            y = sum(components, 2);
        end

        function values = components(p, x, method)
            x = double(x(:)); stride = EELS_ZEST_PeakFit.parameter_stride(method);
            n = numel(p) ./ stride; values = zeros(numel(x), n);
            for ii = 1:n
                q = p((ii - 1) .* stride + (1:stride));
                amp = q(1); center = q(2); fwhm = max(q(3), eps);
                gaussian = exp(-4 .* log(2) .* ((x - center) ./ fwhm) .^ 2);
                lorentz = 1 ./ (1 + 4 .* ((x - center) ./ fwhm) .^ 2);
                if contains(lower(string(method)), 'pseudo')
                    eta = min(max(q(4), 0), 1); values(:, ii) = amp .* ((1 - eta) .* gaussian + eta .* lorentz);
                elseif contains(lower(string(method)), 'lorentz')
                    values(:, ii) = amp .* lorentz;
                else
                    values(:, ii) = amp .* gaussian;
                end
            end
        end

        function params = unpack(p, method)
            stride = EELS_ZEST_PeakFit.parameter_stride(method); n = numel(p) ./ stride;
            q = reshape(p, stride, n).';
            params = struct('amplitude', q(:, 1).', 'center', q(:, 2).', ...
                'fwhm', q(:, 3).', 'eta', NaN(1, n), 'area', NaN(1, n));
            if stride == 4, params.eta = q(:, 4).'; end
        end

        function area = areas(params, method)
            if contains(lower(string(method)), 'pseudo')
                eta = params.eta;
                area = params.amplitude .* params.fwhm .* ...
                    ((1 - eta) .* sqrt(pi ./ (4 .* log(2))) + eta .* pi ./ 2);
            elseif contains(lower(string(method)), 'lorentz')
                area = params.amplitude .* params.fwhm .* pi ./ 2;
            else
                area = params.amplitude .* params.fwhm .* sqrt(pi ./ (4 .* log(2)));
            end
        end

        function result = empty_result(ene, config, npeak)
            result = struct('status', 'not_run', 'success', false, ...
                'amplitude', NaN(1, npeak), 'area', NaN(1, npeak), ...
                'peak_position_meV', NaN(1, npeak), 'fwhm_meV', NaN(1, npeak), ...
                'eta', NaN(1, npeak), 'fitted_spectrum', NaN(size(ene)), ...
                'components', NaN(numel(ene), npeak), 'NRMSE', NaN, ...
                'cost', NaN, 'exitflag', NaN, 'fit_mask', false(size(ene)), ...
                'scale', NaN, 'config', config);
        end

        function batch = empty_batch(map_size, npeak, config, mode)
            batch = struct('mode', mode, 'amplitude', NaN([map_size npeak]), ...
                'area', NaN([map_size npeak]), 'peak_position_meV', NaN([map_size npeak]), ...
                'fwhm_meV', NaN([map_size npeak]), 'eta', NaN([map_size npeak]), ...
                'NRMSE', NaN(map_size), 'success_mask', false(map_size), ...
                'status', {cell(map_size)}, 'config', config, 'cancelled', false);
        end

        function batch = store_batch(batch, fit, yy, xx)
            batch.amplitude(yy, xx, :) = fit.amplitude;
            batch.area(yy, xx, :) = fit.area;
            batch.peak_position_meV(yy, xx, :) = fit.peak_position_meV;
            batch.fwhm_meV(yy, xx, :) = fit.fwhm_meV;
            batch.eta(yy, xx, :) = fit.eta;
            batch.NRMSE(yy, xx) = fit.NRMSE;
            batch.success_mask(yy, xx) = fit.success;
            batch.status{yy, xx} = fit.status;
        end

        function ensure_pool(workers)
            pool = gcp('nocreate');
            workers = max(1, round(workers));
            if isempty(pool)
                parpool('local', workers);
            end
        end
    end
end
