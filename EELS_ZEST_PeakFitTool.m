function fig = EELS_ZEST_PeakFitTool(cube, ene, info)
%EELS_ZEST_PEAKFITTOOL Interactive robust multi-peak fitting GUI.

if nargin < 3 || ~isstruct(info), info = struct(); end
ene = double(ene(:)); [ny, nx, ne] = size(cube);
assert(ne == numel(ene), 'EELS:InvalidDatacube', 'Energy dimension mismatch.');
config0 = EELS_ZEST_PeakFit.defaults();
map = sum(cube, 3, 'omitnan');
fig = figure('Name', 'EELS-ZEST Robust Multi-Peak Fit', 'NumberTitle', 'off', ...
    'Color', 'w', 'Position', [90 70 1320 780], 'CloseRequestFcn', @closeTool);
ax_map = axes(fig, 'Units', 'normalized', 'Position', [0.04 0.41 0.34 0.53]);
imagesc(ax_map, map); axis(ax_map, 'image'); colormap(ax_map, turbo); colorbar(ax_map);
xlabel(ax_map, 'x'); ylabel(ax_map, 'y'); title(ax_map, 'Fit ROI');
ax_spec = axes(fig, 'Units', 'normalized', 'Position', [0.45 0.41 0.51 0.53]);
grid(ax_spec, 'on'); xlabel(ax_spec, 'Energy loss (meV)'); ylabel(ax_spec, 'Intensity');

table_data = [num2cell(true(4, 1)), num2cell(config0.center_init_meV(:)), ...
    num2cell(config0.center_lower_meV(:)), num2cell(config0.center_upper_meV(:)), ...
    num2cell(config0.fwhm_lower_meV(:)), num2cell(config0.fwhm_upper_meV(:)), ...
    num2cell(config0.eta_init(:)), num2cell(config0.eta_lower(:)), num2cell(config0.eta_upper(:))];
parameter_table = uitable(fig, 'Units', 'normalized', 'Position', [0.04 0.055 0.65 0.27], ...
    'Data', table_data, 'ColumnName', {'Use', 'Center', 'Lower', 'Upper', ...
    'FWHM min', 'FWHM max', 'Eta', 'Eta min', 'Eta max'}, ...
    'ColumnEditable', true(1, 9), 'ColumnFormat', {'logical', 'numeric', 'numeric', ...
    'numeric', 'numeric', 'numeric', 'numeric', 'numeric', 'numeric'});
uicontrol(fig, 'Style', 'text', 'String', 'Shape', 'Units', 'normalized', ...
    'Position', [0.72 0.285 0.06 0.03], 'BackgroundColor', 'w');
shape_control = uicontrol(fig, 'Style', 'popupmenu', ...
    'String', {'Gaussian', 'Lorentzian', 'PseudoVoigt'}, 'Units', 'normalized', ...
    'Position', [0.78 0.282 0.15 0.04], 'Callback', @fitCurrent);
uicontrol(fig, 'Style', 'text', 'String', 'Fit range', 'Units', 'normalized', ...
    'Position', [0.72 0.235 0.07 0.03], 'BackgroundColor', 'w');
range_control = uicontrol(fig, 'Style', 'edit', 'String', mat2str(config0.fit_range_meV), ...
    'Units', 'normalized', 'Position', [0.79 0.235 0.14 0.038], 'Callback', @fitCurrent);
parallel_control = uicontrol(fig, 'Style', 'checkbox', 'String', 'Parallel full map', ...
    'Value', 1, 'Units', 'normalized', 'Position', [0.72 0.19 0.14 0.035], 'BackgroundColor', 'w');
workers_control = uicontrol(fig, 'Style', 'edit', 'String', '16', ...
    'Units', 'normalized', 'Position', [0.87 0.19 0.06 0.035]);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Add peak', 'Units', 'normalized', ...
    'Position', [0.72 0.135 0.09 0.042], 'Callback', @addPeak);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Remove peak', 'Units', 'normalized', ...
    'Position', [0.82 0.135 0.11 0.042], 'Callback', @removePeak);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Fit ROI', 'Units', 'normalized', ...
    'Position', [0.72 0.08 0.09 0.045], 'FontWeight', 'bold', 'Callback', @fitCurrent);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Fit x-line', 'Units', 'normalized', ...
    'Position', [0.82 0.08 0.11 0.045], 'Callback', @fitLine);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Fit full map', 'Units', 'normalized', ...
    'Position', [0.72 0.025 0.11 0.045], 'Callback', @fitMap);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Export', 'Units', 'normalized', ...
    'Position', [0.84 0.025 0.09 0.045], 'Callback', @exportResults);

[max_y, max_x] = find(map == max(map(:), [], 'omitnan'), 1);
if isempty(max_y), max_y = 1; max_x = 1; end
roi = drawrectangle(ax_map, 'Position', [max(0.5, max_x - 1.5), max(0.5, max_y - 1.5), 3, 3], ...
    'Color', [0.9 0.15 0.08], 'LineWidth', 1.4);
listeners = {addlistener(roi, 'MovingROI', @fitCurrent), addlistener(roi, 'ROIMoved', @fitCurrent)};
current_result = [];
line_result = [];
map_result = [];
fitCurrent();

    function config = readConfig()
        data = parameter_table.Data;
        use = cell2mat(data(:, 1));
        assert(any(use), 'Enable at least one peak.');
        data = data(use, :);
        config = EELS_ZEST_PeakFit.defaults();
        config.method = shape_control.String{shape_control.Value};
        range = str2num(range_control.String); %#ok<ST2NM>
        assert(numel(range) == 2 && all(isfinite(range)), 'Fit range must contain two numbers.');
        config.fit_range_meV = sort(range(:).');
        config.center_init_meV = cell2mat(data(:, 2)).';
        config.center_lower_meV = cell2mat(data(:, 3)).';
        config.center_upper_meV = cell2mat(data(:, 4)).';
        config.fwhm_lower_meV = cell2mat(data(:, 5)).';
        config.fwhm_upper_meV = cell2mat(data(:, 6)).';
        config.eta_init = cell2mat(data(:, 7)).';
        config.eta_lower = cell2mat(data(:, 8)).';
        config.eta_upper = cell2mat(data(:, 9)).';
    end

    function [spec, yrange, xrange] = roiSpectrum()
        pos = roi.Position;
        x1 = min(nx, max(1, floor(pos(1)))); y1 = min(ny, max(1, floor(pos(2))));
        x2 = min(nx, max(x1, ceil(pos(1) + pos(3)))); y2 = min(ny, max(y1, ceil(pos(2) + pos(4))));
        spec = squeeze(mean(mean(cube(y1:y2, x1:x2, :), 1, 'omitnan'), 2, 'omitnan'));
        yrange = [y1 y2]; xrange = [x1 x2];
    end

    function fitCurrent(varargin) %#ok<INUSD>
        try
            [spec, yrange, xrange] = roiSpectrum(); config = readConfig();
            current_result = EELS_ZEST_PeakFit.fit_spectrum(ene, spec, config);
            current_result.roi_y_range = yrange; current_result.roi_x_range = xrange;
            cla(ax_spec); plot(ax_spec, ene, spec, '.', 'Color', [0.25 0.25 0.25]); hold(ax_spec, 'on');
            if current_result.success
                plot(ax_spec, ene, current_result.fitted_spectrum, 'Color', [0.48 0.10 0.75], 'LineWidth', 1.5);
                plot(ax_spec, ene, current_result.components, 'LineWidth', 0.9);
                title(ax_spec, sprintf('%s | NRMSE %.4g', config.method, current_result.NRMSE));
            else
                title(ax_spec, ['Fit unavailable: ', current_result.status], 'Interpreter', 'none');
            end
            hold(ax_spec, 'off'); grid(ax_spec, 'on'); xlim(ax_spec, config.fit_range_meV);
        catch ME
            title(ax_spec, ME.message, 'Interpreter', 'none');
        end
    end

    function addPeak(~, ~)
        data = parameter_table.Data;
        if size(data, 1) >= 20, return; end
        if isempty(data), center = 20; else, center = cell2mat(data(end, 2)) + 20; end
        data(end + 1, :) = {true, center, center - 10, center + 10, 2, 50, 0.5, 0, 1};
        parameter_table.Data = data; fitCurrent();
    end

    function removePeak(~, ~)
        data = parameter_table.Data;
        if size(data, 1) <= 1, return; end
        data(end, :) = []; parameter_table.Data = data; fitCurrent();
    end

    function fitLine(~, ~)
        try
            [~, yrange] = roiSpectrum(); config = readConfig();
            dlg = waitbar(0, 'Fitting x-line profile...', 'CreateCancelBtn', 'setappdata(gcbf,''cancel'',true)');
            setappdata(dlg, 'cancel', false); cleanup = onCleanup(@()localClose(dlg));
            line_result = EELS_ZEST_PeakFit.fit_line_profile(cube, ene, yrange, config);
            line_result.source_info = info;
            waitbar(1, dlg, sprintf('Complete: %d/%d', nnz(line_result.success_mask), nx));
            clear cleanup; localClose(dlg);
            showBatchSummary(line_result, 'Peak-fit x-line profile');
        catch ME
            errordlg(ME.message, 'Line-profile fit failed');
        end
    end

    function fitMap(~, ~)
        try
            config = readConfig(); workers = max(1, round(str2double(workers_control.String)));
            if ~isfinite(workers), workers = 16; end
            dlg = waitbar(0, 'Fitting datacube...', 'CreateCancelBtn', 'setappdata(gcbf,''cancel'',true)');
            setappdata(dlg, 'cancel', false); cleanup = onCleanup(@()localClose(dlg));
            progress = @(done, total)waitbar(done ./ total, dlg, sprintf('Fitting %d/%d', done, total));
            stop = @()~isgraphics(dlg) || logical(getappdata(dlg, 'cancel'));
            map_result = EELS_ZEST_PeakFit.fit_cube(cube, ene, config, ...
                logical(parallel_control.Value), workers, progress, stop);
            map_result.source_info = info;
            clear cleanup; localClose(dlg);
            showBatchSummary(map_result, 'Peak-fit datacube');
        catch ME
            errordlg(ME.message, 'Datacube fit failed');
        end
    end

    function showBatchSummary(batch, name)
        summary = figure('Name', name, 'Color', 'w');
        tiledlayout(summary, 1, 2, 'TileSpacing', 'compact');
        nexttile; imagesc(batch.NRMSE); axis image; colorbar; title('NRMSE');
        nexttile; imagesc(batch.success_mask); axis image; colorbar; title('Fit success');
    end

    function exportResults(~, ~)
        names = {'Current ROI'}; payloads = {current_result};
        if ~isempty(line_result), names{end + 1} = 'X-line profile'; payloads{end + 1} = line_result; end %#ok<AGROW>
        if ~isempty(map_result), names{end + 1} = 'Full datacube'; payloads{end + 1} = map_result; end %#ok<AGROW>
        choice = EELS_ZEST_ExportDialog('Export multi-peak fit', names, {'MAT', 'CSV', 'FIG', 'SVG', 'PNG', 'TIFF'});
        if choice.cancelled, return; end
        payload = payloads{choice.scope_index}; fmt = lower(choice.format);
        if ismember(fmt, {'svg', 'png', 'tiff'}) && isstruct(payload) && ...
                isfield(payload, 'mode') && strcmp(payload.mode, 'datacube')
            output_dir = uigetdir(pwd, 'Choose folder for parameter-map figures');
            if isequal(output_dir, 0), return; end
            EELS_ZEST_PeakFit.export_parameter_maps(payload, output_dir, fmt);
            if choice.assign_workspace, assignin('base', 'eels_zest_peak_fit', payload); end
            return;
        end
        [file, path] = uiputfile(localFilters(fmt), 'Export multi-peak fit');
        if isequal(file, 0), return; end
        filename = fullfile(path, file);
        if ismember(fmt, {'fig', 'svg', 'png', 'tiff'})
            EELS_ZEST_Export.save_figure(filename, fig, fmt);
        elseif strcmp(fmt, 'csv')
            if isfield(payload, 'mode')
                table_out = EELS_ZEST_PeakFit.long_table(payload);
            else
                peak = (1:numel(payload.amplitude)).';
                table_out = table(peak, payload.amplitude(:), payload.area(:), ...
                    payload.peak_position_meV(:), payload.fwhm_meV(:), payload.eta(:), ...
                    'VariableNames', {'peak', 'amplitude', 'area', 'peak_position_meV', 'fwhm_meV', 'eta'});
            end
            EELS_ZEST_Export.save_payload(filename, table_out, 'csv', 'peak_fit');
        else
            EELS_ZEST_Export.save_payload(filename, payload, 'mat', 'eels_zest_peak_fit');
        end
        if choice.assign_workspace, assignin('base', 'eels_zest_peak_fit', payload); end
    end

    function filters = localFilters(fmt)
        switch fmt
            case 'mat', filters = {'*.mat', 'MAT file'};
            case 'csv', filters = {'*.csv', 'CSV file'};
            case 'fig', filters = {'*.fig', 'MATLAB figure'};
            case 'svg', filters = {'*.svg', 'SVG figure'};
            case 'png', filters = {'*.png', 'PNG figure'};
            otherwise, filters = {'*.tiff', 'TIFF figure'};
        end
    end

    function localClose(handle)
        if isgraphics(handle), delete(handle); end
    end

    function closeTool(varargin) %#ok<INUSD>
        for jj = 1:numel(listeners)
            try delete(listeners{jj}); catch, end
        end
        listeners = {};
        try if isvalid(roi), delete(roi); end; catch, end
        if isgraphics(fig), delete(fig); end
    end
end
