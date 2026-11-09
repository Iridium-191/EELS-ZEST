function fig = EELS_ZEST_StackTool(cube, ene, info, energy_window, completion_callback)
%EELS_ZEST_STACKTOOL Nonblocking point/ROI spectrum stacking window.

if nargin < 3 || ~isstruct(info), info = struct(); end
if nargin < 4 || isempty(energy_window), energy_window = [min(ene), max(ene)]; end
if nargin < 5, completion_callback = []; end
ene = double(ene(:));
[ny, nx, ne] = size(cube);
assert(ne == numel(ene), 'EELS:InvalidDatacube', 'Energy dimension mismatch.');
mask = ene >= min(energy_window) & ene <= max(energy_window);
if ~any(mask), mask = true(size(ene)); end
map = sum(cube(:, :, mask), 3, 'omitnan');

fig = figure('Name', 'EELS-ZEST Stack Spectra', 'NumberTitle', 'off', ...
    'Color', 'w', 'Position', [120 100 1180 700], 'CloseRequestFcn', @closeTool);
ax_map = axes(fig, 'Units', 'normalized', 'Position', [0.05 0.20 0.40 0.72]);
imagesc(ax_map, map); axis(ax_map, 'image'); colormap(ax_map, turbo); colorbar(ax_map);
xlabel(ax_map, 'x'); ylabel(ax_map, 'y'); title(ax_map, 'Select on map');
ax_spec = axes(fig, 'Units', 'normalized', 'Position', [0.53 0.28 0.43 0.64]);
grid(ax_spec, 'on'); xlabel(ax_spec, 'Energy loss (meV)'); ylabel(ax_spec, 'Intensity');
title(ax_spec, 'Stacked spectra');

uicontrol(fig, 'Style', 'text', 'String', 'Mode', 'Units', 'normalized', ...
    'Position', [0.05 0.135 0.05 0.035], 'BackgroundColor', 'w');
mode_control = uicontrol(fig, 'Style', 'popupmenu', 'String', {'Point', 'Rectangle'}, ...
    'Units', 'normalized', 'Position', [0.10 0.135 0.12 0.042]);
name_control = uicontrol(fig, 'Style', 'edit', 'String', 'spectrum_01', ...
    'Units', 'normalized', 'Position', [0.23 0.135 0.14 0.042]);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'New selection', ...
    'Units', 'normalized', 'Position', [0.05 0.075 0.12 0.045], 'Callback', @newSelection);
keep_button = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Keep', ...
    'Units', 'normalized', 'Position', [0.18 0.075 0.08 0.045], ...
    'Enable', 'off', 'Callback', @keepSelection);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Delete', ...
    'Units', 'normalized', 'Position', [0.27 0.075 0.08 0.045], 'Callback', @deleteSelection);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Clear', ...
    'Units', 'normalized', 'Position', [0.36 0.075 0.08 0.045], 'Callback', @clearSelections);
list_control = uicontrol(fig, 'Style', 'listbox', 'String', {}, ...
    'Units', 'normalized', 'Position', [0.53 0.075 0.20 0.15], 'Max', 1);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Up', 'Units', 'normalized', ...
    'Position', [0.74 0.175 0.06 0.045], 'Callback', @(~, ~)moveSelection(-1));
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Down', 'Units', 'normalized', ...
    'Position', [0.74 0.125 0.06 0.045], 'Callback', @(~, ~)moveSelection(1));
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Export', 'Units', 'normalized', ...
    'Position', [0.82 0.125 0.07 0.06], 'FontWeight', 'bold', 'Callback', @exportStack);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Done', 'Units', 'normalized', ...
    'Position', [0.90 0.125 0.06 0.06], 'Callback', @finishTool);

stack = struct('name', {{}}, 'spectrum', {{}}, 'selection', {{}}, ...
    'ene', ene, 'source_info', info, 'aggregation', 'sum');
active_roi = [];
active_listener = {};
mark_handles = {};

    function newSelection(~, ~)
        clearActiveROI();
        hold(ax_map, 'on');
        if mode_control.Value == 1
            active_roi = drawpoint(ax_map, 'Color', [0.55 0.00 0.90], 'LineWidth', 1.5);
        else
            active_roi = drawrectangle(ax_map, 'Color', [0.00 0.50 0.95], 'LineWidth', 1.4);
        end
        hold(ax_map, 'off');
        active_listener = {addlistener(active_roi, 'MovingROI', @previewSelection), ...
            addlistener(active_roi, 'ROIMoved', @previewSelection)};
        keep_button.Enable = 'on';
        previewSelection();
    end

    function previewSelection(varargin) %#ok<INUSD>
        if isempty(active_roi) || ~isvalid(active_roi), return; end
        [spec, ~] = spectrumFromROI(active_roi);
        cla(ax_spec); hold(ax_spec, 'on');
        for jj = 1:numel(stack.spectrum)
            plot(ax_spec, ene, stack.spectrum{jj}, 'LineWidth', 1.0);
        end
        plot(ax_spec, ene, spec, 'k--', 'LineWidth', 1.3);
        hold(ax_spec, 'off'); grid(ax_spec, 'on');
        title(ax_spec, 'Dashed: active selection');
    end

    function keepSelection(~, ~)
        if isempty(active_roi) || ~isvalid(active_roi), return; end
        [spec, selection] = spectrumFromROI(active_roi);
        label = strtrim(name_control.String);
        if isempty(label), label = sprintf('spectrum_%02d', numel(stack.spectrum) + 1); end
        selection.label = label;
        stack.name{end + 1} = label;
        stack.spectrum{end + 1} = spec(:);
        stack.selection{end + 1} = selection;
        mark_handles{end + 1} = makePermanentMark(selection); %#ok<AGROW>
        clearActiveROI();
        name_control.String = sprintf('spectrum_%02d', numel(stack.spectrum) + 1);
        updateListAndPlot();
    end

    function [spec, selection] = spectrumFromROI(roi)
        if isa(roi, 'images.roi.Point')
            xx = min(nx, max(1, round(roi.Position(1))));
            yy = min(ny, max(1, round(roi.Position(2))));
            spec = double(squeeze(cube(yy, xx, :)));
            selection = struct('type', 'point', 'x', xx, 'y', yy);
        else
            pos = roi.Position;
            x1 = min(nx, max(1, floor(pos(1))));
            y1 = min(ny, max(1, floor(pos(2))));
            x2 = min(nx, max(x1, ceil(pos(1) + pos(3))));
            y2 = min(ny, max(y1, ceil(pos(2) + pos(4))));
            spec = double(squeeze(sum(sum(cube(y1:y2, x1:x2, :), 1, 'omitnan'), 2, 'omitnan')));
            selection = struct('type', 'rectangle', 'x_range', [x1 x2], ...
                'y_range', [y1 y2], 'aggregation', 'sum');
        end
    end

    function h = makePermanentMark(selection)
        hold(ax_map, 'on');
        if strcmp(selection.type, 'point')
            h = plot(ax_map, selection.x, selection.y, 'o', 'Color', [0.55 0 0.9], ...
                'MarkerFaceColor', [0.88 0.72 1], 'LineWidth', 1.2);
        else
            h = rectangle(ax_map, 'Position', [selection.x_range(1), selection.y_range(1), ...
                diff(selection.x_range) + 1, diff(selection.y_range) + 1], ...
                'EdgeColor', [0 0.48 0.95], 'LineStyle', '--', 'LineWidth', 1.2);
        end
        hold(ax_map, 'off');
    end

    function updateListAndPlot()
        list_control.String = stack.name;
        if isempty(stack.name), list_control.Value = 1; else, list_control.Value = min(list_control.Value, numel(stack.name)); end
        cla(ax_spec); hold(ax_spec, 'on');
        for jj = 1:numel(stack.spectrum)
            plot(ax_spec, ene, stack.spectrum{jj}, 'LineWidth', 1.2, 'DisplayName', stack.name{jj});
        end
        hold(ax_spec, 'off'); grid(ax_spec, 'on');
        if ~isempty(stack.spectrum), legend(ax_spec, 'show', 'Interpreter', 'none', 'Location', 'best'); end
        title(ax_spec, sprintf('%d stacked spectra', numel(stack.spectrum)));
    end

    function deleteSelection(~, ~)
        if isempty(stack.spectrum), return; end
        idx = min(list_control.Value, numel(stack.spectrum));
        stack.name(idx) = []; stack.spectrum(idx) = []; stack.selection(idx) = [];
        if idx <= numel(mark_handles) && isgraphics(mark_handles{idx}), delete(mark_handles{idx}); end
        mark_handles(idx) = [];
        updateListAndPlot();
    end

    function clearSelections(~, ~)
        for jj = 1:numel(mark_handles)
            if isgraphics(mark_handles{jj}), delete(mark_handles{jj}); end
        end
        mark_handles = {};
        stack.name = {}; stack.spectrum = {}; stack.selection = {};
        updateListAndPlot();
    end

    function moveSelection(direction)
        if numel(stack.spectrum) < 2, return; end
        idx = min(list_control.Value, numel(stack.spectrum));
        dst = min(max(idx + direction, 1), numel(stack.spectrum));
        if dst == idx, return; end
        fields = {'name', 'spectrum', 'selection'};
        for ff = 1:numel(fields)
            temp = stack.(fields{ff}){idx};
            stack.(fields{ff}){idx} = stack.(fields{ff}){dst};
            stack.(fields{ff}){dst} = temp;
        end
        temp = mark_handles{idx}; mark_handles{idx} = mark_handles{dst}; mark_handles{dst} = temp;
        list_control.Value = dst;
        updateListAndPlot();
    end

    function exportStack(~, ~)
        if isempty(stack.spectrum), return; end
        choices = EELS_ZEST_ExportDialog('Export stacked spectra', {'All spectra', 'Selected spectrum'}, ...
            {'MAT', 'CSV', 'FIG', 'SVG', 'PNG', 'TIFF'});
        if choices.cancelled, return; end
        payload = stack;
        if choices.scope_index == 2
            idx = min(list_control.Value, numel(stack.spectrum));
            payload.name = stack.name(idx); payload.spectrum = stack.spectrum(idx); payload.selection = stack.selection(idx);
        end
        fmt = lower(choices.format);
        filters = localFilters(fmt);
        [file, path] = uiputfile(filters, 'Export stacked spectra');
        if isequal(file, 0), return; end
        filename = fullfile(path, file);
        if ismember(fmt, {'fig', 'svg', 'png', 'tiff'})
            EELS_ZEST_Export.save_figure(filename, fig, fmt);
        else
            payload.spectrum = cat(2, payload.spectrum{:});
            EELS_ZEST_Export.save_payload(filename, payload, fmt, 'eels_zest_stacked_spectra');
        end
        if choices.assign_workspace, assignin('base', 'eels_zest_stacked_spectra', payload); end
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

    function clearActiveROI()
        for jj = 1:numel(active_listener)
            try delete(active_listener{jj}); catch, end
        end
        active_listener = {};
        if ~isempty(active_roi)
            try if isvalid(active_roi), delete(active_roi); end; catch, end
        end
        active_roi = [];
        if isgraphics(keep_button), keep_button.Enable = 'off'; end
    end

    function finishTool(~, ~)
        if ~isempty(completion_callback), completion_callback(stack); end
        closeTool();
    end

    function closeTool(varargin) %#ok<INUSD>
        clearActiveROI();
        if isgraphics(fig), delete(fig); end
    end
end
