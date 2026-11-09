function options = EELS_ZEST_CutoffDialog(res, stage)
%EELS_ZEST_CUTOFFDIALOG Select and preview sparse-source truncation.

if nargin < 2, stage = 'stage2'; end
stage = lower(char(string(stage)));
options = struct('cancelled', true, 'mode', '', 'cutoff_meV', NaN, ...
    'tail_meV', 20, 'json_file', '', 'json_meta', struct());
if strcmp(stage, 'stage2')
    modes = {'Manual cutoff', 'Materials Project JSON', 'Attention edge + tail'};
else
    modes = {'Manual cutoff', 'Materials Project JSON'};
end

fig = dialog('Name', 'Sparse phonon source cutoff', 'Position', [220 160 920 660], ...
    'Color', 'w', 'WindowStyle', 'modal', 'CloseRequestFcn', @cancelDialog);
ax_source = axes(fig, 'Units', 'pixels', 'Position', [65 355 805 240]);
ax_signal = axes(fig, 'Units', 'pixels', 'Position', [65 105 805 190]);
uicontrol(fig, 'Style', 'text', 'String', 'Cutoff mode', 'Position', [65 620 100 22], ...
    'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
mode_control = uicontrol(fig, 'Style', 'popupmenu', 'String', modes, ...
    'Position', [165 618 220 28], 'Callback', @refreshPreview);
uicontrol(fig, 'Style', 'text', 'String', 'Cutoff / tail (meV)', ...
    'Position', [410 620 130 22], 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
value_control = uicontrol(fig, 'Style', 'edit', 'String', '20', ...
    'Position', [540 618 90 28], 'Callback', @refreshPreview);
json_button = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Choose JSON', ...
    'Position', [650 618 105 28], 'Callback', @chooseJSON);
apply_button = uicontrol(fig, 'Style', 'pushbutton', 'String', 'Apply', ...
    'Position', [770 618 100 28], 'FontWeight', 'bold', 'Callback', @applyDialog);
status = uicontrol(fig, 'Style', 'text', 'String', '', 'Position', [65 70 805 24], ...
    'BackgroundColor', 'w', 'HorizontalAlignment', 'left', 'ForegroundColor', [0.1 0.25 0.45]);
uicontrol(fig, 'Style', 'pushbutton', 'String', 'Cancel', 'Position', [770 25 100 30], ...
    'Callback', @cancelDialog);

json_file = '';
refreshPreview();
uiwait(fig);
if isgraphics(fig), delete(fig); end

    function chooseJSON(~, ~)
        [file, path] = uigetfile('*.json', 'Choose Materials Project phonon JSON');
        if isequal(file, 0), return; end
        json_file = fullfile(path, file);
        mode_control.Value = find(strcmp(modes, 'Materials Project JSON'), 1);
        refreshPreview();
    end

    function [cutoff, meta] = selectedCutoff()
        mode = modes{mode_control.Value};
        value = str2double(value_control.String);
        meta = struct();
        switch mode
            case 'Manual cutoff'
                assert(isfinite(value) && value > 0, 'Enter a positive cutoff in meV.');
                cutoff = value;
            case 'Materials Project JSON'
                assert(~isempty(json_file), 'Choose a Materials Project JSON file.');
                [cutoff, meta] = EELS_ZEST_DataOps.pdos_cutoff_from_json(json_file);
            otherwise
                if ~isfinite(value), value = 20; end
                cutoff = EELS_ZEST_DataOps.attention_cutoff(res, value);
                assert(isfinite(cutoff), 'This Stage 2 result has no valid hard-attention support.');
        end
    end

    function refreshPreview(varargin) %#ok<INUSD>
        try
            [cutoff, meta] = selectedCutoff();
            [signal_full, full] = EELS_ZEST_DataOps.reconstruct_sparse_signal(res, stage, Inf);
            [signal_cut, cut] = EELS_ZEST_DataOps.reconstruct_sparse_signal(res, stage, cutoff);
            cla(ax_source); hold(ax_source, 'on');
            retained = cut.retained_mask;
            stem(ax_source, full.EM_meV(retained), full.PM_stokes_original(retained), ...
                'Color', [0.90 0.35 0.08], 'Marker', 'none');
            stem(ax_source, -full.EM_meV(retained), full.PM_antistokes_original(retained), ...
                'Color', [0.10 0.45 0.85], 'Marker', 'none');
            removed = ~cut.retained_mask;
            if any(removed)
                stem(ax_source, full.EM_meV(removed), full.PM_stokes_original(removed), ...
                    'Color', [0.72 0.72 0.72], 'Marker', 'x');
                stem(ax_source, -full.EM_meV(removed), full.PM_antistokes_original(removed), ...
                    'Color', [0.72 0.72 0.72], 'Marker', 'x');
            end
            xline(ax_source, cutoff, '--k'); xline(ax_source, -cutoff, '--k');
            hold(ax_source, 'off'); grid(ax_source, 'on');
            xlabel(ax_source, 'Mode energy (meV)'); ylabel(ax_source, 'Sparse source intensity');
            title(ax_source, sprintf('%s source cutoff at %.1f meV', stage, cutoff));
            cla(ax_signal);
            plot(ax_signal, full.energy_meV, signal_full, 'Color', [0.55 0.55 0.55]); hold(ax_signal, 'on');
            plot(ax_signal, cut.energy_meV, signal_cut, 'Color', [0.80 0.12 0.10], 'LineWidth', 1.35);
            xline(ax_signal, cutoff, '--k'); xline(ax_signal, -cutoff, '--k');
            hold(ax_signal, 'off'); grid(ax_signal, 'on');
            xlabel(ax_signal, 'Energy loss (meV)'); ylabel(ax_signal, 'Reconstructed phonon signal');
            legend(ax_signal, {'before cutoff', 'after source cutoff'}, 'Location', 'best');
            retained_max = NaN;
            if any(cut.retained_mask), retained_max = max(cut.EM_meV(cut.retained_mask)); end
            status.String = sprintf(['Cutoff %.1f meV | retained %d/%d source positions | ', ...
                'max retained %.1f meV | output keeps the natural convolution tail'], ...
                cutoff, nnz(cut.retained_mask), numel(cut.retained_mask), retained_max);
            status.ForegroundColor = [0.1 0.35 0.18];
            apply_button.Enable = 'on';
            if isfield(meta, 'formula'), status.String = [status.String, ' | ', char(meta.formula)]; end
        catch ME
            cla(ax_source); cla(ax_signal);
            status.String = ME.message;
            status.ForegroundColor = [0.75 0.08 0.08];
            apply_button.Enable = 'off';
        end
    end

    function applyDialog(~, ~)
        try
            [cutoff, meta] = selectedCutoff();
            options.cancelled = false;
            options.mode = modes{mode_control.Value};
            options.cutoff_meV = cutoff;
            options.tail_meV = str2double(value_control.String);
            options.json_file = json_file;
            options.json_meta = meta;
            uiresume(fig);
        catch ME
            status.String = ME.message;
        end
    end

    function cancelDialog(~, ~)
        options.cancelled = true;
        uiresume(fig);
    end
end
