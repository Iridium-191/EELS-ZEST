classdef EELS_ZEST_App < handle
    %EELS_ZEST_APP GUI for EELS-ZEST.
    %
    % EELS-ZEST = Zero-loss Estimation and Sparse Thermometry.

    properties
        UIFigure
        RootDir
        LayoutMode = 'auto'
        IsCompact = false
        DataFile = ''
        RawCube = []
        Ene = []
        DataInfo = struct()

        Stage1SingleResult = []
        Stage2SingleResult = []
        Stage1CubeResult = []
        Stage2CubeResult = []
        CurrentResult = []
        LastPMResult = []
        RLSourceResult = []
        RLSourceCubeResult = []
        RLDeconvResultCube = []
        RLDiagnosticCache = struct('keys', {{}}, 'values', {{}}, 'capacity', 3)
        FitContext = struct()
        LoadedResultFile = ''
        DataSelectedY = []
        DataSelectedX = []
        StackSpectra = []
        DataCursorHandles = {}
        Stage1CursorHandles = {}
        Stage2CursorHandles = {}
        StackCursorHandles = {}
        OwnedFigures = {}
        DataMapImage = []
        DataMapColorbar = []
        LastLogDrawTime = 0

        StopRequested = false

        HeaderLabel
        StatusArea
        TabGroup
        DataTab
        ReferenceTab
        ReferenceTuner
        ReferenceStrategyActive = false
        Stage1Tab
        Stage2Tab
        DiagnosticTab
        PMTab
        AdvancedTab
        DocsTab

        DataMapAxes
        DataSpectrumAxes
        DataEnergyMinEdit
        DataEnergyMaxEdit
        DataRLDeconvButton
        DataRLDiagnosticButton
        DataRLStatusLabel

        Stage1MapAxes
        Stage1SpectrumAxes
        Stage2MapAxes
        Stage2SpectrumAxes
        Stage1StatusArea
        Stage2StatusArea
        Stage1Controls = struct()
        Stage2Controls = struct()
        DiagnosticDropDown
        ViewAxes1
        ViewAxes2
        AdvancedTable
        AdvancedHelpArea

        PMAxesMap
        PMAxesSpectrum
        PMAxesExtra
        PMControls = struct()
        Settings = struct()
    end

    methods
        function app = EELS_ZEST_App(varargin)
            visible = 'on';
            if nargin >= 2
                for ii = 1:2:nargin
                    if strcmpi(varargin{ii}, 'Visible')
                        visible = varargin{ii + 1};
                    elseif strcmpi(varargin{ii}, 'LayoutMode')
                        app.LayoutMode = lower(char(varargin{ii + 1}));
                    end
                end
            end
            app.RootDir = fileparts(mfilename('fullpath'));
            addpath(app.RootDir);
            app.Settings = EELS_PearsonSparse_FullFit.defaults(struct());
            app.Settings.datacube_sample_strategy = 'all';
            app.Settings.datacube_sample_pixels = Inf;
            app.buildUI(visible);
            app.populateAdvancedTable();
            app.logStatus('EELS-ZEST is ready. Load a .mat or .dm3 datacube to begin.');
        end

        function delete(app)
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                app.closeApp();
            end
        end

        function sendSignalToData(app, stage, options)
            if nargin < 2 || isempty(stage), stage = 'stage2'; end
            if nargin < 3, options = []; end
            stage2 = strcmpi(char(string(stage)), 'stage2');
            app.importFitSignalToData(stage2, options);
        end

        function multiplySTA(app, source, manual_temperature_K)
            if nargin < 2 || isempty(source), source = 'stage2'; end
            if nargin < 3, manual_temperature_K = NaN; end
            app.applySTA(source, manual_temperature_K);
        end

        function settings = fittingSettings(app,stage2)
            % Same settings path used by single-spectrum and full-cube fitting.
            settings=app.makeSettings(stage2,'serial');
        end
    end

    methods (Access = private)
        function buildUI(app, visible)
            [fig_pos, resolved_mode] = app.responsiveFigurePosition(app.LayoutMode);
            app.LayoutMode = resolved_mode;
            app.IsCompact = strcmp(resolved_mode, 'compact');
            app.UIFigure = uifigure('Name', 'EELS-ZEST 2.2 | Zero-loss Estimation and Sparse Thermometry', ...
                'Position', fig_pos, 'Color', [0.93 0.97 1.0], 'Visible', visible);
            app.UIFigure.CloseRequestFcn = @(~, ~)app.closeApp();

            root = uigridlayout(app.UIFigure, [3 1]);
            if app.IsCompact
                root.RowHeight = {72, '1x', 58};
                root.Padding = [8 8 8 8];
            else
                root.RowHeight = {92, '1x', 72};
                root.Padding = [14 12 14 12];
            end
            root.ColumnWidth = {'1x'};
            root.RowSpacing = 10;
            root.BackgroundColor = [0.93 0.97 1.0];

            header = uipanel(root, 'BorderType', 'none', 'BackgroundColor', [0.98 0.995 1.0]);
            header.Layout.Row = 1;
            hg = uigridlayout(header, [1 4]);
            if app.IsCompact
                hg.ColumnWidth = {54, '1x', 120, 92};
                hg.Padding = [10 8 10 8];
            else
                hg.ColumnWidth = {70, '1x', 180, 180};
                hg.Padding = [16 10 16 10];
            end
            hg.BackgroundColor = [0.98 0.995 1.0];
            logo_file = fullfile(app.RootDir, 'logo.png');
            if isfile(logo_file)
                try
                    uiimage(hg, 'ImageSource', logo_file, 'ScaleMethod', 'fit');
                catch
                    uilabel(hg, 'Text', 'Z', 'FontSize', 28, 'FontWeight', 'bold', 'FontColor', [0.04 0.32 0.56]);
                end
            else
                uilabel(hg, 'Text', 'Z', 'FontSize', 28, 'FontWeight', 'bold', 'FontColor', [0.04 0.32 0.56]);
            end
            title_label = uilabel(hg, 'Text', 'EELS-ZEST 2.2', 'FontSize', 30, 'FontWeight', 'bold', ...
                'FontColor', [0.03 0.22 0.38]);
            title_label.Layout.Column = 2;
            subtitle = uilabel(hg, 'Text', 'Zero-loss Estimation and Sparse Thermometry', ...
                'FontSize', 14, 'FontColor', [0.14 0.42 0.62]);
            subtitle.Layout.Column = 3;
            app.HeaderLabel = uilabel(hg, 'Text', 'Zero-loss Estimation and Sparse Thermometry', 'HorizontalAlignment', 'right', ...
                'FontSize', 13, 'FontColor', [0.16 0.38 0.56]);
            app.HeaderLabel.Layout.Column = 4;

            app.TabGroup = uitabgroup(root);
            app.TabGroup.Layout.Row = 2;
            app.DataTab = uitab(app.TabGroup, 'Title', 'Data', 'BackgroundColor', [0.92 0.96 1.0]);
            app.ReferenceTab = uitab(app.TabGroup, 'Title', 'Reference tuning');
            app.Stage1Tab = uitab(app.TabGroup, 'Title', 'Stage 1: Sparse Deconvolution', 'BackgroundColor', [0.92 0.96 1.0]);
            app.Stage2Tab = uitab(app.TabGroup, 'Title', 'Stage 2: Thermometry', 'BackgroundColor', [0.92 0.96 1.0]);
            app.DiagnosticTab = uitab(app.TabGroup, 'Title', 'Diagnostics', 'BackgroundColor', [0.92 0.96 1.0]);
            app.PMTab = uitab(app.TabGroup, 'Title', 'Phonon Intensity', 'BackgroundColor', [0.92 0.96 1.0]);
            app.AdvancedTab = uitab(app.TabGroup, 'Title', 'Advanced', 'BackgroundColor', [0.92 0.96 1.0]);
            app.DocsTab = uitab(app.TabGroup, 'Title', 'Quick Help', 'BackgroundColor', [0.92 0.96 1.0]);

            app.buildDataTab();
            app.buildStageTab(false);
            app.buildStageTab(true);
            app.ReferenceTuner=EELS_ReferenceTunerApp(app.Settings,[],@(s)app.applyTunedSettings(s), ...
                visible,app.ReferenceTab,@()app.resetReferenceSettings(),@()app.referenceInput(),@()app.referenceBaseSettings());
            app.buildDiagnosticTab();
            app.buildPMTab();
            app.buildAdvancedTab();
            app.buildDocsTab();

            app.StatusArea = uitextarea(root, 'Editable', 'off', 'FontName', 'Consolas', ...
                'FontSize', 11, 'BackgroundColor', [0.985 0.995 1.0], 'FontColor', [0.04 0.18 0.30]);
            app.StatusArea.Layout.Row = 3;
            app.updateRLDeconvControl();
        end

        function buildDataTab(app)
            if app.IsCompact
                gl = uigridlayout(app.DataTab, [4 4]);
                gl.RowHeight = {242, '1x', '1x', 46};
            else
                gl = uigridlayout(app.DataTab, [3 4]);
                gl.RowHeight = {194, '1x', 50};
            end
            gl.ColumnWidth = {'1x', '1x', '1x', '1x'};
            gl.Padding = [12 10 12 10];
            gl.BackgroundColor = [0.92 0.96 1.0];

            panel = app.controlPanel(gl, 'Input and basic inspection');
            panel.Layout.Row = 1;
            panel.Layout.Column = [1 4];
            app.enableScrollable(panel);
            if app.IsCompact
                pg = uigridlayout(panel, [6 6]);
                pg.RowHeight = {30, 30, 30, 30, 28, 30};
                pg.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x'};
            else
                pg = uigridlayout(panel, [4 12]);
                pg.RowHeight = {32, 32, 30, 32};
                pg.ColumnWidth = {90, 90, 90, 90, 80, 80, 90, 90, 90, 90, '1x', '1x'};
            end
            pg.RowSpacing = 6;
            pg.ColumnSpacing = 10;
            pg.Padding = [10 8 10 8];
            pg.BackgroundColor = [0.84 0.92 0.99];

            if app.IsCompact
                app.techButton(pg, 'Load data', @(~, ~)app.loadData(), 1, [1 2]);
                app.techButton(pg, 'Load result', @(~, ~)app.loadResult(), 1, [3 4]);
                app.techButton(pg, 'Export data', @(~, ~)app.exportData(), 1, [5 6]);
                app.techButton(pg, 'Energy map', @(~, ~)app.plotEnergyMap(), 2, [1 2]);
                app.techButton(pg, 'Sum spectrum', @(~, ~)app.plotIntegratedSpectrum(), 2, [3 4]);
                app.techButton(pg, 'ROI sum', @(~, ~)app.integrateROI(), 2, [5 6]);
                app.techLabel(pg, 'Map window', 3, 1);
                app.DataEnergyMinEdit = uieditfield(pg, 'numeric', 'Value', -200);
                app.DataEnergyMinEdit.Layout.Row = 3; app.DataEnergyMinEdit.Layout.Column = 2;
                app.DataEnergyMaxEdit = uieditfield(pg, 'numeric', 'Value', 200);
                app.DataEnergyMaxEdit.Layout.Row = 3; app.DataEnergyMaxEdit.Layout.Column = 3;
                app.techLabel(pg, 'meV', 3, 4);
                app.techButton(pg, 'Line profile', @(~, ~)app.plotLineProfile(), 3, 5);
                app.techButton(pg, 'Stack spectra', @(~, ~)app.openStackSpectraTool(), 3, 6);
                app.techButton(pg, 'Apply E^2 loss-only', @(~, ~)app.applyE2LossOnly(), 4, [1 2]);
                app.techButton(pg, 'Basic ZLP align', @(~, ~)app.basicZLPAlign(), 4, [3 4]);
                app.techButton(pg, 'SVD denoise', @(~, ~)app.svdDenoise(), 4, [5 6]);
                app.DataRLDeconvButton = app.techButton(pg, 'RL Conv', @(~, ~)app.openRLDeconvFromData(), 5, [1 2]);
                app.DataRLDeconvButton.Tooltip = 'Batch RL-deconvolve every fitted datacube pixel and replace Data with the latent signal cube.';
                app.DataRLDiagnosticButton = app.techButton(pg, 'RL diag', @(~, ~)app.openRLDiagnosticForDataPixel(), 5, [3 4]);
                app.DataRLDiagnosticButton.Tooltip = 'Click a Data-map pixel, then open its RL diagnostic figure.';
                note = uilabel(pg, 'Text', 'RL Conv replaces Data; RL diag opens only for the selected pixel.', ...
                    'FontColor', [0.05 0.25 0.45]);
                note.Layout.Row = 5; note.Layout.Column = [5 6];
                app.DataRLStatusLabel = note;
                app.techButton(pg, 'Multiply STA', @(~, ~)app.applySTA(), 6, [1 2]);
                app.techButton(pg, 'Multi-peak fit', @(~, ~)app.openPeakFitTool(), 6, [3 4]);
            else
                app.techButton(pg, 'Load data', @(~, ~)app.loadData(), 1, [1 2]);
                app.techButton(pg, 'Load result', @(~, ~)app.loadResult(), 1, [3 4]);
                app.techButton(pg, 'Energy map', @(~, ~)app.plotEnergyMap(), 1, [5 6]);
                app.techButton(pg, 'Sum spectrum', @(~, ~)app.plotIntegratedSpectrum(), 1, [7 8]);
                app.techButton(pg, 'ROI sum', @(~, ~)app.integrateROI(), 1, [9 10]);
                app.techButton(pg, 'Export data', @(~, ~)app.exportData(), 1, [11 12]);
                app.techLabel(pg, 'Map window', 2, 1);
                app.DataEnergyMinEdit = uieditfield(pg, 'numeric', 'Value', -200);
                app.DataEnergyMinEdit.Layout.Row = 2; app.DataEnergyMinEdit.Layout.Column = 2;
                app.DataEnergyMaxEdit = uieditfield(pg, 'numeric', 'Value', 200);
                app.DataEnergyMaxEdit.Layout.Row = 2; app.DataEnergyMaxEdit.Layout.Column = 3;
                app.techLabel(pg, 'meV', 2, 4);
                app.techButton(pg, 'Apply E^2 loss-only', @(~, ~)app.applyE2LossOnly(), 2, [5 6]);
                app.techButton(pg, 'Basic ZLP align', @(~, ~)app.basicZLPAlign(), 2, [7 8]);
                app.techButton(pg, 'SVD denoise', @(~, ~)app.svdDenoise(), 2, [9 10]);
                app.techButton(pg, 'Line profile', @(~, ~)app.plotLineProfile(), 3, [1 2]);
                app.techButton(pg, 'Stack spectra', @(~, ~)app.openStackSpectraTool(), 3, [3 4]);
                app.DataRLDeconvButton = app.techButton(pg, 'RL Conv', @(~, ~)app.openRLDeconvFromData(), 3, [5 6]);
                app.DataRLDeconvButton.Tooltip = 'Batch RL-deconvolve every fitted datacube pixel and replace Data with the latent signal cube.';
                app.DataRLDiagnosticButton = app.techButton(pg, 'RL diag', @(~, ~)app.openRLDiagnosticForDataPixel(), 3, [7 8]);
                app.DataRLDiagnosticButton.Tooltip = 'Click a Data-map pixel, then open its RL diagnostic figure.';
                note = uilabel(pg, 'Text', 'RL Conv replaces Data; RL diag opens only for the selected pixel.');
                note.FontColor = [0.05 0.25 0.45];
                note.Layout.Row = 3; note.Layout.Column = [9 12];
                app.DataRLStatusLabel = note;
                app.techButton(pg, 'Multiply STA', @(~, ~)app.applySTA(), 4, [1 2]);
                app.techButton(pg, 'Multi-peak fit', @(~, ~)app.openPeakFitTool(), 4, [3 4]);
            end

            app.DataMapAxes = uiaxes(gl);
            if app.IsCompact
                app.DataMapAxes.Layout.Row = 2; app.DataMapAxes.Layout.Column = [1 4];
            else
                app.DataMapAxes.Layout.Row = 2; app.DataMapAxes.Layout.Column = [1 2];
            end
            app.styleAxes(app.DataMapAxes);
            title(app.DataMapAxes, 'Energy-filtered map');

            app.DataSpectrumAxes = uiaxes(gl);
            if app.IsCompact
                app.DataSpectrumAxes.Layout.Row = 3; app.DataSpectrumAxes.Layout.Column = [1 4];
            else
                app.DataSpectrumAxes.Layout.Row = 2; app.DataSpectrumAxes.Layout.Column = [3 4];
            end
            app.styleAxes(app.DataSpectrumAxes);
            title(app.DataSpectrumAxes, 'Integrated spectrum');

            help = uilabel(gl, 'Text', ...
                'Tip: use Stage 1 for ZLP + sparse-convolution signal recovery. Use Stage 2 only when the Stokes/anti-Stokes evidence is meaningful.', ...
                'FontColor', [0.05 0.23 0.42], 'FontSize', 12);
            if app.IsCompact
                help.Layout.Row = 4;
            else
                help.Layout.Row = 3;
            end
            help.Layout.Column = [1 4];
        end

        function buildStageTab(app, stage2)
            if stage2
                tab = app.Stage2Tab;
                title_text = 'Stage 2: attention-window thermometry';
            else
                tab = app.Stage1Tab;
                title_text = 'Stage 1: ZLP + sparse convolution';
            end
            if app.IsCompact
                app.buildStageTabCompact(stage2, tab, title_text);
                return;
            end
            gl = uigridlayout(tab, [3 4]);
            gl.RowHeight = {284, '1x', 82};
            gl.ColumnWidth = {'1x', '1x', '1x', '1x'};
            gl.Padding = [12 10 12 10];
            gl.BackgroundColor = [0.92 0.96 1.0];

            panel = app.controlPanel(gl, title_text);
            panel.Layout.Row = 1; panel.Layout.Column = [1 4];
            app.enableScrollable(panel);
            pg = uigridlayout(panel, [6 12]);
            pg.RowHeight = {32, 32, 32, 32, 32, 32};
            pg.ColumnWidth = {64, '1x', 64, '1x', 70, '1x', 70, '1x', 78, '1x', 86, '1x'};
            pg.RowSpacing = 6;
            pg.ColumnSpacing = 10;
            pg.Padding = [10 8 10 8];
            pg.BackgroundColor = [0.84 0.92 0.99];

            app.techLabel(pg, 'Y', 1, 1); py = uispinner(pg, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            py.Layout.Row = 1; py.Layout.Column = 2;
            app.techLabel(pg, 'X', 1, 3); px = uispinner(pg, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            px.Layout.Row = 1; px.Layout.Column = 4;
            app.techLabel(pg, 'Y bin', 1, 5); by = uispinner(pg, 'Limits', [1 64], 'Value', 1, 'Step', 1);
            by.Layout.Row = 1; by.Layout.Column = 6;
            app.techLabel(pg, 'X bin', 1, 7); bx = uispinner(pg, 'Limits', [1 64], 'Value', 1, 'Step', 1);
            bx.Layout.Row = 1; bx.Layout.Column = 8;
            app.techLabel(pg, 'Energy bin', 1, 9); eb = uispinner(pg, 'Limits', [1 32], 'Value', 1, 'Step', 1);
            eb.Layout.Row = 1; eb.Layout.Column = 10;
            app.techLabel(pg, 'Workers', 1, 11); workers = uispinner(pg, 'Limits', [1 128], 'Value', EELS_Parallel.default_workers(), 'Step', 1);
            workers.Layout.Row = 1; workers.Layout.Column = 12;

            app.techLabel(pg, 'Crop min', 2, 1); cmin = uieditfield(pg, 'numeric', 'Value', -250);
            cmin.Layout.Row = 2; cmin.Layout.Column = 2;
            app.techLabel(pg, 'Crop max', 2, 3); cmax = uieditfield(pg, 'numeric', 'Value', 250);
            cmax.Layout.Row = 2; cmax.Layout.Column = 4;
            app.techButton(pg, 'Inherit data window', @(~, ~)app.inheritWindow(stage2), 2, [5 6]);
            resume = uicheckbox(pg, 'Text', 'Resume checkpoint', 'Value', false, 'FontColor', [0.05 0.20 0.35]);
            resume.Layout.Row = 2; resume.Layout.Column = [8 9];
            checkpoint = uicheckbox(pg, 'Text', 'Save checkpoint', 'Value', false, 'FontColor', [0.05 0.20 0.35]);
            checkpoint.Layout.Row = 2; checkpoint.Layout.Column = [10 12];

            app.techLabel(pg, 'Output', 3, 1);
            out = uieditfield(pg, 'text', 'Value', app.defaultOutputDir(stage2));
            out.Layout.Row = 3; out.Layout.Column = [2 10];
            app.techButton(pg, 'Browse', @(~, ~)app.pickOutputDir(out), 3, [11 12]);

            if stage2
                app.techLabel(pg, 'Algorithm', 4, 1);
                alg = uidropdown(pg, 'Items', {'grid_scan'}, 'Value', 'grid_scan');
                alg.Layout.Row = 4; alg.Layout.Column = [2 4];
                spatial = uicheckbox(pg, 'Text', 'Hierarchical spatial model', 'Value', true, 'FontColor', [0.05 0.20 0.35]);
                spatial.Layout.Row = 4; spatial.Layout.Column = [5 7];
                classPrune = uicheckbox(pg, 'Text', 'Class EM pruning', 'Value', true, 'FontColor', [0.05 0.20 0.35]);
                classPrune.Layout.Row = 6; classPrune.Layout.Column = [9 10];
                forceRaw = uicheckbox(pg, 'Text', 'Raw Stage2', 'Value', false, 'FontColor', [0.05 0.20 0.35]);
                forceRaw.Tooltip = 'Disable the Stage 2 shape guard and show the raw hard-BE fit instead of falling back to Stage 1.';
                forceRaw.Layout.Row = 6; forceRaw.Layout.Column = [11 12];
                app.techLabel(pg, 'Map', 4, 8);
            else
                alg = [];
                spatial = [];
                classPrune = [];
                forceRaw = [];
                app.techLabel(pg, 'Map', 4, 8);
            end

            app.techButton(pg, 'Fit pixel', @(~, ~)app.runSingle(stage2), 5, [1 2]);
            app.techButton(pg, 'Run serial', @(~, ~)app.runDatacube(stage2, 'serial'), 5, [3 4]);
            app.techButton(pg, 'Run parallel', @(~, ~)app.runDatacube(stage2, 'parallel'), 5, [5 6]);
            app.techButton(pg, 'Pause', @(~, ~)app.requestPause(), 5, 7);
            app.techButton(pg, 'Pick', @(~, ~)app.pickPixelWithCursor(stage2), 5, 8);
            app.techButton(pg, 'Diagnostic', @(~, ~)app.openDiagnostic(stage2), 5, [9 10]);
            app.techButton(pg, 'Export', @(~, ~)app.exportResult(stage2), 5, [11 12]);

            app.techButton(pg, 'Load result', @(~, ~)app.loadResult(), 6, [1 2]);
            app.techButton(pg, 'Send signal to Data', @(~, ~)app.importFitSignalToData(stage2), 6, [3 5]);
            if stage2
                app.techButton(pg, 'Reference tuner', @(~, ~)app.openReferenceTuner(), 6, [6 8]);
                mapDD = uidropdown(pg, 'Items', app.mapItems(), 'Value', 'T data valid', ...
                    'ValueChangedFcn', @(~, ~)app.plotMap(true));
                mapDD.Layout.Row = 4; mapDD.Layout.Column = [9 12];
            else
                mapDD = uidropdown(pg, 'Items', app.mapItems(), 'Value', 'R2', ...
                    'ValueChangedFcn', @(~, ~)app.plotMap(false));
                mapDD.Layout.Row = 4; mapDD.Layout.Column = [9 12];
            end

            if stage2
                app.Stage2MapAxes = uiaxes(gl); app.Stage2MapAxes.Layout.Row = 2; app.Stage2MapAxes.Layout.Column = [1 2];
                app.Stage2SpectrumAxes = uiaxes(gl); app.Stage2SpectrumAxes.Layout.Row = 2; app.Stage2SpectrumAxes.Layout.Column = [3 4];
                app.Stage2StatusArea = uitextarea(gl, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 10);
                app.Stage2StatusArea.Layout.Row = 3; app.Stage2StatusArea.Layout.Column = [1 4];
                app.styleAxes(app.Stage2MapAxes); app.styleAxes(app.Stage2SpectrumAxes);
                app.Stage2MapAxes.ButtonDownFcn = @(~, ~)app.pickPixelFromMap(true);
                app.Stage2Controls = struct('py', py, 'px', px, 'by', by, 'bx', bx, 'eb', eb, ...
                    'cmin', cmin, 'cmax', cmax, 'workers', workers, ...
                    'resume', resume, 'checkpoint', checkpoint, 'output', out, ...
                    'algorithm', alg, 'spatial', spatial, 'classPrune', classPrune, ...
                    'forceRawStage2', forceRaw, 'map', mapDD);
                title(app.Stage2MapAxes, 'Temperature / reliability map');
                title(app.Stage2SpectrumAxes, 'Stage 2 fit');
            else
                app.Stage1MapAxes = uiaxes(gl); app.Stage1MapAxes.Layout.Row = 2; app.Stage1MapAxes.Layout.Column = [1 2];
                app.Stage1SpectrumAxes = uiaxes(gl); app.Stage1SpectrumAxes.Layout.Row = 2; app.Stage1SpectrumAxes.Layout.Column = [3 4];
                app.Stage1StatusArea = uitextarea(gl, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 10);
                app.Stage1StatusArea.Layout.Row = 3; app.Stage1StatusArea.Layout.Column = [1 4];
                app.styleAxes(app.Stage1MapAxes); app.styleAxes(app.Stage1SpectrumAxes);
                app.Stage1MapAxes.ButtonDownFcn = @(~, ~)app.pickPixelFromMap(false);
                app.Stage1Controls = struct('py', py, 'px', px, 'by', by, 'bx', bx, 'eb', eb, ...
                    'cmin', cmin, 'cmax', cmax, 'workers', workers, ...
                    'resume', resume, 'checkpoint', checkpoint, 'output', out, 'map', mapDD);
                title(app.Stage1MapAxes, 'Fit quality map');
                title(app.Stage1SpectrumAxes, 'Stage 1 fit');
            end
        end

        function buildStageTabCompact(app, stage2, tab, title_text)
            gl = uigridlayout(tab, [4 4]);
            gl.RowHeight = {320, '1x', '1x', 74};
            gl.ColumnWidth = {'1x', '1x', '1x', '1x'};
            gl.Padding = [10 8 10 8];
            gl.BackgroundColor = [0.92 0.96 1.0];

            panel = app.controlPanel(gl, [title_text, ' (compact)']);
            panel.Layout.Row = 1; panel.Layout.Column = [1 4];
            app.enableScrollable(panel);
            pg = uigridlayout(panel, [8 8]);
            pg.RowHeight = {28, 28, 28, 28, 28, 28, 28, 28};
            pg.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            pg.RowSpacing = 6;
            pg.ColumnSpacing = 8;
            pg.Padding = [10 8 10 8];
            pg.BackgroundColor = [0.84 0.92 0.99];

            app.techLabel(pg, 'Y', 1, 1); py = uispinner(pg, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            py.Layout.Row = 1; py.Layout.Column = 2;
            app.techLabel(pg, 'X', 1, 3); px = uispinner(pg, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            px.Layout.Row = 1; px.Layout.Column = 4;
            app.techLabel(pg, 'Y bin', 1, 5); by = uispinner(pg, 'Limits', [1 64], 'Value', 1, 'Step', 1);
            by.Layout.Row = 1; by.Layout.Column = 6;
            app.techLabel(pg, 'X bin', 1, 7); bx = uispinner(pg, 'Limits', [1 64], 'Value', 1, 'Step', 1);
            bx.Layout.Row = 1; bx.Layout.Column = 8;

            app.techLabel(pg, 'Energy bin', 2, 1); eb = uispinner(pg, 'Limits', [1 32], 'Value', 1, 'Step', 1);
            eb.Layout.Row = 2; eb.Layout.Column = 2;
            app.techLabel(pg, 'Workers', 2, 3); workers = uispinner(pg, 'Limits', [1 128], 'Value', EELS_Parallel.default_workers(), 'Step', 1);
            workers.Layout.Row = 2; workers.Layout.Column = 4;
            resume = uicheckbox(pg, 'Text', 'Resume', 'Value', false, 'FontColor', [0.05 0.20 0.35]);
            resume.Layout.Row = 2; resume.Layout.Column = [5 6];
            checkpoint = uicheckbox(pg, 'Text', 'Checkpoint', 'Value', false, 'FontColor', [0.05 0.20 0.35]);
            checkpoint.Layout.Row = 2; checkpoint.Layout.Column = [7 8];

            app.techLabel(pg, 'Crop min', 3, 1); cmin = uieditfield(pg, 'numeric', 'Value', -250);
            cmin.Layout.Row = 3; cmin.Layout.Column = 2;
            app.techLabel(pg, 'Crop max', 3, 3); cmax = uieditfield(pg, 'numeric', 'Value', 250);
            cmax.Layout.Row = 3; cmax.Layout.Column = 4;
            app.techButton(pg, 'Inherit data window', @(~, ~)app.inheritWindow(stage2), 3, [5 8]);

            app.techLabel(pg, 'Output', 4, 1);
            out = uieditfield(pg, 'text', 'Value', app.defaultOutputDir(stage2));
            out.Layout.Row = 4; out.Layout.Column = [2 7];
            app.techButton(pg, 'Browse', @(~, ~)app.pickOutputDir(out), 4, 8);

            if stage2
                app.techLabel(pg, 'Algorithm', 5, 1);
                alg = uidropdown(pg, 'Items', {'grid_scan'}, 'Value', 'grid_scan');
                alg.Layout.Row = 5; alg.Layout.Column = [2 4];
                spatial = uicheckbox(pg, 'Text', 'Hierarchical spatial model', 'Value', true, 'FontColor', [0.05 0.20 0.35]);
                spatial.Layout.Row = 5; spatial.Layout.Column = [5 6];
                classPrune = uicheckbox(pg, 'Text', 'Class EM prune', 'Value', true, 'FontColor', [0.05 0.20 0.35]);
                classPrune.Layout.Row = 5; classPrune.Layout.Column = [7 8];
                forceRaw = uicheckbox(pg, 'Text', 'Raw Stage2', 'Value', false, 'FontColor', [0.05 0.20 0.35]);
                forceRaw.Tooltip = 'Disable the Stage 2 shape guard and show the raw hard-BE fit instead of falling back to Stage 1.';
                forceRaw.Layout.Row = 7; forceRaw.Layout.Column = [4 5];
            else
                alg = [];
                spatial = [];
                classPrune = [];
                forceRaw = [];
            end

            app.techButton(pg, 'Fit pixel', @(~, ~)app.runSingle(stage2), 6, [1 2]);
            app.techButton(pg, 'Run serial', @(~, ~)app.runDatacube(stage2, 'serial'), 6, [3 4]);
            app.techButton(pg, 'Run parallel', @(~, ~)app.runDatacube(stage2, 'parallel'), 6, [5 6]);
            app.techButton(pg, 'Pause', @(~, ~)app.requestPause(), 6, 7);
            app.techButton(pg, 'Pick', @(~, ~)app.pickPixelWithCursor(stage2), 6, 8);
            app.techButton(pg, 'Diagnostic', @(~, ~)app.openDiagnostic(stage2), 7, [1 3]);
            app.techButton(pg, 'Export', @(~, ~)app.exportResult(stage2), 7, [6 8]);

            app.techButton(pg, 'Load result', @(~, ~)app.loadResult(), 8, [1 2]);
            app.techButton(pg, 'Send signal to Data', @(~, ~)app.importFitSignalToData(stage2), 8, [3 4]);
            if stage2
                app.techButton(pg, 'Reference tuner', @(~, ~)app.openReferenceTuner(), 8, [5 6]);
                mapDD = uidropdown(pg, 'Items', app.mapItems(), 'Value', 'T data valid', ...
                    'ValueChangedFcn', @(~, ~)app.plotMap(true));
                mapDD.Layout.Row = 8; mapDD.Layout.Column = [7 8];
            else
                mapDD = uidropdown(pg, 'Items', app.mapItems(), 'Value', 'R2', ...
                    'ValueChangedFcn', @(~, ~)app.plotMap(false));
                mapDD.Layout.Row = 8; mapDD.Layout.Column = [5 8];
            end

            if stage2
                app.Stage2MapAxes = uiaxes(gl); app.Stage2MapAxes.Layout.Row = 2; app.Stage2MapAxes.Layout.Column = [1 4];
                app.Stage2SpectrumAxes = uiaxes(gl); app.Stage2SpectrumAxes.Layout.Row = 3; app.Stage2SpectrumAxes.Layout.Column = [1 4];
                app.Stage2StatusArea = uitextarea(gl, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 10);
                app.Stage2StatusArea.Layout.Row = 4; app.Stage2StatusArea.Layout.Column = [1 4];
                app.styleAxes(app.Stage2MapAxes); app.styleAxes(app.Stage2SpectrumAxes);
                app.Stage2MapAxes.ButtonDownFcn = @(~, ~)app.pickPixelFromMap(true);
                app.Stage2Controls = struct('py', py, 'px', px, 'by', by, 'bx', bx, 'eb', eb, ...
                    'cmin', cmin, 'cmax', cmax, 'workers', workers, ...
                    'resume', resume, 'checkpoint', checkpoint, 'output', out, ...
                    'algorithm', alg, 'spatial', spatial, 'classPrune', classPrune, ...
                    'forceRawStage2', forceRaw, 'map', mapDD);
                title(app.Stage2MapAxes, 'Temperature / reliability map');
                title(app.Stage2SpectrumAxes, 'Stage 2 fit');
            else
                app.Stage1MapAxes = uiaxes(gl); app.Stage1MapAxes.Layout.Row = 2; app.Stage1MapAxes.Layout.Column = [1 4];
                app.Stage1SpectrumAxes = uiaxes(gl); app.Stage1SpectrumAxes.Layout.Row = 3; app.Stage1SpectrumAxes.Layout.Column = [1 4];
                app.Stage1StatusArea = uitextarea(gl, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 10);
                app.Stage1StatusArea.Layout.Row = 4; app.Stage1StatusArea.Layout.Column = [1 4];
                app.styleAxes(app.Stage1MapAxes); app.styleAxes(app.Stage1SpectrumAxes);
                app.Stage1MapAxes.ButtonDownFcn = @(~, ~)app.pickPixelFromMap(false);
                app.Stage1Controls = struct('py', py, 'px', px, 'by', by, 'bx', bx, 'eb', eb, ...
                    'cmin', cmin, 'cmax', cmax, 'workers', workers, ...
                    'resume', resume, 'checkpoint', checkpoint, 'output', out, 'map', mapDD);
                title(app.Stage1MapAxes, 'Fit quality map');
                title(app.Stage1SpectrumAxes, 'Stage 1 fit');
            end
        end

        function buildDiagnosticTab(app)
            gl = uigridlayout(app.DiagnosticTab, [3 4]);
            if app.IsCompact
                gl.RowHeight = {140, '1x', '1x'};
                gl.ColumnWidth = {'1x', '1x', '1x', '1x'};
            else
                gl.RowHeight = {92, '1x', '1x'};
                gl.ColumnWidth = {280, '1x', '1x', 250};
            end
            gl.Padding = [12 10 12 10];
            gl.BackgroundColor = [0.92 0.96 1.0];

            panel = app.controlPanel(gl, 'Selected-spectrum diagnostics');
            panel.Layout.Row = 1; panel.Layout.Column = [1 4];
            app.enableScrollable(panel);
            if app.IsCompact
                pg = uigridlayout(panel, [3 3]);
                pg.RowHeight = {30, 30, 30};
                pg.ColumnWidth = {'1x', '1x', '1x'};
            else
                pg = uigridlayout(panel, [1 6]);
                pg.ColumnWidth = {190, 130, 130, 130, 130, '1x'};
            end
            pg.Padding = [10 8 10 8];
            pg.BackgroundColor = [0.84 0.92 0.99];
            app.DiagnosticDropDown = uidropdown(pg, 'Items', { ...
                'Full diagnostic figure', 'Raw vs total fit', 'ZLP components', ...
                'Sparse convolution / vib', 'Residual and CI', 'Source-intensity sticks', ...
                'Spectrum ln(Iloss/Igain)', 'Hard-BE ln ratio', 'Attention curve', 'Loss budget', 'Log signal vs log ZLP'}, ...
                'Value', 'Raw vs total fit', 'ValueChangedFcn', @(~, ~)app.plotDiagnosticSelection());
            if app.IsCompact
                app.DiagnosticDropDown.Layout.Row = 1; app.DiagnosticDropDown.Layout.Column = [1 3];
                app.techButton(pg, 'Refresh', @(~, ~)app.plotDiagnosticSelection(), 2, 1);
                app.techButton(pg, 'Full figure', @(~, ~)app.openFullDiagnostic(), 2, 2);
                app.techButton(pg, 'Load result', @(~, ~)app.loadResult(), 2, 3);
                app.techButton(pg, 'Export selected', @(~, ~)app.exportCurrentResult(), 3, 1);
                note = uilabel(pg, 'Text', 'Click a map pixel or load a result, then choose a diagnostic item.');
                note.Layout.Row = 3; note.Layout.Column = [2 3];
            else
                app.DiagnosticDropDown.Layout.Row = 1; app.DiagnosticDropDown.Layout.Column = 1;
                app.techButton(pg, 'Refresh', @(~, ~)app.plotDiagnosticSelection(), 1, 2);
                app.techButton(pg, 'Full figure', @(~, ~)app.openFullDiagnostic(), 1, 3);
                app.techButton(pg, 'Load result', @(~, ~)app.loadResult(), 1, 4);
                app.techButton(pg, 'Export selected', @(~, ~)app.exportCurrentResult(), 1, 5);
                note = uilabel(pg, 'Text', 'Click a map pixel or load a result, then choose a diagnostic item.');
                note.Layout.Row = 1; note.Layout.Column = 6;
            end
            note.FontColor = [0.05 0.25 0.45];

            app.ViewAxes1 = uiaxes(gl);
            app.ViewAxes1.Layout.Row = [2 3]; app.ViewAxes1.Layout.Column = [1 2];
            app.ViewAxes2 = uiaxes(gl);
            app.ViewAxes2.Layout.Row = [2 3]; app.ViewAxes2.Layout.Column = [3 4];
            app.styleAxes(app.ViewAxes1); app.styleAxes(app.ViewAxes2);
            title(app.ViewAxes1, 'Diagnostic panel');
            title(app.ViewAxes2, 'Secondary view');
        end

        function buildAdvancedTab(app)
            gl = uigridlayout(app.AdvancedTab, [2 4]);
            if app.IsCompact
                gl.RowHeight = {'1x', 152};
            else
                gl.RowHeight = {'1x', 118};
            end
            gl.ColumnWidth = {'1x', '1x', 150, 150};
            gl.Padding = [12 10 12 10];
            gl.BackgroundColor = [0.92 0.96 1.0];
            app.AdvancedTable = uitable(gl, 'ColumnName', {'Setting', 'Description', 'Value', 'Advice'}, ...
                'ColumnEditable', [false false true false], 'CellEditCallback', @(~, ~)app.applyAdvancedTable());
            app.AdvancedTable.Layout.Row = 1; app.AdvancedTable.Layout.Column = [1 4];
            app.AdvancedHelpArea = uitextarea(gl, 'Editable', 'off', 'FontSize', 11, ...
                'BackgroundColor', [0.98 0.995 1.0]);
            app.AdvancedHelpArea.Layout.Row = 2; app.AdvancedHelpArea.Layout.Column = [1 2];
            app.AdvancedHelpArea.Value = { ...
                'Advanced parameters are applied to the next single-pixel/datacube run.', ...
                'Use defaults first. Stage 2 uses paired attention for hard BE and independent soft source intensities elsewhere.', ...
                'Full datacube Stage 1 can estimate a sparse phonon-mode grid from one representative spectrum, then reuse it for the cube.', ...
                'Change ZLP width/guards for low-frequency leakage, mode-energy density/smoothing for spiky source intensity, and attention parameters for missed peaks.', ...
                'Datacube options in this table are overwritten by the visible GUI controls when you run from the tabs.'};
            b = uibutton(gl, 'Text', 'Apply parameters', 'ButtonPushedFcn', @(~, ~)app.applyAdvancedTable());
            b.Layout.Row = 2; b.Layout.Column = 3;
            b = uibutton(gl, 'Text', 'Reset defaults', 'ButtonPushedFcn', @(~, ~)app.resetAdvanced());
            b.Layout.Row = 2; b.Layout.Column = 4;
        end

        function buildPMTab(app)
            if app.IsCompact
                gl = uigridlayout(app.PMTab, [4 4]);
                gl.RowHeight = {214, '1x', '1x', '1x'};
            else
                gl = uigridlayout(app.PMTab, [3 4]);
                gl.RowHeight = {148, '1x', '1x'};
            end
            gl.ColumnWidth = {'1x', '1x', '1x', '1x'};
            gl.Padding = [12 10 12 10];
            gl.BackgroundColor = [0.92 0.96 1.0];

            panel = app.controlPanel(gl, 'Sparse phonon-mode intensity');
            panel.Layout.Row = 1; panel.Layout.Column = [1 4];
            app.enableScrollable(panel);
            if app.IsCompact
                pg = uigridlayout(panel, [5 6]);
                pg.RowHeight = {28, 28, 28, 28, 28};
                pg.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x'};
            else
                pg = uigridlayout(panel, [3 12]);
                pg.RowHeight = {28, 28, 28};
                pg.ColumnWidth = {90, 90, 76, 76, 76, 76, 86, 86, 86, 86, 86, '1x'};
            end
            pg.Padding = [10 8 10 8];
            pg.BackgroundColor = [0.84 0.92 0.99];
            if app.IsCompact
                app.techLabel(pg, 'Source', 1, 1);
                src = uidropdown(pg, 'Items', {'stage2', 'stage1', 'current'}, 'Value', 'stage2');
                src.Layout.Row = 1; src.Layout.Column = 2;
                app.techLabel(pg, 'EM min', 1, 3); emin = uieditfield(pg, 'numeric', 'Value', -250);
                emin.Layout.Row = 1; emin.Layout.Column = 4;
                app.techLabel(pg, 'EM max', 1, 5); emax = uieditfield(pg, 'numeric', 'Value', 250);
                emax.Layout.Row = 1; emax.Layout.Column = 6;
                app.techLabel(pg, 'Step', 2, 1); step = uieditfield(pg, 'numeric', 'Value', 1);
                step.Layout.Row = 2; step.Layout.Column = 2;
                app.techLabel(pg, 'Window', 2, 3); wmin = uieditfield(pg, 'numeric', 'Value', 40);
                wmin.Layout.Row = 2; wmin.Layout.Column = 4; wmax = uieditfield(pg, 'numeric', 'Value', 120);
                wmax.Layout.Row = 2; wmax.Layout.Column = 5;
                app.techLabel(pg, 'meV', 2, 6);
                app.techLabel(pg, 'NMF k', 3, 1); nk = uispinner(pg, 'Limits', [1 12], 'Value', 3, 'Step', 1);
                nk.Layout.Row = 3; nk.Layout.Column = 2;
                app.techLabel(pg, 'Cluster k', 3, 3); ck = uispinner(pg, 'Limits', [2 16], 'Value', 3, 'Step', 1);
                ck.Layout.Row = 3; ck.Layout.Column = 4;
                smooth = uicheckbox(pg, 'Text', 'Smooth display only', 'Value', true, 'FontColor', [0.05 0.20 0.35]);
                smooth.Layout.Row = 3; smooth.Layout.Column = [5 6];
                app.techButton(pg, 'Build intensity cube', @(~, ~)app.buildPMAnalysis(), 4, [1 2]);
                app.techButton(pg, 'Intensity integral map', @(~, ~)app.plotPMIntegralMap(), 4, [3 4]);
                app.techButton(pg, 'Export intensity cube', @(~, ~)app.exportPMCube(), 4, [5 6]);
                app.techButton(pg, 'Export selected window', @(~, ~)app.exportSelectedPM(), 5, [1 2]);
                app.techButton(pg, 'Intensity NMF', @(~, ~)app.runPMNMF(), 5, [3 4]);
                app.techButton(pg, 'Cluster intensities', @(~, ~)app.runPMCluster(), 5, [5 6]);
            else
                app.techLabel(pg, 'Source', 1, 1);
                src = uidropdown(pg, 'Items', {'stage2', 'stage1', 'current'}, 'Value', 'stage2');
                src.Layout.Row = 1; src.Layout.Column = 2;
                app.techLabel(pg, 'EM min', 1, 3); emin = uieditfield(pg, 'numeric', 'Value', -250);
                emin.Layout.Row = 1; emin.Layout.Column = 4;
                app.techLabel(pg, 'EM max', 1, 5); emax = uieditfield(pg, 'numeric', 'Value', 250);
                emax.Layout.Row = 1; emax.Layout.Column = 6;
                app.techLabel(pg, 'Step', 1, 7); step = uieditfield(pg, 'numeric', 'Value', 1);
                step.Layout.Row = 1; step.Layout.Column = 8;
                app.techLabel(pg, 'Window', 1, 9); wmin = uieditfield(pg, 'numeric', 'Value', 40);
                wmin.Layout.Row = 1; wmin.Layout.Column = 10; wmax = uieditfield(pg, 'numeric', 'Value', 120);
                wmax.Layout.Row = 1; wmax.Layout.Column = 11;

                app.techLabel(pg, 'NMF k', 2, 1); nk = uispinner(pg, 'Limits', [1 12], 'Value', 3, 'Step', 1);
                nk.Layout.Row = 2; nk.Layout.Column = 2;
                app.techLabel(pg, 'Cluster k', 2, 3); ck = uispinner(pg, 'Limits', [2 16], 'Value', 3, 'Step', 1);
                ck.Layout.Row = 2; ck.Layout.Column = 4;
                smooth = uicheckbox(pg, 'Text', 'Smooth display only', 'Value', true, 'FontColor', [0.05 0.20 0.35]);
                smooth.Layout.Row = 2; smooth.Layout.Column = [5 7];
                note = uilabel(pg, 'Text', 'Signed mode energy: gain/anti-Stokes < 0, loss/Stokes > 0. Interpolation gives a DOS-like source-space view, not an absolutely calibrated DOS.', ...
                    'FontColor', [0.05 0.25 0.45]);
                note.Layout.Row = 2; note.Layout.Column = [8 12];

                app.techButton(pg, 'Build intensity cube', @(~, ~)app.buildPMAnalysis(), 3, [1 2]);
                app.techButton(pg, 'Intensity integral map', @(~, ~)app.plotPMIntegralMap(), 3, [3 4]);
                app.techButton(pg, 'Export intensity cube', @(~, ~)app.exportPMCube(), 3, [5 6]);
                app.techButton(pg, 'Export selected window', @(~, ~)app.exportSelectedPM(), 3, [7 8]);
                app.techButton(pg, 'Intensity NMF', @(~, ~)app.runPMNMF(), 3, [9 10]);
                app.techButton(pg, 'Cluster intensities', @(~, ~)app.runPMCluster(), 3, [11 12]);
            end

            app.PMAxesMap = uiaxes(gl);
            app.PMAxesSpectrum = uiaxes(gl);
            app.PMAxesExtra = uiaxes(gl);
            if app.IsCompact
                app.PMAxesMap.Layout.Row = 2; app.PMAxesMap.Layout.Column = [1 4];
                app.PMAxesSpectrum.Layout.Row = 3; app.PMAxesSpectrum.Layout.Column = [1 4];
                app.PMAxesExtra.Layout.Row = 4; app.PMAxesExtra.Layout.Column = [1 4];
            else
                app.PMAxesMap.Layout.Row = [2 3]; app.PMAxesMap.Layout.Column = [1 2];
                app.PMAxesSpectrum.Layout.Row = 2; app.PMAxesSpectrum.Layout.Column = [3 4];
                app.PMAxesExtra.Layout.Row = 3; app.PMAxesExtra.Layout.Column = [3 4];
            end
            app.styleAxes(app.PMAxesMap); app.styleAxes(app.PMAxesSpectrum); app.styleAxes(app.PMAxesExtra);
            app.PMControls = struct('source', src, 'emin', emin, 'emax', emax, 'step', step, ...
                'wmin', wmin, 'wmax', wmax, 'nk', nk, 'ck', ck, 'smooth', smooth);
        end

        function buildDocsTab(app)
            gl = uigridlayout(app.DocsTab, [1 1]);
            gl.BackgroundColor = [0.92 0.96 1.0];
            txt = uitextarea(gl, 'Editable', 'off', 'FontName', 'Consolas', 'FontSize', 12, ...
                'BackgroundColor', [0.97 0.99 1.0]);
            txt.Value = { ...
                'EELS-ZEST', ...
                '', ...
                'Workflow:', ...
                '1. Load a datacube (.mat EELS object / .mat cube / .dm3).', ...
                '2. Stage 1 reconstructs the spectrum using ZLP + sparse convolution without an independent background term.', ...
                '3. Stage 2 reuses Stage 1 source intensities and applies attention-window Bose-Einstein thermometry.', ...
                '4. Phonon Intensity converts fitted coefficients to a signed mode-energy cube for DOS-like maps, NMF, and clustering.', ...
                '', ...
                'Important parameters:', ...
                '- Crop min/max: energy range passed to the fitter.', ...
                '- Y/X/Energy bin: simple averaging before fitting.', ...
                '- Parallel: uses MATLAB Parallel Computing Toolbox when available.', ...
                '- Resume checkpoint: continues from datacube_fit_checkpoint.mat in the output folder.', ...
                '- Stage 2: grid_scan uses attention hard-BE pairs; non-attention mode energies retain independent soft source intensities.', ...
                '', ...
                'See README.md and docs/LOSS_AND_ATTENTION_THEORY_v1.3.md for equations and symbol definitions.'};
        end

        function loadData(app)
            [file, path] = uigetfile({'*.mat;*.dm3', 'EELS datacube (*.mat, *.dm3)'; '*.*', 'All files'}, 'Load EELS data');
            if isequal(file, 0), return; end
            app.DataFile = fullfile(path, file);
            d = uiprogressdlg(app.UIFigure, 'Title', 'Loading data', 'Message', app.DataFile, 'Indeterminate', 'on');
            cleanup = onCleanup(@()app.closeProgressDialog(d));
            [cube, ene, info] = EELS_PearsonSparse_FullFit.load_datacube_mat(app.DataFile);
            app.RawCube = cube;
            app.Ene = ene(:);
            app.DataInfo = info;
            app.FitContext = EELS_ZEST_DataOps.fit_context_from_info(info, size(cube, [1 2]));
            app.LoadedResultFile = '';
            if isfield(app.FitContext, 'result_file') && isfile(app.FitContext.result_file)
                app.LoadedResultFile = app.FitContext.result_file;
            end
            app.Stage1CubeResult = [];
            app.Stage2CubeResult = [];
            app.Stage1SingleResult = [];
            app.Stage2SingleResult = [];
            app.CurrentResult = [];
            app.LastPMResult = [];
            app.RLSourceResult = [];
            app.RLSourceCubeResult = [];
            app.RLDeconvResultCube = [];
            app.RLDiagnosticCache = struct('keys', {{}}, 'values', {{}}, 'capacity', 3);
            app.DataSelectedY = [];
            app.DataSelectedX = [];
            app.updateRLDeconvControl();
            app.updatePixelLimits();
            app.plotEnergyMap();
            app.plotIntegratedSpectrum();
            app.logStatus(sprintf('Loaded %s | cube %s | energy channels %d', file, mat2str(size(cube)), numel(ene)));
        end

        function loadResult(app)
            [file, path] = uigetfile({'*.mat', 'MAT result (*.mat)'}, 'Load EELS-ZEST result');
            if isequal(file, 0), return; end
            app.LoadedResultFile = fullfile(path, file);
            app.Stage1CubeResult = [];
            app.Stage2CubeResult = [];
            app.Stage1SingleResult = [];
            app.Stage2SingleResult = [];
            app.CurrentResult = [];
            app.FitContext = struct();
            app.RLDiagnosticCache = struct('keys', {{}}, 'values', {{}}, 'capacity', 3);
            result_names = {'cuberes', 'stage1_cuberes', 'stage2_cuberes', 'result'};
            listed = whos('-file', app.LoadedResultFile);
            available = intersect(result_names, {listed.name}, 'stable');
            if isempty(available)
                app.logStatus('No supported EELS-ZEST result variable was found in the selected file.');
                return;
            end
            S = load(app.LoadedResultFile, available{:});
            if isfield(S, 'cuberes') || isfield(S, 'stage1_cuberes') || isfield(S, 'stage2_cuberes')
                [cuberes, forced_stage2] = app.extractLoadedCubeResult(S);
                has_stage2 = app.cuberesHasStage2(cuberes);
                if ~isempty(forced_stage2)
                    has_stage2 = forced_stage2;
                end
                if has_stage2
                    app.Stage2CubeResult = cuberes;
                else
                    app.Stage1CubeResult = cuberes;
                end
                app.RLSourceCubeResult = cuberes;
                app.CurrentResult = app.firstValidResult(cuberes.results);
                app.RLSourceResult = app.firstRLDeconvResult(cuberes);
                app.plotMap(has_stage2);
                if ~isempty(app.CurrentResult)
                    app.plotSpectrumResult(app.CurrentResult, has_stage2);
                    app.plotDiagnosticSelection();
                end
                app.updateRLDeconvControl();
                app.logStatus(sprintf('Loaded datacube result from %s. First valid spectrum is selected for diagnostics.', file));
            elseif isfield(S, 'result')
                has_stage2 = app.resultHasStage2(S.result);
                if has_stage2
                    app.Stage2SingleResult = S.result;
                else
                    app.Stage1SingleResult = S.result;
                end
                app.CurrentResult = S.result;
                if app.resultCanRunRLDeconv(S.result)
                    app.RLSourceResult = S.result;
                end
                app.RLSourceCubeResult = [];
                app.RLDeconvResultCube = [];
                app.plotSpectrumResult(S.result, has_stage2);
                app.plotDiagnosticSelection();
                app.updateRLDeconvControl();
                app.logStatus(sprintf('Loaded single-spectrum result from %s.', file));
            else
                app.logStatus('No variable named cuberes, stage1_cuberes, stage2_cuberes, or result was found in the selected file.');
                app.updateRLDeconvControl();
            end
        end

        function runSingle(app, stage2)
            app.requireData();
            ctrl = app.stageControls(stage2);
            yy = round(ctrl.py.Value);
            xx = round(ctrl.px.Value);
            [ny, nx, ~] = size(app.RawCube);
            yy = min(max(yy, 1), ny);
            xx = min(max(xx, 1), nx);
            spec = squeeze(app.RawCube(yy, xx, :));
            settings = app.makeSettings(stage2);
            d = uiprogressdlg(app.UIFigure, 'Title', 'Fitting pixel', 'Message', sprintf('pixel (%d, %d)', yy, xx), 'Indeterminate', 'on');
            cleanup = onCleanup(@()delete(d));
            if stage2 && ~isempty(app.Stage1SingleResult) && isfield(app.Stage1SingleResult, 'pixel_y') && ...
                    app.Stage1SingleResult.pixel_y == yy && app.Stage1SingleResult.pixel_x == xx
                res = EELS_PearsonSparse_FullFit.promote_stage1_to_stage2(app.Stage1SingleResult, settings);
            else
                res = EELS_PearsonSparse_FullFit.solve_stage1(app.Ene, spec, settings);
                if stage2 && (~isfield(res, 'stage2') || isempty(res.stage2))
                    res = EELS_PearsonSparse_FullFit.promote_stage1_to_stage2(res, settings);
                end
            end
            res.pixel_y = yy;
            res.pixel_x = xx;
            app.CurrentResult = res;
            if stage2
                app.Stage2SingleResult = res;
            else
                app.Stage1SingleResult = res;
            end
            app.plotSpectrumResult(res, stage2);
            app.updateRLDeconvControl();
            app.logStage(stage2, app.formatResultLine(res, sprintf('pixel (%d,%d)', yy, xx)));
        end

        function runDatacube(app, stage2, run_mode)
            if nargin < 3 || isempty(run_mode)
                run_mode = 'serial';
            end
            app.requireDataFile();
            settings = app.makeSettings(stage2, run_mode);
            app.StopRequested = false;
            d = uiprogressdlg(app.UIFigure, 'Title', 'EELS-ZEST datacube run', 'Message', 'Starting...', ...
                'Cancelable', 'on', 'Value', 0);
            settings.datacube_progress_callback = @(done, total, elapsed, eta, med, label)app.updateProgress(d, done, total, elapsed, eta, med, label);
            settings.datacube_stop_callback = @()app.datacubeStopRequested(d);
            cleanup = onCleanup(@()app.closeProgressDialog(d));
            try
                if stage2 && ~isempty(app.Stage1CubeResult) && app.isReusableStage1Cube(app.Stage1CubeResult)
                    cuberes = EELS_PearsonSparse_FullFit.promote_datacube_stage1_to_stage2(app.Stage1CubeResult, settings);
                else
                    cuberes = EELS_PearsonSparse_FullFit.run_test_datacube(app.DataFile, settings);
                end
            catch ME
                if strcmp(ME.identifier, 'EELS:DatacubePaused')
                    app.logStatus(ME.message);
                    return;
                end
                rethrow(ME);
            end
            if stage2
                app.Stage2CubeResult = cuberes;
            else
                app.Stage1CubeResult = cuberes;
            end
            app.LoadedResultFile = app.findSavedResultFile(settings.test_output_dir, stage2);
            if isfield(cuberes, 'results') && ~isempty(cuberes.results)
                app.CurrentResult = app.firstValidResult(cuberes.results);
                if ~isempty(app.CurrentResult)
                    app.plotSpectrumResult(app.CurrentResult, stage2);
                end
            end
            app.plotMap(stage2);
            app.updateRLDeconvControl();
            app.logStage(stage2, app.formatCubeLine(cuberes));
            app.closeProgressDialog(d);
            clear cleanup
        end

        function requestPause(app)
            app.StopRequested = true;
            app.logStatus('Pause requested. Fitting will stop at the next checkpoint boundary and save completed pixels.');
        end

        function openReferenceTuner(app)
            app.TabGroup.SelectedTab=app.ReferenceTab;
        end

        function applyTunedSettings(app,settings)
            app.Settings=EELS_PearsonSparse_FullFit.defaults(settings);
            app.ReferenceStrategyActive=true;
            app.syncReferenceControls();
            app.clearFitCaches(true);
            app.populateAdvancedTable();
            app.logStatus('Reference strategy applied to Stage1/2. Previous fit caches cleared; checkpoint resume disabled.');
        end

        function resetReferenceSettings(app)
            app.Settings=EELS_PearsonSparse_FullFit.defaults(struct());
            app.Settings.datacube_sample_strategy='all';app.Settings.datacube_sample_pixels=Inf;
            app.ReferenceStrategyActive=false;app.syncReferenceControls();app.clearFitCaches(true);app.populateAdvancedTable();
            app.logStatus('Stage1/2 fitting parameters restored to defaults.');
        end

        function syncReferenceControls(app)
            for stage2=[false true]
                c=app.stageControls(stage2);c.cmin.Value=min(app.Settings.crop_range);c.cmax.Value=max(app.Settings.crop_range);
                c.resume.Value=false;
                if stage2
                    c.algorithm.Value=app.Settings.stage2_temperature_algorithm;
                    c.spatial.Value=app.Settings.stage1_hierarchical_enabled;
                    if isfield(c,'classPrune')&&~isempty(c.classPrune)&&isfield(app.Settings,'stage2_attention_class_prune_enabled')
                        c.classPrune.Value=app.Settings.stage2_attention_class_prune_enabled;
                    end
                end
            end
        end

        function s=referenceBaseSettings(app)
            app.applyAdvancedTable();s=app.Settings;
            s.crop_range=sort([app.Stage1Controls.cmin.Value,app.Stage1Controls.cmax.Value]);
        end

        function r=referenceInput(app)
            if ~isempty(app.RawCube)
                y=app.DataSelectedY;x=app.DataSelectedX;
                if isempty(y),y=round(app.Stage1Controls.py.Value);end
                if isempty(x),x=round(app.Stage1Controls.px.Value);end
                r=struct('ene',app.Ene,'spectrum',squeeze(app.RawCube(y,x,:)), ...
                    'source',sprintf('%s | pixel Y=%d X=%d',app.DataFile,y,x),'group',app.DataFile);
                r.cube=app.RawCube;
            elseif ~isempty(app.CurrentResult)
                r=struct('ene',app.CurrentResult.ene_raw,'spectrum',app.CurrentResult.raw_data, ...
                    'source','Current fitted spectrum','group',app.LoadedResultFile);
            else,error('EELS:Reference','Import and select a spectrum in Data, or import a reference in the tuning tab.');end
        end

        function settings = makeSettings(app, stage2, run_mode)
            if nargin < 3 || isempty(run_mode)
                run_mode = 'serial';
            end
            app.applyAdvancedTable();
            ctrl = app.stageControls(stage2);
            settings = EELS_PearsonSparse_FullFit.defaults(app.Settings);
            settings.run_stage2 = logical(stage2);
            settings.crop_range = sort([ctrl.cmin.Value, ctrl.cmax.Value]);
            settings.spatial_bin = [max(1, round(ctrl.by.Value)), max(1, round(ctrl.bx.Value))];
            settings.energy_bin = max(1, round(ctrl.eb.Value));
            settings.datacube_sample_strategy = 'all';
            settings.datacube_sample_pixels = Inf;
            settings.test_output_dir = char(ctrl.output.Value);
            % Full-cube GUI runs save numerical results first. Large diagnostic
            % figures are exported on demand to avoid MATLABWindow renderer leaks.
            settings.save_plots = false;
            settings.save_diagnostic_plots = false;
            settings.max_saved_spectrum_plots = 0;
            settings.force_hidden_figures = true;
            settings.close_parallel_pool_after_run = false;
            settings.datacube_progress_every = 1;
            settings.datacube_checkpoint_enabled = logical(ctrl.checkpoint.Value);
            settings.datacube_resume_checkpoint = logical(ctrl.resume.Value);
            settings.datacube_checkpoint_every = 8;
            settings.datacube_parallel_workers = max(1, round(ctrl.workers.Value));
            settings.datacube_parallel_min_pixels = 8;
            settings.use_bg = false;
            switch lower(string(run_mode))
                case "parallel"
                    settings.datacube_parallel_enabled = true;
                    settings.close_parallel_pool_after_run = false;
                otherwise
                    settings.datacube_parallel_enabled = false;
                    % Checkpointing follows the user checkbox in serial mode too.
            end
            if stage2
                settings.stage2_temperature_algorithm = char(ctrl.algorithm.Value);
                settings.stage1_hierarchical_enabled = logical(ctrl.spatial.Value);
                settings.stage1_spatial_zlp_enabled = logical(ctrl.spatial.Value);
                settings.stage2_spatial_map_enabled = logical(ctrl.spatial.Value);
                settings.spatial_refit_enabled = false;
                if isfield(ctrl, 'classPrune') && ~isempty(ctrl.classPrune)
                    settings.stage2_attention_class_prune_enabled = logical(ctrl.classPrune.Value);
                end
                if isfield(ctrl, 'forceRawStage2') && ~isempty(ctrl.forceRawStage2) && ...
                        isvalid(ctrl.forceRawStage2) && logical(ctrl.forceRawStage2.Value)
                    settings.stage2_shape_guard_enabled = false;
                    settings.stage2_force_raw_display = true;
                end
            end
            if ~isfolder(settings.test_output_dir), mkdir(settings.test_output_dir); end
            plot_dir = fullfile(settings.test_output_dir, 'plots');
            if ~isfolder(plot_dir), mkdir(plot_dir); end
        end

        function plotEnergyMap(app)
            app.requireData();
            win = sort([app.DataEnergyMinEdit.Value, app.DataEnergyMaxEdit.Value]);
            mask = app.Ene >= win(1) & app.Ene <= win(2);
            if nnz(mask) < 1, mask = true(size(app.Ene)); end
            img = sum(app.RawCube(:, :, mask), 3, 'omitnan');
            if ~isempty(app.DataMapImage) && isgraphics(app.DataMapImage)
                app.DataMapImage.CData = img;
                app.DataMapImage.XData = [1, size(img, 2)];
                app.DataMapImage.YData = [1, size(img, 1)];
                h = app.DataMapImage;
            else
                cla(app.DataMapAxes);
                h = imagesc(app.DataMapAxes, img);
                app.DataMapImage = h;
            end
            axis(app.DataMapAxes, 'image');
            app.DataMapAxes.XLimMode = 'auto';
            app.DataMapAxes.YLimMode = 'auto';
            app.DataMapAxes.CLimMode = 'auto';
            if isempty(app.DataMapColorbar) || ~isgraphics(app.DataMapColorbar)
                app.DataMapColorbar = colorbar(app.DataMapAxes);
            end
            title(app.DataMapAxes, sprintf('Energy-filtered map [%.1f, %.1f] meV', win(1), win(2)));
            xlabel(app.DataMapAxes, 'x'); ylabel(app.DataMapAxes, 'y');
            h.HitTest = 'on';
            h.ButtonDownFcn = @(~, evt)app.pickRawPixelFromDataMap(evt);
            app.DataMapAxes.ButtonDownFcn = @(~, evt)app.pickRawPixelFromDataMap(evt);
        end

        function pickRawPixelFromDataMap(app, evt)
            if isempty(app.RawCube), return; end
            if nargin >= 2 && ~isempty(evt) && isprop(evt, 'IntersectionPoint')
                cp = evt.IntersectionPoint;
            elseif nargin >= 2 && isstruct(evt) && isfield(evt, 'IntersectionPoint')
                cp = evt.IntersectionPoint;
            else
                cp = app.DataMapAxes.CurrentPoint(1, :);
            end
            xx = round(cp(1));
            yy = round(cp(2));
            [ny, nx, ~] = size(app.RawCube);
            if xx < 1 || xx > nx || yy < 1 || yy > ny, return; end
            spec = squeeze(app.RawCube(yy, xx, :));
            plot(app.DataSpectrumAxes, app.Ene, spec, 'k-', 'LineWidth', 1.2);
            grid(app.DataSpectrumAxes, 'on');
            xlabel(app.DataSpectrumAxes, 'Energy loss (meV)');
            ylabel(app.DataSpectrumAxes, 'Intensity');
            title(app.DataSpectrumAxes, sprintf('Raw spectrum (%d,%d)', yy, xx));
            app.DataSelectedY = yy;
            app.DataSelectedX = xx;
            app.markDataPixel(xx, yy);
            app.updateRLDeconvControl();
        end

        function markDataPixel(app, xx, yy)
            app.deleteGraphicsHandles('DataCursorHandles');
            hold(app.DataMapAxes, 'on');
            app.DataCursorHandles = { ...
                plot(app.DataMapAxes, xx, yy, 'wo', 'MarkerSize', 8, 'LineWidth', 1.5), ...
                plot(app.DataMapAxes, xx, yy, '+', 'Color', [0.85 0.08 0.10], 'MarkerSize', 8, 'LineWidth', 1.3)};
            for ii = 1:numel(app.DataCursorHandles)
                try
                    app.DataCursorHandles{ii}.HitTest = 'off';
                    app.DataCursorHandles{ii}.PickableParts = 'none';
                catch
                end
            end
            hold(app.DataMapAxes, 'off');
        end

        function plotLineProfile(app)
            app.requireData();
            win = sort([app.DataEnergyMinEdit.Value, app.DataEnergyMaxEdit.Value]);
            mask = app.Ene >= win(1) & app.Ene <= win(2);
            if nnz(mask) < 1, mask = true(size(app.Ene)); end
            prof = squeeze(sum(app.RawCube(:, :, mask), 1, 'omitnan'));
            fig = figure('Name', 'EELS-ZEST line profile', 'Color', 'w', 'Position', [160 140 980 620]);
            app.registerOwnedFigure(fig);
            ax = axes(fig);
            imagesc(ax, 1:size(prof, 1), app.Ene(mask), prof.');
            set(ax, 'YDir', 'normal');
            colormap(ax, 'turbo');
            colorbar(ax);
            xlabel(ax, 'x position');
            ylabel(ax, 'Energy loss (meV)');
            title(ax, sprintf('Line profile, summed over y [%.1f, %.1f] meV', win(1), win(2)));
            pbaspect(ax, [1.7 1 1]);
            grid(ax, 'on');
        end

        function plotIntegratedSpectrum(app)
            app.requireData();
            spec = squeeze(sum(sum(app.RawCube, 1, 'omitnan'), 2, 'omitnan'));
            plot(app.DataSpectrumAxes, app.Ene, spec, 'Color', [0.02 0.22 0.40], 'LineWidth', 1.4);
            app.DataSpectrumAxes.XLimMode = 'auto';
            app.DataSpectrumAxes.YLimMode = 'auto';
            grid(app.DataSpectrumAxes, 'on');
            xlabel(app.DataSpectrumAxes, 'Energy loss (meV)'); ylabel(app.DataSpectrumAxes, 'Integrated intensity');
            title(app.DataSpectrumAxes, 'Full-cube integrated spectrum');
        end

        function openRLDeconvFromData(app)
            use_data_context = app.hasDataRLKernelContext();
            cuberes = app.activeRLDeconvCubeResult();
            if ~use_data_context && isempty(cuberes) && app.hasRestorableRLSource()
                app.restoreRLSourceFromContext();
                cuberes = app.activeRLDeconvCubeResult();
            end
            if ~use_data_context && isempty(cuberes)
                app.updateRLDeconvControl();
                app.logStatus('RL Conv needs transferred fitted-ZLP kernels or a Stage 1/2 datacube result.');
                return;
            end
            settings = EELS_PearsonSparse_FullFit.defaults(app.Settings);
            settings.force_hidden_figures = true;
            if use_data_context
                map_size = size(app.RawCube, [1 2]);
                results = [];
                nres = prod(map_size);
            else
                map_size = app.cuberesMapSize(cuberes);
                if isempty(map_size)
                    map_size = [1, numel(cuberes.results)];
                end
                results = cuberes.results;
                nres = numel(results);
            end
            d = uiprogressdlg(app.UIFigure, 'Title', 'RL Conv datacube', ...
                'Message', sprintf('Running RL on %d fitted spectra...', nres), ...
                'Indeterminate', 'off', 'Value', 0, 'Cancelable', 'on');
            cleanup = onCleanup(@()app.closeProgressDialog(d));
            try
                ene_out = [];
                decube = [];
                status_map = zeros(map_size);
                r2_map = NaN(map_size);
                nrmse_map = NaN(map_size);
                fail_reasons = {};
                ok_count = 0;
                fail_count = 0;
                for ii = 1:nres
                    if d.CancelRequested
                        app.logStatus('RL Conv cancelled; Data cube was not replaced.');
                        return;
                    end
                    if use_data_context
                        [yy, xx] = ind2sub(map_size, ii);
                        res = app.currentDataRLResult(yy, xx);
                    else
                        res = app.resultAtIndex(results, ii);
                        [yy, xx] = app.resultPixelPosition(cuberes, ii, res, map_size);
                    end
                    if isempty(res) || yy < 1 || xx < 1 || yy > map_size(1) || xx > map_size(2)
                        fail_count = fail_count + 1;
                        continue;
                    end
                    if ~app.resultCanRunRLDeconv(res)
                        fail_count = fail_count + 1;
                        fail_reasons{end + 1} = 'missing_rl_fields'; %#ok<AGROW>
                        continue;
                    end
                    rl = EELS_PearsonSparse_FullFit.rl_deconvolve_background_subtracted(res, settings);
                    if ~rl.available
                        fail_count = fail_count + 1;
                        fail_reasons{end + 1} = rl.reason; %#ok<AGROW>
                        status_map(yy, xx) = -1;
                        continue;
                    end
                    if isempty(ene_out)
                        ene_out = rl.ene(:);
                        decube = NaN(map_size(1), map_size(2), numel(ene_out));
                    end
                    sig = rl.latent_signal(:);
                    if numel(sig) ~= numel(ene_out) || numel(rl.ene) ~= numel(ene_out) || ...
                            any(abs(rl.ene(:) - ene_out(:)) > 1e-9)
                        sig = interp1(rl.ene(:), sig(:), ene_out(:), 'linear', NaN);
                    end
                    scale = 1;
                    if isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                        scale = res.norm_val;
                    end
                    decube(yy, xx, :) = max(sig(:), 0) .* scale; %#ok<AGROW>
                    status_map(yy, xx) = 1;
                    r2_map(yy, xx) = rl.R2_weighted;
                    nrmse_map(yy, xx) = rl.NRMSE_weighted;
                    ok_count = ok_count + 1;
                    if mod(ii, max(1, floor(nres ./ 50))) == 0 || ii == nres
                        d.Value = ii ./ max(nres, 1);
                        d.Message = sprintf('RL Conv: %d/%d | ok %d | failed %d', ii, nres, ok_count, fail_count);
                        drawnow limitrate;
                    end
                end
                if ok_count < 1 || isempty(decube)
                    app.logStatus('RL Conv did not produce any available spectra; Data cube was not replaced.');
                    return;
                end
                app.RawCube = decube;
                app.Ene = ene_out(:);
                app.DataInfo.source = 'EELS-ZEST RL deconvolved latent signal cube';
                app.DataInfo.rl_conv_status_map = status_map;
                app.DataInfo.rl_conv_R2_map = r2_map;
                app.DataInfo.rl_conv_NRMSE_map = nrmse_map;
                app.DataInfo.rl_conv_ok_count = ok_count;
                app.DataInfo.rl_conv_fail_count = fail_count;
                app.DataInfo.rl_conv_fail_reasons = unique(fail_reasons);
                app.DataInfo.rl_conv_note = ...
                    'Each available fitted pixel was replaced by RL latent_signal * norm_val.';
                app.DataInfo.rl_conv_input_source = 'current Data with transferred fitted-ZLP kernel';
                app.DataInfo.rl_conv_applied = true;
                app.DataInfo.eels_zest_fit_context = app.FitContext;
                app.DataInfo = EELS_ZEST_DataOps.append_history(app.DataInfo, ...
                    'rl_deconvolution', struct('ok_count', ok_count, 'fail_count', fail_count));
                if ~use_data_context
                    app.RLSourceCubeResult = cuberes;
                    app.RLSourceResult = app.firstRLDeconvResult(cuberes);
                end
                app.RLDeconvResultCube = [];
                app.RLDiagnosticCache = struct('keys', {{}}, 'values', {{}}, 'capacity', 3);
                app.DataFile = '';
                app.persistCurrentDataForSolver('rl_deconvolved');
                app.clearFitCaches(true);
                app.updatePixelLimits();
                app.plotEnergyMap();
                app.plotIntegratedSpectrum();
                app.TabGroup.SelectedTab = app.DataTab;
                app.logStatus(sprintf('RL Conv complete: replaced Data with RL latent cube | ok %d | failed %d | output %s', ...
                    ok_count, fail_count, mat2str(size(app.RawCube))));
            catch ME
                app.logStatus(['RL Conv failed: ', ME.message]);
            end
            app.closeProgressDialog(d);
            clear cleanup
            app.updateRLDeconvControl();
        end

        function openRLDiagnosticForDataPixel(app)
            if isempty(app.DataSelectedY) || isempty(app.DataSelectedX)
                app.logStatus('Click a Data-map pixel first, then press RL diag.');
                return;
            end
            yy = app.DataSelectedY;
            xx = app.DataSelectedX;
            cached_rl = app.getRLDiagnosticCache(yy, xx);
            if app.hasDataRLKernelContext()
                res = app.currentDataRLResult(yy, xx);
            else
                [res, cached_rl] = app.rlSourceResultForPixel(yy, xx);
            end
            if isempty(res) && ~app.hasDataRLKernelContext()
                app.restoreRLSourceFromContext();
                [res, cached_rl] = app.rlSourceResultForPixel(yy, xx);
            end
            if isempty(res)
                app.logStatus(sprintf('No RL-capable fitted result is available for Data pixel (%d,%d).', yy, xx));
                return;
            end
            settings = EELS_PearsonSparse_FullFit.defaults(app.Settings);
            settings.force_hidden_figures = true;
            try
                if isempty(cached_rl) || ~isstruct(cached_rl) || ...
                        ~isfield(cached_rl, 'available') || ~cached_rl.available
                    cached_rl = EELS_PearsonSparse_FullFit.rl_deconvolve_background_subtracted(res, settings);
                    if cached_rl.available
                        app.putRLDiagnosticCache(yy, xx, cached_rl);
                    end
                end
                if isstruct(cached_rl) && isfield(cached_rl, 'available') && cached_rl.available
                    app.plotRLDeconvOnDataAxes(cached_rl);
                end
                fig = EELS_PearsonSparse_FullFit.plot_rl_deconvolution_diagnostic( ...
                    res, sprintf('EELS-ZEST RL diagnostic pixel (%d,%d)', yy, xx), settings);
                if ~isempty(fig) && isgraphics(fig)
                    fig.Visible = 'on';
                    app.registerOwnedFigure(fig);
                end
                app.logStatus(sprintf('Opened RL diagnostic for pixel (%d,%d).', yy, xx));
            catch ME
                app.logStatus(['RL diagnostic failed: ', ME.message]);
            end
            app.updateRLDeconvControl();
        end

        function plotRLDeconvOnDataAxes(app, rl)
            ax = app.DataSpectrumAxes;
            cla(ax);
            ene = rl.ene(:);
            plot(ax, ene, rl.target_signal(:), 'k-', 'LineWidth', 0.95); hold(ax, 'on');
            plot(ax, ene, rl.reconvolved_signal(:), 'Color', [0.55 0.10 0.78], 'LineWidth', 1.35);
            plot(ax, ene, rl.latent_signal(:), 'Color', [0.02 0.50 0.72], 'LineWidth', 0.95);
            yline(ax, 0, ':', 'Color', [0.5 0.5 0.5]);
            hold(ax, 'off');
            grid(ax, 'on');
            xlabel(ax, 'Energy loss (meV)');
            ylabel(ax, 'normalized intensity');
            title(ax, sprintf('RL deconv | R2 %.4f | latent BE T %.1f K', ...
                rl.R2_weighted, rl.latent_ratio_T_K));
            legend(ax, {'raw - ZLP - bg', 'RL reconvolved', 'RL latent'}, 'Location', 'best');
        end

        function integrateROI(app)
            app.requireData();
            app.plotEnergyMap();
            app.logStatus('ROI selection: double-click to confirm. Ctrl+Q cancels, Ctrl+H confirms and plots.');
            roi = drawrectangle(app.DataMapAxes, 'Color', [0.0 0.65 1.0], 'LineWidth', 1.2);
            [pos, cancelled] = app.waitForROIWithHotkeys(roi, ...
                'ROI selection: double-click to confirm. Ctrl+Q cancels, Ctrl+H confirms and plots.');
            if cancelled || isempty(pos)
                app.logStatus('ROI selection cancelled.');
                return;
            end
            [ny, nx, ~] = size(app.RawCube);
            x1 = max(1, floor(pos(1)));
            y1 = max(1, floor(pos(2)));
            x2 = min(nx, ceil(pos(1) + pos(3)));
            y2 = min(ny, ceil(pos(2) + pos(4)));
            roi_cube = app.RawCube(y1:y2, x1:x2, :);
            roi_spectrum = squeeze(sum(sum(roi_cube, 1, 'omitnan'), 2, 'omitnan'));
            plot(app.DataSpectrumAxes, app.Ene, roi_spectrum, 'Color', [0.0 0.35 0.75], 'LineWidth', 1.5);
            grid(app.DataSpectrumAxes, 'on');
            title(app.DataSpectrumAxes, sprintf('ROI integrated spectrum: x %d-%d, y %d-%d', x1, x2, y1, y2));
            roi_export = struct('cube', roi_cube, 'spectrum', roi_spectrum, 'ene', app.Ene, ...
                'x_range', [x1 x2], 'y_range', [y1 y2], 'source_file', app.DataFile);
            app.exportPayloadToFile(roi_export, 'Export ROI', 'eels_zest_roi', ...
                {'MAT', 'CSV'}, 'eels_zest_roi');
        end

        function openStackSpectraTool(app)
            app.requireData();
            win = sort([app.DataEnergyMinEdit.Value, app.DataEnergyMaxEdit.Value]);
            fig = EELS_ZEST_StackTool(app.RawCube, app.Ene, app.DataInfo, win, ...
                @(stack)app.receiveStackSpectra(stack));
            app.registerOwnedFigure(fig);
            app.logStatus('Opened the nonblocking Stack Spectra selection window.');
        end

        function receiveStackSpectra(app, stack)
            app.StackSpectra = stack;
            app.logStatus(sprintf('Stack Spectra retained %d selection(s).', numel(stack.spectrum)));
        end

        function [spec, selection] = spectrumFromPointPosition(app, pos)
            [ny, nx, ~] = size(app.RawCube);
            xx = min(nx, max(1, round(pos(1))));
            yy = min(ny, max(1, round(pos(2))));
            spec = squeeze(app.RawCube(yy, xx, :));
            selection = struct('type', 'cursor', 'x', xx, 'y', yy, ...
                'label', sprintf('pixel (%d,%d)', yy, xx));
        end

        function [spec, selection] = spectrumFromROIPosition(app, pos)
            [ny, nx, ~] = size(app.RawCube);
            x1 = max(1, floor(pos(1)));
            y1 = max(1, floor(pos(2)));
            x2 = min(nx, ceil(pos(1) + pos(3)));
            y2 = min(ny, ceil(pos(2) + pos(4)));
            roi_cube = app.RawCube(y1:y2, x1:x2, :);
            spec = squeeze(sum(sum(roi_cube, 1, 'omitnan'), 2, 'omitnan'));
            selection = struct('type', 'roi', 'x_range', [x1 x2], 'y_range', [y1 y2], ...
                'position', [x1 y1 max(1, x2 - x1 + 1) max(1, y2 - y1 + 1)], ...
                'label', sprintf('ROI x%d-%d y%d-%d', x1, x2, y1, y2));
        end

        function plotStackedSpectra(app, stack, fig_name)
            fig = figure('Name', fig_name, 'Color', 'w', 'Position', [120 120 980 640]);
            app.registerOwnedFigure(fig);
            ax = axes(fig); hold(ax, 'on');
            for kk = 1:numel(stack.spectrum)
                y = app.normalizeStackSpectrum(stack.spectrum{kk});
                plot(ax, stack.ene, y, 'LineWidth', 1.35);
            end
            yline(ax, 0, '-', 'Color', [0.72 0.76 0.80], 'HandleVisibility', 'off');
            grid(ax, 'on');
            xlabel(ax, 'Energy loss (meV)');
            ylabel(ax, 'Baseline-aligned, total-intensity normalized');
            title(ax, fig_name, 'Interpreter', 'none');
            legend(ax, stack.name, 'Interpreter', 'none', 'Location', 'best');
            hold(ax, 'off');
        end

        function plotMap(app, stage2)
            ax = app.mapAxes(stage2);
            cuberes = app.cubeResult(stage2);
            cla(ax);
            if isempty(cuberes)
                title(ax, 'No datacube result yet');
                return;
            end
            ctrl = app.stageControls(stage2);
            if isfield(ctrl, 'map') && ~isempty(ctrl.map) && isvalid(ctrl.map)
                map_name = ctrl.map.Value;
            elseif stage2
                map_name = 'T data valid';
            else
                map_name = 'R2';
            end
            [map, ttl] = app.mapData(cuberes, map_name);
            if isempty(map)
                title(ax, 'No map field available');
                return;
            end
            h = imagesc(ax, map);
            axis(ax, 'image');
            colorbar(ax);
            title(ax, ttl);
            xlabel(ax, 'x'); ylabel(ax, 'y');
            h.HitTest = 'on';
            h.ButtonDownFcn = @(~, ~)app.pickPixelFromMap(stage2);
            ax.ButtonDownFcn = @(~, ~)app.pickPixelFromMap(stage2);
        end

        function items = mapItems(~)
            items = {'R2', 'NRMSE', 'T raw', 'T data valid', 'T joint', 'T joint uncertainty', ...
                'Spatial inferred', 'T repaired', 'T spatial refined', 'T uncertainty', ...
                'ZLP cluster', 'Cluster membership', 'Cluster outlier', 'Template deviation', ...
                'ZLP spatial loss', 'ZLP interface score', 'ZLP sharp interface', ...
                'Reliability', 'Attention fraction', 'Warning', ...
                'Bad pixel', 'Interface', 'Neighbor loss', 'Spectrum neighbor loss', ...
                'Temperature neighbor z', 'Attention peak delta', 'Signal peak delta', ...
                'ZLP guard start', 'ZLP guard fraction', 'ZLP core width'};
        end

        function [img, ttl] = mapData(~, cuberes, name)
            img = [];
            ttl = char(name);
            if isempty(cuberes) || ~isstruct(cuberes), return; end
            switch char(name)
                case 'R2', field = 'R2_map'; ttl = 'R2 map';
                case 'NRMSE', field = 'NRMSE_map'; ttl = 'NRMSE map';
                case 'T raw', field = 'T_map_raw'; ttl = 'T raw (K)';
                case 'T data valid', field = 'T_map_data_valid'; ttl = 'Spectrum-constrained T (K)';
                case 'T joint', field = 'T_map_joint'; ttl = 'ZLP-guided T, no temperature smoothing (K)';
                case 'T joint uncertainty', field = 'T_uncertainty_joint_map'; ttl = 'ZLP-guided T uncertainty (K)';
                case 'Spatial inferred', field = 'T_spatial_inferred_mask'; ttl = 'Spatially inferred temperature mask';
                case 'T repaired', field = 'T_map_repaired'; ttl = 'T repaired (K)';
                case 'T spatial refined', field = 'T_map_spatial_refined'; ttl = 'T spatial refined (K)';
                case 'T uncertainty', field = 'T_uncertainty_map'; ttl = 'T uncertainty (K)';
                case 'Reliability', field = 'reliability_score_map'; ttl = 'Reliability';
                case 'Attention fraction', field = 'attention_score_map'; ttl = 'Attention fraction';
                case 'Warning', field = 'warning_map'; ttl = 'Warning map';
                case 'Bad pixel', field = 'bad_pixel_mask'; ttl = 'Bad pixel mask';
                case 'Interface', field = 'interface_mask'; ttl = 'Interface mask';
                case 'Neighbor loss', field = 'neighbor_loss_map'; ttl = 'Neighbor loss';
                case 'Spectrum neighbor loss', field = 'spectrum_neighbor_loss_map'; ttl = 'Spectrum neighbor loss';
                case 'Temperature neighbor z', field = 'temperature_neighbor_z_map'; ttl = 'Temperature neighbor z';
                case 'Attention peak delta', field = 'attention_peak_delta_map'; ttl = 'Attention peak delta';
                case 'Signal peak delta', field = 'signal_peak_delta_map'; ttl = 'Signal peak delta';
                case 'ZLP guard start', field = 'zlp_guard_start_map'; ttl = 'ZLP guard start (meV)';
                case 'ZLP guard fraction', field = 'zlp_guard_fraction_map'; ttl = 'ZLP fraction in guard';
                case 'ZLP core width', field = 'zlp_core_width_map'; ttl = 'ZLP core width (meV)';
                case 'ZLP cluster', field = 'cluster_label_map'; ttl = 'Hierarchical ZLP cluster';
                case 'Cluster membership', field = 'cluster_membership_map'; ttl = 'ZLP cluster membership';
                case 'Cluster outlier', field = 'cluster_outlier_mask'; ttl = 'ZLP cluster outlier';
                case 'Template deviation', field = 'template_deviation_map'; ttl = 'ZLP template deviation';
                case 'ZLP spatial loss', field = 'zlp_spatial_loss_map'; ttl = 'ZLP spatial loss';
                case 'ZLP interface score', field = 'zlp_interface_score_map'; ttl = 'Sharp-interface ZLP score';
                case 'ZLP sharp interface', field = 'zlp_interface_mask'; ttl = 'ZLP sharp-interface mask';
                otherwise, field = 'R2_map'; ttl = 'R2 map';
            end
            if isfield(cuberes, field)
                img = cuberes.(field);
            elseif isfield(cuberes, 'R2_map')
                img = cuberes.R2_map;
                ttl = 'R2 map';
            end
        end

        function plotSpectrumResult(app, res, stage2)
            ax = app.spectrumAxes(stage2);
            cla(ax);
            if isempty(res) || ~isstruct(res)
                title(ax, 'No spectrum result yet');
                return;
            end
            ene = app.fieldOr(res, 'ene_corrected', app.fieldOr(res, 'ene_raw', app.Ene));
            raw = app.fieldOr(res, 'raw_data', []);
            fit = app.fieldOr(res, 'fitted_spectrum', []);
            zlp = app.fieldOr(res, 'zlp_only', []);
            vib = app.fieldOr(res, 'vib_only', []);
            fit_label = 'fit';
            if stage2 && isfield(res, 'stage2') && isstruct(res.stage2) && ...
                    isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                if isfield(res.stage2, 'model_norm') && ~isempty(res.stage2.model_norm)
                    fit = res.stage2.model_norm(:) .* res.norm_val;
                    fit_label = 'stage2 fit';
                end
                if isfield(res.stage2, 'zlp_norm') && ~isempty(res.stage2.zlp_norm)
                    zlp = res.stage2.zlp_norm(:) .* res.norm_val;
                end
                if isfield(res.stage2, 'vib_norm') && ~isempty(res.stage2.vib_norm)
                    vib = res.stage2.vib_norm(:) .* res.norm_val;
                end
            end
            if isempty(raw) || isempty(fit)
                title(ax, 'Result does not contain plottable spectrum fields');
                return;
            end
            h = gobjects(0);
            labels = {};
            h(end + 1) = plot(ax, ene, raw, '.', 'Color', [0.1 0.1 0.1], 'MarkerSize', 4); hold(ax, 'on');
            labels{end + 1} = 'raw';
            if stage2 && strcmp(fit_label, 'stage2 fit')
                fit_color = [0.65 0.05 0.90];
            else
                fit_color = [0 0 0];
            end
            h(end + 1) = plot(ax, ene, fit, '-', 'Color', fit_color, 'LineWidth', 1.4);
            labels{end + 1} = fit_label;
            if ~isempty(zlp)
                h(end + 1) = plot(ax, ene, zlp, '-', 'Color', [0.0 0.45 0.24], 'LineWidth', 1.2);
                labels{end + 1} = 'ZLP';
            end
            if ~isempty(vib)
                h(end + 1) = plot(ax, ene, vib, '-', 'Color', [0.95 0.28 0.10], 'LineWidth', 1.2);
                labels{end + 1} = 'sparse vib';
            end
            legend(ax, h, labels, 'Location', 'best');
            grid(ax, 'on');
            xlabel(ax, 'Energy loss (meV)'); ylabel(ax, 'Intensity');
            title(ax, app.resultTitle(res, stage2));
            hold(ax, 'off');
        end

        function openDiagnostic(app, stage2)
            res = app.activeResult(stage2);
            if isempty(res)
                app.logStatus('No fitted spectrum is available for a diagnostic plot.');
                return;
            end
            app.CurrentResult = res;
            try
                fig = EELS_PearsonSparse_FullFit.plot_result(res, app.resultTitle(res, stage2));
                if ~isempty(fig) && isgraphics(fig)
                    fig.Visible = 'on';
                end
                app.plotDiagnosticSelection();
            catch ME
                app.logStatus(['Diagnostic plotting failed: ', ME.message]);
            end
        end

        function res = activeResult(app, stage2)
            res = app.singleResult(stage2);
            if isempty(res) && ~isempty(app.CurrentResult)
                if ~stage2 || app.resultHasStage2(app.CurrentResult)
                    res = app.CurrentResult;
                end
            end
            if isempty(res)
                cuberes = app.cubeResult(stage2);
                if ~isempty(cuberes) && isfield(cuberes, 'results') && ~isempty(cuberes.results)
                    res = app.firstValidResult(cuberes.results);
                end
            end
        end

        function plotDiagnosticSelection(app)
            if isempty(app.CurrentResult)
                res = app.activeResult(true);
                if isempty(res), res = app.activeResult(false); end
                if isempty(res)
                    app.logStatus('No selected result for diagnostic plotting.');
                    return;
                end
                app.CurrentResult = res;
            end
            if isempty(app.DiagnosticDropDown) || ~isvalid(app.DiagnosticDropDown)
                app.openFullDiagnostic();
                return;
            end
            res = app.CurrentResult;
            if ~(isstruct(res) && isfield(res, 'raw_data'))
                app.logStatus('Selected result does not contain spectrum fields.');
                return;
            end
            cla(app.ViewAxes1); cla(app.ViewAxes2);
            ene = app.fieldOr(res, 'ene_corrected', app.fieldOr(res, 'ene_raw', app.Ene));
            ene = ene(:);
            raw = app.fieldOr(res, 'raw_data', []);
            fit = app.fieldOr(res, 'fitted_spectrum', []);
            zlp = app.fieldOr(res, 'zlp_only', []);
            vib = app.fieldOr(res, 'vib_only', []);
            item = app.DiagnosticDropDown.Value;
            switch char(item)
                case 'Full diagnostic figure'
                    app.openFullDiagnostic();
                case 'Raw vs total fit'
                    plot(app.ViewAxes1, ene, raw(:), 'k.', ene, fit(:), 'r-', 'LineWidth', 1.2);
                    legend(app.ViewAxes1, {'raw', 'fit'}, 'Location', 'best'); grid(app.ViewAxes1, 'on');
                    title(app.ViewAxes1, 'Raw spectrum and total fit');
                    plot(app.ViewAxes2, ene, raw(:) - fit(:), 'b-');
                    yline(app.ViewAxes2, 0, 'k-'); grid(app.ViewAxes2, 'on'); title(app.ViewAxes2, 'Residual');
                case 'ZLP components'
                    plot(app.ViewAxes1, ene, zlp(:), 'Color', [0.05 0.45 0.18], 'LineWidth', 1.4); hold(app.ViewAxes1, 'on');
                    if isfield(res, 'zlp_components') && ~isempty(res.zlp_components)
                        plot(app.ViewAxes1, ene, res.zlp_components, 'LineWidth', 0.8);
                    end
                    hold(app.ViewAxes1, 'off'); grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'ZLP envelope and components');
                    semilogy(app.ViewAxes2, ene, max(raw(:), eps), 'k.', ene, max(zlp(:), eps), 'g-');
                    grid(app.ViewAxes2, 'on'); title(app.ViewAxes2, 'Log ZLP check');
                case 'Sparse convolution / vib'
                    plot(app.ViewAxes1, ene, vib(:), 'Color', [0.9 0.25 0.05], 'LineWidth', 1.2);
                    grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Sparse convolution / vib signal');
                    app.plotPMSticks(app.ViewAxes2, res);
                case 'Residual and CI'
                    residual = raw(:) - fit(:);
                    area(app.ViewAxes1, ene, residual, 'FaceColor', [0.45 0.62 0.95], 'EdgeColor', [0.1 0.35 0.95]);
                    yline(app.ViewAxes1, 0, 'k-'); grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Residual');
                    if isfield(res, 'confidence_interval') && isfield(res.confidence_interval, 'half_width')
                        plot(app.ViewAxes2, ene, res.confidence_interval.half_width(:), 'b-');
                        grid(app.ViewAxes2, 'on'); title(app.ViewAxes2, 'Energy-dependent CI half-width');
                    else
                        title(app.ViewAxes2, 'No confidence interval stored');
                    end
                case 'Source-intensity sticks'
                    app.plotPMSticks(app.ViewAxes1, res);
                case 'Spectrum ln(Iloss/Igain)'
                    app.plotSpectrumLossGainRatio(app.ViewAxes1, res);
                    app.plotHardBERatio(app.ViewAxes2, res);
                case 'Hard-BE ln ratio'
                    app.plotHardBERatio(app.ViewAxes1, res);
                    app.plotPMSticks(app.ViewAxes2, res);
                case 'Attention curve'
                    if isfield(res, 'stage2') && isfield(res.stage2, 'attention') && ...
                            isfield(res.stage2.attention, 'feature_score')
                        plot(app.ViewAxes1, ene, res.stage2.attention.feature_score(:), 'b-');
                        ylim(app.ViewAxes1, [0 1]); grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Attention feature score');
                    else
                        title(app.ViewAxes1, 'No stage2 attention field stored');
                    end
                case 'Loss budget'
                    if isfield(res, 'loss_breakdown') && isfield(res.loss_breakdown, 'fraction')
                        vals = 100 .* res.loss_breakdown.fraction;
                        labels = strrep(res.loss_breakdown.names, '_', ' ');
                        barh(app.ViewAxes1, vals); set(app.ViewAxes1, 'YTick', 1:numel(labels), 'YTickLabel', labels);
                        grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Loss budget (%)');
                    else
                        title(app.ViewAxes1, 'No loss budget stored');
                    end
                case 'Log signal vs log ZLP'
                    plot(app.ViewAxes1, ene, log10(max(raw(:), eps)), 'k.', ...
                        ene, log10(max(zlp(:), eps)), 'g-', ...
                        ene, log10(max(vib(:), eps)), 'r-');
                    legend(app.ViewAxes1, {'log raw', 'log ZLP', 'log vib'}, 'Location', 'best');
                    grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Log signal vs log ZLP');
            end
            xlabel(app.ViewAxes1, 'Energy loss (meV)');
        end

        function plotPMSticks(~, ax, res)
            cla(ax); hold(ax, 'on');
            if isfield(res, 'EM_final') && isfield(res, 'PM_stokes')
                stem(ax, res.EM_final, res.PM_stokes, 'Color', [0.95 0.45 0.12], 'Marker', 'none');
                if isfield(res, 'PM_antistokes')
                    stem(ax, -res.EM_final, res.PM_antistokes, 'Color', [0.1 0.50 0.95], 'Marker', 'none');
                end
            end
            if isfield(res, 'stage2') && isfield(res.stage2, 'EM_hard') && ~isempty(res.stage2.EM_hard)
                yl = ylim(ax);
                plot(ax, res.stage2.EM_hard, 0 .* res.stage2.EM_hard + yl(2) .* 0.94, 'm.', 'MarkerSize', 8);
                plot(ax, -res.stage2.EM_hard, 0 .* res.stage2.EM_hard + yl(2) .* 0.88, 'c.', 'MarkerSize', 8);
            end
            hold(ax, 'off'); grid(ax, 'on');
            xlabel(ax, 'Mode energy (meV)'); ylabel(ax, 'Source intensity');
            title(ax, 'Sparse source intensity and hard-BE markers');
        end

        function plotHardBERatio(~, ax, res)
            cla(ax);
            befit = EELS_PearsonSparse_FullFit.hard_be_ln_ratio_fit(res);
            if befit.available
                scatter(ax, befit.EM_meV, befit.ln_loss_gain, 28, befit.weight, 'filled'); hold(ax, 'on');
                xmax = max(befit.EM_meV, [], 'omitnan');
                xx = linspace(0, max(xmax, eps), 120);
                plot(ax, xx, befit.slope_per_meV .* xx, 'k-', 'LineWidth', 1.35);
                hold(ax, 'off');
                grid(ax, 'on');
                xlabel(ax, 'Mode energy (meV)');
                ylabel(ax, 'ln(I_{source,loss}/I_{source,gain})');
                title(ax, sprintf('Hard-pair BE ratio: T %.1f K, R2 %.3f', befit.T_K, befit.R2_origin));
            else
                text(ax, 0.05, 0.55, 'No valid hard Stokes/anti-Stokes intensity pairs for BE ratio fit', ...
                    'Units', 'normalized');
                title(ax, 'Hard-pair BE ratio unavailable');
                axis(ax, 'off');
            end
        end

        function plotSpectrumLossGainRatio(~, ax, res)
            cla(ax);
            lgfit = EELS_PearsonSparse_FullFit.spectrum_loss_gain_ln_ratio_fit(res);
            if lgfit.available
                scatter(ax, lgfit.E_meV, lgfit.ln_loss_gain, 24, lgfit.weight, 'filled'); hold(ax, 'on');
                if lgfit.fit_available
                    xmax = max(lgfit.E_meV, [], 'omitnan');
                    xx = linspace(0, max(xmax, eps), 120);
                    plot(ax, xx, lgfit.slope_per_meV .* xx, 'k-', 'LineWidth', 1.25);
                    ttl = sprintf('Spectrum ln(Iloss/Igain): T %.1f K, R2 %.3f', lgfit.T_K, lgfit.R2_origin);
                else
                    ttl = sprintf('Spectrum ln(Iloss/Igain): n=%d, no positive BE slope', lgfit.num_points);
                end
                yline(ax, 0, ':', 'Color', [0.45 0.45 0.45]);
                hold(ax, 'off'); grid(ax, 'on');
                xlabel(ax, '|E| (meV)');
                ylabel(ax, 'ln(I_{loss}/I_{gain})');
                title(ax, ttl, 'Interpreter', 'none');
            else
                text(ax, 0.05, 0.55, 'No valid mirrored ZLP/bg-subtracted Iloss/Igain points', ...
                    'Units', 'normalized');
                title(ax, 'Spectrum ln(Iloss/Igain) unavailable');
                axis(ax, 'off');
            end
        end

        function openFullDiagnostic(app)
            if isempty(app.CurrentResult)
                res = app.activeResult(true);
                if isempty(res), res = app.activeResult(false); end
                if isempty(res)
                    app.logStatus('No selected result is available for a full diagnostic figure.');
                    return;
                end
                app.CurrentResult = res;
            end
            try
                fig = EELS_PearsonSparse_FullFit.plot_result(app.CurrentResult, 'EELS-ZEST selected spectrum');
                if ~isempty(fig) && isgraphics(fig)
                    fig.Visible = 'on';
                end
            catch ME
                app.logStatus(['Full diagnostic plot failed: ', ME.message]);
            end
        end

        function buildPMAnalysis(app)
            src = char(app.PMControls.source.Value);
            axis_meV = app.PMControls.emin.Value:app.PMControls.step.Value:app.PMControls.emax.Value;
            if isempty(axis_meV), axis_meV = -250:1:250; end
            [results, map_size] = app.pmSourceResults(src);
            if isempty(results), error('No sparse phonon source-intensity result is available. Run or load Stage 1 / Stage 2 first.'); end
            pm_cube = NaN([map_size, numel(axis_meV)]);
            valid_mask = false(map_size);
            for ii = 1:numel(results)
                res = results{ii};
                if isempty(res) || ~isstruct(res) || isfield(res, 'error'), continue; end
                [em, ps, pa] = app.extractPM(res, src);
                if isempty(em), continue; end
                signed_em = [-flipud(em(:)); em(:)];
                signed_pm = [flipud(pa(:)); ps(:)];
                [yy, xx] = ind2sub(map_size, ii);
                pm_cube(yy, xx, :) = interp1(signed_em, signed_pm, axis_meV, 'linear', 0);
                valid_mask(yy, xx) = true;
            end
            mean_pm = squeeze(mean(mean(pm_cube, 1, 'omitnan'), 2, 'omitnan'));
            app.LastPMResult = struct('axis_meV', axis_meV(:).', 'pm_cube', pm_cube, ...
                'mean_PM', mean_pm(:).', 'valid_mask', valid_mask, 'source', src);
            app.plotPMMean();
            app.plotPMIntegralMap();
            app.logStatus(sprintf('Source-intensity cube built from %s: %d valid spectra.', src, nnz(valid_mask)));
        end

        function plotPMMean(app)
            pmres = app.requirePM();
            y = pmres.mean_PM;
            if app.PMControls.smooth.Value && numel(y) > 7
                y = movmean(y, 7, 'omitnan');
            end
            plot(app.PMAxesSpectrum, pmres.axis_meV, y, 'k-', 'LineWidth', 1.5);
            grid(app.PMAxesSpectrum, 'on');
            xlabel(app.PMAxesSpectrum, 'Signed mode energy (meV)'); ylabel(app.PMAxesSpectrum, 'Mean source intensity');
            title(app.PMAxesSpectrum, sprintf('Mean source intensity, source = %s', pmres.source));
            xline(app.PMAxesSpectrum, 0, '-', 'Color', [0.45 0.45 0.45]);
        end

        function plotPMIntegralMap(app)
            pmres = app.requirePM();
            win = sort([app.PMControls.wmin.Value, app.PMControls.wmax.Value]);
            mask = pmres.axis_meV >= win(1) & pmres.axis_meV <= win(2);
            if nnz(mask) < 1, mask = true(size(pmres.axis_meV)); end
            map = sum(pmres.pm_cube(:, :, mask), 3, 'omitnan');
            imagesc(app.PMAxesMap, map);
            axis(app.PMAxesMap, 'image');
            colorbar(app.PMAxesMap);
            xlabel(app.PMAxesMap, 'x'); ylabel(app.PMAxesMap, 'y');
            title(app.PMAxesMap, sprintf('Source-intensity integral [%.1f, %.1f] meV', min(pmres.axis_meV(mask)), max(pmres.axis_meV(mask))));
        end

        function runPMNMF(app)
            pmres = app.requirePM();
            X = reshape(pmres.pm_cube, [], numel(pmres.axis_meV));
            keep = pmres.valid_mask(:) & all(isfinite(X), 2);
            X = max(X(keep, :), 0);
            k = round(app.PMControls.nk.Value);
            if isempty(X), error('No valid source-intensity rows for NMF.'); end
            if exist('nnmf', 'file') == 2
                [W, H] = nnmf(X, k, 'replicates', 2);
            else
                [U, S, V] = svds(double(X), min(k, min(size(X))));
                W = max(U * S, 0);
                H = max(V', 0);
            end
            maps = NaN([size(pmres.valid_mask), k]);
            for jj = 1:k
                temp = NaN(size(pmres.valid_mask));
                temp(keep) = W(:, jj);
                maps(:, :, jj) = temp;
            end
            imagesc(app.PMAxesExtra, maps(:, :, 1)); axis(app.PMAxesExtra, 'image'); colorbar(app.PMAxesExtra);
            title(app.PMAxesExtra, 'Source-intensity NMF component 1 map');
            fig = figure('Name', 'EELS-ZEST source-intensity NMF spectra');
            app.registerOwnedFigure(fig);
            plot(pmres.axis_meV, H', 'LineWidth', 1.2); grid on;
            xlabel('Signed mode energy (meV)'); ylabel('Component source intensity'); title('Source-intensity NMF component spectra');
            app.LastPMResult.nmf = struct('W', W, 'H', H, 'maps', maps, 'axis_meV', pmres.axis_meV);
            app.logStatus('Source-intensity NMF result is retained with the current intensity cube for file export.');
        end

        function runPMCluster(app)
            pmres = app.requirePM();
            X = reshape(pmres.pm_cube, [], numel(pmres.axis_meV));
            keep = pmres.valid_mask(:) & all(isfinite(X), 2);
            X = X(keep, :);
            k = round(app.PMControls.ck.Value);
            if isempty(X), error('No valid source-intensity rows for clustering.'); end
            if exist('kmeans', 'file') == 2
                labels = kmeans(X, k, 'Replicates', 3, 'MaxIter', 200);
            else
                [~, score] = pca(double(X), 'NumComponents', min(3, size(X, 2)));
                labels = discretize(score(:, 1), k);
            end
            map = NaN(size(pmres.valid_mask));
            map(keep) = labels;
            imagesc(app.PMAxesExtra, map); axis(app.PMAxesExtra, 'image'); colorbar(app.PMAxesExtra);
            title(app.PMAxesExtra, sprintf('Source-intensity clusters, k = %d', k));
            app.LastPMResult.cluster = struct('labels_map', map, 'axis_meV', pmres.axis_meV, 'cluster_count', k);
            app.logStatus('Source-intensity clustering result is retained with the current intensity cube for file export.');
        end

        function exportPMCube(app)
            pmres = app.requirePM();
            app.exportPayloadToFile(pmres, 'Export sparse phonon intensity cube', ...
                'eels_zest_sparse_phonon_intensity_cube', {'MAT', 'CSV'}, 'eels_zest_pm_cube');
        end

        function exportSelectedPM(app)
            pmres = app.requirePM();
            win = sort([app.PMControls.wmin.Value, app.PMControls.wmax.Value]);
            mask = pmres.axis_meV >= win(1) & pmres.axis_meV <= win(2);
            selected = struct('axis_meV', pmres.axis_meV(mask), 'pm_cube', pmres.pm_cube(:, :, mask), ...
                'mean_PM', pmres.mean_PM(mask), 'window_meV', win, 'source', pmres.source);
            app.exportPayloadToFile(selected, 'Export selected sparse phonon intensity window', ...
                'eels_zest_sparse_phonon_window', {'MAT', 'CSV'}, 'eels_zest_selected_pm');
        end

        function populateAdvancedTable(app)
            if isempty(app.AdvancedTable) || ~isvalid(app.AdvancedTable), return; end
            s = EELS_PearsonSparse_FullFit.defaults(app.Settings);
            names = fieldnames(s);
            rows = {};
            for ii = 1:numel(names)
                key = names{ii};
                val = s.(key);
                if isstruct(val) || isa(val, 'function_handle') || iscell(val)
                    continue;
                end
                if ismember(key, {'datacube_sample_pixels', 'datacube_sample_strategy'})
                    continue;
                end
                rows(end + 1, :) = {key, app.settingDisplayName(key), app.valueToText(val), app.settingAdvice(key)}; %#ok<AGROW>
            end
            app.AdvancedTable.Data = rows;
        end

        function applyAdvancedTable(app)
            if isempty(app.AdvancedTable) || ~isvalid(app.AdvancedTable) || isempty(app.AdvancedTable.Data)
                return;
            end
            s = app.Settings;
            if isempty(s), s = EELS_PearsonSparse_FullFit.defaults(struct()); end
            data = app.AdvancedTable.Data;
            for ii = 1:size(data, 1)
                key = data{ii, 1};
                txt = data{ii, 3};
                if isfield(s, key)
                    % Preserve unedited settings exactly across the GUI round trip.
                    if strcmp(char(string(txt)), app.valueToText(s.(key))), continue; end
                    s.(key) = app.textToValue(txt, s.(key));
                end
            end
            app.Settings = EELS_PearsonSparse_FullFit.defaults(s);
        end

        function resetAdvanced(app)
            app.resetReferenceSettings();
            app.ReferenceTuner.Settings=app.Settings;app.ReferenceTuner.refreshParameters();
            app.ReferenceTuner.ActiveLabel.Text='Default settings restored';
        end

        function exportData(app)
            app.requireData();
            [file, path] = uiputfile('*.mat', 'Save current Data as EELS class', 'eels_zest_data.mat');
            if isequal(file, 0), return; end
            filename = fullfile(path, file);
            EELS_ZEST_DataOps.save_eels(filename, app.RawCube, app.Ene, app.DataInfo);
            app.logStatus(sprintf('Current Data saved as EELS class: %s', filename));
        end

        function applySTA(app, source_override, manual_temperature_K)
            app.requireData();
            app.captureFitContext();
            if nargin < 2 || isempty(source_override)
                choices = {'Manual temperature', 'Stage 1 soft T', 'Stage 2 fitted T (fallback to Stage 1)'};
                [idx, ok] = listdlg('PromptString', 'STA temperature source', 'SelectionMode', 'single', ...
                    'ListString', choices, 'InitialValue', 3, 'ListSize', [360 110]);
                if ~ok, return; end
            else
                source_override = lower(char(string(source_override)));
                if strcmp(source_override, 'manual')
                    idx = 1;
                elseif any(strcmp(source_override, {'stage1', 'stage1_soft'}))
                    idx = 2;
                elseif any(strcmp(source_override, {'stage2', 'stage2_fit'}))
                    idx = 3;
                else
                    error('EELS:InvalidTemperatureSource', ...
                        'STA source must be manual, stage1, or stage2.');
                end
            end
            if nargin < 3, manual_temperature_K = NaN; end
            [ny, nx, ~] = size(app.RawCube);
            switch idx
                case 1
                    T = manual_temperature_K;
                    if ~isfinite(T)
                        answer = inputdlg({'Temperature (K)'}, 'Multiply STA', [1 36], {'300'});
                        if isempty(answer), return; end
                        T = str2double(answer{1});
                    end
                    if ~isfinite(T) || T <= 0, error('EELS:InvalidTemperature', 'Enter a positive temperature.'); end
                    [temperature_map, fallback_mask] = EELS_ZEST_DataOps.resolve_sta_temperature( ...
                        'manual', [ny nx], T, [], []);
                    source = 'manual';
                case 2
                    assert(isfield(app.FitContext, 'T_stage1_soft') && ...
                        isequal(size(app.FitContext.T_stage1_soft), [ny nx]), ...
                        'EELS:MissingStage1Temperature', 'No matching Stage 1 soft-temperature map is available.');
                    [temperature_map, fallback_mask] = EELS_ZEST_DataOps.resolve_sta_temperature( ...
                        'stage1', [ny nx], NaN, app.FitContext.T_stage1_soft, []);
                    source = 'stage1_T_soft';
                otherwise
                    assert(isfield(app.FitContext, 'T_stage2_fit') && ...
                        isequal(size(app.FitContext.T_stage2_fit), [ny nx]), ...
                        'EELS:MissingStage2Temperature', 'No matching Stage 2 fitted-temperature map is available.');
                    T1 = [];
                    if isfield(app.FitContext, 'T_stage1_soft'), T1 = app.FitContext.T_stage1_soft; end
                    [temperature_map, fallback_mask] = EELS_ZEST_DataOps.resolve_sta_temperature( ...
                        'stage2_fallback_stage1', [ny nx], NaN, T1, app.FitContext.T_stage2_fit);
                    source = 'stage2_T_fit_with_stage1_T_soft_fallback';
            end
            [cube_out, ene_out, sta_info] = EELS_ZEST_DataOps.apply_sta_cube(app.RawCube, app.Ene, temperature_map);
            app.RawCube = cube_out; app.Ene = ene_out;
            sta_info.temperature_source = source;
            sta_info.stage1_fallback_mask = fallback_mask;
            app.DataInfo.sta = sta_info;
            app.DataInfo = EELS_ZEST_DataOps.append_history(app.DataInfo, 'multiply_STA_loss_only', ...
                struct('temperature_source', source, 'fallback_count', nnz(fallback_mask)));
            app.persistCurrentDataForSolver('sta_loss_only');
            app.clearFitCaches(false);
            app.updatePixelLimits(); app.plotEnergyMap(); app.plotIntegratedSpectrum();
            app.logStatus(sprintf('STA applied to positive-loss Data using %s; valid pixels %d/%d.', ...
                source, nnz(sta_info.valid_temperature_mask), numel(sta_info.valid_temperature_mask)));
        end

        function openPeakFitTool(app)
            app.requireData();
            fig = EELS_ZEST_PeakFitTool(app.RawCube, app.Ene, app.DataInfo);
            app.registerOwnedFigure(fig);
            app.logStatus('Opened Robust Multi-Peak Fit. Full-map fitting retains parameter maps, not a fitted cube.');
        end

        function applyE2LossOnly(app)
            app.requireData();
            mask = app.Ene(:) > 0;
            if nnz(mask) < 2
                app.logStatus('E^2 transform skipped: no positive-loss channels were found.');
                return;
            end
            E = double(app.Ene(mask));
            scale = reshape(E .^ 2, 1, 1, []);
            app.RawCube = double(app.RawCube(:, :, mask)) .* scale;
            app.Ene = E(:);
            app.DataInfo.e2_loss_only = true;
            app.DataInfo.e2_loss_only_note = 'Energy-gain channels were removed; positive-loss intensities were multiplied by E^2.';
            app.persistCurrentDataForSolver('e2_loss_only');
            app.clearFitCaches();
            app.updatePixelLimits();
            app.plotEnergyMap();
            app.plotIntegratedSpectrum();
            app.logStatus('Applied E^2 loss-only transform: all energy-gain channels were removed in memory, and positive-loss data were multiplied by E^2.');
        end

        function basicZLPAlign(app)
            app.requireData();
            settings = EELS_PearsonSparse_FullFit.defaults(app.Settings);
            settings.stage1_zlp_align_enabled = true;
            [ny, nx, ne] = size(app.RawCube);
            aligned = NaN(ny, nx, ne);
            shifts = NaN(ny, nx);
            statuses = cell(ny, nx);
            n_pix = ny * nx;
            d = uiprogressdlg(app.UIFigure, 'Title', 'Basic ZLP align', ...
                'Message', sprintf('Aligning %d spectra to ZLP center...', n_pix), ...
                'Indeterminate', 'off', 'Value', 0, 'Cancelable', 'on');
            cleanup = onCleanup(@()app.closeProgressDialog(d));
            for pix = 1:n_pix
                if d.CancelRequested
                    app.logStatus('Basic ZLP align cancelled; data were not replaced.');
                    return;
                end
                [yy, xx] = ind2sub([ny, nx], pix);
                spec = double(squeeze(app.RawCube(yy, xx, :)));
                norm_val = sum(abs(spec), 'omitnan');
                if ~isfinite(norm_val) || norm_val <= eps
                    norm_val = max(abs(spec), [], 'omitnan');
                end
                if ~isfinite(norm_val) || norm_val <= eps
                    norm_val = 1;
                end
                [~, info] = EELS_PearsonSparse_FullFit.apply_stage1_zlp_alignment(app.Ene, spec ./ norm_val, settings);
                shift = 0;
                if isfield(info, 'applied_shift_meV') && isfinite(info.applied_shift_meV)
                    shift = info.applied_shift_meV;
                end
                shifted_axis = app.Ene(:) - shift;
                ynew = interp1(shifted_axis, spec(:), app.Ene(:), 'linear', NaN);
                bad = ~isfinite(ynew);
                if any(bad)
                    good = find(~bad & isfinite(app.Ene(:)));
                    if isempty(good)
                        ynew = spec(:);
                    else
                        ynew(bad) = interp1(app.Ene(good), ynew(good), app.Ene(bad), 'nearest', 'extrap');
                    end
                end
                aligned(yy, xx, :) = ynew;
                shifts(yy, xx) = shift;
                if isfield(info, 'status')
                    statuses{yy, xx} = info.status;
                else
                    statuses{yy, xx} = 'unknown';
                end
                if mod(pix, max(1, floor(n_pix ./ 40))) == 0 || pix == n_pix
                    d.Value = pix ./ n_pix;
                    d.Message = sprintf('Basic ZLP align: %d/%d | median shift %.3g meV', pix, n_pix, median(shifts(:), 'omitnan'));
                    drawnow limitrate;
                end
            end
            app.RawCube = aligned;
            app.DataInfo.basic_zlp_align = true;
            app.DataInfo.zlp_align_shift_map_meV = shifts;
            app.DataInfo.zlp_align_status_map = statuses;
            app.persistCurrentDataForSolver('basic_zlp_aligned');
            app.clearFitCaches();
            app.updatePixelLimits();
            app.plotEnergyMap();
            app.plotIntegratedSpectrum();
            app.logStatus(sprintf('Basic ZLP align complete. Median shift %.3g meV; max |shift| %.3g meV. Active solver file: %s', ...
                median(shifts(:), 'omitnan'), max(abs(shifts(:)), [], 'omitnan'), app.DataFile));
        end

        function svdDenoise(app)
            app.requireData();
            [ny, nx, ne] = size(app.RawCube);
            npix = ny .* nx;
            default_rank = min([12, max(1, npix - 1), max(1, ne - 1)]);
            answer = inputdlg({ ...
                'SVD/PCA rank (typical: 6-20; larger preserves more weak peaks)', ...
                'Optional energy smoothing after SVD, odd points (1 = off)'}, ...
                'EELS-ZEST SVD denoise', [1 72], {num2str(default_rank), '1'});
            if isempty(answer), return; end
            rk = max(1, round(str2double(answer{1})));
            smooth_pts = max(1, round(str2double(answer{2})));
            if ~isfinite(rk), rk = default_rank; end
            if ~isfinite(smooth_pts), smooth_pts = 1; end
            rk = min(rk, min(npix, ne));
            d = uiprogressdlg(app.UIFigure, 'Title', 'SVD denoise', ...
                'Message', sprintf('Building %d x %d spectrum matrix...', npix, ne), ...
                'Indeterminate', 'off', 'Value', 0);
            cleanup = onCleanup(@()app.closeProgressDialog(d));
            X = reshape(double(app.RawCube), npix, ne);
            finite = isfinite(X);
            col_ref = median(X, 1, 'omitnan');
            col_ref(~isfinite(col_ref)) = 0;
            for jj = 1:ne
                bad = ~finite(:, jj);
                if any(bad)
                    X(bad, jj) = col_ref(jj);
                end
            end
            d.Value = 0.25;
            d.Message = 'Computing low-rank SVD/PCA model...';
            mu = mean(X, 1, 'omitnan');
            mu(~isfinite(mu)) = 0;
            Xc = X - mu;
            try
                if rk < min(size(Xc)) - 1
                    [U, S, V] = svds(Xc, rk);
                else
                    [U, S, V] = svd(Xc, 'econ');
                    U = U(:, 1:rk);
                    S = S(1:rk, 1:rk);
                    V = V(:, 1:rk);
                end
            catch
                [U, S, V] = svd(Xc, 'econ');
                rk = min(rk, size(S, 1));
                U = U(:, 1:rk);
                S = S(1:rk, 1:rk);
                V = V(:, 1:rk);
            end
            d.Value = 0.70;
            d.Message = 'Reconstructing denoised datacube...';
            Xden = U * S * V' + mu;
            if smooth_pts > 1 && ne >= smooth_pts
                if mod(smooth_pts, 2) == 0, smooth_pts = smooth_pts + 1; end
                Xden = movmean(Xden, smooth_pts, 2, 'omitnan');
            end
            app.RawCube = reshape(Xden, ny, nx, ne);
            app.DataInfo.svd_denoise = true;
            app.DataInfo.svd_denoise_rank = rk;
            app.DataInfo.svd_denoise_energy_smooth_points = smooth_pts;
            app.DataInfo.svd_denoise_note = 'Low-rank SVD/PCA denoise.';
            app.persistCurrentDataForSolver('svd_denoised');
            app.clearFitCaches();
            app.updatePixelLimits();
            app.plotEnergyMap();
            app.plotIntegratedSpectrum();
            d.Value = 1;
            app.logStatus(sprintf('SVD denoise complete: rank %d, energy smoothing %d point(s). Active solver file: %s', ...
                rk, smooth_pts, app.DataFile));
        end

        function importFitSignalToData(app, stage2, options)
            if nargin < 2, stage2 = true; end
            if nargin < 3, options = []; end
            stage = 'stage1'; if stage2, stage = 'stage2'; end
            cuberes = app.cubeResult(stage2);
            single_result = app.singleResult(stage2);
            if isempty(single_result) && ~isempty(app.CurrentResult) && (~stage2 || app.resultHasStage2(app.CurrentResult))
                single_result = app.CurrentResult;
            end
            if isempty(single_result) && isempty(cuberes)
                app.logStatus('No fitted result is available to send back to the Data tab.');
                return;
            end
            try
                if ~isempty(cuberes), first = app.firstValidResult(cuberes.results); else, first = single_result; end
                if isempty(first), error('EELS:NoValidResult', 'No valid fitted result is available.'); end
                if isempty(options)
                    options = EELS_ZEST_CutoffDialog(first, stage);
                else
                    options = app.normalizeCutoffOptions(options, stage);
                end
                if options.cancelled, return; end
                if ~isempty(cuberes) && isfield(cuberes, 'results')
                    ene = app.fieldOr(first, 'ene_corrected', app.fieldOr(first, 'ene_raw', app.Ene));
                    map_size = app.cuberesMapSize(cuberes);
                    if isempty(map_size), map_size = [1, numel(cuberes.results)]; end
                    sigcube = NaN(map_size(1), map_size(2), numel(ene));
                    kernelcube = NaN(map_size(1), map_size(2), numel(ene), 'single');
                    cutoff_map = NaN(map_size); valid_map = false(map_size);
                    retained_max_map = NaN(map_size); removed_count_map = NaN(map_size);
                    for ii = 1:numel(cuberes.results)
                        res = app.resultAtIndex(cuberes.results, ii);
                        [yy, xx] = app.resultPixelPosition(cuberes, ii, res, map_size);
                        if isempty(res) || ~isfinite(yy) || ~isfinite(xx) || ...
                                yy < 1 || xx < 1 || yy > map_size(1) || xx > map_size(2)
                            continue;
                        end
                        cutoff = app.cutoffForResult(res, options);
                        if ~isfinite(cutoff), continue; end
                        [sig, detail] = EELS_ZEST_DataOps.reconstruct_sparse_signal(res, stage, cutoff);
                        sigcube(yy, xx, :) = sig(:); cutoff_map(yy, xx) = cutoff; valid_map(yy, xx) = true;
                        kernelcube(yy, xx, :) = single(detail.kernel_norm(:));
                        retained_em = detail.EM_meV(detail.retained_mask);
                        if ~isempty(retained_em), retained_max_map(yy, xx) = max(retained_em); end
                        removed_count_map(yy, xx) = nnz(~detail.retained_mask);
                    end
                    if ~any(valid_map(:))
                        error('EELS:NoValidReconstruction', ...
                            'No fitted pixel produced a valid cutoff/reconvolved spectrum.');
                    end
                    context = app.captureFitContext(map_size);
                    app.RawCube = sigcube;
                    app.Ene = ene(:);
                    label = [stage, ' cutoff/reconvolved phonon cube'];
                else
                    res = single_result;
                    ene = app.fieldOr(res, 'ene_corrected', app.fieldOr(res, 'ene_raw', app.Ene));
                    cutoff = app.cutoffForResult(res, options);
                    if ~isfinite(cutoff), error('EELS:MissingAttentionCutoff', 'No valid attention cutoff is available.'); end
                    [sig, detail] = EELS_ZEST_DataOps.reconstruct_sparse_signal(res, stage, cutoff);
                    app.RawCube = reshape(sig, 1, 1, []);
                    app.Ene = ene(:);
                    cutoff_map = cutoff; valid_map = true;
                    kernelcube = reshape(single(detail.kernel_norm), 1, 1, []);
                    retained_em = detail.EM_meV(detail.retained_mask);
                    retained_max_map = NaN;
                    if ~isempty(retained_em), retained_max_map = max(retained_em); end
                    removed_count_map = nnz(~detail.retained_mask);
                    context = app.captureFitContext([1 1]);
                    label = [stage, ' cutoff/reconvolved phonon spectrum'];
                end
                required_temperature = 'T_stage1_soft';
                if stage2, required_temperature = 'T_stage2_fit'; end
                if ~isfield(context, required_temperature) || ...
                        ~isequal(size(context.(required_temperature)), size(app.RawCube, [1 2]))
                    error('EELS:MissingTransferredTemperature', ...
                        'The fitted-temperature context does not match the reconstructed Data map.');
                end
                app.DataInfo.source = ['EELS-ZEST ', label];
                app.DataInfo.sparse_source_cutoff_mode = options.mode;
                app.DataInfo.sparse_source_cutoff_map_meV = cutoff_map;
                app.DataInfo.sparse_source_cutoff_valid_mask = valid_map;
                app.DataInfo.sparse_source_retained_max_EM_map_meV = retained_max_map;
                app.DataInfo.sparse_source_removed_count_map = removed_count_map;
                app.DataInfo.sparse_source_cutoff_semantics = ...
                    'Cutoff applies to sparse source centers before fitted-ZLP convolution; the output is not hard-clipped.';
                app.DataInfo.sparse_source_stage = stage;
                app.DataInfo.sparse_source_json = options.json_meta;
                app.FitContext.source_stage = stage;
                app.DataInfo.eels_zest_fit_context = app.FitContext;
                app.DataInfo.sparse_source_output_composition = ...
                    'cutoff/reconvolved sparse phonon source only';
                app.DataInfo.sparse_source_background_subtracted = true;
                app.DataInfo.sparse_source_additive_zlp_included = false;
                app.DataInfo.fitted_zlp_kernel_cube = kernelcube;
                app.DataInfo.fitted_zlp_kernel_energy_meV = app.Ene(:);
                app.DataInfo.fitted_zlp_kernel_normalization = 'unit sum per pixel';
                app.DataInfo = EELS_ZEST_DataOps.append_history(app.DataInfo, 'sparse_source_cutoff_reconvolution', ...
                    struct('stage', stage, 'mode', options.mode, 'valid_count', nnz(valid_map)));
                app.DataFile = '';
                app.persistCurrentDataForSolver([stage, '_cutoff_reconvolved']);
                app.clearFitCaches(false);
                app.updatePixelLimits();
                app.plotEnergyMap();
                app.plotIntegratedSpectrum();
                app.TabGroup.SelectedTab = app.DataTab;
                app.logStatus(sprintf(['Sent %s to Data; valid pixels %d/%d; cutoff %.3g-%.3g meV; ', ...
                    'max retained source center %.3g meV.'], label, nnz(valid_map), numel(valid_map), ...
                    min(cutoff_map(valid_map)), max(cutoff_map(valid_map)), ...
                    max(retained_max_map(valid_map), [], 'omitnan')));
            catch ME
                location = '';
                if ~isempty(ME.stack)
                    location = sprintf(' [%s:%d]', ME.stack(1).name, ME.stack(1).line);
                end
                app.logStatus(sprintf('Send signal to Data failed (%s): %s%s', ...
                    ME.identifier, ME.message, location));
            end
        end

        function sig = resultSignalSpectrum(~, res)
            sig = [];
            if ~(isstruct(res) && ~isfield(res, 'error')), return; end
            if isfield(res, 'stage2') && isstruct(res.stage2) && isfield(res.stage2, 'vib_norm') && ...
                    isfield(res, 'norm_val') && ~isempty(res.stage2.vib_norm)
                sig = res.stage2.vib_norm(:) .* res.norm_val;
            elseif isfield(res, 'vib_only') && ~isempty(res.vib_only)
                sig = res.vib_only(:);
            elseif isfield(res, 'raw_data') && isfield(res, 'zlp_only')
                sig = max(res.raw_data(:) - res.zlp_only(:), 0);
            end
        end

        function persistCurrentDataForSolver(app, tag)
            if isempty(app.RawCube) || isempty(app.Ene)
                return;
            end
            if nargin < 2 || isempty(tag)
                tag = 'preprocessed';
            end
            try
                outdir = fullfile(app.RootDir, 'test_outputs', 'gui_preprocessed');
                if ~exist(outdir, 'dir'), mkdir(outdir); end
                safe_tag = regexprep(char(tag), '[^A-Za-z0-9_-]', '_');
                filename = fullfile(outdir, sprintf('eels_zest_%s_%s.mat', ...
                    safe_tag, char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
                cube = app.RawCube;
                ene = app.Ene;
                info = app.DataInfo;
                save(filename, 'cube', 'ene', 'info', '-v7.3');
                app.DataFile = filename;
                app.DataInfo.source_file = filename;
                app.DataInfo.source_variable = 'cube';
                app.DataInfo.preprocessed_for_solver = true;
            catch ME
                app.logStatus(['Could not persist the in-memory dataset for full-datacube solver runs: ', ME.message]);
            end
        end

        function exportResult(app, stage2)
            if stage2
                single = app.Stage2SingleResult;
                cuberes = app.Stage2CubeResult;
                tag = 'stage2';
            else
                single = app.Stage1SingleResult;
                cuberes = app.Stage1CubeResult;
                tag = 'stage1';
            end
            scopes = {};
            payloads = {};
            workspace_names = {};
            if ~isempty(single)
                scopes{end + 1} = 'Selected pixel'; %#ok<AGROW>
                payloads{end + 1} = single; %#ok<AGROW>
                workspace_names{end + 1} = ['eels_zest_', tag, '_selected_result']; %#ok<AGROW>
            end
            if ~isempty(cuberes)
                scopes{end + 1} = 'Complete datacube result'; %#ok<AGROW>
                payloads{end + 1} = cuberes; %#ok<AGROW>
                workspace_names{end + 1} = ['eels_zest_', tag, '_cuberes']; %#ok<AGROW>
            end
            if isempty(scopes)
                app.logStatus('No result is available to export.');
                return;
            end
            options = EELS_ZEST_ExportDialog(['Export ', tag], scopes, {'MAT', 'CSV'});
            if options.cancelled, return; end
            payload = payloads{options.scope_index};
            if strcmpi(options.format, 'CSV')
                if strcmp(scopes{options.scope_index}, 'Selected pixel')
                    payload = app.empmTable(app.makeEMPMExport(payload));
                else
                    payload = app.cubeSummaryTable(payload);
                end
            end
            app.writeChosenExport(payload, options, ['eels_zest_', tag, '_result'], ...
                workspace_names{options.scope_index});
        end

        function exportCurrentResult(app)
            if isempty(app.CurrentResult)
                app.logStatus('No selected result is available to export.');
                return;
            end
            options = EELS_ZEST_ExportDialog('Export selected result', {'Selected pixel'}, {'MAT', 'CSV'});
            if options.cancelled, return; end
            payload = app.CurrentResult;
            if strcmpi(options.format, 'CSV')
                payload = app.empmTable(app.makeEMPMExport(payload));
            end
            app.writeChosenExport(payload, options, 'eels_zest_selected_result', 'eels_zest_current_result');
        end

        function empm = makeEMPMExport(app, res)
            empm = struct();
            [em, ps, pa] = app.extractPM(res, 'stage1');
            if ~isempty(em)
                empm.stage1 = struct('EM_meV', em(:), 'PM_stokes', ps(:), 'PM_antistokes', pa(:), ...
                    'signed_EM_meV', [-flipud(em(:)); em(:)], 'signed_PM', [flipud(pa(:)); ps(:)]);
            end
            if app.resultHasStage2(res)
                [em, ps, pa] = app.extractPM(res, 'stage2');
                empm.stage2 = struct('EM_meV', em(:), 'PM_stokes', ps(:), 'PM_antistokes', pa(:), ...
                    'signed_EM_meV', [-flipud(em(:)); em(:)], 'signed_PM', [flipud(pa(:)); ps(:)]);
            end
        end

        function [results, map_size] = pmSourceResults(app, src)
            switch lower(src)
                case 'stage2'
                    cuberes = app.Stage2CubeResult;
                    single = app.Stage2SingleResult;
                    require_stage2 = true;
                case 'stage1'
                    cuberes = app.Stage1CubeResult;
                    single = app.Stage1SingleResult;
                    require_stage2 = false;
                otherwise
                    cuberes = app.Stage2CubeResult;
                    if isempty(cuberes), cuberes = app.Stage1CubeResult; end
                    single = app.CurrentResult;
                    require_stage2 = false;
            end
            if ~isempty(cuberes) && isfield(cuberes, 'results') && isfield(cuberes, 'cube_info')
                map_size = cuberes.cube_info.size(1:2);
                results = cuberes.results;
                if numel(results) ~= prod(map_size)
                    full_results = cell(1, prod(map_size));
                    if isfield(cuberes, 'sample_pixels')
                        for ii = 1:min(numel(cuberes.sample_pixels), numel(results))
                            full_results{cuberes.sample_pixels(ii)} = results{ii};
                        end
                    end
                    results = full_results;
                end
                if require_stage2 && ~app.cuberesHasStage2(cuberes), results = {}; map_size = [1 1]; end
            elseif ~isempty(single)
                results = {single};
                map_size = [1 1];
            else
                results = {};
                map_size = [1 1];
            end
        end

        function [em, ps, pa] = extractPM(app, res, src)
            em = []; ps = []; pa = [];
            if (strcmpi(src, 'stage2') || app.resultHasStage2(res)) && isfield(res, 'stage2') && ...
                    isfield(res.stage2, 'EM_final') && isfield(res.stage2, 'PM_stokes') && isfield(res.stage2, 'PM_antistokes')
                em = res.stage2.EM_final(:);
                ps = res.stage2.PM_stokes(:);
                pa = res.stage2.PM_antistokes(:);
            elseif isfield(res, 'EM_final') && isfield(res, 'PM_stokes') && isfield(res, 'PM_antistokes')
                em = res.EM_final(:);
                ps = res.PM_stokes(:);
                pa = res.PM_antistokes(:);
            end
        end

        function updateProgress(app, dlg, done, total, elapsed, eta, med, label)
            if ~isvalid(dlg), return; end
            dlg.Value = max(0, min(1, done ./ max(total, 1)));
            if ~isfinite(eta), eta = NaN; end
            dlg.Message = sprintf('%s: %d/%d | elapsed %.1f min | %.2f s/pixel | ETA %.1f min', ...
                char(label), done, total, elapsed ./ 60, med, eta ./ 60);
            drawnow limitrate;
            if mod(done, max(1, round(total ./ 20))) == 0 || done == total
                app.logStatus(dlg.Message);
            end
        end

        function updatePixelLimits(app)
            if isempty(app.RawCube), return; end
            [ny, nx, ~] = size(app.RawCube);
            ctrls = {app.Stage1Controls, app.Stage2Controls};
            for ii = 1:numel(ctrls)
                c = ctrls{ii};
                if isfield(c, 'py') && isvalid(c.py)
                    c.py.Limits = [1 ny]; c.py.Value = min(c.py.Value, ny);
                    c.px.Limits = [1 nx]; c.px.Value = min(c.px.Value, nx);
                end
            end
        end

        function pickPixelFromMap(app, stage2)
            ax = app.mapAxes(stage2);
            cp = ax.CurrentPoint;
            xx = round(cp(1, 1));
            yy = round(cp(1, 2));
            app.selectPixelFromMapCoordinates(stage2, yy, xx);
        end

        function pickPixelWithCursor(app, stage2)
            ax = app.mapAxes(stage2);
            cuberes = app.cubeResult(stage2);
            if isempty(cuberes)
                app.logStatus('No datacube result map is available for cursor picking.');
                return;
            end
            try
                app.logStatus('Map cursor selection: double-click to confirm. Ctrl+Q cancels, Ctrl+H confirms and plots.');
                roi = drawpoint(ax, 'Color', [0.95 0.25 0.05], 'LineWidth', 1.4);
                [pos, cancelled] = app.waitForROIWithHotkeys(roi, ...
                    'Map cursor selection: double-click to confirm. Ctrl+Q cancels, Ctrl+H confirms and plots.');
                if cancelled || isempty(pos)
                    app.logStatus('Map cursor selection cancelled.');
                    return;
                end
            catch ME
                app.logStatus(['Cursor picking requires MATLAB ROI tools: ', ME.message]);
                return;
            end
            app.selectPixelFromMapCoordinates(stage2, round(pos(2)), round(pos(1)));
        end

        function selectPixelFromMapCoordinates(app, stage2, yy, xx)
            cuberes = app.cubeResult(stage2);
            if isempty(cuberes) || ~isfield(cuberes, 'results'), return; end
            map_size = app.cuberesMapSize(cuberes);
            if isempty(map_size) || xx < 1 || xx > map_size(2) || yy < 1 || yy > map_size(1)
                return;
            end
            ctrl = app.stageControls(stage2);
            ctrl.py.Value = yy;
            ctrl.px.Value = xx;
            pix = sub2ind(map_size, yy, xx);
            idx = [];
            if isfield(cuberes, 'sample_pixels') && ~isempty(cuberes.sample_pixels)
                idx = find(cuberes.sample_pixels(:) == pix, 1, 'first');
            elseif numel(cuberes.results) == prod(map_size)
                idx = pix;
            end
            if ~isempty(idx) && idx <= numel(cuberes.results)
                res = cuberes.results{idx};
                if ~(isstruct(res) && isfield(res, 'error'))
                    app.CurrentResult = res;
                    if stage2
                        app.Stage2SingleResult = res;
                    else
                        app.Stage1SingleResult = res;
                    end
                    app.plotSpectrumResult(res, stage2);
                    app.markSelectedPixel(stage2, yy, xx);
                    app.plotDiagnosticSelection();
                    app.updateRLDeconvControl();
                end
            end
        end

        function tf = datacubeStopRequested(app, dlg)
            tf = logical(app.StopRequested);
            if ~tf && ~isempty(dlg) && isvalid(dlg)
                try
                    tf = logical(dlg.CancelRequested);
                catch
                    tf = false;
                end
            end
            if tf
                app.StopRequested = true;
            end
        end

        function markSelectedPixel(app, stage2, yy, xx)
            ax = app.mapAxes(stage2);
            if stage2
                app.deleteGraphicsHandles('Stage2CursorHandles');
            else
                app.deleteGraphicsHandles('Stage1CursorHandles');
            end
            hold(ax, 'on');
            handles = { ...
                plot(ax, xx, yy, 'wo', 'MarkerSize', 9, 'LineWidth', 1.8), ...
                plot(ax, xx, yy, '+', 'Color', [0.85 0.08 0.10], 'MarkerSize', 9, 'LineWidth', 1.4)};
            for ii = 1:numel(handles)
                try
                    handles{ii}.HitTest = 'off';
                    handles{ii}.PickableParts = 'none';
                catch
                end
            end
            if stage2
                app.Stage2CursorHandles = handles;
            else
                app.Stage1CursorHandles = handles;
            end
            hold(ax, 'off');
        end

        function res = firstValidResult(~, results)
            res = [];
            if iscell(results)
                for ii = 1:numel(results)
                    if isstruct(results{ii}) && ~isfield(results{ii}, 'error')
                        res = results{ii};
                        return;
                    end
                end
            elseif isstruct(results)
                res = results(1);
            end
        end

        function map_size = cuberesMapSize(~, cuberes)
            map_size = [];
            if isfield(cuberes, 'R2_map') && ~isempty(cuberes.R2_map)
                map_size = size(cuberes.R2_map);
            elseif isfield(cuberes, 'T_map_raw') && ~isempty(cuberes.T_map_raw)
                map_size = size(cuberes.T_map_raw);
            elseif isfield(cuberes, 'cube_info') && isfield(cuberes.cube_info, 'size')
                map_size = cuberes.cube_info.size(1:2);
            end
        end

        function inheritWindow(app, stage2)
            ctrl = app.stageControls(stage2);
            ctrl.cmin.Value = app.DataEnergyMinEdit.Value;
            ctrl.cmax.Value = app.DataEnergyMaxEdit.Value;
        end

        function pickOutputDir(~, edit_field)
            p = uigetdir(edit_field.Value, 'Choose output folder');
            if isequal(p, 0), return; end
            edit_field.Value = p;
        end

        function closeApp(app)
            if ~isempty(app.ReferenceTuner)&&isvalid(app.ReferenceTuner)&&app.ReferenceTuner.Running
                app.ReferenceTuner.stop();
                app.logStatus('Stopping reference tuning; wait for the current candidate and dataset export before closing.');
                return;
            end
            app.deleteGraphicsHandles('DataCursorHandles');
            app.deleteGraphicsHandles('Stage1CursorHandles');
            app.deleteGraphicsHandles('Stage2CursorHandles');
            app.deleteGraphicsHandles('StackCursorHandles');
            for ii = 1:numel(app.OwnedFigures)
                try
                    if isgraphics(app.OwnedFigures{ii}), delete(app.OwnedFigures{ii}); end
                catch
                end
            end
            app.OwnedFigures = {};
            app.RLDiagnosticCache = struct('keys', {{}}, 'values', {{}}, 'capacity', 3);
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure), delete(app.UIFigure); end
        end

        function clearFitCaches(app, preserve_rl_source)
            if nargin < 2 || isempty(preserve_rl_source)
                preserve_rl_source = false;
            end
            app.Stage1SingleResult = [];
            app.Stage2SingleResult = [];
            app.Stage1CubeResult = [];
            app.Stage2CubeResult = [];
            app.CurrentResult = [];
            app.LastPMResult = [];
            if ~preserve_rl_source
                app.RLSourceResult = [];
                app.RLSourceCubeResult = [];
                app.RLDeconvResultCube = [];
            end
            app.RLDiagnosticCache = struct('keys', {{}}, 'values', {{}}, 'capacity', 3);
            app.deleteGraphicsHandles('Stage1CursorHandles');
            app.deleteGraphicsHandles('Stage2CursorHandles');
            app.updateRLDeconvControl();
        end

        function updateRLDeconvControl(app)
            has_data_context = app.hasDataRLKernelContext();
            can_run = has_data_context || ~isempty(app.activeRLDeconvCubeResult()) || ...
                app.hasRestorableRLSource();
            can_diag = false;
            if ~isempty(app.DataSelectedY) && ~isempty(app.DataSelectedX)
                if has_data_context
                    can_diag = ~isempty(app.currentDataRLResult(app.DataSelectedY, app.DataSelectedX));
                else
                    can_diag = ~isempty(app.rlSourceResultForPixel(app.DataSelectedY, app.DataSelectedX));
                end
            end
            if ~isempty(app.DataRLDeconvButton) && isvalid(app.DataRLDeconvButton)
                if can_run
                    app.DataRLDeconvButton.Enable = 'on';
                    app.DataRLDeconvButton.BackgroundColor = [0.42 0.12 0.68];
                    app.DataRLDeconvButton.Tooltip = 'Batch RL-deconvolve every fitted datacube pixel and replace the Data cube.';
                else
                    app.DataRLDeconvButton.Enable = 'off';
                    app.DataRLDeconvButton.BackgroundColor = [0.55 0.63 0.70];
                    app.DataRLDeconvButton.Tooltip = 'Run or load a Stage 1/2 datacube result first.';
                end
            end
            if ~isempty(app.DataRLDiagnosticButton) && isvalid(app.DataRLDiagnosticButton)
                if can_diag
                    app.DataRLDiagnosticButton.Enable = 'on';
                    app.DataRLDiagnosticButton.BackgroundColor = [0.02 0.34 0.62];
                    app.DataRLDiagnosticButton.Tooltip = 'Open RL diagnostic for the selected Data-map pixel.';
                else
                    app.DataRLDiagnosticButton.Enable = 'off';
                    app.DataRLDiagnosticButton.BackgroundColor = [0.55 0.63 0.70];
                    app.DataRLDiagnosticButton.Tooltip = 'Click a fitted Data-map pixel first.';
                end
            end
            if ~isempty(app.DataRLStatusLabel) && isvalid(app.DataRLStatusLabel)
                if can_run
                    if can_diag
                        app.DataRLStatusLabel.Text = sprintf('RL Conv ready; RL diag ready for pixel (%d,%d).', ...
                            app.DataSelectedY, app.DataSelectedX);
                    else
                        app.DataRLStatusLabel.Text = 'RL Conv ready: click a fitted pixel to enable RL diag.';
                    end
                else
                    app.DataRLStatusLabel.Text = 'RL Conv is locked until a Stage 1/2 datacube result is available.';
                end
            end
            drawnow limitrate;
        end

        function res = activeRLDeconvResult(app)
            candidates = {app.CurrentResult, app.Stage2SingleResult, app.Stage1SingleResult, app.RLSourceResult};
            for ii = 1:numel(candidates)
                if app.resultCanRunRLDeconv(candidates{ii})
                    res = candidates{ii};
                    return;
                end
            end
            res = app.firstRLDeconvResult(app.Stage2CubeResult);
            if ~isempty(res)
                return;
            end
            res = app.firstRLDeconvResult(app.Stage1CubeResult);
        end

        function cuberes = activeRLDeconvCubeResult(app)
            cuberes = [];
            candidates = {app.Stage2CubeResult, app.Stage1CubeResult, app.RLSourceCubeResult};
            for ii = 1:numel(candidates)
                candidate = candidates{ii};
                if isstruct(candidate) && isfield(candidate, 'results') && ...
                        ~isempty(app.firstRLDeconvResult(candidate))
                    cuberes = candidate;
                    return;
                end
            end
        end

        function res = firstRLDeconvResult(app, cuberes)
            res = [];
            if ~(isstruct(cuberes) && isfield(cuberes, 'results') && ~isempty(cuberes.results))
                return;
            end
            results = cuberes.results;
            if iscell(results)
                for ii = 1:numel(results)
                    candidate = results{ii};
                    if app.resultCanRunRLDeconv(candidate)
                        res = candidate;
                        return;
                    end
                end
            elseif isstruct(results)
                for ii = 1:numel(results)
                    candidate = results(ii);
                    if app.resultCanRunRLDeconv(candidate)
                        res = candidate;
                        return;
                    end
                end
            end
        end

        function res = resultAtIndex(~, results, ii)
            res = [];
            if iscell(results)
                if ii <= numel(results)
                    res = results{ii};
                end
            elseif isstruct(results) && ii <= numel(results)
                res = results(ii);
            end
            if ~(isstruct(res) && ~isfield(res, 'error'))
                res = [];
            end
        end

        function [yy, xx] = resultPixelPosition(~, cuberes, ii, res, map_size)
            yy = NaN;
            xx = NaN;
            if isfield(cuberes, 'sample_pixels') && ii <= numel(cuberes.sample_pixels) && ...
                    isfinite(cuberes.sample_pixels(ii)) && cuberes.sample_pixels(ii) >= 1 && ...
                    cuberes.sample_pixels(ii) <= prod(map_size)
                [yy, xx] = ind2sub(map_size, round(cuberes.sample_pixels(ii)));
            elseif isstruct(res) && isfield(res, 'pixel_y') && isfield(res, 'pixel_x') && ...
                    isfinite(res.pixel_y) && isfinite(res.pixel_x)
                yy = round(res.pixel_y);
                xx = round(res.pixel_x);
            elseif isfield(cuberes, 'results') && numel(cuberes.results) == prod(map_size)
                [yy, xx] = ind2sub(map_size, ii);
            end
        end

        function [res, cached_rl] = rlSourceResultForPixel(app, yy, xx)
            res = [];
            cached_rl = app.getRLDiagnosticCache(yy, xx);
            candidates = {app.Stage2CubeResult, app.Stage1CubeResult, app.RLSourceCubeResult};
            for cc = 1:numel(candidates)
                cuberes = candidates{cc};
                if ~(isstruct(cuberes) && isfield(cuberes, 'results')), continue; end
                map_size = app.cuberesMapSize(cuberes);
                if isempty(map_size), continue; end
                results = cuberes.results;
                for ii = 1:numel(results)
                    candidate = app.resultAtIndex(results, ii);
                    if isempty(candidate), continue; end
                    [cy, cx] = app.resultPixelPosition(cuberes, ii, candidate, map_size);
                    if cy == yy && cx == xx && app.resultCanRunRLDeconv(candidate)
                        res = candidate;
                        return;
                    end
                end
            end
        end

        function tf = resultCanRunRLDeconv(~, res)
            direct_input = isstruct(res) && isfield(res, 'rl_input_is_background_subtracted') && ...
                isscalar(res.rl_input_is_background_subtracted) && ...
                logical(res.rl_input_is_background_subtracted) && isfield(res, 'rl_kernel_norm');
            tf = isstruct(res) && isfield(res, 'ene_corrected') && ...
                (isfield(res, 'raw_data_norm') || ...
                (isfield(res, 'raw_data') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0)) && ...
                (direct_input || (isfield(res, 'stage1') && isstruct(res.stage1) && isfield(res.stage1, 'params')) || ...
                isfield(res, 'zlp_params'));
        end

        function [cuberes, forced_stage2] = extractLoadedCubeResult(~, S)
            forced_stage2 = [];
            if isfield(S, 'cuberes')
                cuberes = S.cuberes;
            elseif isfield(S, 'stage2_cuberes')
                cuberes = S.stage2_cuberes;
                forced_stage2 = true;
            else
                cuberes = S.stage1_cuberes;
                forced_stage2 = false;
            end
        end

        function requireData(app)
            if isempty(app.RawCube) || isempty(app.Ene)
                error('Load a datacube first.');
            end
        end

        function requireDataFile(app)
            app.requireData();
            if isempty(app.DataFile) || ~isfile(app.DataFile)
                error('Full datacube processing requires a file-backed .mat or .dm3 input.');
            end
        end

        function pmres = requirePM(app)
            if isempty(app.LastPMResult) || ~isstruct(app.LastPMResult)
                app.buildPMAnalysis();
            end
            pmres = app.LastPMResult;
        end

        function ctrl = stageControls(app, stage2)
            if stage2, ctrl = app.Stage2Controls; else, ctrl = app.Stage1Controls; end
        end

        function ax = mapAxes(app, stage2)
            if stage2, ax = app.Stage2MapAxes; else, ax = app.Stage1MapAxes; end
        end

        function ax = spectrumAxes(app, stage2)
            if stage2, ax = app.Stage2SpectrumAxes; else, ax = app.Stage1SpectrumAxes; end
        end

        function res = singleResult(app, stage2)
            if stage2, res = app.Stage2SingleResult; else, res = app.Stage1SingleResult; end
        end

        function cuberes = cubeResult(app, stage2)
            if stage2, cuberes = app.Stage2CubeResult; else, cuberes = app.Stage1CubeResult; end
        end

        function tf = isReusableStage1Cube(~, cuberes)
            tf = isstruct(cuberes) && isfield(cuberes, 'results') && ~isempty(cuberes.results);
        end

        function tf = resultHasStage2(~, res)
            tf = isstruct(res) && isfield(res, 'stage2') && isstruct(res.stage2) && ...
                isfield(res.stage2, 'enabled') && res.stage2.enabled && ...
                isfield(res.stage2, 'EM_final') && ~isempty(res.stage2.EM_final);
        end

        function tf = cuberesHasStage2(app, cuberes)
            tf = false;
            if ~isstruct(cuberes) || ~isfield(cuberes, 'results'), return; end
            for ii = 1:numel(cuberes.results)
                if app.resultHasStage2(cuberes.results{ii})
                    tf = true;
                    return;
                end
            end
        end

        function ttl = resultTitle(~, res, stage2)
            m = struct();
            if isfield(res, 'metrics'), m = res.metrics; end
            R2 = NaN; NRMSE = NaN;
            if isfield(m, 'R2'), R2 = m.R2; end
            if isfield(m, 'NRMSE'), NRMSE = m.NRMSE; end
            if stage2
                T = NaN; rel = NaN;
                if isfield(res, 'stage2') && isstruct(res.stage2) && isfield(res.stage2, 'T_fit')
                    T = res.stage2.T_fit;
                elseif isfield(res, 'T_fit')
                    T = res.T_fit;
                end
                if isfield(res, 'stage2') && isstruct(res.stage2) && isfield(res.stage2, 'reliability_score')
                    rel = res.stage2.reliability_score;
                elseif isfield(res, 'reliability_score')
                    rel = res.reliability_score;
                end
                ttl = sprintf('Stage 2 | R2 %.6g | NRMSE %.4g | T %.1f K | rel %.2f', R2, NRMSE, T, rel);
            else
                ttl = sprintf('Stage 1 | R2 %.6g | NRMSE %.4g', R2, NRMSE);
            end
        end

        function line = formatResultLine(app, res, label)
            line = sprintf('%s | %s', label, app.resultTitle(res, app.resultHasStage2(res)));
        end

        function line = formatCubeLine(~, cuberes)
            R2 = NaN; nrmse = NaN;
            if isfield(cuberes, 'R2_map'), R2 = median(cuberes.R2_map(:), 'omitnan'); end
            if isfield(cuberes, 'NRMSE_map'), nrmse = median(cuberes.NRMSE_map(:), 'omitnan'); end
            if isfield(cuberes, 'benchmark')
                line = sprintf('Datacube complete | median R2 %.6g | median NRMSE %.4g | %.1f s total', ...
                    R2, nrmse, cuberes.benchmark.total_seconds);
            else
                line = sprintf('Datacube complete | median R2 %.6g | median NRMSE %.4g', R2, nrmse);
            end
        end

        function val = fieldOr(~, s, name, fallback)
            if isstruct(s) && isfield(s, name)
                val = s.(name);
            else
                val = fallback;
            end
        end

        function exportPayloadToFile(app, payload, title_text, base_name, formats, workspace_name)
            options = EELS_ZEST_ExportDialog(title_text, {'Current content'}, formats);
            if options.cancelled, return; end
            if strcmpi(options.format, 'CSV') && isstruct(payload) && isfield(payload, 'pm_cube')
                payload = app.pmLongTable(payload);
            end
            app.writeChosenExport(payload, options, base_name, workspace_name);
        end

        function writeChosenExport(app, payload, options, base_name, workspace_name)
            format = lower(char(options.format));
            switch format
                case 'mat'
                    filter = {'*.mat', 'MAT-file (*.mat)'}; ext = '.mat';
                case 'csv'
                    filter = {'*.csv', 'CSV file (*.csv)'}; ext = '.csv';
                otherwise
                    error('EELS:UnsupportedExportFormat', 'Unsupported export format: %s', format);
            end
            [file, path] = uiputfile(filter, 'Save export', [base_name, ext]);
            if isequal(file, 0), return; end
            filename = fullfile(path, file);
            EELS_ZEST_Export.save_payload(filename, payload, format, workspace_name);
            if options.assign_workspace
                assignin('base', matlab.lang.makeValidName(workspace_name), payload);
            end
            app.logStatus(sprintf('Exported %s.', filename));
        end

        function T = empmTable(~, empm)
            T = table();
            stages = {'stage1', 'stage2'};
            for ii = 1:numel(stages)
                name = stages{ii};
                if ~isfield(empm, name), continue; end
                data = empm.(name);
                n = numel(data.EM_meV);
                Ti = table(repmat(string(name), n, 1), double(data.EM_meV(:)), ...
                    double(data.PM_stokes(:)), double(data.PM_antistokes(:)), ...
                    'VariableNames', {'stage', 'energy_meV', 'loss_intensity', 'gain_intensity'});
                T = [T; Ti]; %#ok<AGROW>
            end
        end

        function T = cubeSummaryTable(app, cuberes)
            map_size = app.cuberesMapSize(cuberes);
            if isempty(map_size)
                error('EELS:MissingMapSize', 'The datacube result has no map dimensions.');
            end
            [yy, xx] = ndgrid(1:map_size(1), 1:map_size(2));
            T = table(yy(:), xx(:), 'VariableNames', {'pixel_y', 'pixel_x'});
            names = {'R2_map', 'NRMSE_map', 'T_map_raw', 'T_map_data_valid', ...
                'T_map_joint', 'T_uncertainty_joint_map', 'reliability_map', 'fit_status_map'};
            for ii = 1:numel(names)
                name = names{ii};
                if ~isfield(cuberes, name), continue; end
                value = cuberes.(name);
                if (isnumeric(value) || islogical(value)) && isequal(size(value), map_size)
                    T.(matlab.lang.makeValidName(name)) = value(:);
                end
            end
        end

        function T = pmLongTable(~, pmres)
            cube = double(pmres.pm_cube);
            axis_meV = double(pmres.axis_meV(:));
            assert(size(cube, 3) == numel(axis_meV), 'EELS:PMShapeMismatch', ...
                'Sparse phonon intensity cube does not match its energy axis.');
            [yy, xx, ee] = ndgrid(1:size(cube, 1), 1:size(cube, 2), 1:size(cube, 3));
            T = table(yy(:), xx(:), axis_meV(ee(:)), cube(:), ...
                'VariableNames', {'pixel_y', 'pixel_x', 'signed_energy_meV', 'source_intensity'});
        end

        function out = defaultOutputDir(app, stage2)
            if stage2, tag = 'stage2'; else, tag = 'stage1'; end
            out = fullfile(app.RootDir, 'test_outputs', ['zest_', tag, '_', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
        end

        function [pos, mode] = responsiveFigurePosition(~, requested_mode)
            try
                ss = get(groot, 'ScreenSize');
            catch
                ss = [1 1 1500 900];
            end
            sw = max(900, ss(3));
            sh = max(700, ss(4));
            mode = lower(char(requested_mode));
            if isempty(mode) || strcmp(mode, 'auto')
                if sw < 1450 || sh < 820
                    mode = 'compact';
                else
                    mode = 'wide';
                end
            end
            if strcmp(mode, 'compact')
                w = min(sw - 40, 1180);
                h = min(sh - 70, 820);
                w = max(w, min(sw - 20, 960));
                h = max(h, min(sh - 40, 700));
            else
                w = min(sw - 80, 1500);
                h = min(sh - 90, 900);
                w = max(w, min(sw - 20, 1200));
                h = max(h, min(sh - 40, 760));
                mode = 'wide';
            end
            x = max(20, ss(1) + (sw - w) / 2);
            y = max(35, ss(2) + (sh - h) / 2);
            pos = [x y w h];
        end

        function panel = controlPanel(~, parent, title_text)
            panel = uipanel(parent, 'Title', title_text, 'FontWeight', 'bold', ...
                'ForegroundColor', [0.03 0.22 0.38], 'BackgroundColor', [0.84 0.92 0.99], ...
                'BorderType', 'line', 'HighlightColor', [0.45 0.72 0.95]);
        end

        function enableScrollable(~, container)
            try
                container.Scrollable = 'on';
            catch
            end
        end

        function closeProgressDialog(~, dlg)
            try
                if ~isempty(dlg) && isvalid(dlg)
                    close(dlg);
                end
            catch
                try
                    if ~isempty(dlg) && isvalid(dlg)
                        delete(dlg);
                    end
                catch
                end
            end
        end

        function b = techButton(~, parent, text, callback, row, col)
            b = uibutton(parent, 'Text', text, 'ButtonPushedFcn', callback, ...
                'BackgroundColor', [0.02 0.34 0.62], 'FontColor', [1 1 1], ...
                'FontWeight', 'bold');
            b.Layout.Row = row; b.Layout.Column = col;
        end

        function lab = techLabel(~, parent, text, row, col)
            lab = uilabel(parent, 'Text', text, 'FontColor', [0.04 0.20 0.36], 'FontWeight', 'bold');
            lab.Layout.Row = row; lab.Layout.Column = col;
        end

        function styleAxes(~, ax)
            ax.Color = [0.985 0.995 1.0];
            ax.XGrid = 'on';
            ax.YGrid = 'on';
            ax.GridColor = [0.72 0.82 0.92];
            ax.FontName = 'Arial';
            ax.FontSize = 11;
        end

        function y = normalizeStackSpectrum(~, spectrum)
            y = double(spectrum(:));
            n = numel(y);
            edge_n = max(3, round(0.08 * n));
            edge_idx = unique([1:edge_n, max(1, n - edge_n + 1):n]);
            base = median(y(edge_idx), 'omitnan');
            if ~isfinite(base), base = median(y, 'omitnan'); end
            if ~isfinite(base), base = 0; end
            y = y - base;
            area = sum(abs(y), 'omitnan');
            if ~isfinite(area) || area <= eps
                area = max(abs(y), [], 'omitnan');
            end
            if ~isfinite(area) || area <= eps, area = 1; end
            y = y ./ area;
        end

        function txt = valueToText(~, val)
            if islogical(val) || isnumeric(val)
                txt = mat2str(val, 17);
            elseif ischar(val)
                txt = val;
            elseif isstring(val)
                txt = char(val);
            else
                txt = '<complex>';
            end
        end

        function val = textToValue(~, txt, old)
            if isstring(txt), txt = char(txt); end
            if islogical(old)
                val = any(strcmpi(strtrim(txt), {'true', '1', 'yes', 'on'}));
            elseif isnumeric(old)
                parsed = str2num(txt); %#ok<ST2NM>
                if isempty(parsed)
                    val = old;
                else
                    val = parsed;
                end
            elseif isstring(old)
                val = string(txt);
            elseif ischar(old)
                val = char(txt);
            else
                val = old;
            end
        end

        function name = settingDisplayName(~, key)
            name = lower(strrep(char(key), '_', ' '));
            pairs = {'stage1','Stage 1'; 'stage2','Stage 2'; 'zlp','ZLP'; ...
                'datacube','datacube'; 'spatial','spatial'; 'neighbor','neighbor'; ...
                'attention','attention'; 'interface','interface'; 'temperature','temperature'; ...
                'highfreq','high-frequency'; 'lowfreq','low-frequency'; 'em','EM'; ...
                'pm','source intensity'; 'be','BE'; 'r2','R2'; 'nrmse','NRMSE'};
            for ii = 1:size(pairs, 1)
                name = strrep(name, pairs{ii, 1}, pairs{ii, 2});
            end
        end

        function tip = settingAdvice(~, key)
            key = char(key);
            if contains(key, 'parallel')
                tip = 'Use parallel for full datacubes; 2-4 workers on modest PCs, 8+ only if memory is comfortable.';
            elseif contains(key, 'checkpoint') || contains(key, 'resume')
                tip = 'Checkpoint saving is off by default. Enable it when needed; resume only with identical data, settings and output folder.';
            elseif contains(key, 'crop_range')
                tip = 'Choose the physically useful energy range. Narrower crops are faster but may miss high-frequency features.';
            elseif contains(key, 'energy_bin')
                tip = 'Use 2 for speed/noise reduction; use 1 when fine spectral structure matters.';
            elseif contains(key, 'zlp_component_width') || contains(key, 'zlp_signal_guard')
                tip = 'Increase guards or lower broad-width limits if the ZLP swallows low-frequency signal.';
            elseif contains(key, 'stage1_auto_sparse_pm') || contains(key, 'stage1_fixed_EM')
                tip = 'Stage 1 sparse source-intensity calibration: estimate a compact mode-energy grid from one representative pixel, then reuse it for the cube.';
            elseif contains(key, 'sigmoid') || contains(key, 'preserve_stage1_outside_attention')
                tip = 'Use sigmoid/tapered transitions so attention-window edges do not create artificial bad pixels.';
            elseif contains(key, 'prune')
                tip = 'Stage 2 default: paired attention selects hard-BE mode energies; all others retain independent soft source intensities.';
            elseif contains(key, 'em_') || contains(key, 'dense')
                tip = 'Controls EM-grid density. More centers improve smoothness but slow the fit.';
            elseif contains(key, 'smooth') || contains(key, 'sparse') || contains(key, 'lambda')
                tip = 'Increase smoothing for spiky source intensities; reduce sparsity/smoothing if real peaks are missed.';
            elseif contains(key, 'attention')
                tip = 'Controls feature-window selection. Lower thresholds or raise gain when obvious peaks are missed.';
            elseif contains(key, 'temperature_algorithm')
                tip = 'Attention hard-BE pairs provide temperature evidence; ZLP-guided inference never smooths data-valid temperatures.';
            elseif contains(key, 'hierarchical') || contains(key, 'spatial_zlp')
                tip = 'Hierarchical Stage 1 shares ZLP shape within spectral clusters while preserving per-pixel amplitudes and phonon source intensities.';
            elseif contains(key, 'spatial_refit')
                tip = 'Use for isolated bad pixels. Keep edge-aware behavior on for interfaces.';
            elseif contains(key, 'reliability') || contains(key, 'uncertainty')
                tip = 'Controls confidence diagnostics and outlier repair.';
            elseif contains(key, 'tail') || contains(key, 'overshoot')
                tip = 'Tune when the high-frequency tail is too high or too low.';
            else
                tip = 'Keep default unless diagnostics point to this parameter.';
            end
        end

        function map = optionalMap(~, s, field_name, template)
            if isstruct(s) && isfield(s, field_name) && isequal(size(s.(field_name)), size(template))
                map = s.(field_name);
            else
                map = NaN(size(template));
            end
        end

        function deleteGraphicsHandles(app, prop_name)
            handles = app.(prop_name);
            if isempty(handles), return; end
            for ii = 1:numel(handles)
                try
                    if iscell(handles), h = handles{ii}; else, h = handles(ii); end
                    if isgraphics(h)
                        delete(h);
                    end
                catch
                end
            end
            app.(prop_name) = {};
        end

        function context = captureFitContext(app, requested_map_size)
            context = struct();
            if nargin >= 2 && ~isempty(requested_map_size)
                map_size = double(requested_map_size(:).');
            else
                stage_cube = app.Stage2CubeResult;
                if isempty(stage_cube), stage_cube = app.Stage1CubeResult; end
                map_size = app.cuberesMapSize(stage_cube);
                if isempty(map_size)
                    if isempty(app.RawCube), return; end
                    map_size = size(app.RawCube, [1 2]);
                end
            end
            if numel(map_size) ~= 2 || any(~isfinite(map_size)) || any(map_size < 1)
                context = struct();
                return;
            end
            map_size = round(map_size);
            context.map_size = map_size;
            context.result_file = app.LoadedResultFile;
            stage1_cube = app.Stage1CubeResult;
            stage2_cube = app.Stage2CubeResult;
            if isempty(stage1_cube) && isstruct(stage2_cube)
                stage1_cube = stage2_cube;
            end
            if ~isempty(stage1_cube) || ~isempty(stage2_cube)
                try
                    [T1, T2, valid2] = EELS_ZEST_DataOps.temperature_maps( ...
                        stage1_cube, stage2_cube, map_size);
                    context.T_stage1_soft = T1;
                    context.T_stage2_fit = T2;
                    context.stage2_valid_mask = valid2;
                catch
                end
            else
                s1 = app.Stage1SingleResult;
                s2 = app.Stage2SingleResult;
                if isempty(s1) && isstruct(app.CurrentResult), s1 = app.CurrentResult; end
                if isempty(s2) && app.resultHasStage2(app.CurrentResult), s2 = app.CurrentResult; end
                T1 = NaN(map_size); T2 = NaN(map_size);
                if prod(map_size) == 1
                    if isstruct(s1) && isfield(s1, 'T_soft') && isfinite(s1.T_soft)
                        T1(1) = s1.T_soft;
                    end
                    if isstruct(s2) && isfield(s2, 'stage2') && isstruct(s2.stage2) && ...
                            isfield(s2.stage2, 'T_fit') && isfinite(s2.stage2.T_fit)
                        T2(1) = s2.stage2.T_fit;
                    elseif isstruct(s2) && isfield(s2, 'T_fit') && isfinite(s2.T_fit)
                        T2(1) = s2.T_fit;
                    end
                    context.T_stage1_soft = T1;
                    context.T_stage2_fit = T2;
                    context.stage2_valid_mask = isfinite(T2) & T2 > 0;
                end
            end
            if isfield(context, 'T_stage1_soft') || isfield(context, 'T_stage2_fit')
                app.FitContext = context;
            end
        end

        function tf = hasRestorableRLSource(app)
            tf = false;
            if isstruct(app.DataInfo) && isfield(app.DataInfo, 'rl_conv_applied') && ...
                    app.DataInfo.rl_conv_applied
                return;
            end
            if isfield(app.FitContext, 'map_size') && prod(app.FitContext.map_size) <= 1
                return;
            end
            filename = app.LoadedResultFile;
            if isempty(filename) && isfield(app.FitContext, 'result_file')
                filename = app.FitContext.result_file;
            end
            if ischar(filename) || isstring(filename)
                tf = isfile(char(filename));
            end
        end

        function tf = hasDataRLKernelContext(app)
            tf = false;
            if isempty(app.RawCube) || isempty(app.Ene) || ~isstruct(app.DataInfo) || ...
                    ~isfield(app.DataInfo, 'fitted_zlp_kernel_cube')
                return;
            end
            if isfield(app.DataInfo, 'rl_conv_applied') && app.DataInfo.rl_conv_applied
                return;
            end
            kernels = app.DataInfo.fitted_zlp_kernel_cube;
            if ~isequal(size(kernels), size(app.RawCube))
                return;
            end
            if isfield(app.DataInfo, 'fitted_zlp_kernel_energy_meV')
                kernel_energy = double(app.DataInfo.fitted_zlp_kernel_energy_meV(:));
                if numel(kernel_energy) ~= numel(app.Ene) || ...
                        any(abs(kernel_energy - app.Ene(:)) > 1e-9)
                    return;
                end
            end
            tf = any(isfinite(kernels(:)) & kernels(:) > 0);
        end

        function res = currentDataRLResult(app, yy, xx)
            res = [];
            if ~app.hasDataRLKernelContext() || yy < 1 || xx < 1 || ...
                    yy > size(app.RawCube, 1) || xx > size(app.RawCube, 2)
                return;
            end
            signal = double(squeeze(app.RawCube(yy, xx, :)));
            kernel = double(squeeze(app.DataInfo.fitted_zlp_kernel_cube(yy, xx, :)));
            if numel(signal) ~= numel(app.Ene) || numel(kernel) ~= numel(app.Ene) || ...
                    ~any(isfinite(signal)) || ~any(isfinite(kernel) & kernel > 0)
                return;
            end
            kernel(~isfinite(kernel) | kernel < 0) = 0;
            kernel = kernel ./ max(sum(kernel, 'omitnan'), eps);
            norm_val = sum(abs(signal), 'omitnan');
            if ~isfinite(norm_val) || norm_val <= eps
                norm_val = max(abs(signal), [], 'omitnan');
            end
            if ~isfinite(norm_val) || norm_val <= eps
                return;
            end
            res = struct('ene_corrected', app.Ene(:), 'ene_raw', app.Ene(:), ...
                'raw_data', signal(:), 'raw_data_norm', signal(:) ./ norm_val, ...
                'norm_val', norm_val, 'rl_input_is_background_subtracted', true, ...
                'rl_kernel_norm', kernel(:));
        end

        function restored = restoreRLSourceFromContext(app)
            restored = false;
            if isstruct(app.DataInfo) && isfield(app.DataInfo, 'rl_conv_applied') && ...
                    app.DataInfo.rl_conv_applied
                return;
            end
            if ~isempty(app.activeRLDeconvCubeResult()) || app.resultCanRunRLDeconv(app.RLSourceResult)
                restored = true;
                return;
            end
            filename = app.LoadedResultFile;
            if isempty(filename) && isfield(app.FitContext, 'result_file')
                filename = app.FitContext.result_file;
            end
            if ~(ischar(filename) || isstring(filename)) || ~isfile(char(filename))
                return;
            end
            filename = char(filename);
            try
                listed = whos('-file', filename);
                names = {listed.name};
                candidates = intersect({'cuberes', 'stage2_cuberes', 'stage1_cuberes', 'result'}, ...
                    names, 'stable');
                if isempty(candidates), return; end
                S = load(filename, candidates{:});
                if isfield(S, 'cuberes') || isfield(S, 'stage2_cuberes') || isfield(S, 'stage1_cuberes')
                    [cuberes, ~] = app.extractLoadedCubeResult(S);
                    if ~isempty(app.firstRLDeconvResult(cuberes))
                        app.RLSourceCubeResult = cuberes;
                        app.RLSourceResult = app.firstRLDeconvResult(cuberes);
                        restored = true;
                    end
                elseif isfield(S, 'result') && app.resultCanRunRLDeconv(S.result)
                    app.RLSourceResult = S.result;
                    restored = true;
                end
                if restored
                    app.LoadedResultFile = filename;
                    app.logStatus('Restored the fitted-result context required by RL deconvolution.');
                end
            catch ME
                app.logStatus(['Could not restore the RL source result: ', ME.message]);
            end
        end

        function cutoff = cutoffForResult(~, res, options)
            cutoff = NaN;
            if contains(lower(string(options.mode)), 'attention')
                cutoff = EELS_ZEST_DataOps.attention_cutoff(res, options.tail_meV);
            elseif isfield(options, 'cutoff_meV') && isfinite(options.cutoff_meV)
                cutoff = options.cutoff_meV;
            end
        end

        function options = normalizeCutoffOptions(~, options, stage)
            assert(isstruct(options), 'EELS:InvalidCutoffOptions', ...
                'Programmatic cutoff options must be a struct.');
            if ~isfield(options, 'cancelled'), options.cancelled = false; end
            if ~isfield(options, 'mode'), options.mode = 'Manual cutoff'; end
            if ~isfield(options, 'cutoff_meV'), options.cutoff_meV = NaN; end
            if ~isfield(options, 'tail_meV'), options.tail_meV = 20; end
            if ~isfield(options, 'json_file'), options.json_file = ''; end
            if ~isfield(options, 'json_meta'), options.json_meta = struct(); end
            if contains(lower(string(options.mode)), 'attention')
                assert(strcmpi(stage, 'stage2'), 'EELS:InvalidCutoffOptions', ...
                    'Attention cutoff is available only for Stage 2.');
            else
                assert(isfinite(options.cutoff_meV) && options.cutoff_meV > 0, ...
                    'EELS:InvalidCutoffOptions', 'Manual/JSON cutoff must be positive.');
            end
        end

        function registerOwnedFigure(app, fig)
            if isempty(fig) || ~isgraphics(fig), return; end
            keep = false(size(app.OwnedFigures));
            for ii = 1:numel(app.OwnedFigures)
                keep(ii) = isgraphics(app.OwnedFigures{ii});
            end
            app.OwnedFigures = app.OwnedFigures(keep);
            app.OwnedFigures{end + 1} = fig;
        end

        function filename = findSavedResultFile(~, output_dir, stage2)
            filename = '';
            if nargin < 3, stage2 = true; end
            if isempty(output_dir) || ~isfolder(output_dir), return; end
            patterns = {'stage2_reused_datacube_result.mat', '*datacube_result.mat'};
            if ~stage2, patterns = {'stage1_datacube_result.mat', '*datacube_result.mat'}; end
            files = [];
            for ii = 1:numel(patterns)
                files = [files; dir(fullfile(output_dir, patterns{ii}))]; %#ok<AGROW>
            end
            if isempty(files), return; end
            [~, order] = sort([files.datenum], 'descend');
            files = files(order);
            filename = fullfile(files(1).folder, files(1).name);
        end

        function cached = getRLDiagnosticCache(app, yy, xx)
            cached = [];
            key = sprintf('%d_%d', yy, xx);
            idx = find(strcmp(app.RLDiagnosticCache.keys, key), 1);
            if isempty(idx), return; end
            cached = app.RLDiagnosticCache.values{idx};
            if idx > 1
                app.RLDiagnosticCache.keys = [{key}, app.RLDiagnosticCache.keys([1:idx-1, idx+1:end])];
                app.RLDiagnosticCache.values = [{cached}, app.RLDiagnosticCache.values([1:idx-1, idx+1:end])];
            end
        end

        function putRLDiagnosticCache(app, yy, xx, value)
            key = sprintf('%d_%d', yy, xx);
            idx = find(strcmp(app.RLDiagnosticCache.keys, key), 1);
            if ~isempty(idx)
                app.RLDiagnosticCache.keys(idx) = [];
                app.RLDiagnosticCache.values(idx) = [];
            end
            app.RLDiagnosticCache.keys = [{key}, app.RLDiagnosticCache.keys];
            app.RLDiagnosticCache.values = [{value}, app.RLDiagnosticCache.values];
            cap = max(1, app.RLDiagnosticCache.capacity);
            app.RLDiagnosticCache.keys = app.RLDiagnosticCache.keys(1:min(cap, end));
            app.RLDiagnosticCache.values = app.RLDiagnosticCache.values(1:min(cap, end));
        end

        function appendGraphicsHandle(app, prop_name, handle)
            handles = app.(prop_name);
            if ~iscell(handles), handles = num2cell(handles); end
            handles{end + 1} = handle;
            app.(prop_name) = handles;
        end

        function [pos, cancelled] = waitForROIWithHotkeys(app, roi, prompt)
            pos = [];
            cancelled = false;
            fig = app.UIFigure;
            token = ['EELS_ZEST_ROI_', char(java.util.UUID.randomUUID())];
            token = matlab.lang.makeValidName(strrep(token, '-', '_'));
            status_key = [token, '_status'];
            pos_key = [token, '_position'];
            old_key_fcn = [];
            try
                old_key_fcn = fig.WindowKeyPressFcn;
            catch
            end
            setappdata(fig, status_key, 'waiting');
            setappdata(fig, pos_key, []);
            try
                fig.WindowKeyPressFcn = @(src, evt)app.roiHotkey(src, evt, roi, status_key, pos_key);
            catch
            end
            app.logStatus(prompt);
            try
                wait(roi);
            catch
            end
            try
                forced_pos = getappdata(fig, pos_key);
                status = getappdata(fig, status_key);
            catch
                forced_pos = [];
                status = 'waiting';
            end
            if strcmp(status, 'cancelled')
                cancelled = true;
            elseif ~isempty(forced_pos)
                pos = forced_pos;
            else
                try
                    if isvalid(roi)
                        pos = roi.Position;
                    end
                catch
                    pos = [];
                end
            end
            try
                if isvalid(roi)
                    delete(roi);
                end
            catch
            end
            try
                fig.WindowKeyPressFcn = old_key_fcn;
            catch
            end
            try
                rmappdata(fig, status_key);
                rmappdata(fig, pos_key);
            catch
            end
            if cancelled
                pos = [];
            end
        end

        function roiHotkey(~, src, evt, roi, status_key, pos_key)
            key = lower(string(evt.Key));
            mods = lower(string(evt.Modifier));
            has_ctrl = any(mods == "control" | mods == "ctrl");
            if ~has_ctrl, return; end
            if key == "q" || key == "h"
                pos = [];
                try
                    if isvalid(roi)
                        pos = roi.Position;
                    end
                catch
                end
                if key == "q"
                    setappdata(src, status_key, 'cancelled');
                else
                    setappdata(src, status_key, 'confirmed');
                    setappdata(src, pos_key, pos);
                end
                try
                    if isvalid(roi)
                        resume(roi);
                    end
                catch
                end
            end
        end

        function logStage(app, stage2, msg)
            app.logStatus(msg);
            if stage2 && ~isempty(app.Stage2StatusArea) && isvalid(app.Stage2StatusArea)
                values = [cellstr(app.Stage2StatusArea.Value); {msg}];
                app.Stage2StatusArea.Value = values(max(1, end - 249):end);
            elseif ~stage2 && ~isempty(app.Stage1StatusArea) && isvalid(app.Stage1StatusArea)
                values = [cellstr(app.Stage1StatusArea.Value); {msg}];
                app.Stage1StatusArea.Value = values(max(1, end - 249):end);
            end
        end

        function logStatus(app, msg)
            stamp = char(datetime('now', 'Format', 'HH:mm:ss'));
            line = sprintf('[%s] %s', stamp, msg);
            if isempty(app.StatusArea) || ~isvalid(app.StatusArea)
                fprintf('%s\n', line);
            else
                values = [cellstr(app.StatusArea.Value); {line}];
                app.StatusArea.Value = values(max(1, end - 399):end);
                stamp_now = now;
                if app.LastLogDrawTime == 0 || (stamp_now - app.LastLogDrawTime) .* 86400 >= 0.15
                    drawnow limitrate;
                    app.LastLogDrawTime = stamp_now;
                end
            end
        end
    end
end
