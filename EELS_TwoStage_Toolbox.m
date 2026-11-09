classdef EELS_TwoStage_Toolbox < handle
    % GUI wrapper for the Pearson sparse full-fit workflow.
    %
    % Usage:
    %   addpath('pearson_sparse_fullfit');
    %   app = EELS_TwoStage_Toolbox;

    properties
        UIFigure
        TabGroup
        DataTab
        ReferenceTab
        ReferenceTuner
        Stage1Tab
        Stage2Tab
        ViewTab
        AdvancedTab

        DataAxes
        SpectrumAxes
        Stage1MapAxes
        Stage1SpectrumAxes
        Stage2MapAxes
        Stage2SpectrumAxes
        ViewAxes1
        ViewAxes2

        StatusArea
        Stage1StatusArea
        Stage2StatusArea
        AdvancedHelpArea

        FileEdit
        OutputDirEdit
        EnergyMinEdit
        EnergyMaxEdit
        PixelYSpinner
        PixelXSpinner
        Stage1PixelYSpinner
        Stage1PixelXSpinner
        Stage2PixelYSpinner
        Stage2PixelXSpinner
        SpatialYBinSpinner
        SpatialXBinSpinner
        EnergyBinSpinner
        WorkerSpinner
        ParallelCheckBox
        CheckpointCheckBox
        ResumeCheckBox
        SavePlotsCheckBox
        MapDropDown
        Stage2MapDropDown
        Stage2AlgorithmDropDown
        DiagnosticDropDown
        AlignMethodDropDown
        BM3DFactorEdit
        BGMethodDropDown
        BGWin1MinEdit
        BGWin1MaxEdit
        BGWin2MinEdit
        BGWin2MaxEdit
        DataIOActionDropDown
        DataToolDropDown
        Stage1Controls = struct()
        Stage2Controls = struct()

        AdvancedTable

        RootDir
        LegacyToolboxDir
        LayoutMode = 'auto'
        IsCompact = false
        DataFile = ''
        Ene = []
        Data = []
        DataLabel = ''
        Settings
        CurrentResult = []
        CubeResult = []
        Stage1CubeResult = []
        Stage2CubeResult = []
        Stage1SingleResult = []
        Stage2SingleResult = []
        LastOutputDir = ''
        SelectedROICube = []
        SelectedROISpectrum = []
        SelectedROIRange = []
        StackSpectra = []
        DataCursorHandles = []
        Stage1CursorHandles = []
        Stage2CursorHandles = []
        StackCursorHandles = []
        DatacubePauseRequested = false
        ReferenceTunerFigures = {}
    end

    methods
        function app = EELS_TwoStage_Toolbox(varargin)
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
            app.LegacyToolboxDir = fullfile(app.RootDir, 'external_legacy_tools');
            addpath(app.RootDir);
            if isfolder(app.LegacyToolboxDir)
                addpath(genpath(app.LegacyToolboxDir));
            end
            app.Settings = EELS_PearsonSparse_FullFit.defaults(struct());
            app.Settings.datacube_sample_strategy = 'all';
            app.Settings.datacube_sample_pixels = Inf;
            app.LastOutputDir = fullfile(app.RootDir, 'test_outputs', ...
                ['gui_run_', char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))]);
            app.buildUI(visible);
            app.populateAdvancedTable();
            app.logStatus('GUI ready. Load a datacube or result. Stage1 fits sparse spectra; Stage2 estimates temperature and diagnostics.');
        end

        function delete(app)
            for k=1:numel(app.ReferenceTunerFigures)
                if isgraphics(app.ReferenceTunerFigures{k}),delete(app.ReferenceTunerFigures{k});end
            end
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end
    end

    methods (Access = private)
        function openReferenceTuner(app)
            app.TabGroup.SelectedTab=app.ReferenceTab;
        end

        function applyTunedSettings(app,settings)
            app.Settings=settings;
            app.clearFitCaches();
            for stage2=[false true]
                c=app.controlsForStage(stage2);c.resume.Value=false;
                c.cropMin.Value=min(settings.crop_range);c.cropMax.Value=max(settings.crop_range);
                if stage2&&isfield(c,'hierarchical'),c.hierarchical.Value=settings.stage1_hierarchical_enabled;end
            end
            app.Stage2AlgorithmDropDown.Value=settings.stage2_temperature_algorithm;
            app.populateAdvancedTable();
            app.logStatus('Reference strategy applied. Validate independent spectra before full-map processing.');
        end

        function s=referenceBaseSettings(app)
            app.applyAdvancedTable();s=app.Settings;s.crop_range=app.stageEnergyWindow(false);
        end

        function r=referenceInput(app)
            if ~isempty(app.Data)
                y=round(app.PixelYSpinner.Value);x=round(app.PixelXSpinner.Value);
                r=struct('ene',app.Ene,'spectrum',squeeze(app.Data(y,x,:)), ...
                    'source',sprintf('%s | Y=%d X=%d',app.DataFile,y,x),'group',app.DataFile);
                r.cube=app.Data;
            elseif ~isempty(app.CurrentResult)
                r=struct('ene',app.CurrentResult.ene_raw,'spectrum',app.CurrentResult.raw_data,'source','Current fitted spectrum');
            else,error('EELS:Reference','Import a reference spectrum first.');end
        end

        function closeWithTuner(app)
            if ~isempty(app.ReferenceTuner)&&app.ReferenceTuner.Running
                app.ReferenceTuner.stop();return;
            end
            delete(app);
        end

        function buildUI(app, visible)
            [fig_pos, resolved_mode] = app.responsiveFigurePosition(app.LayoutMode);
            app.LayoutMode = resolved_mode;
            app.IsCompact = strcmp(resolved_mode, 'compact');
            app.UIFigure = uifigure('Name', 'EELS-ZEST 2.2 Two-Stage Toolbox', ...
                'Position', fig_pos, 'Color', [0.94 0.97 1.00], 'Visible', visible);
            app.UIFigure.CloseRequestFcn=@(~,~)app.closeWithTuner();
            rootLayout = uigridlayout(app.UIFigure, [1 1]);
            rootLayout.Padding = [0 0 0 0];
            rootLayout.RowSpacing = 0;
            rootLayout.ColumnSpacing = 0;
            rootLayout.BackgroundColor = [0.94 0.97 1.00];
            app.TabGroup = uitabgroup(rootLayout);
            app.TabGroup.Layout.Row = 1;
            app.TabGroup.Layout.Column = 1;
            if strcmpi(visible, 'on') && ~app.IsCompact
                try
                    app.UIFigure.WindowState = 'maximized';
                catch
                end
            end
            app.DataTab = uitab(app.TabGroup, 'Title', 'Data and preprocessing', 'BackgroundColor', [0.94 0.97 1.00]);
            app.ReferenceTab = uitab(app.TabGroup, 'Title', 'Reference tuning');
            app.Stage1Tab = uitab(app.TabGroup, 'Title', 'Stage1: Sparse fit', 'BackgroundColor', [0.94 0.97 1.00]);
            app.Stage2Tab = uitab(app.TabGroup, 'Title', 'Stage2: Temperature', 'BackgroundColor', [0.94 0.97 1.00]);
            app.ViewTab = uitab(app.TabGroup, 'Title', 'Diagnostics', 'BackgroundColor', [0.94 0.97 1.00]);
            app.AdvancedTab = uitab(app.TabGroup, 'Title', 'Advanced', 'BackgroundColor', [0.94 0.97 1.00]);
            app.buildDataTab();
            app.buildStage1Tab();
            app.buildStage2Tab();
            app.buildViewTab();
            app.buildAdvancedTab();
            app.ReferenceTuner=EELS_ReferenceTunerApp(app.Settings,[],@(s)app.applyTunedSettings(s),visible, ...
                app.ReferenceTab,@()app.resetAdvanced(),@()app.referenceInput(),@()app.referenceBaseSettings());
        end

        function buildDataTab(app)
            if app.IsCompact
                app.buildDataTabCompact();
                return;
            end
            gl = uigridlayout(app.DataTab, [3 4]);
            gl.RowHeight = {210, '1x', 135};
            gl.ColumnWidth = {360, '1x', '1x', 360};
            gl.Padding = [14 14 14 14];
            gl.BackgroundColor = [0.94 0.97 1.00];

            p = uipanel(gl, 'Title', 'Data / Preprocessing / ROI / Background', ...
                'FontWeight', 'bold', 'ForegroundColor', [0.05 0.18 0.32], ...
                'BackgroundColor', [0.88 0.94 1.00], 'BorderType', 'none');
            p.Layout.Row = 1; p.Layout.Column = [1 4];
            app.enableScrollable(p);
            g = uigridlayout(p, [4 12]);
            g.RowHeight = {30, 30, 30, 30};
            g.ColumnWidth = {92, 96, 96, 96, 78, 78, 78, 78, 96, 110, 110, 150};
            g.Padding = [10 8 10 8];
            g.BackgroundColor = [0.88 0.94 1.00];

            b = app.techButton(g, 'Load MAT', @(~, ~)app.onLoadData()); b.Layout.Row = 1; b.Layout.Column = 1;
            b = app.techButton(g, 'Load result', @(~, ~)app.onLoadResult()); b.Layout.Row = 1; b.Layout.Column = 2;
            b = app.techButton(g, 'Align', @(~, ~)app.onAlign()); b.Layout.Row = 1; b.Layout.Column = 3;
            b = app.techButton(g, 'BM3D denoise', @(~, ~)app.onBM3D()); b.Layout.Row = 1; b.Layout.Column = 4;
            lab = app.techLabel(g, 'Align method'); lab.Layout.Row = 1; lab.Layout.Column = 5;
            app.AlignMethodDropDown = uidropdown(g, 'Items', {'xcorr', 'gaussian', 'COM'}, 'Value', 'xcorr');
            app.AlignMethodDropDown.Layout.Row = 1; app.AlignMethodDropDown.Layout.Column = 6;
            lab = app.techLabel(g, 'BM factor'); lab.Layout.Row = 1; lab.Layout.Column = 7;
            app.BM3DFactorEdit = uieditfield(g, 'numeric', 'Value', 1, 'Limits', [0.1 Inf]);
            app.BM3DFactorEdit.Layout.Row = 1; app.BM3DFactorEdit.Layout.Column = 8;
            b = app.techButton(g, 'Energy map', @(~, ~)app.plotEnergyMap()); b.Layout.Row = 1; b.Layout.Column = 9;
            b = app.techButton(g, 'LineProfile', @(~, ~)app.plotLineProfile()); b.Layout.Row = 1; b.Layout.Column = 10;
            b = app.techButton(g, 'Export data', @(~, ~)app.exportData()); b.Layout.Row = 1; b.Layout.Column = 11;
            app.FileEdit = uieditfield(g, 'text', 'Editable', 'off');
            app.FileEdit.Layout.Row = 2; app.FileEdit.Layout.Column = [1 8];
            lab = app.techLabel(g, 'Energy window'); lab.Layout.Row = 2; lab.Layout.Column = 9;
            app.EnergyMinEdit = uieditfield(g, 'numeric', 'Value', -250);
            app.EnergyMinEdit.Layout.Row = 2; app.EnergyMinEdit.Layout.Column = 10;
            app.EnergyMaxEdit = uieditfield(g, 'numeric', 'Value', 250);
            app.EnergyMaxEdit.Layout.Row = 2; app.EnergyMaxEdit.Layout.Column = 11;
            b = app.techButton(g, 'Refresh', @(~, ~)app.refreshDataPreview()); b.Layout.Row = 2; b.Layout.Column = 12;
            b = app.techButton(g, 'Integrate ROI', @(~, ~)app.selectROIAndIntegrate()); b.Layout.Row = 3; b.Layout.Column = [1 2];
            b = app.techButton(g, 'Export ROI', @(~, ~)app.exportROI()); b.Layout.Row = 3; b.Layout.Column = 3;
            b = app.techButton(g, 'Import fit signal', @(~, ~)app.importFitSignalToData()); b.Layout.Row = 3; b.Layout.Column = [4 5];

            lab = app.techLabel(g, 'Background'); lab.Layout.Row = 4; lab.Layout.Column = 1;
            app.BGMethodDropDown = uidropdown(g, 'Items', {'Pearson', 'ExpPoly3', 'ExpPoly2', 'ExpPoly4', ...
                'Power', 'Exponential_2', 'PseudoVoigt'}, 'Value', 'ExpPoly3');
            app.BGMethodDropDown.Layout.Row = 4; app.BGMethodDropDown.Layout.Column = [2 3];
            app.BGWin1MinEdit = uieditfield(g, 'numeric', 'Value', -120);
            app.BGWin1MinEdit.Layout.Row = 4; app.BGWin1MinEdit.Layout.Column = 4;
            app.BGWin1MaxEdit = uieditfield(g, 'numeric', 'Value', -90);
            app.BGWin1MaxEdit.Layout.Row = 4; app.BGWin1MaxEdit.Layout.Column = 5;
            app.BGWin2MinEdit = uieditfield(g, 'numeric', 'Value', 120);
            app.BGWin2MinEdit.Layout.Row = 4; app.BGWin2MinEdit.Layout.Column = 6;
            app.BGWin2MaxEdit = uieditfield(g, 'numeric', 'Value', 180);
            app.BGWin2MaxEdit.Layout.Row = 4; app.BGWin2MaxEdit.Layout.Column = 7;
            app.DataIOActionDropDown = uidropdown(g, 'Items', {'Import MAT', 'Import result', 'Export data', 'Export ROI', 'Export current result'}, 'Value', 'Import MAT');
            app.DataIOActionDropDown.Layout.Row = 3; app.DataIOActionDropDown.Layout.Column = [6 7];
            b = app.techButton(g, 'Run import/export', @(~, ~)app.runDataIOAction()); b.Layout.Row = 3; b.Layout.Column = [8 9];
            app.DataToolDropDown = uidropdown(g, 'Items', {'NMF', 'STA correction', 'Background (interactive)', 'Separate IF mode', 'MLLS residual', 'Peak fitting', 'Drift correction', 'Stack spectra'}, 'Value', 'NMF');
            app.DataToolDropDown.Layout.Row = 3; app.DataToolDropDown.Layout.Column = [10 11];
            app.DataToolDropDown.Items = [app.DataToolDropDown.Items, {'Phonon intensity'}];
            b = app.techButton(g, 'Open tool', @(~, ~)app.runDataToolAction()); b.Layout.Row = 3; b.Layout.Column = 12;
            b = app.techButton(g, 'Subtract background', @(~, ~)app.applyTraditionalBackground()); b.Layout.Row = 4; b.Layout.Column = [8 10];
            b = app.techButton(g, 'Stack spectra', @(~, ~)app.openStackSpectraTool()); b.Layout.Row = 4; b.Layout.Column = [11 12];

            app.DataAxes = uiaxes(gl);
            app.DataAxes.Layout.Row = 2; app.DataAxes.Layout.Column = [1 2];
            app.styleAxes(app.DataAxes);
            title(app.DataAxes, 'Energy-filtered image / ROI map'); axis(app.DataAxes, 'image');
            app.SpectrumAxes = uiaxes(gl);
            app.SpectrumAxes.Layout.Row = 2; app.SpectrumAxes.Layout.Column = [3 4];
            app.styleAxes(app.SpectrumAxes);
            title(app.SpectrumAxes, 'Selected spectrum / ROI sum');

            app.StatusArea = uitextarea(gl, 'Editable', 'off');
            app.StatusArea.Layout.Row = 3; app.StatusArea.Layout.Column = [1 4];
        end

        function buildDataTabCompact(app)
            gl = uigridlayout(app.DataTab, [4 4]);
            gl.RowHeight = {245, '1x', '1x', 105};
            gl.ColumnWidth = {'1x', '1x', '1x', '1x'};
            gl.Padding = [10 8 10 8];
            gl.BackgroundColor = [0.94 0.97 1.00];

            p = uipanel(gl, 'Title', 'Data / Preprocessing / ROI / Tools (compact)', ...
                'FontWeight', 'bold', 'ForegroundColor', [0.05 0.18 0.32], ...
                'BackgroundColor', [0.88 0.94 1.00], 'BorderType', 'none');
            p.Layout.Row = 1; p.Layout.Column = [1 4];
            app.enableScrollable(p);
            g = uigridlayout(p, [6 8]);
            g.RowHeight = {29, 29, 29, 29, 29, 29};
            g.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            g.Padding = [10 8 10 8];
            g.RowSpacing = 7;
            g.ColumnSpacing = 8;
            g.BackgroundColor = [0.88 0.94 1.00];

            b = app.techButton(g, 'Load MAT', @(~, ~)app.onLoadData()); b.Layout.Row = 1; b.Layout.Column = [1 2];
            b = app.techButton(g, 'Load result', @(~, ~)app.onLoadResult()); b.Layout.Row = 1; b.Layout.Column = [3 4];
            b = app.techButton(g, 'Export data', @(~, ~)app.exportData()); b.Layout.Row = 1; b.Layout.Column = [5 6];
            b = app.techButton(g, 'Refresh', @(~, ~)app.refreshDataPreview()); b.Layout.Row = 1; b.Layout.Column = [7 8];

            app.FileEdit = uieditfield(g, 'text', 'Editable', 'off');
            app.FileEdit.Layout.Row = 2; app.FileEdit.Layout.Column = [1 8];

            lab = app.techLabel(g, 'Energy window'); lab.Layout.Row = 3; lab.Layout.Column = 1;
            app.EnergyMinEdit = uieditfield(g, 'numeric', 'Value', -250);
            app.EnergyMinEdit.Layout.Row = 3; app.EnergyMinEdit.Layout.Column = 2;
            app.EnergyMaxEdit = uieditfield(g, 'numeric', 'Value', 250);
            app.EnergyMaxEdit.Layout.Row = 3; app.EnergyMaxEdit.Layout.Column = 3;
            b = app.techButton(g, 'Energy map', @(~, ~)app.plotEnergyMap()); b.Layout.Row = 3; b.Layout.Column = 4;
            b = app.techButton(g, 'LineProfile', @(~, ~)app.plotLineProfile()); b.Layout.Row = 3; b.Layout.Column = [5 6];
            b = app.techButton(g, 'Stack spectra', @(~, ~)app.openStackSpectraTool()); b.Layout.Row = 3; b.Layout.Column = [7 8];

            b = app.techButton(g, 'Integrate ROI', @(~, ~)app.selectROIAndIntegrate()); b.Layout.Row = 4; b.Layout.Column = [1 2];
            b = app.techButton(g, 'Export ROI', @(~, ~)app.exportROI()); b.Layout.Row = 4; b.Layout.Column = [3 4];
            b = app.techButton(g, 'Import fit signal', @(~, ~)app.importFitSignalToData()); b.Layout.Row = 4; b.Layout.Column = [5 6];
            b = app.techButton(g, 'Align', @(~, ~)app.onAlign()); b.Layout.Row = 4; b.Layout.Column = 7;
            b = app.techButton(g, 'BM3D', @(~, ~)app.onBM3D()); b.Layout.Row = 4; b.Layout.Column = 8;

            lab = app.techLabel(g, 'Align method'); lab.Layout.Row = 5; lab.Layout.Column = 1;
            app.AlignMethodDropDown = uidropdown(g, 'Items', {'xcorr', 'gaussian', 'COM'}, 'Value', 'xcorr');
            app.AlignMethodDropDown.Layout.Row = 5; app.AlignMethodDropDown.Layout.Column = 2;
            lab = app.techLabel(g, 'BM factor'); lab.Layout.Row = 5; lab.Layout.Column = 3;
            app.BM3DFactorEdit = uieditfield(g, 'numeric', 'Value', 1, 'Limits', [0.1 Inf]);
            app.BM3DFactorEdit.Layout.Row = 5; app.BM3DFactorEdit.Layout.Column = 4;
            app.DataIOActionDropDown = uidropdown(g, 'Items', {'Import MAT', 'Import result', 'Export data', 'Export ROI', 'Export current result'}, 'Value', 'Import MAT');
            app.DataIOActionDropDown.Layout.Row = 5; app.DataIOActionDropDown.Layout.Column = [5 6];
            b = app.techButton(g, 'Run import/export', @(~, ~)app.runDataIOAction()); b.Layout.Row = 5; b.Layout.Column = [7 8];

            app.DataToolDropDown = uidropdown(g, 'Items', {'NMF', 'STA correction', 'Background (interactive)', 'Separate IF mode', 'MLLS residual', 'Peak fitting', 'Drift correction', 'Stack spectra', 'Phonon intensity'}, 'Value', 'NMF');
            app.DataToolDropDown.Layout.Row = 6; app.DataToolDropDown.Layout.Column = [1 3];
            b = app.techButton(g, 'Open tool', @(~, ~)app.runDataToolAction()); b.Layout.Row = 6; b.Layout.Column = [4 5];
            app.BGMethodDropDown = uidropdown(g, 'Items', {'Pearson', 'ExpPoly3', 'ExpPoly2', 'ExpPoly4', ...
                'Power', 'Exponential_2', 'PseudoVoigt'}, 'Value', 'ExpPoly3');
            app.BGMethodDropDown.Layout.Row = 6; app.BGMethodDropDown.Layout.Column = [6 7];
            b = app.techButton(g, 'Background tool', @(~, ~)app.applyTraditionalBackground()); b.Layout.Row = 6; b.Layout.Column = 8;
            app.BGWin1MinEdit = uieditfield(app.UIFigure, 'numeric', 'Value', -120, 'Visible', 'off');
            app.BGWin1MaxEdit = uieditfield(app.UIFigure, 'numeric', 'Value', -90, 'Visible', 'off');
            app.BGWin2MinEdit = uieditfield(app.UIFigure, 'numeric', 'Value', 120, 'Visible', 'off');
            app.BGWin2MaxEdit = uieditfield(app.UIFigure, 'numeric', 'Value', 180, 'Visible', 'off');

            app.DataAxes = uiaxes(gl);
            app.DataAxes.Layout.Row = 2; app.DataAxes.Layout.Column = [1 4];
            app.styleAxes(app.DataAxes);
            title(app.DataAxes, 'Energy-filtered image / ROI map'); axis(app.DataAxes, 'image');
            app.SpectrumAxes = uiaxes(gl);
            app.SpectrumAxes.Layout.Row = 3; app.SpectrumAxes.Layout.Column = [1 4];
            app.styleAxes(app.SpectrumAxes);
            title(app.SpectrumAxes, 'Selected spectrum / ROI sum');

            app.StatusArea = uitextarea(gl, 'Editable', 'off');
            app.StatusArea.Layout.Row = 4; app.StatusArea.Layout.Column = [1 4];
        end

        function buildStage1Tab(app)
            if app.IsCompact
                gl = uigridlayout(app.Stage1Tab, [4 4]);
                gl.RowHeight = {280, '1x', '1x', 105};
                gl.ColumnWidth = {'1x', '1x', '1x', '1x'};
            else
                gl = uigridlayout(app.Stage1Tab, [3 4]);
                gl.RowHeight = {320, '1x', 128};
                gl.ColumnWidth = {380, '1x', '1x', 360};
            end
            gl.Padding = [14 14 14 14];
            gl.BackgroundColor = [0.94 0.97 1.00];

            p = uipanel(gl, 'Title', 'Stage1: ZLP + sparse convolution', ...
                'FontWeight', 'bold', 'ForegroundColor', [0.05 0.18 0.32], ...
                'BackgroundColor', [0.88 0.94 1.00], 'BorderType', 'none');
            p.Layout.Row = 1; p.Layout.Column = [1 4];
            app.enableScrollable(p);
            app.buildRunControls(p, false);

            app.Stage1MapAxes = uiaxes(gl);
            if app.IsCompact
                app.Stage1MapAxes.Layout.Row = 2; app.Stage1MapAxes.Layout.Column = [1 4];
            else
                app.Stage1MapAxes.Layout.Row = 2; app.Stage1MapAxes.Layout.Column = [1 2];
            end
            app.styleAxes(app.Stage1MapAxes);
            title(app.Stage1MapAxes, 'Stage1 / datacube map'); axis(app.Stage1MapAxes, 'image');
            app.Stage1SpectrumAxes = uiaxes(gl);
            if app.IsCompact
                app.Stage1SpectrumAxes.Layout.Row = 3; app.Stage1SpectrumAxes.Layout.Column = [1 4];
            else
                app.Stage1SpectrumAxes.Layout.Row = 2; app.Stage1SpectrumAxes.Layout.Column = [3 4];
            end
            app.styleAxes(app.Stage1SpectrumAxes);
            title(app.Stage1SpectrumAxes, 'Current Stage1 spectrum');

            app.Stage1StatusArea = uitextarea(gl, 'Editable', 'off');
            if app.IsCompact
                app.Stage1StatusArea.Layout.Row = 4;
            else
                app.Stage1StatusArea.Layout.Row = 3;
            end
            app.Stage1StatusArea.Layout.Column = [1 4];
        end

        function buildStage2Tab(app)
            if app.IsCompact
                gl = uigridlayout(app.Stage2Tab, [4 4]);
                gl.RowHeight = {280, '1x', '1x', 105};
                gl.ColumnWidth = {'1x', '1x', '1x', '1x'};
            else
                gl = uigridlayout(app.Stage2Tab, [3 4]);
                gl.RowHeight = {320, '1x', 128};
                gl.ColumnWidth = {380, '1x', '1x', 360};
            end
            gl.Padding = [14 14 14 14];
            gl.BackgroundColor = [0.94 0.97 1.00];

            p = uipanel(gl, 'Title', 'Stage2: Attention window + BE thermometry', ...
                'FontWeight', 'bold', 'ForegroundColor', [0.05 0.18 0.32], ...
                'BackgroundColor', [0.88 0.94 1.00], 'BorderType', 'none');
            p.Layout.Row = 1; p.Layout.Column = [1 4];
            app.enableScrollable(p);
            app.buildRunControls(p, true);

            app.Stage2MapAxes = uiaxes(gl);
            if app.IsCompact
                app.Stage2MapAxes.Layout.Row = 2; app.Stage2MapAxes.Layout.Column = [1 4];
            else
                app.Stage2MapAxes.Layout.Row = 2; app.Stage2MapAxes.Layout.Column = [1 2];
            end
            app.styleAxes(app.Stage2MapAxes);
            title(app.Stage2MapAxes, 'Temperature / reliability / interface map'); axis(app.Stage2MapAxes, 'image');
            app.Stage2SpectrumAxes = uiaxes(gl);
            if app.IsCompact
                app.Stage2SpectrumAxes.Layout.Row = 3; app.Stage2SpectrumAxes.Layout.Column = [1 4];
            else
                app.Stage2SpectrumAxes.Layout.Row = 2; app.Stage2SpectrumAxes.Layout.Column = [3 4];
            end
            app.styleAxes(app.Stage2SpectrumAxes);
            title(app.Stage2SpectrumAxes, 'Current Stage2 spectrum');

            app.Stage2StatusArea = uitextarea(gl, 'Editable', 'off');
            if app.IsCompact
                app.Stage2StatusArea.Layout.Row = 4;
            else
                app.Stage2StatusArea.Layout.Row = 3;
            end
            app.Stage2StatusArea.Layout.Column = [1 4];
        end

        function buildRunControls(app, parent, stage2)
            if app.IsCompact
                app.buildRunControlsCompact(parent, stage2);
                return;
            end
            g = uigridlayout(parent, [6 12]);
            g.RowHeight = {31, 31, 31, 31, 31, 28};
            g.ColumnWidth = {92, 82, 92, 82, 92, 82, 92, 82, 120, 125, 125, '1x'};
            g.Padding = [16 12 16 12];
            g.RowSpacing = 8;
            g.ColumnSpacing = 12;
            g.BackgroundColor = [0.88 0.94 1.00];

            lab = app.techLabel(g, 'Pixel Y'); lab.Layout.Row = 1; lab.Layout.Column = 1;
            py = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            py.Layout.Row = 1; py.Layout.Column = 2;
            lab = app.techLabel(g, 'Pixel X'); lab.Layout.Row = 1; lab.Layout.Column = 3;
            px = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            px.Layout.Row = 1; px.Layout.Column = 4;
            lab = app.techLabel(g, 'Y bin'); lab.Layout.Row = 1; lab.Layout.Column = 5;
            sy = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            sy.Layout.Row = 1; sy.Layout.Column = 6;
            lab = app.techLabel(g, 'X bin'); lab.Layout.Row = 1; lab.Layout.Column = 7;
            sx = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            sx.Layout.Row = 1; sx.Layout.Column = 8;
            lab = app.techLabel(g, 'Energy bin'); lab.Layout.Row = 1; lab.Layout.Column = 9;
            eb = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            eb.Layout.Row = 1; eb.Layout.Column = 10;
            fullLabel = uilabel(g, 'Text', 'Full map by default');
            fullLabel.Layout.Row = 1; fullLabel.Layout.Column = [11 12];
            fullLabel.FontColor = [0.03 0.30 0.48];
            fullLabel.FontWeight = 'bold';

            lab = app.techLabel(g, 'Crop range');
            lab.Layout.Row = 2; lab.Layout.Column = 1;
            cropMin = uieditfield(g, 'numeric', 'Value', app.EnergyMinEdit.Value);
            cropMin.Layout.Row = 2; cropMin.Layout.Column = 2;
            cropMax = uieditfield(g, 'numeric', 'Value', app.EnergyMaxEdit.Value);
            cropMax.Layout.Row = 2; cropMax.Layout.Column = 3;
            inheritBtn = app.techButton(g, 'Use Data range', @(~, ~)app.inheritStageEnergyWindow(stage2));
            inheritBtn.Layout.Row = 2; inheritBtn.Layout.Column = [4 5];
            cropNote = uilabel(g, 'Text', 'Units: meV. A narrower range can reduce fitting time.');
            cropNote.Layout.Row = 2; cropNote.Layout.Column = [6 12];
            cropNote.FontColor = [0.05 0.30 0.48];

            if stage2
                app.Stage2PixelYSpinner = py; app.Stage2PixelXSpinner = px;
            else
                app.Stage1PixelYSpinner = py; app.Stage1PixelXSpinner = px;
            end

            workers = uispinner(g, 'Limits', [1 64], 'Value', EELS_Parallel.default_workers(), 'Step', 1);
            checkpoint = uicheckbox(g, 'Text', 'Save checkpoint', 'Value', false);
            resume = uicheckbox(g, 'Text', 'Resume checkpoint', 'Value', false);
            saveplots = uicheckbox(g, 'Text', 'Save diagnostics', 'Value', false);
            lab = app.techLabel(g, 'Output folder');
            outdir = uieditfield(g, 'text', 'Value', app.LastOutputDir);
            lab.Layout.Row = 3; lab.Layout.Column = 1;
            outdir.Layout.Row = 3; outdir.Layout.Column = [2 12];
            lab = app.techLabel(g, 'Workers');
            lab.Layout.Row = 4; lab.Layout.Column = 1;
            workers.Layout.Row = 4; workers.Layout.Column = 2;
            checkpoint.Layout.Row = 4; checkpoint.Layout.Column = [3 4];
            resume.Layout.Row = 4; resume.Layout.Column = [5 6];
            saveplots.Layout.Row = 4; saveplots.Layout.Column = [7 8];
            if stage2
                lab = app.techLabel(g, 'Algorithm');
                lab.Layout.Row = 4; lab.Layout.Column = 9;
                app.Stage2AlgorithmDropDown = uidropdown(g, ...
                    'Items', {'grid_scan'}, 'Value', 'grid_scan');
                app.Stage2AlgorithmDropDown.Layout.Row = 4; app.Stage2AlgorithmDropDown.Layout.Column = [10 12];
            end

            ctrl = struct('pixelY', py, 'pixelX', px, 'ybin', sy, 'xbin', sx, ...
                'ebin', eb, 'cropMin', cropMin, 'cropMax', cropMax, ...
                'workers', workers, 'checkpoint', checkpoint, 'resume', resume, ...
                'saveplots', saveplots, 'outdir', outdir);

            app.SpatialYBinSpinner = sy;
            app.SpatialXBinSpinner = sx;
            app.EnergyBinSpinner = eb;
            app.ParallelCheckBox = [];
            app.WorkerSpinner = workers;
            app.CheckpointCheckBox = checkpoint;
            app.ResumeCheckBox = resume;
            app.SavePlotsCheckBox = saveplots;
            app.OutputDirEdit = outdir;

            if stage2
                b = app.techButton(g, 'Fit pixel', @(~, ~)app.runSingle(true)); b.Layout.Row = 5; b.Layout.Column = 1;
                b = app.techButton(g, 'Pick pixel', @(~, ~)app.pickMapPixelWithCursor(true)); b.Layout.Row = 5; b.Layout.Column = 2;
                b = app.techButton(g, 'Full map: serial', @(~, ~)app.runDatacube(true, 'serial')); b.Layout.Row = 5; b.Layout.Column = 3;
                b = app.techButton(g, 'Full map: parallel', @(~, ~)app.runDatacube(true, 'parallel')); b.Layout.Row = 5; b.Layout.Column = 4;
                b = app.techButton(g, 'Export result', @(~, ~)app.exportResult('stage2')); b.Layout.Row = 5; b.Layout.Column = [5 6];
                b = app.techButton(g, 'Signal to Data', @(~, ~)app.importFitSignalToData()); b.Layout.Row = 5; b.Layout.Column = [7 8];
                lab = app.techLabel(g, 'Show map');
                lab.Layout.Row = 5; lab.Layout.Column = 9;
                app.Stage2MapDropDown = uidropdown(g, 'Items', app.mapItems(), 'Value', 'T data valid', ...
                    'ValueChangedFcn', @(~, ~)app.refreshResultMap(true));
                app.Stage2MapDropDown.Layout.Row = 5; app.Stage2MapDropDown.Layout.Column = [10 12];
                b = app.techButton(g, 'Reference tuning', @(~, ~)app.openReferenceTuner());
                b.Layout.Row = 6; b.Layout.Column = [1 3];
                outlierNote = uilabel(g, 'Text', '2-sigma refit removed; calibrate with independent reference spectra.');
                outlierNote.Layout.Row = 6; outlierNote.Layout.Column = [4 6];
                outlierNote.FontColor = [0.05 0.30 0.48];
                hierarchical = uicheckbox(g, 'Text', 'Hierarchical ZLP / interface', 'Value', true);
                hierarchical.Tooltip = 'Enable ZLP clustering, class-template constraints, Stage1 spatial smoothing and joint temperature MAP.';
                hierarchical.Layout.Row = 6; hierarchical.Layout.Column = [7 9];
                forceRawStage2 = uicheckbox(g, 'Text', 'Show raw Stage2', 'Value', false);
                forceRawStage2.Tooltip = 'Disable the Stage2 shape guard and Stage1 fallback to inspect the hard-BE fitted shape.';
                forceRawStage2.Layout.Row = 6; forceRawStage2.Layout.Column = [10 12];
                ctrl.hierarchical = hierarchical;
                ctrl.forceRawStage2 = forceRawStage2;
                ctrl.mapDropDown = app.Stage2MapDropDown;
                app.Stage2Controls = ctrl;
            else
                b = app.techButton(g, 'Fit pixel', @(~, ~)app.runSingle(false)); b.Layout.Row = 5; b.Layout.Column = 1;
                b = app.techButton(g, 'Pick pixel', @(~, ~)app.pickMapPixelWithCursor(false)); b.Layout.Row = 5; b.Layout.Column = 2;
                b = app.techButton(g, 'Full map: serial', @(~, ~)app.runDatacube(false, 'serial')); b.Layout.Row = 5; b.Layout.Column = 3;
                b = app.techButton(g, 'Full map: parallel', @(~, ~)app.runDatacube(false, 'parallel')); b.Layout.Row = 5; b.Layout.Column = 4;
                b = app.techButton(g, 'Export result', @(~, ~)app.exportResult('stage1')); b.Layout.Row = 5; b.Layout.Column = [5 6];
                b = app.techButton(g, 'Signal to Data', @(~, ~)app.importFitSignalToData()); b.Layout.Row = 5; b.Layout.Column = [7 8];
                lab = app.techLabel(g, 'Show map');
                lab.Layout.Row = 5; lab.Layout.Column = 9;
                app.MapDropDown = uidropdown(g, 'Items', app.mapItems(), 'Value', 'R2');
                app.MapDropDown.ValueChangedFcn = @(~, ~)app.refreshResultMap(false);
                app.MapDropDown.Layout.Row = 5; app.MapDropDown.Layout.Column = [10 12];
                ctrl.mapDropDown = app.MapDropDown;
                app.Stage1Controls = ctrl;
            end
            note = uilabel(g, 'Text', 'Checkpoints are optional for serial and parallel runs. Resume requires identical data/settings. Export diagnostic figures on demand.');
            note.Layout.Row = 6; note.Layout.Column = [1 12];
            note.FontColor = [0.18 0.33 0.50];
        end

        function buildRunControlsCompact(app, parent, stage2)
            g = uigridlayout(parent, [8 8]);
            g.RowHeight = {29, 29, 29, 29, 29, 29, 29, 29};
            g.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x', '1x', '1x'};
            g.Padding = [10 8 10 8];
            g.RowSpacing = 7;
            g.ColumnSpacing = 8;
            g.BackgroundColor = [0.88 0.94 1.00];

            lab = app.techLabel(g, 'Pixel Y'); lab.Layout.Row = 1; lab.Layout.Column = 1;
            py = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            py.Layout.Row = 1; py.Layout.Column = 2;
            lab = app.techLabel(g, 'Pixel X'); lab.Layout.Row = 1; lab.Layout.Column = 3;
            px = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            px.Layout.Row = 1; px.Layout.Column = 4;
            lab = app.techLabel(g, 'Y bin'); lab.Layout.Row = 1; lab.Layout.Column = 5;
            sy = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            sy.Layout.Row = 1; sy.Layout.Column = 6;
            lab = app.techLabel(g, 'X bin'); lab.Layout.Row = 1; lab.Layout.Column = 7;
            sx = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            sx.Layout.Row = 1; sx.Layout.Column = 8;

            lab = app.techLabel(g, 'Energy bin'); lab.Layout.Row = 2; lab.Layout.Column = 1;
            eb = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            eb.Layout.Row = 2; eb.Layout.Column = 2;
            lab = app.techLabel(g, 'Workers'); lab.Layout.Row = 2; lab.Layout.Column = 3;
            workers = uispinner(g, 'Limits', [1 64], 'Value', EELS_Parallel.default_workers(), 'Step', 1);
            workers.Layout.Row = 2; workers.Layout.Column = 4;
            checkpoint = uicheckbox(g, 'Text', 'Save checkpoint', 'Value', false);
            checkpoint.Layout.Row = 2; checkpoint.Layout.Column = [5 6];
            resume = uicheckbox(g, 'Text', 'Resume checkpoint', 'Value', false);
            resume.Layout.Row = 2; resume.Layout.Column = [7 8];

            lab = app.techLabel(g, 'Crop range'); lab.Layout.Row = 3; lab.Layout.Column = 1;
            cropMin = uieditfield(g, 'numeric', 'Value', app.EnergyMinEdit.Value);
            cropMin.Layout.Row = 3; cropMin.Layout.Column = 2;
            cropMax = uieditfield(g, 'numeric', 'Value', app.EnergyMaxEdit.Value);
            cropMax.Layout.Row = 3; cropMax.Layout.Column = 3;
            inheritBtn = app.techButton(g, 'Use Data range', @(~, ~)app.inheritStageEnergyWindow(stage2));
            inheritBtn.Layout.Row = 3; inheritBtn.Layout.Column = [4 5];
            saveplots = uicheckbox(g, 'Text', 'Save diagnostics', 'Value', false);
            saveplots.Layout.Row = 3; saveplots.Layout.Column = [6 8];

            lab = app.techLabel(g, 'Output folder'); lab.Layout.Row = 4; lab.Layout.Column = 1;
            outdir = uieditfield(g, 'text', 'Value', app.LastOutputDir);
            outdir.Layout.Row = 4; outdir.Layout.Column = [2 8];

            if stage2
                lab = app.techLabel(g, 'Algorithm'); lab.Layout.Row = 5; lab.Layout.Column = 1;
                app.Stage2AlgorithmDropDown = uidropdown(g, ...
                    'Items', {'grid_scan'}, 'Value', 'grid_scan');
                app.Stage2AlgorithmDropDown.Layout.Row = 5;
                app.Stage2AlgorithmDropDown.Layout.Column = [2 4];
                hierarchical = uicheckbox(g, 'Text', 'Hierarchical ZLP / interface', 'Value', true);
                hierarchical.Tooltip = 'Enable ZLP clustering, class-template constraints, Stage1 spatial smoothing and joint temperature MAP.';
                hierarchical.Layout.Row = 5; hierarchical.Layout.Column = [5 6];
                forceRawStage2 = uicheckbox(g, 'Text', 'Show raw Stage2', 'Value', false);
                forceRawStage2.Tooltip = 'Disable the Stage2 shape guard and Stage1 fallback to inspect the hard-BE fitted shape.';
                forceRawStage2.Layout.Row = 5; forceRawStage2.Layout.Column = [7 8];
                app.Stage2PixelYSpinner = py; app.Stage2PixelXSpinner = px;
            else
                fullLabel = uilabel(g, 'Text', 'Datacubes use the full map; pixel coordinates affect only Fit pixel.');
                fullLabel.Layout.Row = 5; fullLabel.Layout.Column = [1 8];
                fullLabel.FontColor = [0.03 0.30 0.48];
                fullLabel.FontWeight = 'bold';
                app.Stage1PixelYSpinner = py; app.Stage1PixelXSpinner = px;
            end

            ctrl = struct('pixelY', py, 'pixelX', px, 'ybin', sy, 'xbin', sx, ...
                'ebin', eb, 'cropMin', cropMin, 'cropMax', cropMax, ...
                'workers', workers, 'checkpoint', checkpoint, 'resume', resume, ...
                'saveplots', saveplots, 'outdir', outdir);

            app.SpatialYBinSpinner = sy;
            app.SpatialXBinSpinner = sx;
            app.EnergyBinSpinner = eb;
            app.ParallelCheckBox = [];
            app.WorkerSpinner = workers;
            app.CheckpointCheckBox = checkpoint;
            app.ResumeCheckBox = resume;
            app.SavePlotsCheckBox = saveplots;
            app.OutputDirEdit = outdir;

            if stage2
                b = app.techButton(g, 'Fit pixel', @(~, ~)app.runSingle(true)); b.Layout.Row = 6; b.Layout.Column = [1 2];
                b = app.techButton(g, 'Pick pixel', @(~, ~)app.pickMapPixelWithCursor(true)); b.Layout.Row = 6; b.Layout.Column = [3 4];
                b = app.techButton(g, 'Full map: serial', @(~, ~)app.runDatacube(true, 'serial')); b.Layout.Row = 6; b.Layout.Column = [5 6];
                b = app.techButton(g, 'Full map: parallel', @(~, ~)app.runDatacube(true, 'parallel')); b.Layout.Row = 6; b.Layout.Column = [7 8];
                b = app.techButton(g, 'Export result', @(~, ~)app.exportResult('stage2')); b.Layout.Row = 7; b.Layout.Column = [1 2];
                b = app.techButton(g, 'Signal to Data', @(~, ~)app.importFitSignalToData()); b.Layout.Row = 7; b.Layout.Column = [3 4];
                b = app.techButton(g, 'Reference tuning', @(~, ~)app.openReferenceTuner()); b.Layout.Row = 7; b.Layout.Column = [5 6];
                lab = app.techLabel(g, 'Show map'); lab.Layout.Row = 7; lab.Layout.Column = 7;
                app.Stage2MapDropDown = uidropdown(g, 'Items', app.mapItems(), 'Value', 'T data valid', ...
                    'ValueChangedFcn', @(~, ~)app.refreshResultMap(true));
                app.Stage2MapDropDown.Layout.Row = 7; app.Stage2MapDropDown.Layout.Column = 8;
                ctrl.hierarchical = hierarchical;
                ctrl.forceRawStage2 = forceRawStage2;
                ctrl.mapDropDown = app.Stage2MapDropDown;
                app.Stage2Controls = ctrl;
            else
                b = app.techButton(g, 'Fit pixel', @(~, ~)app.runSingle(false)); b.Layout.Row = 6; b.Layout.Column = [1 2];
                b = app.techButton(g, 'Pick pixel', @(~, ~)app.pickMapPixelWithCursor(false)); b.Layout.Row = 6; b.Layout.Column = [3 4];
                b = app.techButton(g, 'Full map: serial', @(~, ~)app.runDatacube(false, 'serial')); b.Layout.Row = 6; b.Layout.Column = [5 6];
                b = app.techButton(g, 'Full map: parallel', @(~, ~)app.runDatacube(false, 'parallel')); b.Layout.Row = 6; b.Layout.Column = [7 8];
                b = app.techButton(g, 'Export result', @(~, ~)app.exportResult('stage1')); b.Layout.Row = 7; b.Layout.Column = [1 2];
                b = app.techButton(g, 'Signal to Data', @(~, ~)app.importFitSignalToData()); b.Layout.Row = 7; b.Layout.Column = [3 4];
                lab = app.techLabel(g, 'Show map'); lab.Layout.Row = 7; lab.Layout.Column = 5;
                app.MapDropDown = uidropdown(g, 'Items', app.mapItems(), 'Value', 'R2', ...
                    'ValueChangedFcn', @(~, ~)app.refreshResultMap(false));
                app.MapDropDown.Layout.Row = 7; app.MapDropDown.Layout.Column = [6 8];
                ctrl.mapDropDown = app.MapDropDown;
                app.Stage1Controls = ctrl;
            end

            note = uilabel(g, 'Text', 'Checkpoints are optional for serial and parallel runs. Resume requires identical data/settings. Export diagnostic figures on demand.');
            note.Layout.Row = 8; note.Layout.Column = [1 8];
            note.FontColor = [0.18 0.33 0.50];
        end

        function buildViewTab(app)
            gl = uigridlayout(app.ViewTab, [3 4]);
            if app.IsCompact
                gl.RowHeight = {140, '1x', '1x'};
                gl.ColumnWidth = {'1x', '1x', '1x', '1x'};
            else
                gl.RowHeight = {92, '1x', '1x'};
                gl.ColumnWidth = {280, '1x', '1x', 250};
            end
            gl.Padding = [12 12 12 12];

            p = uipanel(gl, 'Title', 'Diagnostic', 'FontWeight', 'bold');
            p.Layout.Row = 1; p.Layout.Column = [1 4];
            if app.IsCompact
                g = uigridlayout(p, [3 3]);
                g.RowHeight = {30, 30, 30};
                g.ColumnWidth = {'1x', '1x', '1x'};
            else
                g = uigridlayout(p, [1 6]);
                g.ColumnWidth = {170, 150, 120, 120, 140, '1x'};
            end
            app.DiagnosticDropDown = uidropdown(g, 'Items', { ...
                'Full diagnostic figure', 'Raw / total fit', 'ZLP components', 'Sparse convolution / vibration', ...
                'Residual / confidence interval', 'Sparse intensity sticks', 'log(Iloss/Igain)', 'Attention curve', 'Loss components', ...
                'Log signal / log ZLP'}, 'Value', 'Full diagnostic figure');
            if app.IsCompact
                app.DiagnosticDropDown.Layout.Row = 1; app.DiagnosticDropDown.Layout.Column = [1 3];
                b = uibutton(g, 'Text', 'Refresh view', 'ButtonPushedFcn', @(~, ~)app.plotDiagnosticSelection()); b.Layout.Row = 2; b.Layout.Column = 1;
                b = uibutton(g, 'Text', 'Open full figure', 'ButtonPushedFcn', @(~, ~)app.openFullDiagnostic()); b.Layout.Row = 2; b.Layout.Column = 2;
                b = uibutton(g, 'Text', 'Export current result', 'ButtonPushedFcn', @(~, ~)app.exportResult('current')); b.Layout.Row = 2; b.Layout.Column = 3;
                b = uibutton(g, 'Text', 'Load result MAT', 'ButtonPushedFcn', @(~, ~)app.onLoadResult()); b.Layout.Row = 3; b.Layout.Column = 1;
                lab = uilabel(g, 'Text', 'Select a section of the main diagnostic figure.'); lab.Layout.Row = 3; lab.Layout.Column = [2 3];
            else
                uibutton(g, 'Text', 'Refresh view', 'ButtonPushedFcn', @(~, ~)app.plotDiagnosticSelection());
                uibutton(g, 'Text', 'Open full figure', 'ButtonPushedFcn', @(~, ~)app.openFullDiagnostic());
                uibutton(g, 'Text', 'Export current result', 'ButtonPushedFcn', @(~, ~)app.exportResult('current'));
                uibutton(g, 'Text', 'Load result MAT', 'ButtonPushedFcn', @(~, ~)app.onLoadResult());
                uilabel(g, 'Text', 'Select a section of the main diagnostic figure.');
            end

            app.ViewAxes1 = uiaxes(gl);
            app.ViewAxes1.Layout.Row = [2 3]; app.ViewAxes1.Layout.Column = [1 2];
            app.ViewAxes2 = uiaxes(gl);
            app.ViewAxes2.Layout.Row = [2 3]; app.ViewAxes2.Layout.Column = [3 4];
        end

        function buildAdvancedTab(app)
            gl = uigridlayout(app.AdvancedTab, [2 4]);
            gl.RowHeight = {'1x', 120};
            gl.ColumnWidth = {'1x', '1x', 160, 160};
            gl.Padding = [12 12 12 12];
            app.AdvancedTable = uitable(gl, 'ColumnName', {'Parameter name', 'Description', 'Current value', 'Advice'}, ...
                'ColumnEditable', [false false true false], 'CellEditCallback', @(~, ~)app.applyAdvancedTable());
            app.AdvancedTable.Layout.Row = 1; app.AdvancedTable.Layout.Column = [1 4];
            app.AdvancedHelpArea = uitextarea(gl, 'Editable', 'off');
            app.AdvancedHelpArea.Layout.Row = 2; app.AdvancedHelpArea.Layout.Column = [1 2];
            app.AdvancedHelpArea.Value = { ...
                'Advanced settings:', ...
                '1. Stage1 fits the spectrum. Start with EM spacing, ZLP width and tail/feature weights.', ...
                '2. Stage2 estimates temperature from hard-BE pairs in the attention window. Inspect attention and reliability if unstable.', ...
                '3. Stage1/2 process the full map. Pixel coordinates apply only to Fit pixel; checkpoints are optional.'};
            b1 = uibutton(gl, 'Text', 'Apply settings', 'ButtonPushedFcn', @(~, ~)app.applyAdvancedTable());
            b1.Layout.Row = 2; b1.Layout.Column = 3;
            b2 = uibutton(gl, 'Text', 'Reset defaults', 'ButtonPushedFcn', @(~, ~)app.resetAdvanced());
            b2.Layout.Row = 2; b2.Layout.Column = 4;
        end

        function items = mapItems(~)
            items = {'R2', 'NRMSE', 'T raw', 'T data valid', 'T joint', ...
                'T joint uncertainty', 'Spatial inferred', 'T repaired', 'T spatial refined', ...
                'T uncertainty', 'Reliability', 'Attention fraction', 'Warning', ...
                'Bad pixel', 'Interface', 'Neighbor loss', 'Spectrum neighbor loss', ...
                'Temperature neighbor z', 'Attention peak delta', 'Signal peak delta', ...
                'ZLP guard start', 'ZLP guard fraction', 'ZLP core width', ...
                'ZLP cluster', 'Cluster membership', 'Cluster outlier', ...
                'Template deviation', 'ZLP spatial loss', 'ZLP interface score', ...
                'ZLP sharp interface'};
        end

        function tf = hasExternalTool(app, function_name, display_name)
            if nargin < 3 || isempty(display_name)
                display_name = function_name;
            end
            tf = exist(function_name, 'file') == 2 || exist(function_name, 'file') == 6;
            if ~tf
                app.logStatus(sprintf('%s was not found on the MATLAB path.', display_name));
            end
        end

        function onLoadData(app)
            [file, path] = uigetfile({'*.mat;*.dm3', 'EELS data (*.mat, *.dm3)'; '*.mat', 'MAT files (*.mat)'; '*.dm3', 'DM3 files (*.dm3)'}, 'Select EELS datacube MAT/DM3');
            if isequal(file, 0), return; end
            filename = fullfile(path, file);
            d = uiprogressdlg(app.UIFigure, 'Title', 'Load data', 'Message', 'Reading datacube...', 'Indeterminate', 'on');
            try
                [cube, ene, info] = EELS_PearsonSparse_FullFit.load_datacube_mat(filename);
                app.Data = cube;
                app.Ene = ene(:);
                app.DataFile = filename;
                app.DataLabel = sprintf('%s (%s)', file, info.source_variable);
                app.FileEdit.Value = filename;
                app.clearFitCaches();
                app.setPixelLimits();
                app.refreshDataPreview();
                app.logStatus(sprintf('Loaded %s, size = [%s], energy channels = %d', app.DataLabel, num2str(size(cube)), numel(ene)));
            catch ME
                app.showError(ME);
            end
            close(d);
        end

        function onLoadResult(app)
            [file, path] = uigetfile({'*.mat', 'MAT files (*.mat)'}, 'Select result MAT');
            if isequal(file, 0), return; end
            filename = fullfile(path, file);
            try
                S = load(filename);
                if isfield(S, 'cuberes')
                    app.CubeResult = S.cuberes;
                    if app.cuberesHasStage2(S.cuberes)
                        app.Stage2CubeResult = S.cuberes;
                    else
                        app.Stage1CubeResult = S.cuberes;
                    end
                    if isfield(app.CubeResult, 'results') && ~isempty(app.CubeResult.results)
                        app.CurrentResult = app.firstValidResult(app.CubeResult.results);
                    end
                    if isempty(app.Data) && isfield(app.CubeResult, 'R2_map') && ~isempty(app.CurrentResult) && ...
                            isfield(app.CurrentResult, 'ene_corrected')
                        [ny0, nx0] = size(app.CubeResult.R2_map);
                        app.Ene = app.CurrentResult.ene_corrected(:);
                        app.Data = NaN(ny0, nx0, numel(app.Ene));
                        app.DataLabel = 'loaded result placeholder';
                        app.FileEdit.Value = filename;
                        app.setPixelLimits();
                    end
                    app.refreshResultMap(true);
                    app.logStatus(['Loaded datacube result: ', filename]);
                elseif isfield(S, 'result')
                    app.CurrentResult = S.result;
                    if app.resultHasStage2(S.result)
                        app.Stage2SingleResult = S.result;
                    else
                        app.Stage1SingleResult = S.result;
                    end
                    if isempty(app.Data) && isfield(S.result, 'ene_corrected')
                        app.Ene = S.result.ene_corrected(:);
                        app.Data = reshape(S.result.raw_data(:), 1, 1, []);
                        app.DataLabel = 'loaded single result';
                        app.FileEdit.Value = filename;
                        app.setPixelLimits();
                    end
                    app.plotDiagnosticSelection();
                    app.logStatus(['Loaded single-spectrum result: ', filename]);
                elseif isfield(S, 'testres')
                    app.CurrentResult = app.firstValidResult(S.testres.results);
                    app.plotDiagnosticSelection();
                    app.logStatus(['Loaded simulation result: ', filename]);
                else
                    names = fieldnames(S);
                    app.CurrentResult = S.(names{1});
                    app.plotDiagnosticSelection();
                end
            catch ME
                app.showError(ME);
            end
        end

        function onAlign(app)
            if ~app.hasExternalTool('Align', 'Align')
                return;
            end
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            d = uiprogressdlg(app.UIFigure, 'Title', 'Align', 'Message', 'Running 3D-EELS alignment...', 'Indeterminate', 'on');
            try
                obj = app.makeEELSObject(app.Data, app.Ene);
                out = Align(obj, app.AlignMethodDropDown.Value, 0);
                app.Data = double(out.data);
                app.Ene = double(out.ene(:));
                app.clearFitCaches();
                app.setPixelLimits();
                app.refreshDataPreview();
                app.logStatus(['Alignment complete. Method: ', app.AlignMethodDropDown.Value]);
            catch ME
                app.showError(ME);
            end
            close(d);
        end

        function onBM3D(app)
            if ~app.hasExternalTool('BM3Ddenoise', 'BM3D denoise')
                return;
            end
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            d = uiprogressdlg(app.UIFigure, 'Title', 'BM3D', 'Message', 'Running BM3D denoising...', 'Indeterminate', 'on');
            try
                obj = app.makeEELSObject(app.Data, app.Ene);
                win = app.energyWindow();
                out = BM3Ddenoise(obj, win, app.BM3DFactorEdit.Value);
                app.Data = double(out.data);
                app.Ene = double(out.ene(:));
                app.clearFitCaches();
                app.refreshDataPreview();
                app.logStatus(sprintf('BM3D complete, window [%.1f, %.1f] meV.', win(1), win(2)));
            catch ME
                app.showError(ME);
            end
            close(d);
        end

        function runSingle(app, stage2)
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            settings = app.collectSettings(stage2);
            if stage2
                settings = app.applyStage2SingleFastSettings(settings);
            end
            [yy, xx] = app.currentPixel(stage2);
            spec = squeeze(app.Data(yy, xx, :));
            d = uiprogressdlg(app.UIFigure, 'Title', 'Single-pixel fit', 'Message', 'Fitting the selected pixel...', 'Indeterminate', 'on');
            try
                reuse = stage2 && ~isempty(app.Stage1SingleResult) && ...
                    isfield(app.Stage1SingleResult, 'pixel_y') && isfield(app.Stage1SingleResult, 'pixel_x') && ...
                    app.Stage1SingleResult.pixel_y == yy && app.Stage1SingleResult.pixel_x == xx;
                if reuse
                    res = EELS_PearsonSparse_FullFit.promote_stage1_to_stage2(app.Stage1SingleResult, settings);
                    d.Message = 'Reusing Stage1; estimating Stage2 temperature...';
                else
                    res = EELS_PearsonSparse_FullFit.solve_spectrum(app.Ene, spec, settings);
                end
                [ny, nx, ~] = size(app.Data);
                res.pixel_y = yy; res.pixel_x = xx; res.pixel_index = sub2ind([ny, nx], yy, xx);
                app.CurrentResult = res;
                if stage2
                    app.Stage2SingleResult = res;
                else
                    app.Stage1SingleResult = res;
                end
                app.updateSingleFitAxes(stage2);
                app.plotDiagnosticSelection();
                msg = sprintf('Pixel (%d,%d) complete: R2 %.6g, NRMSE %.4g, T %.1f K, reliability %.2f', ...
                    yy, xx, res.metrics.R2, res.metrics.NRMSE, res.T_fit, res.reliability_score);
                app.logStage(stage2, msg);
            catch ME
                app.showError(ME);
            end
            close(d);
        end

        function runDatacube(app, stage2, run_mode)
            if nargin < 3 || isempty(run_mode)
                run_mode = 'serial';
            end
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            app.DatacubePauseRequested = false;
            settings = app.collectSettings(stage2, run_mode);
            if ~isfolder(settings.test_output_dir), mkdir(settings.test_output_dir); end
            input_file = app.writeGuiCube(settings.test_output_dir);
            d = uiprogressdlg(app.UIFigure, 'Title', 'Datacube processing', ...
                'Message', sprintf('Running datacube fitting in %s mode.', app.runModeLabel(run_mode)), ...
                'Indeterminate', 'off', 'Value', 0, 'Cancelable', 'on');
            settings.datacube_progress_callback = @(done,total,elapsed,eta,med,label)app.updateProgressDialog(d, done, total, elapsed, eta, med, label);
            settings.datacube_stop_callback = @()app.shouldPauseDatacube(d);
            cleanup = onCleanup(@()app.closeProgressDialog(d));
            try
                if stage2 && ~isempty(app.Stage1CubeResult) && app.isReusableStage1Cube(app.Stage1CubeResult)
                    d.Message = 'Reusing Stage1 datacube results for Stage2 temperature...';
                    cuberes = EELS_PearsonSparse_FullFit.promote_datacube_stage1_to_stage2(app.Stage1CubeResult, settings);
                else
                    cuberes = EELS_PearsonSparse_FullFit.run_test_datacube(input_file, settings);
                end
                app.CubeResult = cuberes;
                if stage2
                    app.Stage2CubeResult = cuberes;
                else
                    app.Stage1CubeResult = cuberes;
                end
                app.CurrentResult = app.firstValidResult(cuberes.results);
                app.refreshResultMap(stage2);
                app.plotDiagnosticSelection();
                msg = sprintf('Datacube complete: success %d, median R2 %.6g, median NRMSE %.4g, output=%s', ...
                    nnz(cuberes.fit_status_map > 0), median(cuberes.R2_map(:), 'omitnan'), ...
                    median(cuberes.NRMSE_map(:), 'omitnan'), settings.test_output_dir);
                app.logStage(stage2, msg);
            catch ME
                if strcmp(ME.identifier, 'EELS:DatacubePaused')
                    app.logStage(stage2, sprintf('Datacube paused; checkpoint: %s. Enable Resume checkpoint and keep the output folder to continue.', settings.test_output_dir));
                else
                    app.showError(ME);
                end
            end
            app.closeProgressDialog(d);
            clear cleanup
        end

        function settings = collectSettings(app, stage2, run_mode)
            if nargin < 3 || isempty(run_mode)
                run_mode = 'serial';
            end
            app.applyAdvancedTable();
            ctrl = app.controlsForStage(stage2);
            settings = app.Settings;
            settings.run_stage2 = logical(stage2);
            settings.test_output_dir = ctrl.outdir.Value;
            app.LastOutputDir = settings.test_output_dir;
            settings.crop_range = app.stageEnergyWindow(stage2);
            settings.datacube_y_bin = round(ctrl.ybin.Value);
            settings.datacube_x_bin = round(ctrl.xbin.Value);
            settings.datacube_energy_bin = round(ctrl.ebin.Value);
            settings.datacube_sample_strategy = 'all';
            settings.datacube_sample_pixels = Inf;
            settings.datacube_parallel_workers = max(1, round(ctrl.workers.Value));
            settings.datacube_checkpoint_enabled = logical(ctrl.checkpoint.Value);
            settings.datacube_resume_checkpoint = logical(ctrl.resume.Value);
            settings.save_plots = logical(ctrl.saveplots.Value);
            settings.save_diagnostic_plots = logical(ctrl.saveplots.Value);
            settings.force_hidden_figures = true;
            settings.close_parallel_pool_after_run = false;
            switch lower(string(run_mode))
                case "parallel"
                    settings.datacube_parallel_enabled = true;
                    settings.close_parallel_pool_after_run = false;
                otherwise
                    settings.datacube_parallel_enabled = false;
                    % Checkpointing follows the user checkbox in serial mode too.
            end
            settings.max_saved_spectrum_plots = 8;
            settings.datacube_progress_every = 1;
            if stage2 && ~isempty(app.Stage2AlgorithmDropDown) && isvalid(app.Stage2AlgorithmDropDown)
                settings.stage2_temperature_algorithm = app.Stage2AlgorithmDropDown.Value;
            end
            if stage2 && isfield(ctrl, 'hierarchical') && ~isempty(ctrl.hierarchical) && ...
                    isvalid(ctrl.hierarchical)
                use_hierarchy = logical(ctrl.hierarchical.Value);
                settings.stage1_hierarchical_enabled = use_hierarchy;
                settings.stage1_spatial_zlp_enabled = use_hierarchy;
                settings.stage2_spatial_map_enabled = use_hierarchy;
                settings.spatial_refit_enabled = false;
            end
            if stage2 && isfield(ctrl, 'forceRawStage2') && ~isempty(ctrl.forceRawStage2) && ...
                    isvalid(ctrl.forceRawStage2) && logical(ctrl.forceRawStage2.Value)
                settings.stage2_shape_guard_enabled = false;
                settings.stage2_force_raw_display = true;
            end
        end

        function label = runModeLabel(~, run_mode)
            if strcmpi(char(run_mode), 'parallel')
                label = 'parallel';
            else
                label = 'serial';
            end
        end

        function settings = applyStage2SingleFastSettings(~, settings)
            if isfield(settings, 'stage2_single_fast_enabled') && ~settings.stage2_single_fast_enabled
                return;
            end
            if isfield(settings, 'stage2_single_fast_T_grid')
                settings.stage2_T_grid = settings.stage2_single_fast_T_grid;
            else
                settings.stage2_T_grid = 80:60:1800;
            end
            if isfield(settings, 'stage2_single_fast_fine_step_K')
                settings.stage2_fine_step_K = settings.stage2_single_fast_fine_step_K;
            else
                settings.stage2_fine_step_K = 12;
            end
            if isfield(settings, 'stage2_single_fast_fine_half_width_K')
                settings.stage2_fine_half_width_K = settings.stage2_single_fast_fine_half_width_K;
            else
                settings.stage2_fine_half_width_K = 90;
            end
            if isfield(settings, 'stage2_single_fast_hard_em_max_centers')
                settings.stage2_hard_em_max_centers = min(settings.stage2_hard_em_max_centers, ...
                    settings.stage2_single_fast_hard_em_max_centers);
            end
        end

        function refreshDataPreview(app)
            if isempty(app.Data), return; end
            app.plotEnergyMap();
            [yy, xx] = app.currentPixel(false);
            app.plotRawSpectrum(yy, xx, app.SpectrumAxes);
        end

        function plotEnergyMap(app)
            if isempty(app.Data), return; end
            win = app.energyWindow();
            mask = app.Ene >= win(1) & app.Ene <= win(2);
            if ~any(mask), mask = true(size(app.Ene)); end
            img = sum(app.Data(:, :, mask), 3, 'omitnan');
            imagesc(app.DataAxes, img);
            axis(app.DataAxes, 'image'); colorbar(app.DataAxes);
            title(app.DataAxes, sprintf('Energy-filtered image [%.1f, %.1f] meV', win(1), win(2)));
            h = app.DataAxes.Children(1);
            h.HitTest = 'on';
            h.ButtonDownFcn = @(~, evt)app.onDataMapClick(evt);
            app.DataAxes.ButtonDownFcn = @(~, evt)app.onDataMapClick(evt);
        end

        function onDataMapClick(app, evt)
            if isempty(app.Data), return; end
            if isstruct(evt) && isfield(evt, 'IntersectionPoint')
                cp = evt.IntersectionPoint;
            elseif isprop(evt, 'IntersectionPoint')
                cp = evt.IntersectionPoint;
            else
                cp = app.DataAxes.CurrentPoint(1, :);
            end
            xx = max(1, round(cp(1)));
            yy = max(1, round(cp(2)));
            [ny, nx, ~] = size(app.Data);
            if xx > nx || yy > ny, return; end
            if ~isempty(app.Stage1PixelYSpinner) && isvalid(app.Stage1PixelYSpinner)
                app.Stage1PixelYSpinner.Value = min(max(yy, 1), app.Stage1PixelYSpinner.Limits(2));
                app.Stage1PixelXSpinner.Value = min(max(xx, 1), app.Stage1PixelXSpinner.Limits(2));
            end
            app.plotRawSpectrum(yy, xx, app.SpectrumAxes);
            app.markDataPixel(xx, yy);
        end

        function markDataPixel(app, xx, yy)
            app.deleteGraphicsHandles('DataCursorHandles');
            hold(app.DataAxes, 'on');
            app.DataCursorHandles = [ ...
                plot(app.DataAxes, xx, yy, 'wo', 'MarkerSize', 8, 'LineWidth', 1.5), ...
                plot(app.DataAxes, xx, yy, '+', 'Color', [0.85 0.08 0.10], 'MarkerSize', 8, 'LineWidth', 1.3)];
            for ii = 1:numel(app.DataCursorHandles)
                try
                    app.DataCursorHandles(ii).HitTest = 'off';
                    app.DataCursorHandles(ii).PickableParts = 'none';
                catch
                end
            end
            hold(app.DataAxes, 'off');
        end

        function plotLineProfile(app)
            if isempty(app.Data), return; end
            win = app.energyWindow();
            mask = app.Ene >= win(1) & app.Ene <= win(2);
            if ~any(mask), mask = true(size(app.Ene)); end
            prof = squeeze(sum(app.Data(:, :, mask), 1, 'omitnan'));
            imagesc(app.DataAxes, 1:size(prof, 1), app.Ene(mask), prof.');
            set(app.DataAxes, 'YDir', 'normal');
            pbaspect(app.DataAxes, [1.8 1 1]);
            colorbar(app.DataAxes);
            xlabel(app.DataAxes, 'x position'); ylabel(app.DataAxes, 'Energy loss (meV)');
            title(app.DataAxes, 'Line profile: sum over y');
        end

        function selectROIAndIntegrate(app)
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            app.plotEnergyMap();
            app.logStatus('Draw a rectangular ROI on the energy map; double-click or press Enter to confirm.');
            try
                roi = drawrectangle(app.DataAxes, 'Color', [0.05 0.65 1.00], 'LineWidth', 1.5);
                [pos, cancelled] = app.waitForROIWithHotkeys(roi, ...
                    'ROI: double-click to confirm; Ctrl+Q cancels; Ctrl+H confirms and plots.');
                if cancelled || isempty(pos)
                    app.logStatus('ROI selection cancelled.');
                    return;
                end
                if isempty(pos)
                    error('ROI position is unavailable after rectangle selection.');
                end
                if isempty(pos) || any(~isfinite(pos)), return; end
            catch ME
                app.showError(ME);
                return;
            end
            [ny, nx, ~] = size(app.Data);
            x1 = max(1, floor(pos(1)));
            y1 = max(1, floor(pos(2)));
            x2 = min(nx, ceil(pos(1) + pos(3)));
            y2 = min(ny, ceil(pos(2) + pos(4)));
            if x2 < x1 || y2 < y1
                app.logStatus('Invalid ROI. Select it again.');
                return;
            end
            app.SelectedROIRange = [y1 y2 x1 x2];
            app.SelectedROICube = app.Data(y1:y2, x1:x2, :);
            app.SelectedROISpectrum = squeeze(sum(app.SelectedROICube, [1 2], 'omitnan'));
            cla(app.SpectrumAxes);
            plot(app.SpectrumAxes, app.Ene, app.SelectedROISpectrum, 'Color', [0.05 0.20 0.38], 'LineWidth', 1.4);
            grid(app.SpectrumAxes, 'on');
            xlabel(app.SpectrumAxes, 'Energy loss (meV)'); ylabel(app.SpectrumAxes, 'Integrated intensity');
            title(app.SpectrumAxes, sprintf('ROI integrated spectrum: y %d-%d, x %d-%d', y1, y2, x1, x2));
            app.logStatus(sprintf('ROI integrated: y %d-%d, x %d-%d. Use Export ROI to save.', y1, y2, x1, x2));
        end

        function exportROI(app)
            if isempty(app.SelectedROICube)
                app.logStatus('No ROI result. Use Integrate ROI first.');
                return;
            end
            assignin('base', 'eels_roi_cube', app.SelectedROICube);
            assignin('base', 'eels_roi_spectrum', app.SelectedROISpectrum);
            assignin('base', 'eels_roi_energy', app.Ene);
            assignin('base', 'eels_roi_range_y1y2x1x2', app.SelectedROIRange);
            app.logStatus('Exported eels_roi_cube / eels_roi_spectrum / eels_roi_energy / eels_roi_range_y1y2x1x2 to the workspace.');
        end

        function openStackSpectraTool(app)
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            app.plotEnergyMap();
            answer = inputdlg({'Number of spectra to stack', 'Selection mode: cursor or roi'}, ...
                'Stack spectra', [1 42], {'3', 'cursor'});
            if isempty(answer), return; end
            n_pick = max(1, round(str2double(answer{1})));
            if ~isfinite(n_pick), n_pick = 3; end
            mode = lower(strtrim(answer{2}));
            if ~ismember(mode, {'cursor', 'roi'}), mode = 'cursor'; end
            app.deleteGraphicsHandles('StackCursorHandles');
            app.logStatus('Stack selection: double-click confirms; Ctrl+Q cancels; Ctrl+H confirms and plots. Stack cursors are purple/blue.');
            stack = struct('name', {{}}, 'spectrum', {{}}, 'selection', {{}}, ...
                'mode', mode, 'ene', app.Ene(:), 'source_file', app.FileEdit.Value);
            hold(app.DataAxes, 'on');
            for kk = 1:n_pick
                try
                    if strcmp(mode, 'roi')
                        roi = drawrectangle(app.DataAxes, 'Color', [0.00 0.55 1.00], 'LineWidth', 1.4);
                        [pos, cancelled] = app.waitForROIWithHotkeys(roi, ...
                            sprintf('Spectrum %d/%d: double-click confirms; Ctrl+Q cancels; Ctrl+H confirms and plots.', kk, n_pick));
                        if cancelled || isempty(pos), break; end
                        [spec, selection] = app.spectrumFromROIPosition(pos);
                        hmark = rectangle(app.DataAxes, 'Position', selection.position, ...
                            'EdgeColor', [0.00 0.48 0.95], 'LineStyle', '--', 'LineWidth', 1.2);
                        app.appendGraphicsHandle('StackCursorHandles', hmark);
                    else
                        roi = drawpoint(app.DataAxes, 'Color', [0.55 0.00 0.95], 'LineWidth', 1.6);
                        [pos, cancelled] = app.waitForROIWithHotkeys(roi, ...
                            sprintf('Spectrum %d/%d: double-click confirms; Ctrl+Q cancels; Ctrl+H confirms and plots.', kk, n_pick));
                        if cancelled || isempty(pos), break; end
                        [spec, selection] = app.spectrumFromPointPosition(pos);
                        hmark = plot(app.DataAxes, selection.x, selection.y, 'o', ...
                            'Color', [0.50 0.00 0.95], 'MarkerFaceColor', [0.90 0.75 1.00], ...
                            'MarkerSize', 7, 'LineWidth', 1.2);
                        app.appendGraphicsHandle('StackCursorHandles', hmark);
                    end
                catch ME
                    hold(app.DataAxes, 'off');
                    app.showError(ME);
                    break;
                end
                default_name = sprintf('%s_%02d', mode, kk);
                name_answer = inputdlg({'Spectrum name'}, 'Name spectrum', [1 42], {default_name});
                if isempty(name_answer)
                    label = default_name;
                else
                    label = strtrim(name_answer{1});
                    if isempty(label), label = default_name; end
                end
                stack.name{end + 1} = label;
                stack.spectrum{end + 1} = spec(:);
                stack.selection{end + 1} = selection;
            end
            hold(app.DataAxes, 'off');
            if isempty(stack.spectrum), return; end
            app.StackSpectra = stack;
            app.plotStackedSpectra(stack, 'EELS stacked spectra');
            assignin('base', 'eels_stacked_spectra', stack);
            app.logStatus('Stack exported to workspace variable eels_stacked_spectra.');
        end

        function [spec, selection] = spectrumFromPointPosition(app, pos)
            [ny, nx, ~] = size(app.Data);
            xx = min(nx, max(1, round(pos(1))));
            yy = min(ny, max(1, round(pos(2))));
            spec = squeeze(app.Data(yy, xx, :));
            selection = struct('type', 'cursor', 'x', xx, 'y', yy, ...
                'label', sprintf('pixel (%d,%d)', yy, xx));
        end

        function [spec, selection] = spectrumFromROIPosition(app, pos)
            [ny, nx, ~] = size(app.Data);
            x1 = max(1, floor(pos(1)));
            y1 = max(1, floor(pos(2)));
            x2 = min(nx, ceil(pos(1) + pos(3)));
            y2 = min(ny, ceil(pos(2) + pos(4)));
            roi_cube = app.Data(y1:y2, x1:x2, :);
            spec = squeeze(sum(roi_cube, [1 2], 'omitnan'));
            selection = struct('type', 'roi', 'x_range', [x1 x2], 'y_range', [y1 y2], ...
                'position', [x1 y1 max(1, x2 - x1 + 1) max(1, y2 - y1 + 1)], ...
                'label', sprintf('ROI x%d-%d y%d-%d', x1, x2, y1, y2));
        end

        function plotStackedSpectra(app, stack, fig_name)
            fig = figure('Name', fig_name, 'Color', 'w', 'Position', [120 120 980 640]);
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

        function applyTraditionalBackground(app, interactive)
            if ~app.hasExternalTool('BackgroundSubtraction', 'Traditional background subtraction')
                return;
            end
            if nargin < 2, interactive = true; end
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            method = app.BGMethodDropDown.Value;
            win = [app.BGWin1MinEdit.Value, app.BGWin1MaxEdit.Value, ...
                app.BGWin2MinEdit.Value, app.BGWin2MaxEdit.Value];
            if any(~isfinite(win)) || win(1) >= win(2) || win(3) >= win(4)
                app.logStatus('Invalid background windows. Enter [w1min w1max w2min w2max] as in 3D-EELS Toolbox.');
                return;
            end
            d = uiprogressdlg(app.UIFigure, 'Title', 'Background', ...
                'Message', sprintf('Running BackgroundSubtraction(..., windows, %s)...', method), 'Indeterminate', 'on');
            if interactive && ~isempty(d) && isvalid(d)
                close(d);
                d = [];
            end
            try
                obj = app.makeEELSObject(app.Data, app.Ene);
                out = BackgroundSubtraction(obj, win, method, logical(interactive));
                app.Data = double(out.data);
                app.Ene = double(out.ene(:));
                app.clearFitCaches();
                app.refreshDataPreview();
                app.logStatus(sprintf('Background subtraction complete: %s, windows [%.1f %.1f %.1f %.1f] meV.', ...
                    method, win(1), win(2), win(3), win(4)));
            catch ME
                app.showError(ME);
            end
            if ~isempty(d) && isvalid(d), close(d); end
        end

        function importFitSignalToData(app)
            if isempty(app.CurrentResult) && isempty(app.CubeResult)
                app.logStatus('No fitted result to import.');
                return;
            end
            try
                if ~isempty(app.CubeResult) && isfield(app.CubeResult, 'results')
                    first = app.firstValidResult(app.CubeResult.results);
                    if isempty(first), app.logStatus('No valid pixel results in cuberes.'); return; end
                    ene = first.ene_corrected(:);
                    [ny, nx] = size(app.CubeResult.R2_map);
                    sigcube = NaN(ny, nx, numel(ene));
                    for ii = 1:numel(app.CubeResult.results)
                        res = app.CubeResult.results{ii};
                        if isstruct(res) && ~isfield(res, 'error') && isfield(res, 'pixel_y') && isfield(res, 'pixel_x')
                            sigcube(res.pixel_y, res.pixel_x, :) = app.resultSignalSpectrum(res);
                        end
                    end
                    app.Data = sigcube;
                    app.Ene = ene;
                    app.DataLabel = 'fit signal cube';
                else
                    ene = app.CurrentResult.ene_corrected(:);
                    app.Data = reshape(app.resultSignalSpectrum(app.CurrentResult), 1, 1, []);
                    app.Ene = ene;
                    app.DataLabel = 'fit signal single spectrum';
                end
                app.FileEdit.Value = ['Imported from fit: ', app.DataLabel];
                app.clearFitCaches();
                app.setPixelLimits();
                app.refreshDataPreview();
                app.logStatus('Imported fitted signal/vib_only into Data for ROI analysis, processing or export.');
            catch ME
                app.showError(ME);
            end
        end

        function clearFitCaches(app)
            app.CurrentResult = [];
            app.CubeResult = [];
            app.Stage1CubeResult = [];
            app.Stage2CubeResult = [];
            app.Stage1SingleResult = [];
            app.Stage2SingleResult = [];
            app.SelectedROICube = [];
            app.SelectedROISpectrum = [];
            app.SelectedROIRange = [];
        end

        function runDataIOAction(app)
            switch app.DataIOActionDropDown.Value
                case 'Import MAT'
                    app.onLoadData();
                case 'Import result'
                    app.onLoadResult();
                case 'Export data'
                    app.exportData();
                case 'Export ROI'
                    app.exportROI();
                case 'Export current result'
                    app.exportResult('current');
            end
        end

        function runDataToolAction(app)
            switch app.DataToolDropDown.Value
                case 'NMF'
                    app.openNMFWindow();
                case 'STA correction'
                    app.applySTACorrection();
                case 'Background (interactive)'
                    app.applyTraditionalBackground(true);
                case 'Separate IF mode'
                    app.openSeparateIFMode();
                case 'MLLS residual'
                    app.openMLLSWindow();
                case 'Peak fitting'
                    app.openPeakFitting();
                case 'Drift correction'
                    app.openDriftCorrection();
                case 'Phonon intensity'
                    app.openPMAnalysisWindow();
                case 'Stack spectra'
                    app.openStackSpectraTool();
            end
        end

        function applySTACorrection(app)
            if ~app.hasExternalTool('CorrectSTA', 'STA correction')
                return;
            end
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            try
                out = CorrectSTA(app.makeEELSObject(app.Data, app.Ene));
                app.Data = double(out.data);
                app.Ene = double(out.ene(:));
                app.clearFitCaches();
                app.setPixelLimits();
                app.refreshDataPreview();
                app.logStatus('STA correction complete; CorrectSTA updated the current data.');
            catch ME
                app.showError(ME);
            end
        end

        function openSeparateIFMode(app)
            if ~app.hasExternalTool('SeparateIFmode_Interactive', 'Separate IF mode')
                return;
            end
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            try
                SeparateIFmode_Interactive(app.makeEELSObject(app.Data, app.Ene));
                app.logStatus('Opened the Separate IF mode tool.');
            catch ME
                app.showError(ME);
            end
        end

        function openPeakFitting(app)
            if ~app.hasExternalTool('PeakFitting', 'Peak fitting')
                return;
            end
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            try
                PeakFitting(app.makeEELSObject(app.Data, app.Ene));
                app.logStatus('Opened the PeakFitting tool.');
            catch ME
                app.showError(ME);
            end
        end

        function openDriftCorrection(app)
            if ~app.hasExternalTool('CorrDrift_Custom', 'Interactive drift correction')
                return;
            end
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            try
                out = CorrDrift_Custom(app.makeEELSObject(app.Data, app.Ene));
                if exist('out', 'var') && ((isstruct(out) && isfield(out, 'data') && isfield(out, 'ene')) || ...
                        (isobject(out) && isprop(out, 'data') && isprop(out, 'ene')))
                    app.Data = double(out.data);
                    app.Ene = double(out.ene(:));
                    app.clearFitCaches();
                    app.setPixelLimits();
                    app.refreshDataPreview();
                end
                app.logStatus('Drift tool completed. Data updated if a drift file was selected.');
            catch ME
                app.showError(ME);
            end
        end

        function openNMFWindow(app)
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            fig = uifigure('Name', 'NMF decomposition', 'Position', [120 120 1060 680], 'Color', [0.94 0.97 1.00]);
            gl = uigridlayout(fig, [3 4]);
            gl.RowHeight = {48, '1x', '1x'};
            gl.ColumnWidth = {180, 180, 160, '1x'};
            gl.Padding = [12 12 12 12];
            uilabel(gl, 'Text', 'Components', 'FontWeight', 'bold');
            rankEdit = uieditfield(gl, 'numeric', 'Value', 3, 'Limits', [1 20]);
            uilabel(gl, 'Text', 'Uses the energy window from the main GUI', 'FontColor', [0.1 0.25 0.42]);
            axMap = uiaxes(gl); axMap.Layout.Row = 2; axMap.Layout.Column = [1 2]; app.styleAxes(axMap);
            axSpec = uiaxes(gl); axSpec.Layout.Row = [2 3]; axSpec.Layout.Column = [3 4]; app.styleAxes(axSpec);
            axCoef = uiaxes(gl); axCoef.Layout.Row = 3; axCoef.Layout.Column = [1 2]; app.styleAxes(axCoef);
            result = struct();
            b = app.techButton(gl, 'Run NMF', @(~, ~)runNMF()); b.Layout.Row = 1; b.Layout.Column = 4;
            b2 = app.techButton(gl, 'Export NMF', @(~, ~)exportNMF()); b2.Layout.Row = 1; b2.Layout.Column = 3;

            function runNMF()
                win = app.energyWindow();
                mask = app.Ene >= win(1) & app.Ene <= win(2);
                if nnz(mask) < 3, mask = true(size(app.Ene)); end
                [ny, nx, ~] = size(app.Data);
                X = reshape(app.Data(:, :, mask), ny * nx, nnz(mask));
                X(~isfinite(X)) = 0;
                X = X - min(X(:));
                k = max(1, round(rankEdit.Value));
                if exist('nnmf', 'file') == 2
                    [W, H] = nnmf(max(X, 0), k, 'replicates', 3, 'algorithm', 'mult');
                else
                    [U, S, V] = svds(double(X), k);
                    W = max(U * S, 0);
                    H = max(V', 0);
                end
                result = struct('W', W, 'H', H, 'energy', app.Ene(mask), 'rank', k, ...
                    'maps', reshape(W, ny, nx, k));
                imagesc(axMap, result.maps(:, :, 1)); axis(axMap, 'image'); colorbar(axMap);
                title(axMap, 'NMF component 1 map');
                plot(axSpec, result.energy, H', 'LineWidth', 1.2); grid(axSpec, 'on');
                xlabel(axSpec, 'Energy loss (meV)'); ylabel(axSpec, 'Component spectrum');
                title(axSpec, 'NMF spectra');
                bar(axCoef, sum(W, 1)); grid(axCoef, 'on'); title(axCoef, 'Component weights');
            end
            function exportNMF()
                if isempty(fieldnames(result)), runNMF(); end
                assignin('base', 'eels_nmf_result', result);
                app.logStatus('Exported eels_nmf_result to the base workspace.');
            end
        end

        function openMLLSWindow(app)
            if isempty(app.Data), app.logStatus('Load a datacube first.'); return; end
            fig = uifigure('Name', 'ROI / Reference MLLS residual', 'Position', [160 150 1120 700], 'Color', [0.94 0.97 1.00]);
            gl = uigridlayout(fig, [3 4]);
            gl.RowHeight = {48, '1x', '1x'};
            gl.ColumnWidth = {180, 220, 180, '1x'};
            vars = evalin('base', 'who');
            uilabel(gl, 'Text', 'Reference matrix variable (energy x nRef)', 'FontWeight', 'bold');
            dd = uidropdown(gl, 'Items', [{'Current ROI spectrum'}; vars(:)]);
            uilabel(gl, 'Text', 'Uses the energy window from the main GUI');
            axMap = uiaxes(gl); axMap.Layout.Row = [2 3]; axMap.Layout.Column = [1 2]; app.styleAxes(axMap);
            axSpec = uiaxes(gl); axSpec.Layout.Row = [2 3]; axSpec.Layout.Column = [3 4]; app.styleAxes(axSpec);
            result = struct();
            b = app.techButton(gl, 'Run MLLS', @(~, ~)runMLLS()); b.Layout.Row = 1; b.Layout.Column = 4;
            b2 = app.techButton(gl, 'Export MLLS', @(~, ~)exportMLLS()); b2.Layout.Row = 1; b2.Layout.Column = 3;
            function runMLLS()
                if isempty(dd.Value), return; end
                if strcmp(dd.Value, 'Current ROI spectrum')
                    if isempty(app.SelectedROISpectrum)
                        app.logStatus('Select an ROI on the Data tab or choose a reference matrix from the workspace.');
                        return;
                    end
                    ref = app.SelectedROISpectrum(:);
                else
                    ref = evalin('base', dd.Value);
                end
                ref = double(ref);
                if size(ref, 1) ~= numel(app.Ene) && size(ref, 2) == numel(app.Ene)
                    ref = ref.';
                end
                win = app.energyWindow();
                mask = app.Ene >= win(1) & app.Ene <= win(2);
                if nnz(mask) < 3, mask = true(size(app.Ene)); end
                A = ref(mask, :);
                [ny, nx, ~] = size(app.Data);
                residual_map = NaN(ny, nx);
                coeff = NaN(ny, nx, size(A, 2));
                residual_cube = NaN(ny, nx, nnz(mask));
                for yy = 1:ny
                    for xx = 1:nx
                        y = squeeze(app.Data(yy, xx, mask));
                        if exist('lsqnonneg', 'file') == 2
                            c = lsqnonneg(A, y(:));
                        else
                            c = max(A \ y(:), 0);
                        end
                        r = y(:) - A * c;
                        residual_map(yy, xx) = sqrt(mean(r .^ 2, 'omitnan'));
                        coeff(yy, xx, :) = c;
                        residual_cube(yy, xx, :) = r;
                    end
                end
                result = struct('residual_map', residual_map, 'coeff', coeff, ...
                    'residual_cube', residual_cube, 'energy', app.Ene(mask), 'reference_name', dd.Value);
                imagesc(axMap, residual_map); axis(axMap, 'image'); colorbar(axMap); title(axMap, 'MLLS residual RMS');
                [~, idx] = max(residual_map(:)); [yy, xx] = ind2sub([ny, nx], idx);
                plot(axSpec, result.energy, squeeze(residual_cube(yy, xx, :)), 'b-'); grid(axSpec, 'on');
                title(axSpec, sprintf('Max residual spectrum (%d,%d)', yy, xx));
            end
            function exportMLLS()
                if isempty(fieldnames(result)), runMLLS(); end
                assignin('base', 'eels_mlls_residual_result', result);
                app.logStatus('Exported eels_mlls_residual_result to the base workspace.');
            end
        end

        function openPMAnalysisWindow(app)
            defaultSource = 'current';
            if ~isempty(app.Stage2CubeResult) || ~isempty(app.Stage2SingleResult) || ...
                    app.cuberesHasStage2(app.CubeResult) || app.resultHasStage2(app.CurrentResult)
                defaultSource = 'stage2';
            elseif ~isempty(app.Stage1CubeResult) || ~isempty(app.Stage1SingleResult)
                defaultSource = 'stage1';
            end

            pm_pos = app.responsiveFigurePosition(app.LayoutMode);
            fig = uifigure('Name', 'Sparse phonon intensity / Loss-gain DOS-like analysis', ...
                'Position', pm_pos, 'Color', [0.93 0.97 1.00]);
            gl = uigridlayout(fig, [3 1]);
            if app.IsCompact
                gl.RowHeight = {220, '1x', 118};
            else
                gl.RowHeight = {150, '1x', 64};
            end
            gl.ColumnWidth = {'1x'};
            gl.Padding = [12 12 12 12];
            gl.RowSpacing = 10;
            gl.BackgroundColor = [0.93 0.97 1.00];

            panel = uipanel(gl, 'Title', 'Intensity interpolation, selection, NMF and clustering', ...
                'FontWeight', 'bold', 'ForegroundColor', [0.05 0.18 0.32], ...
                'BackgroundColor', [0.88 0.94 1.00], 'BorderType', 'none');
            panel.Layout.Row = 1; panel.Layout.Column = 1;
            app.enableScrollable(panel);
            if app.IsCompact
                pg = uigridlayout(panel, [5 6]);
                pg.RowHeight = {28, 28, 28, 28, 28};
                pg.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x'};
            else
                pg = uigridlayout(panel, [3 12]);
                pg.RowHeight = {28, 28, 28};
                pg.ColumnWidth = {72, 105, 92, 86, 86, 70, 92, 86, 86, 70, 86, '1x'};
            end
            pg.Padding = [10 8 10 8];
            pg.BackgroundColor = [0.88 0.94 1.00];

            lab = app.techLabel(pg, 'Intensity source'); lab.Layout.Row = 1; lab.Layout.Column = 1;
            sourceDD = uidropdown(pg, 'Items', {'stage2', 'stage1', 'current'}, 'Value', defaultSource);
            sourceDD.Layout.Row = 1; sourceDD.Layout.Column = 2;
            lab = app.techLabel(pg, 'Interpolation range'); lab.Layout.Row = 1; lab.Layout.Column = 3;
            emMinEdit = uieditfield(pg, 'numeric', 'Value', -250);
            emMinEdit.Layout.Row = 1; emMinEdit.Layout.Column = 4;
            emMaxEdit = uieditfield(pg, 'numeric', 'Value', 250);
            emMaxEdit.Layout.Row = 1; emMaxEdit.Layout.Column = 5;
            if app.IsCompact
                lab = app.techLabel(pg, 'meV'); lab.Layout.Row = 1; lab.Layout.Column = 6;
                lab = app.techLabel(pg, 'Step'); lab.Layout.Row = 2; lab.Layout.Column = 1;
            else
                lab = app.techLabel(pg, 'Step'); lab.Layout.Row = 1; lab.Layout.Column = 6;
            end
            stepEdit = uieditfield(pg, 'numeric', 'Value', 1, 'Limits', [0.1 20]);
            if app.IsCompact
                stepEdit.Layout.Row = 2; stepEdit.Layout.Column = 2;
                lab = app.techLabel(pg, 'Select window'); lab.Layout.Row = 2; lab.Layout.Column = 3;
            else
                stepEdit.Layout.Row = 1; stepEdit.Layout.Column = 7;
                lab = app.techLabel(pg, 'Select window'); lab.Layout.Row = 1; lab.Layout.Column = 8;
            end
            selMinEdit = uieditfield(pg, 'numeric', 'Value', 40);
            if app.IsCompact
                selMinEdit.Layout.Row = 2; selMinEdit.Layout.Column = 4;
            else
                selMinEdit.Layout.Row = 1; selMinEdit.Layout.Column = 9;
            end
            selMaxEdit = uieditfield(pg, 'numeric', 'Value', 120);
            if app.IsCompact
                selMaxEdit.Layout.Row = 2; selMaxEdit.Layout.Column = 5;
                lab = app.techLabel(pg, 'NMF'); lab.Layout.Row = 2; lab.Layout.Column = 6;
            else
                selMaxEdit.Layout.Row = 1; selMaxEdit.Layout.Column = 10;
                lab = app.techLabel(pg, 'NMF'); lab.Layout.Row = 1; lab.Layout.Column = 11;
            end
            rankEdit = uieditfield(pg, 'numeric', 'Value', 3, 'Limits', [1 20]);
            if app.IsCompact
                rankEdit.Layout.Row = 3; rankEdit.Layout.Column = 1;
                lab = app.techLabel(pg, 'Clusters'); lab.Layout.Row = 3; lab.Layout.Column = 2;
            else
                rankEdit.Layout.Row = 1; rankEdit.Layout.Column = 12;
                lab = app.techLabel(pg, 'Clusters'); lab.Layout.Row = 2; lab.Layout.Column = 1;
            end
            clusterEdit = uieditfield(pg, 'numeric', 'Value', 3, 'Limits', [2 20]);
            if app.IsCompact
                clusterEdit.Layout.Row = 3; clusterEdit.Layout.Column = 3;
                lab = app.techLabel(pg, 'Display smoothing'); lab.Layout.Row = 3; lab.Layout.Column = 4;
            else
                clusterEdit.Layout.Row = 2; clusterEdit.Layout.Column = 2;
                lab = app.techLabel(pg, 'Display smoothing'); lab.Layout.Row = 2; lab.Layout.Column = 3;
            end
            smoothEdit = uieditfield(pg, 'numeric', 'Value', 3, 'Limits', [1 51]);
            if app.IsCompact
                smoothEdit.Layout.Row = 3; smoothEdit.Layout.Column = 5;
            else
                smoothEdit.Layout.Row = 2; smoothEdit.Layout.Column = 4;
            end
            smoothCheck = uicheckbox(pg, 'Text', 'Smooth display/export only; preserve source intensity', 'Value', true);
            if app.IsCompact
                smoothCheck.Layout.Row = 3; smoothCheck.Layout.Column = 6;
            else
                smoothCheck.Layout.Row = 2; smoothCheck.Layout.Column = [5 8];
            end
            sparseNote = uilabel(pg, 'Text', 'Attention-selected pairs constrain temperature; soft phonon intensities reconstruct other energies. Hierarchical ZLP protects shape.', ...
                'FontColor', [0.55 0.16 0.02]);
            if app.IsCompact
                sparseNote.Layout.Row = 4; sparseNote.Layout.Column = [1 6];
            else
                sparseNote.Layout.Row = 2; sparseNote.Layout.Column = [9 12];
            end
            note = uilabel(pg, 'Text', 'Signed source intensity: negative energy for gain/anti-Stokes, positive for loss/Stokes. This DOS-like representation is not an absolute phonon DOS calibration.', ...
                'FontColor', [0.08 0.25 0.42]);
            if app.IsCompact
                note.Layout.Row = 5; note.Layout.Column = [1 6];
            else
                note.Layout.Row = 3; note.Layout.Column = [1 12];
            end

            tabs = uitabgroup(gl);
            tabs.Layout.Row = 2; tabs.Layout.Column = 1;
            tabOverview = uitab(tabs, 'Title', 'Overview / selection');
            tabNMF = uitab(tabs, 'Title', 'Intensity NMF');
            tabCluster = uitab(tabs, 'Title', 'Clusters');

            og = uigridlayout(tabOverview, [2 2]);
            og.RowHeight = {'1x', 170};
            og.ColumnWidth = {'1x', '1.25x'};
            og.Padding = [10 10 10 10];
            og.RowSpacing = 10; og.ColumnSpacing = 12;
            axMap = uiaxes(og); axMap.Layout.Row = [1 2]; axMap.Layout.Column = 1; app.styleAxes(axMap);
            axSpec = uiaxes(og); axSpec.Layout.Row = 1; axSpec.Layout.Column = 2; app.styleAxes(axSpec);
            axSummary = uiaxes(og); axSummary.Layout.Row = 2; axSummary.Layout.Column = 2; app.styleAxes(axSummary);

            ng = uigridlayout(tabNMF, [1 2]);
            ng.ColumnWidth = {'1x', '1.2x'}; ng.Padding = [10 10 10 10]; ng.ColumnSpacing = 12;
            axNMFMap = uiaxes(ng); axNMFMap.Layout.Row = 1; axNMFMap.Layout.Column = 1; app.styleAxes(axNMFMap);
            axNMFSpectrum = uiaxes(ng); axNMFSpectrum.Layout.Row = 1; axNMFSpectrum.Layout.Column = 2; app.styleAxes(axNMFSpectrum);

            cg = uigridlayout(tabCluster, [1 2]);
            cg.ColumnWidth = {'1x', '1.2x'}; cg.Padding = [10 10 10 10]; cg.ColumnSpacing = 12;
            axCluster = uiaxes(cg); axCluster.Layout.Row = 1; axCluster.Layout.Column = 1; app.styleAxes(axCluster);
            axClusterSummary = uiaxes(cg); axClusterSummary.Layout.Row = 1; axClusterSummary.Layout.Column = 2; app.styleAxes(axClusterSummary);

            pmres = struct();
            nmfres = struct(); %#ok<NASGU>
            clusterres = struct(); %#ok<NASGU>
            modeMap = []; %#ok<NASGU>

            bp = uipanel(gl, 'BorderType', 'none', 'BackgroundColor', [0.88 0.94 1.00]);
            bp.Layout.Row = 3; bp.Layout.Column = 1;
            if app.IsCompact
                bg = uigridlayout(bp, [2 3]);
                bg.RowHeight = {40, 40};
                bg.ColumnWidth = {'1x', '1x', '1x'};
            else
                bg = uigridlayout(bp, [1 6]);
                bg.ColumnWidth = {'1x', '1x', '1x', '1x', '1x', '1x'};
            end
            bg.Padding = [8 7 8 7]; bg.ColumnSpacing = 10;
            if app.IsCompact
                b = app.techButton(bg, 'Build interpolation', @(~, ~)buildPM()); b.Layout.Row = 1; b.Layout.Column = 1;
                b = app.techButton(bg, 'Export interpolation', @(~, ~)exportPM()); b.Layout.Row = 1; b.Layout.Column = 2;
                b = app.techButton(bg, 'Export selection', @(~, ~)exportSelectedPM()); b.Layout.Row = 1; b.Layout.Column = 3;
                b = app.techButton(bg, 'Integrated mode map', @(~, ~)plotModeMap()); b.Layout.Row = 2; b.Layout.Column = 1;
                b = app.techButton(bg, 'Intensity NMF', @(~, ~)runPMNMF()); b.Layout.Row = 2; b.Layout.Column = 2;
                b = app.techButton(bg, 'Run clustering', @(~, ~)runPMCluster()); b.Layout.Row = 2; b.Layout.Column = 3;
            else
                b = app.techButton(bg, 'Build interpolation', @(~, ~)buildPM()); b.Layout.Row = 1; b.Layout.Column = 1;
                b = app.techButton(bg, 'Export interpolation', @(~, ~)exportPM()); b.Layout.Row = 1; b.Layout.Column = 2;
                b = app.techButton(bg, 'Export selection', @(~, ~)exportSelectedPM()); b.Layout.Row = 1; b.Layout.Column = 3;
                b = app.techButton(bg, 'Integrated mode map', @(~, ~)plotModeMap()); b.Layout.Row = 1; b.Layout.Column = 4;
                b = app.techButton(bg, 'Intensity NMF', @(~, ~)runPMNMF()); b.Layout.Row = 1; b.Layout.Column = 5;
                b = app.techButton(bg, 'Run clustering', @(~, ~)runPMCluster()); b.Layout.Row = 1; b.Layout.Column = 6;
            end

            buildPM();

            function buildPM()
                try
                    range = sort([emMinEdit.Value, emMaxEdit.Value]);
                    pmres = app.buildPMInterpolatedCube(sourceDD.Value, range, stepEdit.Value);
                    cla(axSpec);
                    meanPM = displayMeanPM(pmres.mean_PM);
                    plot(axSpec, pmres.axis_meV, meanPM, 'k-', 'LineWidth', 1.5); hold(axSpec, 'on');
                    xline(axSpec, 0, 'Color', [0.35 0.35 0.35]);
                    xline(axSpec, selMinEdit.Value, 'Color', [0.1 0.45 0.95], 'LineStyle', '--');
                    xline(axSpec, selMaxEdit.Value, 'Color', [0.1 0.45 0.95], 'LineStyle', '--');
                    hold(axSpec, 'off'); grid(axSpec, 'on');
                    xlabel(axSpec, 'Signed mode energy (meV), gain < 0, loss > 0');
                    ylabel(axSpec, 'Mean source intensity');
                    title(axSpec, sprintf('Mean intensity, source=%s, valid pixels=%d', pmres.source, nnz(pmres.valid_mask)));
                    plotModeMap();
                    cla(axSummary);
                    bar(axSummary, [nnz(pmres.axis_meV < 0), nnz(pmres.axis_meV == 0), nnz(pmres.axis_meV > 0)]);
                    set(axSummary, 'XTickLabel', {'gain', 'zero', 'loss'});
                    ylabel(axSummary, 'Interpolation bins'); title(axSummary, 'Signed intensity coverage'); grid(axSummary, 'on');
                catch ME
                    app.showError(ME);
                end
            end

            function plotModeMap()
                if isempty(fieldnames(pmres)), buildPM(); end
                if isempty(fieldnames(pmres)), return; end
                [modeMap, useAxis] = selectedPMIntegral(pmres);
                imagesc(axMap, modeMap); axis(axMap, 'image'); colorbar(axMap);
                title(axMap, sprintf('Integrated intensity [%g, %g] meV', min(useAxis), max(useAxis)));
                xlabel(axMap, 'x'); ylabel(axMap, 'y');
                assignin('base', 'eels_pm_mode_map', struct('map', modeMap, ...
                    'axis_meV', useAxis, 'source', pmres.source));
                app.logStatus('Integrated intensity map created and exported.');
            end

            function exportPM()
                if isempty(fieldnames(pmres)), buildPM(); end
                assignin('base', 'eels_pm_interpolated', pmres);
                app.logStatus('Exported signed mode energies, intensity cube, mean spectrum and valid-pixel mask.');
            end

            function exportSelectedPM()
                if isempty(fieldnames(pmres)), buildPM(); end
                range = sort([selMinEdit.Value, selMaxEdit.Value]);
                mask = pmres.axis_meV >= range(1) & pmres.axis_meV <= range(2);
                if ~any(mask)
                    app.logStatus('Empty intensity energy window. Adjust the selection.');
                    return;
                end
                selected = pmres;
                selected.axis_meV = pmres.axis_meV(mask);
                selected.pm_cube = pmres.pm_cube(:, :, mask);
                selected.mean_PM = pmres.mean_PM(mask);
                selected.selection_range = range;
                assignin('base', 'eels_pm_selected', selected);
                app.logStatus('Selected intensity window exported for spectral or spatial analysis.');
            end

            function runPMNMF()
                if isempty(fieldnames(pmres)), buildPM(); end
                X = reshape(pmres.pm_cube, [], numel(pmres.axis_meV));
                valid = pmres.valid_mask(:) & any(isfinite(X) & X > 0, 2);
                X = X(valid, :);
                X(~isfinite(X)) = 0;
                X = max(X, 0);
                if isempty(X)
                    app.logStatus('No valid intensity data for NMF.');
                    return;
                end
                k = min(max(1, round(rankEdit.Value)), min(size(X)));
                if exist('nnmf', 'file') == 2
                    [W, H] = nnmf(X, k, 'replicates', 3, 'algorithm', 'mult');
                else
                    [U, S, V] = svds(double(X), k);
                    W = max(U * S, 0);
                    H = max(V', 0);
                end
                maps_flat = NaN(numel(pmres.valid_mask), k);
                maps_flat(valid, :) = W;
                maps = reshape(maps_flat, pmres.map_size(1), pmres.map_size(2), k);
                nmfres = struct('W', W, 'H', H, 'maps', maps, 'axis_meV', pmres.axis_meV, ...
                    'source', pmres.source, 'rank', k, 'valid_mask', pmres.valid_mask);
                imagesc(axNMFMap, maps(:, :, 1)); axis(axNMFMap, 'image'); colorbar(axNMFMap);
                title(axNMFMap, 'Intensity NMF component 1 map');
                plot(axNMFSpectrum, pmres.axis_meV, H', 'LineWidth', 1.2); grid(axNMFSpectrum, 'on');
                xlabel(axNMFSpectrum, 'Signed mode energy (meV)'); ylabel(axNMFSpectrum, 'Component source intensity');
                title(axNMFSpectrum, 'Intensity NMF component spectra');
                assignin('base', 'eels_pm_nmf_result', nmfres);
                app.logStatus('Intensity NMF complete and exported.');
            end

            function runPMCluster()
                if isempty(fieldnames(pmres)), buildPM(); end
                X = reshape(pmres.pm_cube, [], numel(pmres.axis_meV));
                valid = pmres.valid_mask(:) & any(isfinite(X) & X > 0, 2);
                X = X(valid, :);
                X(~isfinite(X)) = 0;
                X = max(X, 0);
                if isempty(X)
                    app.logStatus('No valid intensity data for clustering.');
                    return;
                end
                row_sum = sum(X, 2);
                Xn = X ./ max(row_sum, eps);
                Xn = Xn - mean(Xn, 1);
                ncomp = min([6, size(Xn, 1) - 1, size(Xn, 2)]);
                if ncomp >= 1
                    [U, S, ~] = svds(double(Xn), ncomp);
                    feat = U * S;
                else
                    feat = Xn;
                end
                k = min(max(2, round(clusterEdit.Value)), size(feat, 1));
                if exist('kmeans', 'file') == 2
                    labels = kmeans(feat, k, 'Replicates', 5, 'MaxIter', 300, 'Display', 'off');
                else
                    labels = localKMeans(feat, k, 80);
                end
                label_flat = NaN(numel(pmres.valid_mask), 1);
                label_flat(valid) = labels;
                label_map = reshape(label_flat, pmres.map_size(1), pmres.map_size(2));
                clusterres = struct('label_map', label_map, 'labels', labels, 'features', feat, ...
                    'axis_meV', pmres.axis_meV, 'source', pmres.source, 'k', k, ...
                    'valid_mask', pmres.valid_mask);
                imagesc(axCluster, label_map); axis(axCluster, 'image'); colorbar(axCluster);
                title(axCluster, sprintf('Intensity clusters, k=%d', k));
                cla(axClusterSummary); hold(axClusterSummary, 'on');
                for cc = 1:k
                    plot(axClusterSummary, pmres.axis_meV, mean(X(labels == cc, :), 1, 'omitnan'), 'LineWidth', 1.2);
                end
                hold(axClusterSummary, 'off'); grid(axClusterSummary, 'on');
                xlabel(axClusterSummary, 'Signed mode energy (meV)'); ylabel(axClusterSummary, 'Cluster mean source intensity');
                title(axClusterSummary, 'Mean intensity spectrum per cluster');
                assignin('base', 'eels_pm_cluster_result', clusterres);
                app.logStatus('Intensity clustering complete and exported.');
            end

            function y = displayMeanPM(y)
                y = y(:).';
                if ~logical(smoothCheck.Value)
                    return;
                end
                w = max(1, round(smoothEdit.Value));
                if mod(w, 2) == 0
                    w = w + 1;
                end
                if w > 1 && numel(y) >= w
                    y = movmean(y, w, 'omitnan');
                end
            end

            function [map, useAxis] = selectedPMIntegral(src)
                range = sort([selMinEdit.Value, selMaxEdit.Value]);
                mask = src.axis_meV >= range(1) & src.axis_meV <= range(2);
                if ~any(mask)
                    mask = true(size(src.axis_meV));
                end
                useAxis = src.axis_meV(mask);
                block = src.pm_cube(:, :, mask);
                block(~isfinite(block)) = 0;
                if numel(useAxis) > 1
                    scale = median(abs(diff(useAxis)), 'omitnan');
                else
                    scale = 1;
                end
                map = sum(block, 3) .* scale;
                map(~src.valid_mask) = NaN;
            end

            function labels = localKMeans(feat, k, maxIter)
                n = size(feat, 1);
                if n == 0
                    labels = [];
                    return;
                end
                [~, order] = sort(feat(:, 1), 'ascend');
                seed_idx = order(unique(round(linspace(1, n, k))));
                while numel(seed_idx) < k
                    seed_idx(end + 1) = order(end); %#ok<AGROW>
                end
                centers = feat(seed_idx(1:k), :);
                labels = ones(n, 1);
                for iter = 1:maxIter
                    dist = zeros(n, k);
                    for cc = 1:k
                        delta = feat - centers(cc, :);
                        dist(:, cc) = sum(delta .^ 2, 2);
                    end
                    [~, newLabels] = min(dist, [], 2);
                    if iter > 1 && isequal(newLabels, labels)
                        break;
                    end
                    labels = newLabels;
                    for cc = 1:k
                        if any(labels == cc)
                            centers(cc, :) = mean(feat(labels == cc, :), 1);
                        end
                    end
                end
            end
        end

        function pmres = buildPMInterpolatedCube(app, source, em_range, em_step)
            cuberes = app.pmSourceCuberes(source);
            if isempty(cuberes) || ~isfield(cuberes, 'results') || isempty(cuberes.results)
                error('No %s datacube/spectrum result. Run or load Stage1/Stage2 first.', char(source));
            end
            em_range = sort(double(em_range(:).'));
            em_step = max(double(em_step), eps);
            signed_axis = em_range(1):em_step:em_range(2);
            if isempty(signed_axis)
                signed_axis = linspace(em_range(1), em_range(2), 101);
            end
            if isfield(cuberes, 'R2_map') && ~isempty(cuberes.R2_map)
                [ny, nx] = size(cuberes.R2_map);
            elseif isfield(cuberes, 'cube_info') && isfield(cuberes.cube_info, 'size') && numel(cuberes.cube_info.size) >= 2
                ny = cuberes.cube_info.size(1);
                nx = cuberes.cube_info.size(2);
            else
                ny = 1;
                nx = numel(cuberes.results);
            end
            nE = numel(signed_axis);
            pm_cube = NaN(ny, nx, nE);
            valid_mask = false(ny, nx);
            stage_used = strings(ny, nx);
            for ii = 1:numel(cuberes.results)
                res = cuberes.results{ii};
                if ~(isstruct(res) && ~isfield(res, 'error'))
                    continue;
                end
                [em, pm_s, pm_a, stageLabel] = app.extractPMVectors(res, source);
                if isempty(em)
                    continue;
                end
                [yy, xx] = app.resultPixelPosition(cuberes, res, ii, ny, nx);
                if yy < 1 || yy > ny || xx < 1 || xx > nx
                    continue;
                end
                em = abs(em(:));
                pm_s = max(pm_s(:), 0);
                pm_a = max(pm_a(:), 0);
                n = min([numel(em), numel(pm_s), numel(pm_a)]);
                em = em(1:n);
                pm_s = pm_s(1:n);
                pm_a = pm_a(1:n);
                signed_em = [-flipud(em); em];
                signed_pm = [flipud(pm_a); pm_s];
                good = isfinite(signed_em) & isfinite(signed_pm);
                signed_em = signed_em(good);
                signed_pm = signed_pm(good);
                if numel(signed_em) < 2
                    continue;
                end
                [signed_em, ord] = sort(signed_em);
                signed_pm = signed_pm(ord);
                key = round(signed_em .* 1e6) ./ 1e6;
                [uniq_em, ~, ic] = unique(key);
                uniq_pm = accumarray(ic, signed_pm, [], @sum);
                vals = interp1(uniq_em, uniq_pm, signed_axis, 'linear', 0);
                pm_cube(yy, xx, :) = vals;
                valid_mask(yy, xx) = true;
                stage_used(yy, xx) = string(stageLabel);
            end
            X = reshape(pm_cube, [], nE);
            X(~valid_mask(:), :) = NaN;
            mean_pm = mean(X, 1, 'omitnan');
            pmres = struct('source', char(source), 'axis_meV', signed_axis(:).', ...
                'pm_cube', pm_cube, 'mean_PM', mean_pm(:).', 'valid_mask', valid_mask, ...
                'stage_used', stage_used, 'map_size', [ny, nx], ...
                'description', ['Signed phonon-mode source-intensity cube: negative mode energy is gain/anti-Stokes, ', ...
                'positive mode energy is loss/Stokes. This DOS-like interpolation is not a new fit or an absolutely calibrated DOS.']);
        end

        function cuberes = pmSourceCuberes(app, source)
            source = lower(char(source));
            cuberes = [];
            single = [];
            switch source
                case 'stage2'
                    cuberes = app.firstUsableCuberes({app.Stage2CubeResult, app.CubeResult}, true);
                    single = app.firstUsableSingle({app.Stage2SingleResult, app.CurrentResult}, true);
                case 'stage1'
                    cuberes = app.firstUsableCuberes({app.Stage1CubeResult, app.CubeResult, app.Stage2CubeResult}, false);
                    single = app.firstUsableSingle({app.Stage1SingleResult, app.CurrentResult, app.Stage2SingleResult}, false);
                otherwise
                    cuberes = app.firstUsableCuberes({app.CubeResult, app.Stage2CubeResult, app.Stage1CubeResult}, false);
                    single = app.firstUsableSingle({app.CurrentResult, app.Stage2SingleResult, app.Stage1SingleResult}, false);
            end
            if isempty(cuberes) && ~isempty(single) && isstruct(single) && ~isfield(single, 'error')
                if ~isfield(single, 'pixel_y'), single.pixel_y = 1; end
                if ~isfield(single, 'pixel_x'), single.pixel_x = 1; end
                cuberes = struct('results', {{single}}, 'sample_pixels', 1, ...
                    'R2_map', NaN(1, 1), 'cube_info', struct('size', [1, 1]));
            end
        end

        function cuberes = firstUsableCuberes(app, candidates, requireStage2)
            cuberes = [];
            for ii = 1:numel(candidates)
                c = candidates{ii};
                if ~(isstruct(c) && isfield(c, 'results') && ~isempty(c.results))
                    continue;
                end
                if requireStage2 && ~app.cuberesHasStage2(c)
                    continue;
                end
                cuberes = c;
                return;
            end
        end

        function single = firstUsableSingle(app, candidates, requireStage2)
            single = [];
            for ii = 1:numel(candidates)
                c = candidates{ii};
                if ~(isstruct(c) && ~isfield(c, 'error'))
                    continue;
                end
                if requireStage2 && ~app.resultHasStage2(c)
                    continue;
                end
                single = c;
                return;
            end
        end

        function [em, pm_s, pm_a, stageLabel] = extractPMVectors(app, res, source)
            em = [];
            pm_s = [];
            pm_a = [];
            stageLabel = 'stage1';
            wantStage2 = strcmpi(char(source), 'stage2') || ...
                (strcmpi(char(source), 'current') && app.resultHasStage2(res));
            if wantStage2 && isfield(res, 'stage2') && isstruct(res.stage2) && ...
                    isfield(res.stage2, 'EM_final') && isfield(res.stage2, 'PM_stokes') && ...
                    isfield(res.stage2, 'PM_antistokes') && ~isempty(res.stage2.EM_final)
                em = res.stage2.EM_final(:);
                pm_s = res.stage2.PM_stokes(:);
                pm_a = res.stage2.PM_antistokes(:);
                stageLabel = 'stage2';
            elseif isfield(res, 'EM_final') && isfield(res, 'PM_stokes') && ...
                    isfield(res, 'PM_antistokes') && ~isempty(res.EM_final)
                em = res.EM_final(:);
                pm_s = res.PM_stokes(:);
                pm_a = res.PM_antistokes(:);
                stageLabel = 'stage1';
            end
        end

        function [yy, xx] = resultPixelPosition(~, cuberes, res, ii, ny, nx)
            if isfield(res, 'pixel_y') && isfield(res, 'pixel_x') && ...
                    isfinite(res.pixel_y) && isfinite(res.pixel_x)
                yy = round(res.pixel_y);
                xx = round(res.pixel_x);
            elseif isfield(cuberes, 'sample_pixels') && numel(cuberes.sample_pixels) >= ii
                [yy, xx] = ind2sub([ny, nx], cuberes.sample_pixels(ii));
            else
                yy = 1;
                xx = ii;
            end
        end

        function tf = resultHasStage2(~, res)
            tf = isstruct(res) && isfield(res, 'stage2') && isstruct(res.stage2) && ...
                isfield(res.stage2, 'status') && strcmp(res.stage2.status, 'ok') && ...
                isfield(res.stage2, 'EM_final') && isfield(res.stage2, 'PM_stokes') && ...
                isfield(res.stage2, 'PM_antistokes') && ~isempty(res.stage2.EM_final);
        end

        function tf = cuberesHasStage2(app, cuberes)
            tf = false;
            if ~(isstruct(cuberes) && isfield(cuberes, 'results') && ~isempty(cuberes.results)), return; end
            for ii = 1:numel(cuberes.results)
                if app.resultHasStage2(cuberes.results{ii})
                    tf = true;
                    return;
                end
            end
        end

        function tf = isReusableStage1Cube(app, cuberes)
            tf = isstruct(cuberes) && isfield(cuberes, 'results') && isfield(cuberes, 'sample_pixels') && ...
                ~isempty(cuberes.results) && ~app.cuberesHasStage2(cuberes);
        end

        function updateProgressDialog(app, dlg, done, total, elapsed, eta, med, label)
            if isempty(dlg) || ~isvalid(dlg), return; end
            if total > 0
                dlg.Value = min(1, max(0, done ./ total));
            end
            if isfinite(eta)
                eta_txt = sprintf('ETA %.1f min', eta ./ 60);
            else
                eta_txt = 'Estimating ETA';
            end
            if ~isfinite(med), med = NaN; end
            if app.shouldPauseDatacube(dlg)
                pause_txt = ' | Pause requested; save after the current pixel if checkpointing is enabled';
            else
                pause_txt = '';
            end
            dlg.Message = sprintf('%s: %d/%d | elapsed %.1f min | %.2f s/pixel | %s', ...
                char(label), done, total, elapsed ./ 60, med, eta_txt);
            dlg.Message = [dlg.Message, pause_txt];
            drawnow limitrate;
        end

        function plotRawSpectrum(app, yy, xx, ax)
            cla(ax);
            spec = squeeze(app.Data(yy, xx, :));
            plot(ax, app.Ene, spec, 'k-');
            grid(ax, 'on'); xlabel(ax, 'Energy loss (meV)'); ylabel(ax, 'Intensity');
            title(ax, sprintf('Raw spectrum (%d,%d)', yy, xx));
        end

        function refreshResultMap(app, stage2)
            if isempty(app.CubeResult), return; end
            if stage2
                ax = app.Stage2MapAxes; dd = app.Stage2Controls.mapDropDown;
            else
                ax = app.Stage1MapAxes; dd = app.Stage1Controls.mapDropDown;
            end
            [img, ttl] = app.mapData(dd.Value);
            imagesc(ax, img); axis(ax, 'image'); colorbar(ax); title(ax, ttl);
            h = ax.Children(1);
            h.HitTest = 'on';
            h.ButtonDownFcn = @(~, evt)app.onMapClick(evt, stage2);
            ax.ButtonDownFcn = @(~, evt)app.onMapClick(evt, stage2);
        end

        function [img, ttl] = mapData(app, name)
            cr = app.CubeResult;
            switch name
                case 'R2', img = cr.R2_map; ttl = 'R2 map';
                case 'NRMSE', img = cr.NRMSE_map; ttl = 'NRMSE map';
                case 'T raw', img = cr.T_map_raw; ttl = 'T raw (K)';
                case 'T data valid', img = app.optionalMap(cr, 'T_map_data_valid', cr.R2_map); ttl = 'T data-valid (K)';
                case 'T joint', img = app.optionalMap(cr, 'T_map_joint', cr.R2_map); ttl = 'ZLP-guided T, no T smoothing (K)';
                case 'T joint uncertainty', img = app.optionalMap(cr, 'T_uncertainty_joint_map', cr.R2_map); ttl = 'ZLP-guided T uncertainty (K)';
                case 'Spatial inferred', img = app.optionalMap(cr, 'T_spatial_inferred_mask', cr.R2_map); ttl = 'Spatial-inferred mask';
                case 'T repaired', img = cr.T_map_repaired; ttl = 'T repaired (K)';
                case 'T spatial refined', img = cr.T_map_spatial_refined; ttl = 'T spatial refined (K)';
                case 'T uncertainty', img = cr.T_uncertainty_map; ttl = 'T uncertainty (K)';
                case 'Reliability', img = cr.reliability_score_map; ttl = 'Reliability';
                case 'Attention fraction', img = cr.attention_score_map; ttl = 'Attention fraction';
                case 'Warning', img = cr.warning_map; ttl = 'Warning map';
                case 'Bad pixel', img = cr.bad_pixel_mask; ttl = 'Bad pixel mask';
                case 'Interface', img = cr.interface_mask; ttl = 'Interface mask';
                case 'Neighbor loss', img = cr.neighbor_loss_map; ttl = 'Neighbor loss';
                case 'Spectrum neighbor loss', img = cr.spectrum_neighbor_loss_map; ttl = 'Spectrum neighbor loss';
                case 'Temperature neighbor z', img = cr.temperature_neighbor_z_map; ttl = 'Temperature neighbor z';
                case 'Attention peak delta', img = cr.attention_peak_delta_map; ttl = 'Attention peak delta';
                case 'Signal peak delta', img = cr.signal_peak_delta_map; ttl = 'Signal peak delta';
                case 'ZLP guard start', img = cr.zlp_guard_start_map; ttl = 'ZLP guard start (meV)';
                case 'ZLP guard fraction', img = cr.zlp_guard_fraction_map; ttl = 'ZLP fraction in guard';
                case 'ZLP core width', img = cr.zlp_core_width_map; ttl = 'ZLP core width (meV)';
                case 'ZLP cluster', img = app.optionalMap(cr, 'cluster_label_map', cr.R2_map); ttl = 'ZLP cluster';
                case 'Cluster membership', img = app.optionalMap(cr, 'cluster_membership_map', cr.R2_map); ttl = 'Cluster membership';
                case 'Cluster outlier', img = app.optionalMap(cr, 'cluster_outlier_mask', cr.R2_map); ttl = 'Cluster outlier mask';
                case 'Template deviation', img = app.optionalMap(cr, 'template_deviation_map', cr.R2_map); ttl = 'Template deviation';
                case 'ZLP spatial loss', img = app.optionalMap(cr, 'zlp_spatial_loss_map', cr.R2_map); ttl = 'ZLP spatial loss';
                case 'ZLP interface score', img = app.optionalMap(cr, 'zlp_interface_score_map', cr.R2_map); ttl = 'Sharp-interface ZLP score';
                case 'ZLP sharp interface', img = app.optionalMap(cr, 'zlp_interface_mask', cr.R2_map); ttl = 'ZLP sharp-interface mask';
                otherwise, img = cr.R2_map; ttl = 'R2 map';
            end
        end

        function onMapClick(app, evt, stage2)
            if isstruct(evt) && isfield(evt, 'IntersectionPoint')
                cp = evt.IntersectionPoint;
            elseif isprop(evt, 'IntersectionPoint')
                cp = evt.IntersectionPoint;
            else
                if stage2, ax = app.Stage2MapAxes; else, ax = app.Stage1MapAxes; end
                cp = ax.CurrentPoint(1, :);
            end
            xx = max(1, round(cp(1)));
            yy = max(1, round(cp(2)));
            if stage2
                app.Stage2PixelYSpinner.Value = yy; app.Stage2PixelXSpinner.Value = xx;
            else
                app.Stage1PixelYSpinner.Value = yy; app.Stage1PixelXSpinner.Value = xx;
            end
            app.selectResultByPixel(yy, xx);
            app.updateSingleFitAxes(stage2);
            app.markSelectedMapPixel(stage2, yy, xx);
            app.plotDiagnosticSelection();
        end

        function markSelectedMapPixel(app, stage2, yy, xx)
            if stage2
                ax = app.Stage2MapAxes;
                app.deleteGraphicsHandles('Stage2CursorHandles');
            else
                ax = app.Stage1MapAxes;
                app.deleteGraphicsHandles('Stage1CursorHandles');
            end
            hold(ax, 'on');
            handles = [ ...
                plot(ax, xx, yy, 'wo', 'MarkerSize', 9, 'LineWidth', 1.8), ...
                plot(ax, xx, yy, '+', 'Color', [0.85 0.08 0.10], 'MarkerSize', 9, 'LineWidth', 1.4)];
            for ii = 1:numel(handles)
                try
                    handles(ii).HitTest = 'off';
                    handles(ii).PickableParts = 'none';
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

        function pickMapPixelWithCursor(app, stage2)
            if isempty(app.CubeResult)
                app.logStage(stage2, 'No datacube map available for pixel selection.');
                return;
            end
            if stage2, ax = app.Stage2MapAxes; else, ax = app.Stage1MapAxes; end
            try
                app.logStage(stage2, 'Pixel: double-click confirms; Ctrl+Q cancels; Ctrl+H confirms and plots.');
                roi = drawpoint(ax, 'Color', [0.95 0.20 0.05], 'LineWidth', 1.4);
                [pos, cancelled] = app.waitForROIWithHotkeys(roi, ...
                    'Pixel: double-click confirms; Ctrl+Q cancels; Ctrl+H confirms and plots.');
                if cancelled || isempty(pos)
                    app.logStage(stage2, 'Pixel selection cancelled.');
                    return;
                end
            catch ME
                app.logStage(stage2, ['Pixel selection requires MATLAB ROI tools: ', ME.message]);
                return;
            end
            evt = struct('IntersectionPoint', [pos(1), pos(2), 0]);
            app.onMapClick(evt, stage2);
        end

        function selectResultByPixel(app, yy, xx)
            if isempty(app.CubeResult) || ~isfield(app.CubeResult, 'sample_pixels'), return; end
            sz = size(app.CubeResult.R2_map);
            if yy > sz(1) || xx > sz(2), return; end
            pix = sub2ind(sz, yy, xx);
            idx = find(app.CubeResult.sample_pixels == pix, 1);
            if ~isempty(idx)
                res = app.CubeResult.results{idx};
                if isstruct(res) && ~isfield(res, 'error')
                    app.CurrentResult = res;
                end
            end
        end

        function updateSingleFitAxes(app, stage2)
            if isempty(app.CurrentResult), return; end
            if stage2, ax = app.Stage2SpectrumAxes; else, ax = app.Stage1SpectrumAxes; end
            cla(ax);
            res = app.CurrentResult;
            plot(ax, res.ene_corrected, res.raw_data, 'k.', 'MarkerSize', 4); hold(ax, 'on');
            plot(ax, res.ene_corrected, res.fitted_spectrum, 'k-', 'LineWidth', 1.2);
            plot(ax, res.ene_corrected, res.zlp_only, 'Color', [0.05 0.45 0.18], 'LineWidth', 1.0);
            plot(ax, res.ene_corrected, res.vib_only, 'Color', [0.9 0.25 0.05], 'LineWidth', 1.0);
            hold(ax, 'off'); grid(ax, 'on');
            legend(ax, {'raw', 'fit', 'ZLP', 'vib'}, 'Location', 'best');
            title(ax, sprintf('R2 %.6g | NRMSE %.4g | T %.1f K | rel %.2f', ...
                res.metrics.R2, res.metrics.NRMSE, res.T_fit, res.reliability_score));
            xlabel(ax, 'Energy loss (meV)'); ylabel(ax, 'Intensity');
        end

        function plotDiagnosticSelection(app)
            if isempty(app.CurrentResult), return; end
            res = app.CurrentResult;
            item = app.DiagnosticDropDown.Value;
            cla(app.ViewAxes1); cla(app.ViewAxes2);
            ene = res.ene_corrected(:);
            switch item
                case 'Full diagnostic figure'
                    app.openFullDiagnostic();
                case 'Raw / total fit'
                    plot(app.ViewAxes1, ene, res.raw_data, 'k.', ene, res.fitted_spectrum, 'r-', 'LineWidth', 1.2);
                    legend(app.ViewAxes1, {'raw', 'fit'}); grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Raw vs total fit');
                    plot(app.ViewAxes2, ene, res.raw_data(:) - res.fitted_spectrum(:), 'b-');
                    yline(app.ViewAxes2, 0, 'k-'); grid(app.ViewAxes2, 'on'); title(app.ViewAxes2, 'Residual');
                case 'ZLP components'
                    plot(app.ViewAxes1, ene, res.zlp_only, 'Color', [0.05 0.45 0.18], 'LineWidth', 1.4); hold(app.ViewAxes1, 'on');
                    if isfield(res, 'zlp_components') && ~isempty(res.zlp_components)
                        plot(app.ViewAxes1, ene, res.zlp_components, 'LineWidth', 0.8);
                    end
                    hold(app.ViewAxes1, 'off'); grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'ZLP envelope and components');
                    semilogy(app.ViewAxes2, ene, max(res.raw_data, eps), 'k.', ene, max(res.zlp_only, eps), 'g-');
                    grid(app.ViewAxes2, 'on'); title(app.ViewAxes2, 'Log ZLP check');
                case 'Sparse convolution / vibration'
                    plot(app.ViewAxes1, ene, res.vib_only, 'Color', [0.9 0.25 0.05], 'LineWidth', 1.2);
                    grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Sparse convolution / vib');
                    app.plotPM(app.ViewAxes2, res);
                case 'Residual / confidence interval'
                    residual = res.raw_data(:) - res.fitted_spectrum(:);
                    area(app.ViewAxes1, ene, residual, 'FaceColor', [0.45 0.62 0.95], 'EdgeColor', [0.1 0.35 0.95]);
                    yline(app.ViewAxes1, 0, 'k-'); grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Residual');
                    if isfield(res, 'confidence_interval') && isfield(res.confidence_interval, 'half_width')
                        plot(app.ViewAxes2, ene, res.confidence_interval.half_width, 'b-');
                        grid(app.ViewAxes2, 'on'); title(app.ViewAxes2, 'Energy-dependent CI half-width');
                    end
                case 'Sparse intensity sticks'
                    app.plotPM(app.ViewAxes1, res);
                case 'log(Iloss/Igain)'
                    app.plotSpectrumLossGainRatio(app.ViewAxes1, res);
                    app.plotPM(app.ViewAxes2, res);
                case 'Attention curve'
                    if isfield(res, 'stage2') && isfield(res.stage2, 'attention')
                        plot(app.ViewAxes1, ene, res.stage2.attention.feature_score, 'b-');
                        ylim(app.ViewAxes1, [0 1]); grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Attention feature score');
                    end
                case 'Loss components'
                    if isfield(res, 'loss_breakdown')
                        vals = 100 .* res.loss_breakdown.fraction;
                        labels = strrep(res.loss_breakdown.names, '_', ' ');
                        barh(app.ViewAxes1, vals); set(app.ViewAxes1, 'YTick', 1:numel(labels), 'YTickLabel', labels);
                        grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Loss budget (%)');
                    end
                case 'Log signal / log ZLP'
                    plot(app.ViewAxes1, ene, log10(max(res.raw_data, eps)), 'k.', ...
                        ene, log10(max(res.zlp_only, eps)), 'g-', ...
                        ene, log10(max(res.vib_only, eps)), 'r-');
                    legend(app.ViewAxes1, {'log raw', 'log ZLP', 'log vib'}, 'Location', 'best');
                    grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Log signal vs log ZLP');
            end
        end

        function plotPM(~, ax, res)
            cla(ax); hold(ax, 'on');
            if isfield(res, 'EM_final') && isfield(res, 'PM_stokes')
                stem(ax, res.EM_final, res.PM_stokes, 'Color', [0.95 0.45 0.12], 'Marker', 'none');
                stem(ax, -res.EM_final, res.PM_antistokes, 'Color', [0.1 0.50 0.95], 'Marker', 'none');
            end
            if isfield(res, 'stage2') && isfield(res.stage2, 'EM_hard') && ~isempty(res.stage2.EM_hard)
                yl = ylim(ax);
                plot(ax, res.stage2.EM_hard, 0 .* res.stage2.EM_hard + yl(2) .* 0.94, 'm.', 'MarkerSize', 8);
                plot(ax, -res.stage2.EM_hard, 0 .* res.stage2.EM_hard + yl(2) .* 0.88, 'c.', 'MarkerSize', 8);
            end
            hold(ax, 'off'); grid(ax, 'on');
            xlabel(ax, 'Mode energy (meV)'); ylabel(ax, 'Source intensity');
            title(ax, 'Sparse source intensity and hard-BE markers');
            legend(ax, {'Stokes intensity', 'Anti-Stokes intensity'}, 'Location', 'best');
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
            if isempty(app.CurrentResult), return; end
            fig = EELS_PearsonSparse_FullFit.plot_result(app.CurrentResult, 'GUI selected spectrum');
            fig.Visible = 'on';
        end

        function populateAdvancedTable(app)
            s = EELS_PearsonSparse_FullFit.defaults(struct());
            names = fieldnames(s);
            rows = {};
            for ii = 1:numel(names)
                val = s.(names{ii});
                if isstruct(val) || isa(val, 'function_handle') || iscell(val)
                    continue;
                end
                if ismember(names{ii}, {'datacube_sample_pixels', 'datacube_sample_strategy'})
                    continue;
                end
                rows(end + 1, :) = {names{ii}, app.zhName(names{ii}), app.valueToText(val), app.settingAdvice(names{ii})}; %#ok<AGROW>
            end
            app.AdvancedTable.Data = rows;
        end

        function applyAdvancedTable(app)
            if isempty(app.AdvancedTable) || isempty(app.AdvancedTable.Data), return; end
            s = app.Settings;
            data = app.AdvancedTable.Data;
            for ii = 1:size(data, 1)
                key = data{ii, 1};
                txt = data{ii, 3};
                if isfield(s, key)
                    s.(key) = app.textToValue(txt, s.(key));
                end
            end
            app.Settings = EELS_PearsonSparse_FullFit.defaults(s);
        end

        function resetAdvanced(app)
            app.Settings = EELS_PearsonSparse_FullFit.defaults(struct());
            app.Settings.datacube_sample_strategy = 'all';
            app.Settings.datacube_sample_pixels = Inf;
            app.applyTunedSettings(app.Settings);
            if ~isempty(app.ReferenceTuner)
                app.ReferenceTuner.Settings=app.Settings;app.ReferenceTuner.refreshParameters();
                app.ReferenceTuner.ActiveLabel.Text='Default settings restored';
            end
            app.logStatus('Advanced settings reset to defaults.');
        end

        function exportData(app)
            if isempty(app.Data), return; end
            assignin('base', 'eels_gui_data', app.Data);
            assignin('base', 'eels_gui_energy', app.Ene);
            app.logStatus('Exported eels_gui_data / eels_gui_energy to the base workspace.');
        end

        function exportResult(app, tag)
            selected = [];
            cuberes = [];
            switch lower(string(tag))
                case "stage1"
                    selected = app.Stage1SingleResult;
                    cuberes = app.Stage1CubeResult;
                    if isempty(selected) && ~isempty(cuberes) && isfield(cuberes, 'results')
                        selected = app.firstValidResult(cuberes.results);
                    end
                    if isempty(selected) && ~isempty(app.CurrentResult) && ~app.resultHasStage2(app.CurrentResult)
                        selected = app.CurrentResult;
                    end
                case "stage2"
                    selected = app.Stage2SingleResult;
                    cuberes = app.Stage2CubeResult;
                    if isempty(selected) && ~isempty(cuberes) && isfield(cuberes, 'results')
                        selected = app.firstValidResult(cuberes.results);
                    end
                    if isempty(selected) && ~isempty(app.CurrentResult) && app.resultHasStage2(app.CurrentResult)
                        selected = app.CurrentResult;
                    end
                otherwise
                    selected = app.CurrentResult;
                    cuberes = app.CubeResult;
            end
            exported = {};
            if ~isempty(selected)
                assignin('base', ['eels_', char(tag), '_selected_result'], selected);
                exported{end + 1} = ['eels_', char(tag), '_selected_result']; %#ok<AGROW>
                empm = app.makeEMPMExport(selected);
                assignin('base', ['eels_', char(tag), '_em_pm'], empm);
                exported{end + 1} = ['eels_', char(tag), '_em_pm']; %#ok<AGROW>
                try
                    empm_table = app.makeEMPMTable(empm);
                    assignin('base', ['eels_', char(tag), '_em_pm_table'], empm_table);
                    exported{end + 1} = ['eels_', char(tag), '_em_pm_table']; %#ok<AGROW>
                catch
                end
            end
            if ~isempty(cuberes)
                assignin('base', ['eels_', char(tag), '_cuberes'], cuberes);
                exported{end + 1} = ['eels_', char(tag), '_cuberes']; %#ok<AGROW>
            end
            if isempty(exported)
                app.logStatus(['No exportable ', char(tag), ' result. Run or load the corresponding stage first.']);
            else
                app.logStatus(['Exported to the base workspace: ', strjoin(exported, ' / ')]);
            end
            app.logStatus(['Result exported to the base workspace; tag: ', tag]);
        end

        function empm = makeEMPMExport(~, res)
            empm = struct();
            if isempty(res) || ~isstruct(res)
                return;
            end
            if isfield(res, 'EM_final') && isfield(res, 'PM_stokes') && isfield(res, 'PM_antistokes')
                em = res.EM_final(:);
                ps = res.PM_stokes(:);
                pa = res.PM_antistokes(:);
                n = min([numel(em), numel(ps), numel(pa)]);
                empm.stage1 = struct('EM_meV', em(1:n), 'PM_stokes', ps(1:n), ...
                    'PM_antistokes', pa(1:n), 'signed_EM_meV', [-flipud(em(1:n)); em(1:n)], ...
                    'signed_PM', [flipud(pa(1:n)); ps(1:n)]);
            end
            if isfield(res, 'stage2') && isstruct(res.stage2) && ...
                    isfield(res.stage2, 'EM_final') && isfield(res.stage2, 'PM_stokes') && ...
                    isfield(res.stage2, 'PM_antistokes')
                em = res.stage2.EM_final(:);
                ps = res.stage2.PM_stokes(:);
                pa = res.stage2.PM_antistokes(:);
                n = min([numel(em), numel(ps), numel(pa)]);
                hard = false(n, 1);
                if isfield(res.stage2, 'EM_hard') && ~isempty(res.stage2.EM_hard)
                    hard = ismembertol(em(1:n), res.stage2.EM_hard(:), 1e-6);
                end
                empm.stage2 = struct('EM_meV', em(1:n), 'PM_stokes', ps(1:n), ...
                    'PM_antistokes', pa(1:n), 'is_hard_BE', hard, ...
                    'signed_EM_meV', [-flipud(em(1:n)); em(1:n)], ...
                    'signed_PM', [flipud(pa(1:n)); ps(1:n)]);
            end
        end

        function T = makeEMPMTable(~, empm)
            rows_stage = {};
            em_all = [];
            pm_s_all = [];
            pm_a_all = [];
            hard_all = [];
            if isfield(empm, 'stage1')
                n = numel(empm.stage1.EM_meV);
                rows_stage = [rows_stage; repmat({'stage1'}, n, 1)]; %#ok<AGROW>
                em_all = [em_all; empm.stage1.EM_meV(:)]; %#ok<AGROW>
                pm_s_all = [pm_s_all; empm.stage1.PM_stokes(:)]; %#ok<AGROW>
                pm_a_all = [pm_a_all; empm.stage1.PM_antistokes(:)]; %#ok<AGROW>
                hard_all = [hard_all; false(n, 1)]; %#ok<AGROW>
            end
            if isfield(empm, 'stage2')
                n = numel(empm.stage2.EM_meV);
                rows_stage = [rows_stage; repmat({'stage2'}, n, 1)]; %#ok<AGROW>
                em_all = [em_all; empm.stage2.EM_meV(:)]; %#ok<AGROW>
                pm_s_all = [pm_s_all; empm.stage2.PM_stokes(:)]; %#ok<AGROW>
                pm_a_all = [pm_a_all; empm.stage2.PM_antistokes(:)]; %#ok<AGROW>
                hard_all = [hard_all; empm.stage2.is_hard_BE(:)]; %#ok<AGROW>
            end
            T = table(string(rows_stage), em_all, pm_s_all, pm_a_all, hard_all, ...
                'VariableNames', {'stage', 'EM_meV', 'PM_stokes', 'PM_antistokes', 'is_hard_BE'});
        end

        function obj = makeEELSObject(~, cube, ene)
            obj = EELS;
            obj.data = cube;
            obj.ene = ene(:).';
            escale = median(diff(ene), 'omitnan');
            obj.info = struct('xscale', 1, 'yscale', 1, 'escale', escale, ...
                'isaligned', 0, 'isnormalized', 0, 'issaturated', 0);
        end

        function sig = resultSignalSpectrum(~, res)
            if isfield(res, 'stage2') && isfield(res.stage2, 'vib_norm') && ...
                    isfield(res, 'norm_val') && ~isempty(res.stage2.vib_norm)
                sig = res.stage2.vib_norm(:) .* res.norm_val;
            elseif isfield(res, 'vib_only')
                sig = res.vib_only(:);
            elseif isfield(res, 'raw_data') && isfield(res, 'zlp_only')
                sig = max(res.raw_data(:) - res.zlp_only(:), 0);
            else
                sig = [];
            end
        end

        function filename = writeGuiCube(app, outdir)
            if ~isfolder(outdir), mkdir(outdir); end
            filename = fullfile(outdir, 'gui_current_datacube.mat');
            gui_cube = app.Data;
            ene = app.Ene;
            save(filename, 'gui_cube', 'ene', '-v7.3');
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

        function b = techButton(~, parent, txt, callback)
            b = uibutton(parent, 'Text', txt, 'ButtonPushedFcn', callback);
            b.FontWeight = 'bold';
            b.FontColor = [1 1 1];
            b.BackgroundColor = [0.05 0.33 0.62];
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

        function lab = techLabel(~, parent, txt)
            lab = uilabel(parent, 'Text', txt);
            lab.FontWeight = 'bold';
            lab.FontColor = [0.05 0.18 0.32];
        end

        function styleAxes(~, ax)
            ax.Color = [0.985 0.992 1.000];
            ax.XColor = [0.08 0.16 0.28];
            ax.YColor = [0.08 0.16 0.28];
            ax.GridColor = [0.55 0.68 0.82];
            ax.FontName = 'Microsoft YaHei UI';
            grid(ax, 'on');
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

        function deleteGraphicsHandles(app, prop_name)
            handles = app.(prop_name);
            if isempty(handles), return; end
            for ii = 1:numel(handles)
                try
                    if isgraphics(handles(ii))
                        delete(handles(ii));
                    end
                catch
                end
            end
            app.(prop_name) = [];
        end

        function appendGraphicsHandle(app, prop_name, handle)
            if isempty(app.(prop_name))
                app.(prop_name) = handle;
            else
                app.(prop_name)(end + 1) = handle;
            end
        end

        function [pos, cancelled] = waitForROIWithHotkeys(app, roi, prompt)
            pos = [];
            cancelled = false;
            fig = app.UIFigure;
            token = ['EELS_ROI_', char(java.util.UUID.randomUUID())];
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
                        delete(roi);
                    end
                catch
                end
            end
        end

        function setPixelLimits(app)
            if isempty(app.Data), return; end
            [ny, nx, ~] = size(app.Data);
            spinners = {app.Stage1PixelYSpinner, app.Stage1PixelXSpinner, ...
                app.Stage2PixelYSpinner, app.Stage2PixelXSpinner};
            for ii = 1:numel(spinners)
                if isempty(spinners{ii}) || ~isvalid(spinners{ii}), continue; end
                if mod(ii, 2) == 1
                    spinners{ii}.Limits = [1 ny];
                    spinners{ii}.Value = min(max(round(spinners{ii}.Value), 1), ny);
                else
                    spinners{ii}.Limits = [1 nx];
                    spinners{ii}.Value = min(max(round(spinners{ii}.Value), 1), nx);
                end
            end
        end

        function [yy, xx] = currentPixel(app, stage2)
            ctrl = app.controlsForStage(stage2);
            yy = round(ctrl.pixelY.Value);
            xx = round(ctrl.pixelX.Value);
            if ~isempty(app.Data)
                [ny, nx, ~] = size(app.Data);
                yy = min(max(yy, 1), ny);
                xx = min(max(xx, 1), nx);
            end
        end

        function win = energyWindow(app)
            win = sort([app.EnergyMinEdit.Value, app.EnergyMaxEdit.Value]);
            if any(~isfinite(win)) || win(1) == win(2)
                win = [-Inf, Inf];
            end
        end

        function win = stageEnergyWindow(app, stage2)
            ctrl = app.controlsForStage(stage2);
            if isfield(ctrl, 'cropMin') && isfield(ctrl, 'cropMax') && ...
                    isvalid(ctrl.cropMin) && isvalid(ctrl.cropMax)
                win = sort([ctrl.cropMin.Value, ctrl.cropMax.Value]);
                if all(isfinite(win)) && win(1) ~= win(2)
                    return;
                end
            end
            win = app.energyWindow();
        end

        function inheritStageEnergyWindow(app, stage2)
            ctrl = app.controlsForStage(stage2);
            win = app.energyWindow();
            if isfield(ctrl, 'cropMin') && isvalid(ctrl.cropMin)
                ctrl.cropMin.Value = win(1);
            end
            if isfield(ctrl, 'cropMax') && isvalid(ctrl.cropMax)
                ctrl.cropMax.Value = win(2);
            end
            app.logStage(stage2, sprintf('Using Data energy window: [%.1f, %.1f] meV.', win(1), win(2)));
        end

        function tf = shouldPauseDatacube(app, dlg)
            tf = logical(app.DatacubePauseRequested);
            if ~tf && ~isempty(dlg) && isvalid(dlg)
                try
                    tf = logical(dlg.CancelRequested);
                catch
                    tf = false;
                end
            end
            if tf
                app.DatacubePauseRequested = true;
            end
        end

        function map = optionalMap(~, s, field_name, template)
            if isstruct(s) && isfield(s, field_name) && isequal(size(s.(field_name)), size(template))
                map = s.(field_name);
            else
                map = NaN(size(template));
            end
        end

        function ctrl = controlsForStage(app, stage2)
            if stage2
                ctrl = app.Stage2Controls;
            else
                ctrl = app.Stage1Controls;
            end
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

        function logStage(app, stage2, msg)
            if stage2
                app.Stage2StatusArea.Value = [app.Stage2StatusArea.Value; {msg}];
            else
                app.Stage1StatusArea.Value = [app.Stage1StatusArea.Value; {msg}];
            end
            app.logStatus(msg);
        end

        function logStatus(app, msg)
            if isempty(app.StatusArea), return; end
            stamp = char(datetime('now', 'Format', 'HH:mm:ss'));
            app.StatusArea.Value = [app.StatusArea.Value; {[stamp, '  ', msg]}];
            drawnow limitrate;
        end

        function showError(app, ME)
            msg = sprintf('%s\n%s', ME.message, app.shortStack(ME));
            app.logStatus(['Error: ', ME.message]);
            uialert(app.UIFigure, msg, 'EELS Toolbox Error');
        end

        function txt = shortStack(~, ME)
            parts = cell(1, min(4, numel(ME.stack)));
            for ii = 1:numel(parts)
                parts{ii} = sprintf('%s:%d', ME.stack(ii).name, ME.stack(ii).line);
            end
            txt = strjoin(parts, newline);
        end

        function txt = valueToText(~, val)
            if islogical(val)
                txt = mat2str(val);
            elseif isnumeric(val)
                txt = mat2str(val);
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

        function name = zhName(~, key)
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

    end
end

