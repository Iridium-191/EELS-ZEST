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
            app.logStatus('GUI 已启动。先加载 datacube 或已有 result；Stage1 用于稀疏反卷扣背地，Stage2 用于温度/可靠性/界面诊断。');
        end

        function delete(app)
            if ~isempty(app.UIFigure) && isvalid(app.UIFigure)
                delete(app.UIFigure);
            end
        end
    end

    methods (Access = private)
        function buildUI(app, visible)
            [fig_pos, resolved_mode] = app.responsiveFigurePosition(app.LayoutMode);
            app.LayoutMode = resolved_mode;
            app.IsCompact = strcmp(resolved_mode, 'compact');
            app.UIFigure = uifigure('Name', 'EELS-ZEST 2.1 Two-Stage Toolbox', ...
                'Position', fig_pos, 'Color', [0.94 0.97 1.00], 'Visible', visible);
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
            app.DataTab = uitab(app.TabGroup, 'Title', '数据与预处理', 'BackgroundColor', [0.94 0.97 1.00]);
            app.Stage1Tab = uitab(app.TabGroup, 'Title', 'Stage1 稀疏反卷', 'BackgroundColor', [0.94 0.97 1.00]);
            app.Stage2Tab = uitab(app.TabGroup, 'Title', 'Stage2 温度拟合', 'BackgroundColor', [0.94 0.97 1.00]);
            app.ViewTab = uitab(app.TabGroup, 'Title', '诊断图查看', 'BackgroundColor', [0.94 0.97 1.00]);
            app.AdvancedTab = uitab(app.TabGroup, 'Title', '高级参数', 'BackgroundColor', [0.94 0.97 1.00]);
            app.buildDataTab();
            app.buildStage1Tab();
            app.buildStage2Tab();
            app.buildViewTab();
            app.buildAdvancedTab();
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

            p = uipanel(gl, 'Title', '数据导入 / 预处理 / ROI / 传统扣背地', ...
                'FontWeight', 'bold', 'ForegroundColor', [0.05 0.18 0.32], ...
                'BackgroundColor', [0.88 0.94 1.00], 'BorderType', 'none');
            p.Layout.Row = 1; p.Layout.Column = [1 4];
            app.enableScrollable(p);
            g = uigridlayout(p, [4 12]);
            g.RowHeight = {30, 30, 30, 30};
            g.ColumnWidth = {92, 96, 96, 96, 78, 78, 78, 78, 96, 110, 110, 150};
            g.Padding = [10 8 10 8];
            g.BackgroundColor = [0.88 0.94 1.00];

            b = app.techButton(g, '加载MAT', @(~, ~)app.onLoadData()); b.Layout.Row = 1; b.Layout.Column = 1;
            b = app.techButton(g, '加载结果', @(~, ~)app.onLoadResult()); b.Layout.Row = 1; b.Layout.Column = 2;
            b = app.techButton(g, 'Align', @(~, ~)app.onAlign()); b.Layout.Row = 1; b.Layout.Column = 3;
            b = app.techButton(g, 'BM3D降噪', @(~, ~)app.onBM3D()); b.Layout.Row = 1; b.Layout.Column = 4;
            lab = app.techLabel(g, 'Align方法'); lab.Layout.Row = 1; lab.Layout.Column = 5;
            app.AlignMethodDropDown = uidropdown(g, 'Items', {'xcorr', 'gaussian', 'COM'}, 'Value', 'xcorr');
            app.AlignMethodDropDown.Layout.Row = 1; app.AlignMethodDropDown.Layout.Column = 6;
            lab = app.techLabel(g, 'BM因子'); lab.Layout.Row = 1; lab.Layout.Column = 7;
            app.BM3DFactorEdit = uieditfield(g, 'numeric', 'Value', 1, 'Limits', [0.1 Inf]);
            app.BM3DFactorEdit.Layout.Row = 1; app.BM3DFactorEdit.Layout.Column = 8;
            b = app.techButton(g, '能量图', @(~, ~)app.plotEnergyMap()); b.Layout.Row = 1; b.Layout.Column = 9;
            b = app.techButton(g, 'LineProfile', @(~, ~)app.plotLineProfile()); b.Layout.Row = 1; b.Layout.Column = 10;
            b = app.techButton(g, '导出数据', @(~, ~)app.exportData()); b.Layout.Row = 1; b.Layout.Column = 11;
            app.FileEdit = uieditfield(g, 'text', 'Editable', 'off');
            app.FileEdit.Layout.Row = 2; app.FileEdit.Layout.Column = [1 8];
            lab = app.techLabel(g, '能量窗'); lab.Layout.Row = 2; lab.Layout.Column = 9;
            app.EnergyMinEdit = uieditfield(g, 'numeric', 'Value', -250);
            app.EnergyMinEdit.Layout.Row = 2; app.EnergyMinEdit.Layout.Column = 10;
            app.EnergyMaxEdit = uieditfield(g, 'numeric', 'Value', 250);
            app.EnergyMaxEdit.Layout.Row = 2; app.EnergyMaxEdit.Layout.Column = 11;
            b = app.techButton(g, '刷新预览', @(~, ~)app.refreshDataPreview()); b.Layout.Row = 2; b.Layout.Column = 12;
            b = app.techButton(g, '框选ROI积分', @(~, ~)app.selectROIAndIntegrate()); b.Layout.Row = 3; b.Layout.Column = [1 2];
            b = app.techButton(g, '导出ROI', @(~, ~)app.exportROI()); b.Layout.Row = 3; b.Layout.Column = 3;
            b = app.techButton(g, '导入拟合信号', @(~, ~)app.importFitSignalToData()); b.Layout.Row = 3; b.Layout.Column = [4 5];

            lab = app.techLabel(g, '传统扣背地'); lab.Layout.Row = 4; lab.Layout.Column = 1;
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
            app.DataIOActionDropDown = uidropdown(g, 'Items', {'导入MAT', '导入结果', '导出数据', '导出ROI', '导出当前结果'}, 'Value', '导入MAT');
            app.DataIOActionDropDown.Layout.Row = 3; app.DataIOActionDropDown.Layout.Column = [6 7];
            b = app.techButton(g, '执行导入/导出', @(~, ~)app.runDataIOAction()); b.Layout.Row = 3; b.Layout.Column = [8 9];
            app.DataToolDropDown = uidropdown(g, 'Items', {'NMF分解', '扣STA', '传统扣背地(交互)', 'Separate IF mode', 'MLLS残差', 'Peak fitting', '交互漂移矫正', '叠谱工具'}, 'Value', 'NMF分解');
            app.DataToolDropDown.Layout.Row = 3; app.DataToolDropDown.Layout.Column = [10 11];
            app.DataToolDropDown.Items = [app.DataToolDropDown.Items, {'声子源强度分析'}];
            b = app.techButton(g, '打开工具', @(~, ~)app.runDataToolAction()); b.Layout.Row = 3; b.Layout.Column = 12;
            b = app.techButton(g, '调用原版扣背地', @(~, ~)app.applyTraditionalBackground()); b.Layout.Row = 4; b.Layout.Column = [8 10];
            b = app.techButton(g, '叠谱工具', @(~, ~)app.openStackSpectraTool()); b.Layout.Row = 4; b.Layout.Column = [11 12];

            app.DataAxes = uiaxes(gl);
            app.DataAxes.Layout.Row = 2; app.DataAxes.Layout.Column = [1 2];
            app.styleAxes(app.DataAxes);
            title(app.DataAxes, '能量过滤像 / ROI map'); axis(app.DataAxes, 'image');
            app.SpectrumAxes = uiaxes(gl);
            app.SpectrumAxes.Layout.Row = 2; app.SpectrumAxes.Layout.Column = [3 4];
            app.styleAxes(app.SpectrumAxes);
            title(app.SpectrumAxes, '选中像素谱线 / ROI积分谱线');

            app.StatusArea = uitextarea(gl, 'Editable', 'off');
            app.StatusArea.Layout.Row = 3; app.StatusArea.Layout.Column = [1 4];
        end

        function buildDataTabCompact(app)
            gl = uigridlayout(app.DataTab, [4 4]);
            gl.RowHeight = {245, '1x', '1x', 105};
            gl.ColumnWidth = {'1x', '1x', '1x', '1x'};
            gl.Padding = [10 8 10 8];
            gl.BackgroundColor = [0.94 0.97 1.00];

            p = uipanel(gl, 'Title', '数据导入 / 预处理 / ROI / 工具 (compact)', ...
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

            b = app.techButton(g, '加载MAT', @(~, ~)app.onLoadData()); b.Layout.Row = 1; b.Layout.Column = [1 2];
            b = app.techButton(g, '加载结果', @(~, ~)app.onLoadResult()); b.Layout.Row = 1; b.Layout.Column = [3 4];
            b = app.techButton(g, '导出数据', @(~, ~)app.exportData()); b.Layout.Row = 1; b.Layout.Column = [5 6];
            b = app.techButton(g, '刷新预览', @(~, ~)app.refreshDataPreview()); b.Layout.Row = 1; b.Layout.Column = [7 8];

            app.FileEdit = uieditfield(g, 'text', 'Editable', 'off');
            app.FileEdit.Layout.Row = 2; app.FileEdit.Layout.Column = [1 8];

            lab = app.techLabel(g, '能量窗'); lab.Layout.Row = 3; lab.Layout.Column = 1;
            app.EnergyMinEdit = uieditfield(g, 'numeric', 'Value', -250);
            app.EnergyMinEdit.Layout.Row = 3; app.EnergyMinEdit.Layout.Column = 2;
            app.EnergyMaxEdit = uieditfield(g, 'numeric', 'Value', 250);
            app.EnergyMaxEdit.Layout.Row = 3; app.EnergyMaxEdit.Layout.Column = 3;
            b = app.techButton(g, '能量图', @(~, ~)app.plotEnergyMap()); b.Layout.Row = 3; b.Layout.Column = 4;
            b = app.techButton(g, 'LineProfile', @(~, ~)app.plotLineProfile()); b.Layout.Row = 3; b.Layout.Column = [5 6];
            b = app.techButton(g, '叠谱工具', @(~, ~)app.openStackSpectraTool()); b.Layout.Row = 3; b.Layout.Column = [7 8];

            b = app.techButton(g, '框选ROI积分', @(~, ~)app.selectROIAndIntegrate()); b.Layout.Row = 4; b.Layout.Column = [1 2];
            b = app.techButton(g, '导出ROI', @(~, ~)app.exportROI()); b.Layout.Row = 4; b.Layout.Column = [3 4];
            b = app.techButton(g, '导入拟合信号', @(~, ~)app.importFitSignalToData()); b.Layout.Row = 4; b.Layout.Column = [5 6];
            b = app.techButton(g, 'Align', @(~, ~)app.onAlign()); b.Layout.Row = 4; b.Layout.Column = 7;
            b = app.techButton(g, 'BM3D', @(~, ~)app.onBM3D()); b.Layout.Row = 4; b.Layout.Column = 8;

            lab = app.techLabel(g, 'Align方法'); lab.Layout.Row = 5; lab.Layout.Column = 1;
            app.AlignMethodDropDown = uidropdown(g, 'Items', {'xcorr', 'gaussian', 'COM'}, 'Value', 'xcorr');
            app.AlignMethodDropDown.Layout.Row = 5; app.AlignMethodDropDown.Layout.Column = 2;
            lab = app.techLabel(g, 'BM因子'); lab.Layout.Row = 5; lab.Layout.Column = 3;
            app.BM3DFactorEdit = uieditfield(g, 'numeric', 'Value', 1, 'Limits', [0.1 Inf]);
            app.BM3DFactorEdit.Layout.Row = 5; app.BM3DFactorEdit.Layout.Column = 4;
            app.DataIOActionDropDown = uidropdown(g, 'Items', {'导入MAT', '导入结果', '导出数据', '导出ROI', '导出当前结果'}, 'Value', '导入MAT');
            app.DataIOActionDropDown.Layout.Row = 5; app.DataIOActionDropDown.Layout.Column = [5 6];
            b = app.techButton(g, '执行导入/导出', @(~, ~)app.runDataIOAction()); b.Layout.Row = 5; b.Layout.Column = [7 8];

            app.DataToolDropDown = uidropdown(g, 'Items', {'NMF分解', '扣STA', '传统扣背地(交互)', 'Separate IF mode', 'MLLS残差', 'Peak fitting', '交互漂移矫正', '叠谱工具', '声子源强度分析'}, 'Value', 'NMF分解');
            app.DataToolDropDown.Layout.Row = 6; app.DataToolDropDown.Layout.Column = [1 3];
            b = app.techButton(g, '打开工具', @(~, ~)app.runDataToolAction()); b.Layout.Row = 6; b.Layout.Column = [4 5];
            app.BGMethodDropDown = uidropdown(g, 'Items', {'Pearson', 'ExpPoly3', 'ExpPoly2', 'ExpPoly4', ...
                'Power', 'Exponential_2', 'PseudoVoigt'}, 'Value', 'ExpPoly3');
            app.BGMethodDropDown.Layout.Row = 6; app.BGMethodDropDown.Layout.Column = [6 7];
            b = app.techButton(g, '原版扣背地', @(~, ~)app.applyTraditionalBackground()); b.Layout.Row = 6; b.Layout.Column = 8;
            app.BGWin1MinEdit = uieditfield(app.UIFigure, 'numeric', 'Value', -120, 'Visible', 'off');
            app.BGWin1MaxEdit = uieditfield(app.UIFigure, 'numeric', 'Value', -90, 'Visible', 'off');
            app.BGWin2MinEdit = uieditfield(app.UIFigure, 'numeric', 'Value', 120, 'Visible', 'off');
            app.BGWin2MaxEdit = uieditfield(app.UIFigure, 'numeric', 'Value', 180, 'Visible', 'off');

            app.DataAxes = uiaxes(gl);
            app.DataAxes.Layout.Row = 2; app.DataAxes.Layout.Column = [1 4];
            app.styleAxes(app.DataAxes);
            title(app.DataAxes, '能量过滤像 / ROI map'); axis(app.DataAxes, 'image');
            app.SpectrumAxes = uiaxes(gl);
            app.SpectrumAxes.Layout.Row = 3; app.SpectrumAxes.Layout.Column = [1 4];
            app.styleAxes(app.SpectrumAxes);
            title(app.SpectrumAxes, '选中像素谱线 / ROI积分谱线');

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

            p = uipanel(gl, 'Title', 'Stage1: ZLP + sparse convolution 扣背地', ...
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
            title(app.Stage1SpectrumAxes, '当前单点 Stage1');

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

            p = uipanel(gl, 'Title', 'Stage2: 注意力窗口 + BE 温度拟合', ...
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
            title(app.Stage2MapAxes, '温度 / 可靠性 / 界面 map'); axis(app.Stage2MapAxes, 'image');
            app.Stage2SpectrumAxes = uiaxes(gl);
            if app.IsCompact
                app.Stage2SpectrumAxes.Layout.Row = 3; app.Stage2SpectrumAxes.Layout.Column = [1 4];
            else
                app.Stage2SpectrumAxes.Layout.Row = 2; app.Stage2SpectrumAxes.Layout.Column = [3 4];
            end
            app.styleAxes(app.Stage2SpectrumAxes);
            title(app.Stage2SpectrumAxes, '当前单点 Stage2');

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

            lab = app.techLabel(g, '单点Y'); lab.Layout.Row = 1; lab.Layout.Column = 1;
            py = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            py.Layout.Row = 1; py.Layout.Column = 2;
            lab = app.techLabel(g, '单点X'); lab.Layout.Row = 1; lab.Layout.Column = 3;
            px = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            px.Layout.Row = 1; px.Layout.Column = 4;
            lab = app.techLabel(g, 'Y bin'); lab.Layout.Row = 1; lab.Layout.Column = 5;
            sy = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            sy.Layout.Row = 1; sy.Layout.Column = 6;
            lab = app.techLabel(g, 'X bin'); lab.Layout.Row = 1; lab.Layout.Column = 7;
            sx = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            sx.Layout.Row = 1; sx.Layout.Column = 8;
            lab = app.techLabel(g, '能量bin'); lab.Layout.Row = 1; lab.Layout.Column = 9;
            eb = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            eb.Layout.Row = 1; eb.Layout.Column = 10;
            fullLabel = uilabel(g, 'Text', 'datacube默认全空间');
            fullLabel.Layout.Row = 1; fullLabel.Layout.Column = [11 12];
            fullLabel.FontColor = [0.03 0.30 0.48];
            fullLabel.FontWeight = 'bold';

            lab = app.techLabel(g, '截取范围');
            lab.Layout.Row = 2; lab.Layout.Column = 1;
            cropMin = uieditfield(g, 'numeric', 'Value', app.EnergyMinEdit.Value);
            cropMin.Layout.Row = 2; cropMin.Layout.Column = 2;
            cropMax = uieditfield(g, 'numeric', 'Value', app.EnergyMaxEdit.Value);
            cropMax.Layout.Row = 2; cropMax.Layout.Column = 3;
            inheritBtn = app.techButton(g, '继承数据页', @(~, ~)app.inheritStageEnergyWindow(stage2));
            inheritBtn.Layout.Row = 2; inheritBtn.Layout.Column = [4 5];
            cropNote = uilabel(g, 'Text', '单位 meV；如 -190 到 190 可明显提速');
            cropNote.Layout.Row = 2; cropNote.Layout.Column = [6 12];
            cropNote.FontColor = [0.05 0.30 0.48];

            if stage2
                app.Stage2PixelYSpinner = py; app.Stage2PixelXSpinner = px;
            else
                app.Stage1PixelYSpinner = py; app.Stage1PixelXSpinner = px;
            end

            workers = uispinner(g, 'Limits', [1 64], 'Value', 8, 'Step', 1);
            checkpoint = uicheckbox(g, 'Text', '断点保存', 'Value', true);
            resume = uicheckbox(g, 'Text', '断点续跑', 'Value', false);
            saveplots = uicheckbox(g, 'Text', '保存诊断图', 'Value', false);
            lab = app.techLabel(g, '输出目录');
            outdir = uieditfield(g, 'text', 'Value', app.LastOutputDir);
            lab.Layout.Row = 3; lab.Layout.Column = 1;
            outdir.Layout.Row = 3; outdir.Layout.Column = [2 12];
            lab = app.techLabel(g, '并行worker');
            lab.Layout.Row = 4; lab.Layout.Column = 1;
            workers.Layout.Row = 4; workers.Layout.Column = 2;
            checkpoint.Layout.Row = 4; checkpoint.Layout.Column = [3 4];
            resume.Layout.Row = 4; resume.Layout.Column = [5 6];
            saveplots.Layout.Row = 4; saveplots.Layout.Column = [7 8];
            if stage2
                lab = app.techLabel(g, '算法');
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
                b = app.techButton(g, '拟合单点', @(~, ~)app.runSingle(true)); b.Layout.Row = 5; b.Layout.Column = 1;
                b = app.techButton(g, '光标选点', @(~, ~)app.pickMapPixelWithCursor(true)); b.Layout.Row = 5; b.Layout.Column = 2;
                b = app.techButton(g, '串行全量', @(~, ~)app.runDatacube(true, 'serial')); b.Layout.Row = 5; b.Layout.Column = 3;
                b = app.techButton(g, '并行全量', @(~, ~)app.runDatacube(true, 'parallel')); b.Layout.Row = 5; b.Layout.Column = 4;
                b = app.techButton(g, '导出结果', @(~, ~)app.exportResult('stage2')); b.Layout.Row = 5; b.Layout.Column = [5 6];
                b = app.techButton(g, '送信号回数据页', @(~, ~)app.importFitSignalToData()); b.Layout.Row = 5; b.Layout.Column = [7 8];
                lab = app.techLabel(g, '显示map');
                lab.Layout.Row = 5; lab.Layout.Column = 9;
                app.Stage2MapDropDown = uidropdown(g, 'Items', app.mapItems(), 'Value', 'T data valid', ...
                    'ValueChangedFcn', @(~, ~)app.refreshResultMap(true));
                app.Stage2MapDropDown.Layout.Row = 5; app.Stage2MapDropDown.Layout.Column = [10 12];
                b = app.techButton(g, '2σ联合修正离群', @(~, ~)app.refitReliabilityOutliers());
                b.Layout.Row = 6; b.Layout.Column = [1 3];
                outlierNote = uilabel(g, 'Text', '复用Stage1，联合检查Reliability/R2/NRMSE/T不确定度的2σ离群点');
                outlierNote.Layout.Row = 6; outlierNote.Layout.Column = [4 6];
                outlierNote.FontColor = [0.05 0.30 0.48];
                hierarchical = uicheckbox(g, 'Text', '分层ZLP + ZLP界面引导', 'Value', true);
                hierarchical.Tooltip = '启用ZLP谱形聚类、类模板约束、Stage1空间平滑和联合温度MAP。';
                hierarchical.Layout.Row = 6; hierarchical.Layout.Column = [7 9];
                forceRawStage2 = uicheckbox(g, 'Text', '强制显示Stage2', 'Value', false);
                forceRawStage2.Tooltip = '勾选后关闭Stage2谱形保护，不回退到Stage1，用于查看硬BE强行拟合的真实形状。';
                forceRawStage2.Layout.Row = 6; forceRawStage2.Layout.Column = [10 12];
                ctrl.hierarchical = hierarchical;
                ctrl.forceRawStage2 = forceRawStage2;
                ctrl.mapDropDown = app.Stage2MapDropDown;
                app.Stage2Controls = ctrl;
            else
                b = app.techButton(g, '拟合单点', @(~, ~)app.runSingle(false)); b.Layout.Row = 5; b.Layout.Column = 1;
                b = app.techButton(g, '光标选点', @(~, ~)app.pickMapPixelWithCursor(false)); b.Layout.Row = 5; b.Layout.Column = 2;
                b = app.techButton(g, '串行全量', @(~, ~)app.runDatacube(false, 'serial')); b.Layout.Row = 5; b.Layout.Column = 3;
                b = app.techButton(g, '并行全量', @(~, ~)app.runDatacube(false, 'parallel')); b.Layout.Row = 5; b.Layout.Column = 4;
                b = app.techButton(g, '导出结果', @(~, ~)app.exportResult('stage1')); b.Layout.Row = 5; b.Layout.Column = [5 6];
                b = app.techButton(g, '送信号回数据页', @(~, ~)app.importFitSignalToData()); b.Layout.Row = 5; b.Layout.Column = [7 8];
                lab = app.techLabel(g, '显示map');
                lab.Layout.Row = 5; lab.Layout.Column = 9;
                app.MapDropDown = uidropdown(g, 'Items', app.mapItems(), 'Value', 'R2');
                app.MapDropDown.ValueChangedFcn = @(~, ~)app.refreshResultMap(false);
                app.MapDropDown.Layout.Row = 5; app.MapDropDown.Layout.Column = [10 12];
                ctrl.mapDropDown = app.MapDropDown;
                app.Stage1Controls = ctrl;
            end
            note = uilabel(g, 'Text', '串行默认保存断点；并行先保存数值结果，诊断图建议按需导出。');
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

            lab = app.techLabel(g, '单点Y'); lab.Layout.Row = 1; lab.Layout.Column = 1;
            py = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            py.Layout.Row = 1; py.Layout.Column = 2;
            lab = app.techLabel(g, '单点X'); lab.Layout.Row = 1; lab.Layout.Column = 3;
            px = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            px.Layout.Row = 1; px.Layout.Column = 4;
            lab = app.techLabel(g, 'Y bin'); lab.Layout.Row = 1; lab.Layout.Column = 5;
            sy = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            sy.Layout.Row = 1; sy.Layout.Column = 6;
            lab = app.techLabel(g, 'X bin'); lab.Layout.Row = 1; lab.Layout.Column = 7;
            sx = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            sx.Layout.Row = 1; sx.Layout.Column = 8;

            lab = app.techLabel(g, '能量bin'); lab.Layout.Row = 2; lab.Layout.Column = 1;
            eb = uispinner(g, 'Limits', [1 Inf], 'Value', 1, 'Step', 1);
            eb.Layout.Row = 2; eb.Layout.Column = 2;
            lab = app.techLabel(g, '并行worker'); lab.Layout.Row = 2; lab.Layout.Column = 3;
            workers = uispinner(g, 'Limits', [1 64], 'Value', 8, 'Step', 1);
            workers.Layout.Row = 2; workers.Layout.Column = 4;
            checkpoint = uicheckbox(g, 'Text', '断点保存', 'Value', true);
            checkpoint.Layout.Row = 2; checkpoint.Layout.Column = [5 6];
            resume = uicheckbox(g, 'Text', '断点续跑', 'Value', false);
            resume.Layout.Row = 2; resume.Layout.Column = [7 8];

            lab = app.techLabel(g, '截取范围'); lab.Layout.Row = 3; lab.Layout.Column = 1;
            cropMin = uieditfield(g, 'numeric', 'Value', app.EnergyMinEdit.Value);
            cropMin.Layout.Row = 3; cropMin.Layout.Column = 2;
            cropMax = uieditfield(g, 'numeric', 'Value', app.EnergyMaxEdit.Value);
            cropMax.Layout.Row = 3; cropMax.Layout.Column = 3;
            inheritBtn = app.techButton(g, '继承数据页', @(~, ~)app.inheritStageEnergyWindow(stage2));
            inheritBtn.Layout.Row = 3; inheritBtn.Layout.Column = [4 5];
            saveplots = uicheckbox(g, 'Text', '保存诊断图', 'Value', false);
            saveplots.Layout.Row = 3; saveplots.Layout.Column = [6 8];

            lab = app.techLabel(g, '输出目录'); lab.Layout.Row = 4; lab.Layout.Column = 1;
            outdir = uieditfield(g, 'text', 'Value', app.LastOutputDir);
            outdir.Layout.Row = 4; outdir.Layout.Column = [2 8];

            if stage2
                lab = app.techLabel(g, '算法'); lab.Layout.Row = 5; lab.Layout.Column = 1;
                app.Stage2AlgorithmDropDown = uidropdown(g, ...
                    'Items', {'grid_scan'}, 'Value', 'grid_scan');
                app.Stage2AlgorithmDropDown.Layout.Row = 5;
                app.Stage2AlgorithmDropDown.Layout.Column = [2 4];
                hierarchical = uicheckbox(g, 'Text', '分层ZLP + ZLP界面引导', 'Value', true);
                hierarchical.Tooltip = '启用ZLP谱形聚类、类模板约束、Stage1空间平滑和联合温度MAP。';
                hierarchical.Layout.Row = 5; hierarchical.Layout.Column = [5 6];
                forceRawStage2 = uicheckbox(g, 'Text', '强制显示Stage2', 'Value', false);
                forceRawStage2.Tooltip = '勾选后关闭Stage2谱形保护，不回退到Stage1，用于查看硬BE强行拟合的真实形状。';
                forceRawStage2.Layout.Row = 5; forceRawStage2.Layout.Column = [7 8];
                app.Stage2PixelYSpinner = py; app.Stage2PixelXSpinner = px;
            else
                fullLabel = uilabel(g, 'Text', 'datacube 默认全空间；单点坐标只影响“拟合单点”。');
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
                b = app.techButton(g, '拟合单点', @(~, ~)app.runSingle(true)); b.Layout.Row = 6; b.Layout.Column = [1 2];
                b = app.techButton(g, '光标选点', @(~, ~)app.pickMapPixelWithCursor(true)); b.Layout.Row = 6; b.Layout.Column = [3 4];
                b = app.techButton(g, '串行全量', @(~, ~)app.runDatacube(true, 'serial')); b.Layout.Row = 6; b.Layout.Column = [5 6];
                b = app.techButton(g, '并行全量', @(~, ~)app.runDatacube(true, 'parallel')); b.Layout.Row = 6; b.Layout.Column = [7 8];
                b = app.techButton(g, '导出结果', @(~, ~)app.exportResult('stage2')); b.Layout.Row = 7; b.Layout.Column = [1 2];
                b = app.techButton(g, '送信号回数据页', @(~, ~)app.importFitSignalToData()); b.Layout.Row = 7; b.Layout.Column = [3 4];
                b = app.techButton(g, '2σ联合修正离群', @(~, ~)app.refitReliabilityOutliers()); b.Layout.Row = 7; b.Layout.Column = [5 6];
                lab = app.techLabel(g, '显示map'); lab.Layout.Row = 7; lab.Layout.Column = 7;
                app.Stage2MapDropDown = uidropdown(g, 'Items', app.mapItems(), 'Value', 'T data valid', ...
                    'ValueChangedFcn', @(~, ~)app.refreshResultMap(true));
                app.Stage2MapDropDown.Layout.Row = 7; app.Stage2MapDropDown.Layout.Column = 8;
                ctrl.hierarchical = hierarchical;
                ctrl.forceRawStage2 = forceRawStage2;
                ctrl.mapDropDown = app.Stage2MapDropDown;
                app.Stage2Controls = ctrl;
            else
                b = app.techButton(g, '拟合单点', @(~, ~)app.runSingle(false)); b.Layout.Row = 6; b.Layout.Column = [1 2];
                b = app.techButton(g, '光标选点', @(~, ~)app.pickMapPixelWithCursor(false)); b.Layout.Row = 6; b.Layout.Column = [3 4];
                b = app.techButton(g, '串行全量', @(~, ~)app.runDatacube(false, 'serial')); b.Layout.Row = 6; b.Layout.Column = [5 6];
                b = app.techButton(g, '并行全量', @(~, ~)app.runDatacube(false, 'parallel')); b.Layout.Row = 6; b.Layout.Column = [7 8];
                b = app.techButton(g, '导出结果', @(~, ~)app.exportResult('stage1')); b.Layout.Row = 7; b.Layout.Column = [1 2];
                b = app.techButton(g, '送信号回数据页', @(~, ~)app.importFitSignalToData()); b.Layout.Row = 7; b.Layout.Column = [3 4];
                lab = app.techLabel(g, '显示map'); lab.Layout.Row = 7; lab.Layout.Column = 5;
                app.MapDropDown = uidropdown(g, 'Items', app.mapItems(), 'Value', 'R2', ...
                    'ValueChangedFcn', @(~, ~)app.refreshResultMap(false));
                app.MapDropDown.Layout.Row = 7; app.MapDropDown.Layout.Column = [6 8];
                ctrl.mapDropDown = app.MapDropDown;
                app.Stage1Controls = ctrl;
            end

            note = uilabel(g, 'Text', '串行默认保存断点；并行先保存数值结果，诊断图建议按需导出。');
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

            p = uipanel(gl, 'Title', '诊断项', 'FontWeight', 'bold');
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
                '完整诊断图', '原始-总拟合', 'ZLP与分量', '稀疏卷积/振动', ...
                '残差与置信区间', '稀疏源强度棍棍图', 'log(Iloss/Igain)', '注意力曲线', '损失占比', ...
                'log信号 vs log ZLP'}, 'Value', '完整诊断图');
            if app.IsCompact
                app.DiagnosticDropDown.Layout.Row = 1; app.DiagnosticDropDown.Layout.Column = [1 3];
                b = uibutton(g, 'Text', '刷新当前项', 'ButtonPushedFcn', @(~, ~)app.plotDiagnosticSelection()); b.Layout.Row = 2; b.Layout.Column = 1;
                b = uibutton(g, 'Text', '完整弹窗', 'ButtonPushedFcn', @(~, ~)app.openFullDiagnostic()); b.Layout.Row = 2; b.Layout.Column = 2;
                b = uibutton(g, 'Text', '导出当前结果', 'ButtonPushedFcn', @(~, ~)app.exportResult('current')); b.Layout.Row = 2; b.Layout.Column = 3;
                b = uibutton(g, 'Text', '加载结果MAT', 'ButtonPushedFcn', @(~, ~)app.onLoadResult()); b.Layout.Row = 3; b.Layout.Column = 1;
                lab = uilabel(g, 'Text', '下拉菜单对应主诊断图中的分区内容。'); lab.Layout.Row = 3; lab.Layout.Column = [2 3];
            else
                uibutton(g, 'Text', '刷新当前项', 'ButtonPushedFcn', @(~, ~)app.plotDiagnosticSelection());
                uibutton(g, 'Text', '完整弹窗', 'ButtonPushedFcn', @(~, ~)app.openFullDiagnostic());
                uibutton(g, 'Text', '导出当前结果', 'ButtonPushedFcn', @(~, ~)app.exportResult('current'));
                uibutton(g, 'Text', '加载结果MAT', 'ButtonPushedFcn', @(~, ~)app.onLoadResult());
                uilabel(g, 'Text', '下拉菜单对应主诊断图中的分区内容。');
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
            app.AdvancedTable = uitable(gl, 'ColumnName', {'参数名', '中文说明', '当前值', '调整建议'}, ...
                'ColumnEditable', [false false true false], 'CellEditCallback', @(~, ~)app.applyAdvancedTable());
            app.AdvancedTable.Layout.Row = 1; app.AdvancedTable.Layout.Column = [1 4];
            app.AdvancedHelpArea = uitextarea(gl, 'Editable', 'off');
            app.AdvancedHelpArea.Layout.Row = 2; app.AdvancedHelpArea.Layout.Column = [1 2];
            app.AdvancedHelpArea.Value = { ...
                '高级参数说明：', ...
                '1. Stage1 主要管谱线复现和扣背地；优先调 EM 步长、ZLP宽度、tail/feature权重。', ...
                '2. Stage2 只在注意力窗口内用硬BE估温；温度不稳时先看注意力窗口和可靠性。', ...
                '3. Stage1/Stage2 datacube 默认全空间；单点坐标只用于“拟合单点”，断点续跑建议关闭并行，稳定优先。'};
            b1 = uibutton(gl, 'Text', '应用参数', 'ButtonPushedFcn', @(~, ~)app.applyAdvancedTable());
            b1.Layout.Row = 2; b1.Layout.Column = 3;
            b2 = uibutton(gl, 'Text', '重置默认', 'ButtonPushedFcn', @(~, ~)app.resetAdvanced());
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
            [file, path] = uigetfile({'*.mat;*.dm3', 'EELS data (*.mat, *.dm3)'; '*.mat', 'MAT files (*.mat)'; '*.dm3', 'DM3 files (*.dm3)'}, '选择 EELS datacube MAT/DM3');
            if isequal(file, 0), return; end
            filename = fullfile(path, file);
            d = uiprogressdlg(app.UIFigure, 'Title', '加载数据', 'Message', '正在读取 datacube...', 'Indeterminate', 'on');
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
                app.logStatus(sprintf('已加载 %s, size = [%s], energy points = %d', app.DataLabel, num2str(size(cube)), numel(ene)));
            catch ME
                app.showError(ME);
            end
            close(d);
        end

        function onLoadResult(app)
            [file, path] = uigetfile({'*.mat', 'MAT files (*.mat)'}, '选择 result MAT');
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
                    app.logStatus(['已加载 datacube result: ', filename]);
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
                    app.logStatus(['已加载 single result: ', filename]);
                elseif isfield(S, 'testres')
                    app.CurrentResult = app.firstValidResult(S.testres.results);
                    app.plotDiagnosticSelection();
                    app.logStatus(['已加载 simulation result: ', filename]);
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
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
            d = uiprogressdlg(app.UIFigure, 'Title', 'Align', 'Message', '正在调用 3D-EELS Align...', 'Indeterminate', 'on');
            try
                obj = app.makeEELSObject(app.Data, app.Ene);
                out = Align(obj, app.AlignMethodDropDown.Value, 0);
                app.Data = double(out.data);
                app.Ene = double(out.ene(:));
                app.clearFitCaches();
                app.setPixelLimits();
                app.refreshDataPreview();
                app.logStatus(['Align 完成，方法：', app.AlignMethodDropDown.Value]);
            catch ME
                app.showError(ME);
            end
            close(d);
        end

        function onBM3D(app)
            if ~app.hasExternalTool('BM3Ddenoise', 'BM3D denoise')
                return;
            end
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
            d = uiprogressdlg(app.UIFigure, 'Title', 'BM3D', 'Message', '正在调用 BM3D 降噪...', 'Indeterminate', 'on');
            try
                obj = app.makeEELSObject(app.Data, app.Ene);
                win = app.energyWindow();
                out = BM3Ddenoise(obj, win, app.BM3DFactorEdit.Value);
                app.Data = double(out.data);
                app.Ene = double(out.ene(:));
                app.clearFitCaches();
                app.refreshDataPreview();
                app.logStatus(sprintf('BM3D 完成，窗口 [%.1f, %.1f] meV。', win(1), win(2)));
            catch ME
                app.showError(ME);
            end
            close(d);
        end

        function runSingle(app, stage2)
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
            settings = app.collectSettings(stage2);
            if stage2
                settings = app.applyStage2SingleFastSettings(settings);
            end
            [yy, xx] = app.currentPixel(stage2);
            spec = squeeze(app.Data(yy, xx, :));
            d = uiprogressdlg(app.UIFigure, 'Title', '单点拟合', 'Message', '正在拟合当前像素...', 'Indeterminate', 'on');
            try
                reuse = stage2 && ~isempty(app.Stage1SingleResult) && ...
                    isfield(app.Stage1SingleResult, 'pixel_y') && isfield(app.Stage1SingleResult, 'pixel_x') && ...
                    app.Stage1SingleResult.pixel_y == yy && app.Stage1SingleResult.pixel_x == xx;
                if reuse
                    res = EELS_PearsonSparse_FullFit.promote_stage1_to_stage2(app.Stage1SingleResult, settings);
                    d.Message = '已复用当前单点 Stage1，正在计算 Stage2 温度层...';
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
                msg = sprintf('单点完成 (%d,%d): R2 %.6g, NRMSE %.4g, T %.1f K, rel %.2f', ...
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
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
            app.DatacubePauseRequested = false;
            settings = app.collectSettings(stage2, run_mode);
            if ~isfolder(settings.test_output_dir), mkdir(settings.test_output_dir); end
            input_file = app.writeGuiCube(settings.test_output_dir);
            d = uiprogressdlg(app.UIFigure, 'Title', 'Datacube处理', ...
                'Message', sprintf('正在以%s模式运行 datacube 拟合。', app.runModeLabel(run_mode)), ...
                'Indeterminate', 'off', 'Value', 0, 'Cancelable', 'on');
            settings.datacube_progress_callback = @(done,total,elapsed,eta,med,label)app.updateProgressDialog(d, done, total, elapsed, eta, med, label);
            settings.datacube_stop_callback = @()app.shouldPauseDatacube(d);
            cleanup = onCleanup(@()app.closeProgressDialog(d));
            try
                if stage2 && ~isempty(app.Stage1CubeResult) && app.isReusableStage1Cube(app.Stage1CubeResult)
                    d.Message = '已检测到 Stage1 datacube 结果，正在复用并计算 Stage2 温度层...';
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
                msg = sprintf('datacube 完成: success %d, median R2 %.6g, median NRMSE %.4g, output=%s', ...
                    nnz(cuberes.fit_status_map > 0), median(cuberes.R2_map(:), 'omitnan'), ...
                    median(cuberes.NRMSE_map(:), 'omitnan'), settings.test_output_dir);
                app.logStage(stage2, msg);
            catch ME
                if strcmp(ME.identifier, 'EELS:DatacubePaused')
                    app.logStage(stage2, sprintf('datacube 已暂停并保存断点：%s。勾选“断点续跑”并保持输出目录即可继续。', settings.test_output_dir));
                else
                    app.showError(ME);
                end
            end
            app.closeProgressDialog(d);
            clear cleanup
        end

        function refitReliabilityOutliers(app)
            if isempty(app.Stage2CubeResult) || ~isfield(app.Stage2CubeResult, 'reliability_score_map')
                app.logStatus('请先获得 Stage2 datacube 结果，再做2σ联合离群修正。');
                return;
            end
            cuberes = app.Stage2CubeResult;
            rel = cuberes.reliability_score_map;
            valid = isfinite(rel);
            if isfield(cuberes, 'fit_status_map') && isequal(size(cuberes.fit_status_map), size(rel))
                valid = valid & cuberes.fit_status_map > 0;
            end
            if nnz(valid) < 8
                app.logStatus('Stage2有效像素太少，无法做2σ联合离群判定。');
                return;
            end

            sigmaN = 2;
            [rel_mask, rel_stats, rel_z] = app.metricOutlierMask(rel, valid, sigmaN, 'lower');
            [r2_mask, r2_stats, r2_z] = app.metricOutlierMask(app.optionalMap(cuberes, 'R2_map', rel), valid, sigmaN, 'lower');
            [nrmse_mask, nrmse_stats, nrmse_z] = app.metricOutlierMask(app.optionalMap(cuberes, 'NRMSE_map', rel), valid, sigmaN, 'upper');
            [unc_mask, unc_stats, unc_z] = app.metricOutlierMask(app.optionalMap(cuberes, 'T_uncertainty_map', rel), valid, sigmaN, 'upper');
            mask = valid & (rel_mask | r2_mask | nrmse_mask | unc_mask);
            severity = rel_z + r2_z + nrmse_z + unc_z;
            metric_stats = struct('sigma_N', sigmaN, ...
                'reliability_low', rel_stats, 'R2_low', r2_stats, ...
                'NRMSE_high', nrmse_stats, 'T_uncertainty_high', unc_stats);

            stat_msg = sprintf('rel低:%d, R2低:%d, NRMSE高:%d, T不确定度高:%d', ...
                rel_stats.count, r2_stats.count, nrmse_stats.count, unc_stats.count);
            if ~any(mask(:))
                app.logStatus(['2σ联合检查完成：' stat_msg '，未发现需要重拟合的离群点。']);
                return;
            end

            settings = app.collectSettings(true);
            settings.run_stage2 = true;
            settings.stage2_temperature_algorithm = 'grid_scan';
            settings.save_plots = false;
            settings.stage2_pm_smooth_weight = max(settings.stage2_pm_smooth_weight, 7e-4);

            pix_list = find(mask(:)).';
            if isfield(cuberes, 'sample_pixels')
                [~, keep] = intersect(cuberes.sample_pixels(:).', pix_list, 'stable');
                pix_list = cuberes.sample_pixels(keep);
                result_idx = keep;
            else
                result_idx = pix_list;
            end
            if isempty(pix_list)
                app.logStatus('2σ联合检查发现的离群点不在当前结果列表中，未执行重拟合。');
                return;
            end
            pix_score = severity(pix_list);
            pix_score(~isfinite(pix_score)) = 0;
            [~, ord] = sort(pix_score, 'descend');
            pix_list = pix_list(ord);
            result_idx = result_idx(ord);
            max_refit = min(numel(pix_list), 80);
            pix_list = pix_list(1:max_refit);
            result_idx = result_idx(1:max_refit);
            d = uiprogressdlg(app.UIFigure, 'Title', '2σ联合离群修正', ...
                'Message', sprintf('正在重拟合 %d 个2σ联合离群像素...', numel(pix_list)), ...
                'Indeterminate', 'off', 'Value', 0, 'Cancelable', 'on');
            accepted = 0;
            attempted = 0;
            run_timer = tic;
            try
                for ii = 1:numel(pix_list)
                    if d.CancelRequested
                        break;
                    end
                    pix = pix_list(ii);
                    idx = result_idx(ii);
                    [yy, xx] = ind2sub(size(rel), pix);
                    if idx < 1 || idx > numel(cuberes.results)
                        continue;
                    end
                    old_res = cuberes.results{idx};
                    if ~(isstruct(old_res) && ~isfield(old_res, 'error'))
                        continue;
                    end
                    base_res = old_res;
                    if ~isempty(app.Stage1CubeResult) && isfield(app.Stage1CubeResult, 'results') && ...
                            isfield(app.Stage1CubeResult, 'sample_pixels')
                        sidx = find(app.Stage1CubeResult.sample_pixels(:).' == pix, 1, 'first');
                        if ~isempty(sidx) && sidx <= numel(app.Stage1CubeResult.results)
                            candidate_base = app.Stage1CubeResult.results{sidx};
                            if isstruct(candidate_base) && ~isfield(candidate_base, 'error')
                                base_res = candidate_base;
                            end
                        end
                    end
                    pixel_settings = settings;
                    attempted = attempted + 1;
                    try
                        new_res = EELS_PearsonSparse_FullFit.promote_stage1_to_stage2(base_res, pixel_settings);
                    catch
                        continue;
                    end
                    old_rel = app.resultMetricValue(old_res, 'reliability');
                    new_rel = app.resultMetricValue(new_res, 'reliability');
                    old_R2 = app.resultMetricValue(old_res, 'R2');
                    new_R2 = app.resultMetricValue(new_res, 'R2');
                    old_N = app.resultMetricValue(old_res, 'NRMSE');
                    new_N = app.resultMetricValue(new_res, 'NRMSE');
                    old_unc = app.resultMetricValue(old_res, 'T_uncertainty');
                    new_unc = app.resultMetricValue(new_res, 'T_uncertainty');

                    improves_rel = isfinite(new_rel) && (~isfinite(old_rel) || new_rel >= old_rel + 0.025);
                    improves_R2 = isfinite(new_R2) && (~isfinite(old_R2) || new_R2 >= old_R2 + 1e-4);
                    improves_N = isfinite(new_N) && (~isfinite(old_N) || new_N <= old_N - 5e-4);
                    improves_unc = isfinite(new_unc) && (~isfinite(old_unc) || new_unc <= old_unc - max(10, 0.08 .* abs(old_unc)));
                    preserves_R2 = ~isfinite(old_R2) || (isfinite(new_R2) && new_R2 >= old_R2 - 8e-4);
                    preserves_N = ~isfinite(old_N) || (isfinite(new_N) && new_N <= old_N + 0.003);
                    preserves_rel = ~isfinite(old_rel) || (isfinite(new_rel) && new_rel >= old_rel - 0.05);
                    preserves_unc = ~isfinite(old_unc) || ~isfinite(new_unc) || new_unc <= old_unc + max(25, 0.15 .* abs(old_unc));
                    preserves_fit = preserves_R2 && preserves_N && preserves_rel && preserves_unc;
                    if preserves_fit && (improves_rel || improves_R2 || improves_N || improves_unc)
                        outlier_flags = struct('reliability_low', rel_mask(pix), ...
                            'R2_low', r2_mask(pix), 'NRMSE_high', nrmse_mask(pix), ...
                            'T_uncertainty_high', unc_mask(pix));
                        new_res.metric_outlier_refit = struct('accepted', true, ...
                            'sigma_N', sigmaN, 'outlier_score', severity(pix), ...
                            'flags', outlier_flags, ...
                            'old_reliability', old_rel, 'new_reliability', new_rel, ...
                            'old_R2', old_R2, 'new_R2', new_R2, ...
                            'old_NRMSE', old_N, 'new_NRMSE', new_N, ...
                            'old_T_uncertainty', old_unc, 'new_T_uncertainty', new_unc, ...
                            'target_T_K', local_T, 'metric_stats', metric_stats);
                        new_res.reliability_outlier_refit = new_res.metric_outlier_refit;
                        cuberes.results{idx} = new_res;
                        accepted = accepted + 1;
                    end
                    d.Value = ii ./ max(numel(pix_list), 1);
                    d.Message = sprintf('2σ联合离群修正: %d/%d, accepted %d', ii, numel(pix_list), accepted);
                    drawnow limitrate;
                end
                if isfield(cuberes, 'sample_pixels')
                    pixel_time = NaN(1, numel(cuberes.sample_pixels));
                else
                    pixel_time = NaN(1, numel(cuberes.results));
                end
                cuberes = EELS_PearsonSparse_FullFit.rebuild_datacube_maps_from_results(cuberes, settings, pixel_time, run_timer, true);
                cuberes.metric_outlier_refit = struct('sigma_N', sigmaN, ...
                    'attempted', attempted, 'accepted', accepted, ...
                    'candidate_count', nnz(mask), 'metric_stats', metric_stats, ...
                    'candidate_mask', mask, 'severity_map', severity);
                cuberes.reliability_outlier_refit = cuberes.metric_outlier_refit;
                app.Stage2CubeResult = cuberes;
                app.CubeResult = cuberes;
                app.CurrentResult = app.firstValidResult(cuberes.results);
                app.refreshResultMap(true);
                app.plotDiagnosticSelection();
                if ~isfolder(settings.test_output_dir), mkdir(settings.test_output_dir); end
                save(fullfile(settings.test_output_dir, 'stage2_metric_outlier_refit_result.mat'), 'cuberes', '-v7.3');
                app.logStage(true, sprintf('2σ联合离群修正完成：%s；候选 %d，尝试 %d，接受 %d。', stat_msg, nnz(mask), attempted, accepted));
            catch ME
                app.showError(ME);
            end
            close(d);
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
                    settings.datacube_resume_checkpoint = false;
                    settings.close_parallel_pool_after_run = true;
                otherwise
                    settings.datacube_parallel_enabled = false;
                    settings.datacube_checkpoint_enabled = true;
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
                label = '并行';
            else
                label = '串行';
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
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
            app.plotEnergyMap();
            app.logStatus('请在能量过滤像上拖拽矩形 ROI，双击或回车确认。');
            try
                roi = drawrectangle(app.DataAxes, 'Color', [0.05 0.65 1.00], 'LineWidth', 1.5);
                [pos, cancelled] = app.waitForROIWithHotkeys(roi, ...
                    'ROI框选：双击确认；Ctrl+Q取消；Ctrl+H直接确认并作图。');
                if cancelled || isempty(pos)
                    app.logStatus('ROI框选已取消。');
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
                app.logStatus('ROI无效，请重新框选。');
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
            app.logStatus(sprintf('ROI积分完成：y %d-%d, x %d-%d；可点击“导出ROI”。', y1, y2, x1, x2));
        end

        function exportROI(app)
            if isempty(app.SelectedROICube)
                app.logStatus('还没有ROI结果。请先点击“框选ROI积分”。');
                return;
            end
            assignin('base', 'eels_roi_cube', app.SelectedROICube);
            assignin('base', 'eels_roi_spectrum', app.SelectedROISpectrum);
            assignin('base', 'eels_roi_energy', app.Ene);
            assignin('base', 'eels_roi_range_y1y2x1x2', app.SelectedROIRange);
            app.logStatus('已导出 eels_roi_cube / eels_roi_spectrum / eels_roi_energy / eels_roi_range_y1y2x1x2 到工作区。');
        end

        function openStackSpectraTool(app)
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
            app.plotEnergyMap();
            answer = inputdlg({'需要叠谱的数量', '选择方式：cursor 或 roi'}, ...
                '叠谱工具', [1 42], {'3', 'cursor'});
            if isempty(answer), return; end
            n_pick = max(1, round(str2double(answer{1})));
            if ~isfinite(n_pick), n_pick = 3; end
            mode = lower(strtrim(answer{2}));
            if ~ismember(mode, {'cursor', 'roi'}), mode = 'cursor'; end
            app.deleteGraphicsHandles('StackCursorHandles');
            app.logStatus('叠谱选择：双击确认每个点/ROI；Ctrl+Q取消；Ctrl+H直接确认并作图。叠谱专用光标为紫/蓝色。');
            stack = struct('name', {{}}, 'spectrum', {{}}, 'selection', {{}}, ...
                'mode', mode, 'ene', app.Ene(:), 'source_file', app.FileEdit.Value);
            hold(app.DataAxes, 'on');
            for kk = 1:n_pick
                try
                    if strcmp(mode, 'roi')
                        roi = drawrectangle(app.DataAxes, 'Color', [0.00 0.55 1.00], 'LineWidth', 1.4);
                        [pos, cancelled] = app.waitForROIWithHotkeys(roi, ...
                            sprintf('叠谱 %d/%d：双击确认；Ctrl+Q取消；Ctrl+H直接确认并作图。', kk, n_pick));
                        if cancelled || isempty(pos), break; end
                        [spec, selection] = app.spectrumFromROIPosition(pos);
                        hmark = rectangle(app.DataAxes, 'Position', selection.position, ...
                            'EdgeColor', [0.00 0.48 0.95], 'LineStyle', '--', 'LineWidth', 1.2);
                        app.appendGraphicsHandle('StackCursorHandles', hmark);
                    else
                        roi = drawpoint(app.DataAxes, 'Color', [0.55 0.00 0.95], 'LineWidth', 1.6);
                        [pos, cancelled] = app.waitForROIWithHotkeys(roi, ...
                            sprintf('叠谱 %d/%d：双击确认；Ctrl+Q取消；Ctrl+H直接确认并作图。', kk, n_pick));
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
                name_answer = inputdlg({'谱线名称'}, '命名谱线', [1 42], {default_name});
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
            app.logStatus('叠谱结果已导出到工作区 eels_stacked_spectra。');
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
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
            method = app.BGMethodDropDown.Value;
            win = [app.BGWin1MinEdit.Value, app.BGWin1MaxEdit.Value, ...
                app.BGWin2MinEdit.Value, app.BGWin2MaxEdit.Value];
            if any(~isfinite(win)) || win(1) >= win(2) || win(3) >= win(4)
                app.logStatus('传统扣背地窗口无效。请按 3D-EELS Toolbox 习惯输入 [w1min w1max w2min w2max]。');
                return;
            end
            d = uiprogressdlg(app.UIFigure, 'Title', '传统扣背地', ...
                'Message', sprintf('正在调用原版 BackgroundSubtraction(..., windows, %s)...', method), 'Indeterminate', 'on');
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
                app.logStatus(sprintf('原版 BackgroundSubtraction 扣背地完成：%s, windows [%.1f %.1f %.1f %.1f] meV。', ...
                    method, win(1), win(2), win(3), win(4)));
            catch ME
                app.showError(ME);
            end
            if ~isempty(d) && isvalid(d), close(d); end
        end

        function importFitSignalToData(app)
            if isempty(app.CurrentResult) && isempty(app.CubeResult)
                app.logStatus('还没有拟合结果可导入。');
                return;
            end
            try
                if ~isempty(app.CubeResult) && isfield(app.CubeResult, 'results')
                    first = app.firstValidResult(app.CubeResult.results);
                    if isempty(first), app.logStatus('cuberes 中没有有效单点结果。'); return; end
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
                app.logStatus('已将拟合后的 signal/vib_only 导入到第一个窗口，可继续做ROI、传统处理或导出。');
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
                case '导入MAT'
                    app.onLoadData();
                case '导入结果'
                    app.onLoadResult();
                case '导出数据'
                    app.exportData();
                case '导出ROI'
                    app.exportROI();
                case '导出当前结果'
                    app.exportResult('current');
            end
        end

        function runDataToolAction(app)
            switch app.DataToolDropDown.Value
                case 'NMF分解'
                    app.openNMFWindow();
                case '扣STA'
                    app.applySTACorrection();
                case '传统扣背地(交互)'
                    app.applyTraditionalBackground(true);
                case 'Separate IF mode'
                    app.openSeparateIFMode();
                case 'MLLS残差'
                    app.openMLLSWindow();
                case 'Peak fitting'
                    app.openPeakFitting();
                case '交互漂移矫正'
                    app.openDriftCorrection();
                case '声子源强度分析'
                    app.openPMAnalysisWindow();
                case '叠谱工具'
                    app.openStackSpectraTool();
            end
        end

        function applySTACorrection(app)
            if ~app.hasExternalTool('CorrectSTA', 'STA correction')
                return;
            end
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
            try
                out = CorrectSTA(app.makeEELSObject(app.Data, app.Ene));
                app.Data = double(out.data);
                app.Ene = double(out.ene(:));
                app.clearFitCaches();
                app.setPixelLimits();
                app.refreshDataPreview();
                app.logStatus('STA 校正完成：已调用原版 CorrectSTA 并更新当前数据。');
            catch ME
                app.showError(ME);
            end
        end

        function openSeparateIFMode(app)
            if ~app.hasExternalTool('SeparateIFmode_Interactive', 'Separate IF mode')
                return;
            end
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
            try
                SeparateIFmode_Interactive(app.makeEELSObject(app.Data, app.Ene));
                app.logStatus('已打开 Separate IF mode 原版交互窗口。');
            catch ME
                app.showError(ME);
            end
        end

        function openPeakFitting(app)
            if ~app.hasExternalTool('PeakFitting', 'Peak fitting')
                return;
            end
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
            try
                PeakFitting(app.makeEELSObject(app.Data, app.Ene));
                app.logStatus('已打开原版 PeakFitting 交互窗口。');
            catch ME
                app.showError(ME);
            end
        end

        function openDriftCorrection(app)
            if ~app.hasExternalTool('CorrDrift_Custom', 'Interactive drift correction')
                return;
            end
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
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
                app.logStatus('交互漂移矫正窗口已运行；如果选择了 drift 文件，当前数据已更新。');
            catch ME
                app.showError(ME);
            end
        end

        function openNMFWindow(app)
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
            fig = uifigure('Name', 'NMF 分解', 'Position', [120 120 1060 680], 'Color', [0.94 0.97 1.00]);
            gl = uigridlayout(fig, [3 4]);
            gl.RowHeight = {48, '1x', '1x'};
            gl.ColumnWidth = {180, 180, 160, '1x'};
            gl.Padding = [12 12 12 12];
            uilabel(gl, 'Text', '成分数', 'FontWeight', 'bold');
            rankEdit = uieditfield(gl, 'numeric', 'Value', 3, 'Limits', [1 20]);
            uilabel(gl, 'Text', '能量窗使用主界面设置', 'FontColor', [0.1 0.25 0.42]);
            axMap = uiaxes(gl); axMap.Layout.Row = 2; axMap.Layout.Column = [1 2]; app.styleAxes(axMap);
            axSpec = uiaxes(gl); axSpec.Layout.Row = [2 3]; axSpec.Layout.Column = [3 4]; app.styleAxes(axSpec);
            axCoef = uiaxes(gl); axCoef.Layout.Row = 3; axCoef.Layout.Column = [1 2]; app.styleAxes(axCoef);
            result = struct();
            b = app.techButton(gl, '运行NMF', @(~, ~)runNMF()); b.Layout.Row = 1; b.Layout.Column = 4;
            b2 = app.techButton(gl, '导出NMF结果', @(~, ~)exportNMF()); b2.Layout.Row = 1; b2.Layout.Column = 3;

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
                app.logStatus('已导出 eels_nmf_result 到 base workspace。');
            end
        end

        function openMLLSWindow(app)
            if isempty(app.Data), app.logStatus('请先加载 datacube。'); return; end
            fig = uifigure('Name', 'ROI / Reference MLLS residual', 'Position', [160 150 1120 700], 'Color', [0.94 0.97 1.00]);
            gl = uigridlayout(fig, [3 4]);
            gl.RowHeight = {48, '1x', '1x'};
            gl.ColumnWidth = {180, 220, 180, '1x'};
            vars = evalin('base', 'who');
            uilabel(gl, 'Text', '参考谱矩阵变量(energy x nRef)', 'FontWeight', 'bold');
            dd = uidropdown(gl, 'Items', [{'当前ROI积分谱'}; vars(:)]);
            uilabel(gl, 'Text', '能量窗使用主界面设置');
            axMap = uiaxes(gl); axMap.Layout.Row = [2 3]; axMap.Layout.Column = [1 2]; app.styleAxes(axMap);
            axSpec = uiaxes(gl); axSpec.Layout.Row = [2 3]; axSpec.Layout.Column = [3 4]; app.styleAxes(axSpec);
            result = struct();
            b = app.techButton(gl, '运行MLLS', @(~, ~)runMLLS()); b.Layout.Row = 1; b.Layout.Column = 4;
            b2 = app.techButton(gl, '导出MLLS结果', @(~, ~)exportMLLS()); b2.Layout.Row = 1; b2.Layout.Column = 3;
            function runMLLS()
                if isempty(dd.Value), return; end
                if strcmp(dd.Value, '当前ROI积分谱')
                    if isempty(app.SelectedROISpectrum)
                        app.logStatus('MLLS 需要先在第一个页签框选 ROI，或在 workspace 中选择参考谱矩阵。');
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
                app.logStatus('已导出 eels_mlls_residual_result 到 base workspace。');
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
            fig = uifigure('Name', '稀疏声子模态源强度 / Loss-gain DOS-like 工具', ...
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

            panel = uipanel(gl, 'Title', '源强度插值、选区、NMF与聚类', ...
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

            lab = app.techLabel(pg, '强度来源'); lab.Layout.Row = 1; lab.Layout.Column = 1;
            sourceDD = uidropdown(pg, 'Items', {'stage2', 'stage1', 'current'}, 'Value', defaultSource);
            sourceDD.Layout.Row = 1; sourceDD.Layout.Column = 2;
            lab = app.techLabel(pg, '插值范围'); lab.Layout.Row = 1; lab.Layout.Column = 3;
            emMinEdit = uieditfield(pg, 'numeric', 'Value', -250);
            emMinEdit.Layout.Row = 1; emMinEdit.Layout.Column = 4;
            emMaxEdit = uieditfield(pg, 'numeric', 'Value', 250);
            emMaxEdit.Layout.Row = 1; emMaxEdit.Layout.Column = 5;
            if app.IsCompact
                lab = app.techLabel(pg, 'meV'); lab.Layout.Row = 1; lab.Layout.Column = 6;
                lab = app.techLabel(pg, '步长'); lab.Layout.Row = 2; lab.Layout.Column = 1;
            else
                lab = app.techLabel(pg, '步长'); lab.Layout.Row = 1; lab.Layout.Column = 6;
            end
            stepEdit = uieditfield(pg, 'numeric', 'Value', 1, 'Limits', [0.1 20]);
            if app.IsCompact
                stepEdit.Layout.Row = 2; stepEdit.Layout.Column = 2;
                lab = app.techLabel(pg, '选择窗口'); lab.Layout.Row = 2; lab.Layout.Column = 3;
            else
                stepEdit.Layout.Row = 1; stepEdit.Layout.Column = 7;
                lab = app.techLabel(pg, '选择窗口'); lab.Layout.Row = 1; lab.Layout.Column = 8;
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
                lab = app.techLabel(pg, '聚类数'); lab.Layout.Row = 3; lab.Layout.Column = 2;
            else
                rankEdit.Layout.Row = 1; rankEdit.Layout.Column = 12;
                lab = app.techLabel(pg, '聚类数'); lab.Layout.Row = 2; lab.Layout.Column = 1;
            end
            clusterEdit = uieditfield(pg, 'numeric', 'Value', 3, 'Limits', [2 20]);
            if app.IsCompact
                clusterEdit.Layout.Row = 3; clusterEdit.Layout.Column = 3;
                lab = app.techLabel(pg, '显示平滑'); lab.Layout.Row = 3; lab.Layout.Column = 4;
            else
                clusterEdit.Layout.Row = 2; clusterEdit.Layout.Column = 2;
                lab = app.techLabel(pg, '显示平滑'); lab.Layout.Row = 2; lab.Layout.Column = 3;
            end
            smoothEdit = uieditfield(pg, 'numeric', 'Value', 3, 'Limits', [1 51]);
            if app.IsCompact
                smoothEdit.Layout.Row = 3; smoothEdit.Layout.Column = 5;
            else
                smoothEdit.Layout.Row = 2; smoothEdit.Layout.Column = 4;
            end
            smoothCheck = uicheckbox(pg, 'Text', '仅平滑显示/导出不改原源强度', 'Value', true);
            if app.IsCompact
                smoothCheck.Layout.Row = 3; smoothCheck.Layout.Column = 6;
            else
                smoothCheck.Layout.Row = 2; smoothCheck.Layout.Column = [5 8];
            end
            sparseNote = uilabel(pg, 'Text', 'v2.1：attention 可信对约束温度，其余能区由软声子源强度重建；分层ZLP保护谱形。', ...
                'FontColor', [0.55 0.16 0.02]);
            if app.IsCompact
                sparseNote.Layout.Row = 4; sparseNote.Layout.Column = [1 6];
            else
                sparseNote.Layout.Row = 2; sparseNote.Layout.Column = [9 12];
            end
            note = uilabel(pg, 'Text', '有符号源强度：gain/anti-Stokes 为负能量，loss/Stokes 为正能量；它是 DOS-like 源空间表示，不是绝对标定的声子 DOS。', ...
                'FontColor', [0.08 0.25 0.42]);
            if app.IsCompact
                note.Layout.Row = 5; note.Layout.Column = [1 6];
            else
                note.Layout.Row = 3; note.Layout.Column = [1 12];
            end

            tabs = uitabgroup(gl);
            tabs.Layout.Row = 2; tabs.Layout.Column = 1;
            tabOverview = uitab(tabs, 'Title', '总览 / 选区积分');
            tabNMF = uitab(tabs, 'Title', '源强度 NMF');
            tabCluster = uitab(tabs, 'Title', '聚类');

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
                b = app.techButton(bg, '构建源强度插值', @(~, ~)buildPM()); b.Layout.Row = 1; b.Layout.Column = 1;
                b = app.techButton(bg, '导出插值强度', @(~, ~)exportPM()); b.Layout.Row = 1; b.Layout.Column = 2;
                b = app.techButton(bg, '导出所选强度', @(~, ~)exportSelectedPM()); b.Layout.Row = 1; b.Layout.Column = 3;
                b = app.techButton(bg, '积分模式map', @(~, ~)plotModeMap()); b.Layout.Row = 2; b.Layout.Column = 1;
                b = app.techButton(bg, '源强度 NMF', @(~, ~)runPMNMF()); b.Layout.Row = 2; b.Layout.Column = 2;
                b = app.techButton(bg, '聚类分析', @(~, ~)runPMCluster()); b.Layout.Row = 2; b.Layout.Column = 3;
            else
                b = app.techButton(bg, '构建源强度插值', @(~, ~)buildPM()); b.Layout.Row = 1; b.Layout.Column = 1;
                b = app.techButton(bg, '导出插值强度', @(~, ~)exportPM()); b.Layout.Row = 1; b.Layout.Column = 2;
                b = app.techButton(bg, '导出所选强度', @(~, ~)exportSelectedPM()); b.Layout.Row = 1; b.Layout.Column = 3;
                b = app.techButton(bg, '积分模式map', @(~, ~)plotModeMap()); b.Layout.Row = 1; b.Layout.Column = 4;
                b = app.techButton(bg, '源强度 NMF', @(~, ~)runPMNMF()); b.Layout.Row = 1; b.Layout.Column = 5;
                b = app.techButton(bg, '聚类分析', @(~, ~)runPMCluster()); b.Layout.Row = 1; b.Layout.Column = 6;
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
                    title(axSpec, sprintf('平均源强度，来源=%s，有效像素=%d', pmres.source, nnz(pmres.valid_mask)));
                    plotModeMap();
                    cla(axSummary);
                    bar(axSummary, [nnz(pmres.axis_meV < 0), nnz(pmres.axis_meV == 0), nnz(pmres.axis_meV > 0)]);
                    set(axSummary, 'XTickLabel', {'gain', 'zero', 'loss'});
                    ylabel(axSummary, '插值bin数'); title(axSummary, '有符号源强度轴覆盖'); grid(axSummary, 'on');
                catch ME
                    app.showError(ME);
                end
            end

            function plotModeMap()
                if isempty(fieldnames(pmres)), buildPM(); end
                if isempty(fieldnames(pmres)), return; end
                [modeMap, useAxis] = selectedPMIntegral(pmres);
                imagesc(axMap, modeMap); axis(axMap, 'image'); colorbar(axMap);
                title(axMap, sprintf('源强度积分图 [%g, %g] meV', min(useAxis), max(useAxis)));
                xlabel(axMap, 'x'); ylabel(axMap, 'y');
                assignin('base', 'eels_pm_mode_map', struct('map', modeMap, ...
                    'axis_meV', useAxis, 'source', pmres.source));
                app.logStatus('源强度积分图已生成并导出。');
            end

            function exportPM()
                if isempty(fieldnames(pmres)), buildPM(); end
                assignin('base', 'eels_pm_interpolated', pmres);
                app.logStatus('已导出源强度插值结果：包含有符号模态能量轴、强度 cube、均值谱和有效像素 mask。');
            end

            function exportSelectedPM()
                if isempty(fieldnames(pmres)), buildPM(); end
                range = sort([selMinEdit.Value, selMaxEdit.Value]);
                mask = pmres.axis_meV >= range(1) & pmres.axis_meV <= range(2);
                if ~any(mask)
                    app.logStatus('所选源强度能量窗口为空，请调整选择范围。');
                    return;
                end
                selected = pmres;
                selected.axis_meV = pmres.axis_meV(mask);
                selected.pm_cube = pmres.pm_cube(:, :, mask);
                selected.mean_PM = pmres.mean_PM(mask);
                selected.selection_range = range;
                assignin('base', 'eels_pm_selected', selected);
                app.logStatus('已导出所选源强度窗口，可用于手动谱学或空间特征分析。');
            end

            function runPMNMF()
                if isempty(fieldnames(pmres)), buildPM(); end
                X = reshape(pmres.pm_cube, [], numel(pmres.axis_meV));
                valid = pmres.valid_mask(:) & any(isfinite(X) & X > 0, 2);
                X = X(valid, :);
                X(~isfinite(X)) = 0;
                X = max(X, 0);
                if isempty(X)
                    app.logStatus('没有可用于 NMF 的有效源强度数据。');
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
                title(axNMFMap, '源强度 NMF component 1 map');
                plot(axNMFSpectrum, pmres.axis_meV, H', 'LineWidth', 1.2); grid(axNMFSpectrum, 'on');
                xlabel(axNMFSpectrum, 'Signed mode energy (meV)'); ylabel(axNMFSpectrum, 'Component source intensity');
                title(axNMFSpectrum, '源强度 NMF component spectra');
                assignin('base', 'eels_pm_nmf_result', nmfres);
                app.logStatus('源强度 NMF 已完成并导出。');
            end

            function runPMCluster()
                if isempty(fieldnames(pmres)), buildPM(); end
                X = reshape(pmres.pm_cube, [], numel(pmres.axis_meV));
                valid = pmres.valid_mask(:) & any(isfinite(X) & X > 0, 2);
                X = X(valid, :);
                X(~isfinite(X)) = 0;
                X = max(X, 0);
                if isempty(X)
                    app.logStatus('没有可用于聚类的有效源强度数据。');
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
                title(axCluster, sprintf('源强度聚类，k=%d', k));
                cla(axClusterSummary); hold(axClusterSummary, 'on');
                for cc = 1:k
                    plot(axClusterSummary, pmres.axis_meV, mean(X(labels == cc, :), 1, 'omitnan'), 'LineWidth', 1.2);
                end
                hold(axClusterSummary, 'off'); grid(axClusterSummary, 'on');
                xlabel(axClusterSummary, 'Signed mode energy (meV)'); ylabel(axClusterSummary, 'Cluster mean source intensity');
                title(axClusterSummary, '各簇平均源强度谱');
                assignin('base', 'eels_pm_cluster_result', clusterres);
                app.logStatus('源强度聚类已完成并导出。');
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
                error('没有可分析的 %s datacube/单点拟合结果。请先运行或加载Stage1/Stage2。', char(source));
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
                eta_txt = 'ETA 估算中';
            end
            if ~isfinite(med), med = NaN; end
            if app.shouldPauseDatacube(dlg)
                pause_txt = ' | 已请求暂停，当前像素结束后保存断点';
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
                app.logStage(stage2, '还没有 datacube map 可选点。');
                return;
            end
            if stage2, ax = app.Stage2MapAxes; else, ax = app.Stage1MapAxes; end
            try
                app.logStage(stage2, '光标选点：双击确认；Ctrl+Q取消；Ctrl+H直接确认并作图。');
                roi = drawpoint(ax, 'Color', [0.95 0.20 0.05], 'LineWidth', 1.4);
                [pos, cancelled] = app.waitForROIWithHotkeys(roi, ...
                    '光标选点：双击确认；Ctrl+Q取消；Ctrl+H直接确认并作图。');
                if cancelled || isempty(pos)
                    app.logStage(stage2, '光标选点已取消。');
                    return;
                end
            catch ME
                app.logStage(stage2, ['光标选点需要 MATLAB ROI 工具：', ME.message]);
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
                case '完整诊断图'
                    app.openFullDiagnostic();
                case '原始-总拟合'
                    plot(app.ViewAxes1, ene, res.raw_data, 'k.', ene, res.fitted_spectrum, 'r-', 'LineWidth', 1.2);
                    legend(app.ViewAxes1, {'raw', 'fit'}); grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Raw vs total fit');
                    plot(app.ViewAxes2, ene, res.raw_data(:) - res.fitted_spectrum(:), 'b-');
                    yline(app.ViewAxes2, 0, 'k-'); grid(app.ViewAxes2, 'on'); title(app.ViewAxes2, 'Residual');
                case 'ZLP与分量'
                    plot(app.ViewAxes1, ene, res.zlp_only, 'Color', [0.05 0.45 0.18], 'LineWidth', 1.4); hold(app.ViewAxes1, 'on');
                    if isfield(res, 'zlp_components') && ~isempty(res.zlp_components)
                        plot(app.ViewAxes1, ene, res.zlp_components, 'LineWidth', 0.8);
                    end
                    hold(app.ViewAxes1, 'off'); grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'ZLP envelope and components');
                    semilogy(app.ViewAxes2, ene, max(res.raw_data, eps), 'k.', ene, max(res.zlp_only, eps), 'g-');
                    grid(app.ViewAxes2, 'on'); title(app.ViewAxes2, 'Log ZLP check');
                case '稀疏卷积/振动'
                    plot(app.ViewAxes1, ene, res.vib_only, 'Color', [0.9 0.25 0.05], 'LineWidth', 1.2);
                    grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Sparse convolution / vib');
                    app.plotPM(app.ViewAxes2, res);
                case '残差与置信区间'
                    residual = res.raw_data(:) - res.fitted_spectrum(:);
                    area(app.ViewAxes1, ene, residual, 'FaceColor', [0.45 0.62 0.95], 'EdgeColor', [0.1 0.35 0.95]);
                    yline(app.ViewAxes1, 0, 'k-'); grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Residual');
                    if isfield(res, 'confidence_interval') && isfield(res.confidence_interval, 'half_width')
                        plot(app.ViewAxes2, ene, res.confidence_interval.half_width, 'b-');
                        grid(app.ViewAxes2, 'on'); title(app.ViewAxes2, 'Energy-dependent CI half-width');
                    end
                case '稀疏源强度棍棍图'
                    app.plotPM(app.ViewAxes1, res);
                case 'log(Iloss/Igain)'
                    app.plotSpectrumLossGainRatio(app.ViewAxes1, res);
                    app.plotPM(app.ViewAxes2, res);
                case '注意力曲线'
                    if isfield(res, 'stage2') && isfield(res.stage2, 'attention')
                        plot(app.ViewAxes1, ene, res.stage2.attention.feature_score, 'b-');
                        ylim(app.ViewAxes1, [0 1]); grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Attention feature score');
                    end
                case '损失占比'
                    if isfield(res, 'loss_breakdown')
                        vals = 100 .* res.loss_breakdown.fraction;
                        labels = strrep(res.loss_breakdown.names, '_', ' ');
                        barh(app.ViewAxes1, vals); set(app.ViewAxes1, 'YTick', 1:numel(labels), 'YTickLabel', labels);
                        grid(app.ViewAxes1, 'on'); title(app.ViewAxes1, 'Loss budget (%)');
                    end
                case 'log信号 vs log ZLP'
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
            app.populateAdvancedTable();
            app.logStatus('高级参数已重置为默认值。');
        end

        function exportData(app)
            if isempty(app.Data), return; end
            assignin('base', 'eels_gui_data', app.Data);
            assignin('base', 'eels_gui_energy', app.Ene);
            app.logStatus('已导出 eels_gui_data / eels_gui_energy 到 base workspace。');
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
                app.logStatus(['没有可导出的 ', char(tag), ' 结果。请先运行或加载对应阶段结果。']);
            else
                app.logStatus(['已导出到 base workspace：', strjoin(exported, ' / ')]);
            end
            app.logStatus(['已导出结果到 base workspace，标签：', tag]);
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
            app.logStage(stage2, sprintf('已继承数据页能量窗：[%.1f, %.1f] meV。', win(1), win(2)));
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

        function val = localMedianMapValue(~, map, yy, xx, exclude_mask)
            val = NaN;
            if isempty(map) || yy < 1 || xx < 1 || yy > size(map, 1) || xx > size(map, 2)
                return;
            end
            y1 = max(1, yy - 1); y2 = min(size(map, 1), yy + 1);
            x1 = max(1, xx - 1); x2 = min(size(map, 2), xx + 1);
            block = map(y1:y2, x1:x2);
            if nargin >= 5 && ~isempty(exclude_mask) && isequal(size(exclude_mask), size(map))
                ex = exclude_mask(y1:y2, x1:x2);
                block(ex) = NaN;
            end
            vals = block(isfinite(block));
            if isempty(vals)
                vals = map(isfinite(map));
            end
            if ~isempty(vals)
                val = median(vals, 'omitnan');
            end
        end

        function [mask, stats, zmap] = metricOutlierMask(~, map, valid, sigmaN, side)
            mask = false(size(valid));
            zmap = zeros(size(valid));
            stats = struct('available', false, 'side', side, 'mean', NaN, 'std', NaN, ...
                'threshold', NaN, 'count', 0);
            if isempty(map) || ~isequal(size(map), size(valid))
                return;
            end
            map = double(map);
            valid_map = valid & isfinite(map);
            vals = map(valid_map);
            if numel(vals) < 8
                return;
            end
            mu = mean(vals, 'omitnan');
            sig = std(vals, 0, 'omitnan');
            if ~(isfinite(mu) && isfinite(sig) && sig > 0)
                return;
            end
            stats.available = true;
            stats.mean = mu;
            stats.std = sig;
            switch lower(side)
                case 'lower'
                    thr = mu - sigmaN .* sig;
                    mask = valid_map & map < thr;
                    zmap(valid_map) = max(0, (mu - map(valid_map)) ./ sig);
                otherwise
                    thr = mu + sigmaN .* sig;
                    mask = valid_map & map > thr;
                    zmap(valid_map) = max(0, (map(valid_map) - mu) ./ sig);
            end
            stats.threshold = thr;
            stats.count = nnz(mask);
        end

        function map = optionalMap(~, s, field_name, template)
            if isstruct(s) && isfield(s, field_name) && isequal(size(s.(field_name)), size(template))
                map = s.(field_name);
            else
                map = NaN(size(template));
            end
        end

        function val = resultMetricValue(~, res, name)
            val = NaN;
            if ~(isstruct(res) && ~isfield(res, 'error'))
                return;
            end
            switch lower(name)
                case 'reliability'
                    if isfield(res, 'reliability_score')
                        val = res.reliability_score;
                    end
                case 'r2'
                    if isfield(res, 'metrics') && isfield(res.metrics, 'R2')
                        val = res.metrics.R2;
                    end
                case 'nrmse'
                    if isfield(res, 'metrics') && isfield(res.metrics, 'NRMSE')
                        val = res.metrics.NRMSE;
                    end
                case 't_uncertainty'
                    if isfield(res, 'T_uncertainty')
                        val = res.T_uncertainty;
                    elseif isfield(res, 'T_uncertainty_K')
                        val = res.T_uncertainty_K;
                    elseif isfield(res, 'error_budget') && isfield(res.error_budget, 'total_T_uncertainty_K')
                        val = res.error_budget.total_T_uncertainty_K;
                    end
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
            app.logStatus(['错误：', ME.message]);
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
            switch key
                case 'run_stage2', name = '启用第二轮温度拟合';
                case 'crop_range', name = '参与拟合的能量截取区间';
                case 'zlp_model_type', name = '零峰模型类型';
                case 'num_zlp_components', name = '零峰分量数量';
                case 'pearson_m_bounds', name = 'Pearson形状指数上下限';
                case 'pvoigt_eta_bounds', name = 'pseudo-Voigt混合比例上下限';
                case 'zlp_component_width_max_meV', name = '每个零峰分量的宽度上限';
                case 'zlp_component_width_min_meV', name = '每个零峰分量的宽度下限';
                case 'stage2_temperature_algorithm', name = '第二轮温度优化算法';
                case 'stage2_T_grid', name = '第二轮网格扫温的温度列表';
                case 'stage2_window_meV', name = '第二轮默认温度证据窗口';
                case 'stage2_attention_window_meV', name = '第二轮注意力候选窗口';
                case 'stage2_shape_guard_enabled', name = '第二轮谱形保护开关';
                case 'nonlinear_inner_solver', name = '第一轮线性内层求解器';
                case 'stage2_scan_inner_solver', name = '第二轮扫描内层求解器';
                case 'datacube_y_bin', name = 'datacube的Y方向空间合并';
                case 'datacube_x_bin', name = 'datacube的X方向空间合并';
                case 'datacube_energy_bin', name = 'datacube的能量轴合并';
                case 'datacube_parallel_enabled', name = 'datacube并行处理开关';
                case 'datacube_parallel_workers', name = '并行工作进程数量';
                case 'datacube_checkpoint_enabled', name = '断点保存开关';
                case 'datacube_resume_checkpoint', name = '从断点继续运行';
                case 'spatial_refit_enabled', name = '空间邻域二次拟合开关';
                case 'spatial_refit_iterations', name = '空间邻域二次拟合轮数';
                otherwise
                    name = lower(strrep(char(key), '_', ' '));
                    name = strrep(name, 'stage1', '第一轮');
                    name = strrep(name, 'stage2', '第二轮');
                    name = strrep(name, 'zlp', '零峰');
                    name = strrep(name, 'datacube', '数据立方');
                    name = strrep(name, 'spatial', '空间');
                    name = strrep(name, 'neighbor', '邻域');
                    name = strrep(name, 'region', '区域');
                    name = strrep(name, 'attention', '注意力');
                    name = strrep(name, 'interface', '界面');
                    name = strrep(name, 'bad', '坏点');
                    name = strrep(name, 'r2', 'R2拟合优度');
                    name = strrep(name, 'abs', '绝对');
                    name = strrep(name, 'local', '局部');
                    name = strrep(name, 'drop', '下降');
                    name = strrep(name, 'priority', '优先级');
                    name = strrep(name, 'bonus', '加分');
                    name = strrep(name, 'signal', '信号');
                    name = strrep(name, 'spectrum', '谱线');
                    name = strrep(name, 'temperature', '温度');
                    name = strrep(name, 'highfreq', '高频');
                    name = strrep(name, 'lowfreq', '低频');
                    name = strrep(name, 'wide', '宽尾');
                    name = strrep(name, 'overshoot', '过冲');
                    name = strrep(name, 'transition', '过渡区');
                    name = strrep(name, 'continuum', '连续谱');
                    name = strrep(name, 'reliability', '可靠性');
                    name = strrep(name, 'uncertainty', '不确定度');
                    name = strrep(name, 'component', '分量');
                    name = strrep(name, 'components', '分量');
                    name = strrep(name, 'width', '宽度');
                    name = strrep(name, 'bounds', '上下限');
                    name = strrep(name, 'center', '中心');
                    name = strrep(name, 'index', '序号');
                    name = strrep(name, 'taper', '过渡宽度');
                    name = strrep(name, 'limit', '限制');
                    name = strrep(name, 'fraction', '比例');
                    name = strrep(name, 'range', '范围');
                    name = strrep(name, 'guard', '保护');
                    name = strrep(name, 'window', '窗口');
                    name = strrep(name, 'weight', '权重');
                    name = strrep(name, 'threshold', '阈值');
                    name = strrep(name, 'enabled', '开关');
                    name = strrep(name, 'min', '最小');
                    name = strrep(name, 'max', '最大');
                    name = strrep(name, 'start', '起点');
                    name = strrep(name, 'end', '终点');
                    name = strrep(name, 'floor', '下限');
                    name = strrep(name, 'sigma', '尺度');
                    name = strrep(name, 'smooth', '平滑');
                    name = strrep(name, 'sparse', '稀疏');
                    name = strrep(name, 'dense', '加密');
                    name = strrep(name, 'step', '步长');
                    name = strrep(name, 'iterations', '迭代轮数');
                    name = strrep(name, 'iter', '迭代');
                    name = strrep(name, 'solver', '求解器');
                    name = strrep(name, 'inner', '内层');
                    name = strrep(name, 'scan', '扫描');
                    name = strrep(name, 'continuous', '连续优化');
                    name = strrep(name, 'coarse', '粗定位');
                    name = strrep(name, 'fine', '精修');
                    name = strrep(name, 'single', '单点');
                    name = strrep(name, 'fast', '快速');
                    name = strrep(name, 'parallel', '并行');
                    name = strrep(name, 'workers', '工作进程');
                    name = strrep(name, 'checkpoint', '断点');
                    name = strrep(name, 'resume', '续跑');
                    name = strrep(name, 'sample', '采样');
                    name = strrep(name, 'strategy', '策略');
                    name = strrep(name, 'pixels', '像素数');
                    name = strrep(name, 'energy', '能量');
                    name = strrep(name, 'crop', '截取');
                    name = strrep(name, 'bin', '合并');
                    name = strrep(name, 'shape', '谱形');
                    name = strrep(name, 'preserve', '保持');
                    name = strrep(name, 'fallback', '回退');
                    name = strrep(name, 'ratio', '比例温度');
                    name = strrep(name, 'prior', '先验');
                    name = strrep(name, 'noise', '噪声');
                    name = strrep(name, 'tail', '高频拖尾');
                    name = strrep(name, 'hard', '硬约束');
                    name = strrep(name, 'soft', '软约束');
                    name = strrep(name, 'core', '核心');
                    name = strrep(name, 'lambda', '正则');
                    name = strrep(name, 'l1', 'L1正则');
                    name = strrep(name, 'l2', 'L2正则');
                    name = strrep(name, 'mev', 'meV');
                    name = strrep(name, 'model', '模型');
                    name = strrep(name, 'type', '类型');
                    name = strrep(name, 'em ', '能量模 ');
                    name = strrep(name, 'pm ', '峰强 ');
                    name = strrep(name, ' be ', ' 细致平衡 ');
                    name = strrep(strtrim(name), '  ', ' ');
            end
        end

        function tip = settingAdvice(~, key)
            key = char(key);
            if contains(key, 'parallel')
                tip = '全量数据可打开并行；普通电脑建议2-4个worker，高配机器可逐步加到8-16。';
            elseif contains(key, 'checkpoint') || contains(key, 'resume')
                tip = '长时间全量运行建议开启断点保存；续跑时保持输出目录一致。';
            elseif contains(key, 'crop_range')
                tip = '按实验信息设定能量范围；不要强制固定，常用[-200 200]或[-250 250]。';
            elseif contains(key, 'energy_bin')
                tip = '速度慢或噪声大可设为2；需要看细峰形时设为1。';
            elseif contains(key, 'zlp_component_width') || contains(key, 'zlp_width')
                tip = '控制零峰宽度，太宽会吃低频光学/声学支，太窄会漏高频拖尾。';
            elseif contains(key, 'zlp_signal_guard') || contains(key, 'zlp_lowfreq_guard')
                tip = '保护真实信号不被零峰吞掉；低频被吃时提高保护权重或提前保护起点。';
            elseif contains(key, 'em_') || contains(key, 'dense_em')
                tip = '控制源强度/模态能量网格密度；峰形毛刺多时可适当加密，但会变慢。';
            elseif contains(key, 'lambda') || contains(key, 'smooth') || contains(key, 'sparse') || contains(key, 'l1') || contains(key, 'l2')
                tip = '正则项控制平滑和稀疏；毛刺大提高平滑，丢峰则降低稀疏惩罚。';
            elseif contains(key, 'attention')
                tip = '控制注意力选区；明显峰被漏掉时提高增益或降低阈值。';
            elseif contains(key, 'stage2_temperature_algorithm')
                tip = 'v2.1 使用 grid_scan；Data-valid 温度不做空间平滑，仅对无效点做 ZLP 界面引导的单向推断。';
            elseif contains(key, 'stage2_T_grid') || contains(key, 'fine_step') || contains(key, 'coarse_step')
                tip = '控制温度搜索精度和速度；步长越小越准但越慢。';
            elseif contains(key, 'shape_guard')
                tip = '建议开启，防止第二轮为了温度牺牲谱线拟合形状。';
            elseif contains(key, 'tail') || contains(key, 'overshoot')
                tip = '控制高频拖尾和过冲；高频抬尾时提高过冲约束或提前拖尾起点。';
            elseif contains(key, 'spatial_bad_R2') || contains(key, 'low_R2')
                tip = '用于识别局部R2明显偏低的坏点，并提高这些像素在空间二次拟合队列中的优先级。阈值太高会增加耗时，太低会漏掉坏点。';
            elseif contains(key, 'spatial_refit')
                tip = '温度图有孤立坏点时开启；界面样品要结合interface/peak shift避免抹平真实突变。';
            elseif contains(key, 'spatial_post') || contains(key, 'region_post')
                tip = '后处理温度图毛刺；权重太大会过度平滑，界面样品需谨慎。';
            elseif contains(key, 'reliability') || contains(key, 'uncertainty')
                tip = '控制可靠性/不确定度阈值；错得很自信时需调高惩罚或降低可靠性。';
            elseif contains(key, 'save')
                tip = '保存图和结果会增加运行时间；调试时建议开，批量摸参数可关。';
            else
                tip = '建议先保持默认；只有诊断图、R2/NRMSE或温度图明确指向该项时再调整。';
            end
        end
    end
end

