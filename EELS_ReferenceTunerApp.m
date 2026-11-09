classdef EELS_ReferenceTunerApp < handle
    properties
        UIFigure
        Container
        OwnsFigure=true
        ParameterTable
        InventoryTable
        GroupTable
        ModelChoice
        ModelFlags
        ConfigurationOverrides=struct()
        ParameterPreset='all'
        Stage1WeightAxes
        Stage2WeightAxes
        ConstraintAxes
        CostComponentAxes
        Temperature
        Budget
        Algorithm
        SwarmSize
        AnnealingTemperature
        FitRatioLimit
        SearchRadius
        DirectionTable
        TrainCount
        ValidationCount
        SampleSeed
        ParallelEnabled
        Workers
        JointAxes
        JointTab
        JointTable
        ValidationSpectra=[]
        SampleInfo=struct()
        Columns
        SpectrumAxes
        ValleyAxes
        LossAxes
        ZLPAxes
        ResidualAxes
        ProjectionAxes
        HeatmapAxes
        FeatureTable
        HistoryTable
        Status
        ActiveLabel
        RunButton
        ApplyButton
        ExportButton
        OutputEdit
        GroupEdit
        PlotTabs
        OverviewTab
        ReferenceViewTab
        Report=[]
        Energy=[]
        Spectrum=[]
        Settings=struct()
        PendingSettings=[]
        ApplyCallback=[]
        ResetCallback=[]
        ReferenceCallback=[]
        SettingsCallback=[]
        StopRequested=false
        Running=false
        Source=''
        LastExportDirectory=''
        LastError=''
        LastErrorReport=''
        ProgressHistory=[]
        ProgressBaseline=[]
        LastDraw=0
    end
    methods
        function app=EELS_ReferenceTunerApp(settings,result,apply_callback,visible,parent,reset_callback,reference_callback,settings_callback)
            if nargin<1,settings=struct();end
            if nargin<2,result=[];end
            if nargin<3,apply_callback=[];end
            if nargin<4,visible='on';end
            if nargin<5,parent=[];end
            if nargin<6,reset_callback=[];end
            if nargin<7,reference_callback=[];end
            if nargin<8,settings_callback=[];end
            app.Settings=EELS_PearsonSparse_FullFit.defaults(settings);app.ApplyCallback=apply_callback;
            app.ResetCallback=reset_callback;app.ReferenceCallback=reference_callback;app.SettingsCallback=settings_callback;
            if isempty(parent)
                app.UIFigure=uifigure('Name','2.2 | Reference spectrum tuning','Position',[60 60 1360 850],'Visible',visible);
                app.UIFigure.CloseRequestFcn=@(~,~)app.close();parent=app.UIFigure;
            else
                app.OwnsFigure=false;app.UIFigure=ancestor(parent,'figure');
            end
            app.Container=uigridlayout(parent,[1 2]);app.Container.ColumnWidth={340,'1x'};app.Container.Padding=[8 8 8 8];
            left=uigridlayout(app.Container,[14 2]);left.ColumnWidth={'1x','1x'};
            left.RowHeight={30,30,30,30,30,30,26,62,30,30,26,30,30,75};left.Scrollable='on';left.RowSpacing=5;
            b=uibutton(left,'Text','Import reference','ButtonPushedFcn',@(~,~)app.importSpectrum());b.Layout.Row=1;b.Layout.Column=1;
            b=uibutton(left,'Text','Use Data selection','ButtonPushedFcn',@(~,~)app.useCurrent());b.Layout.Row=1;b.Layout.Column=2;
            if isempty(reference_callback),b.Enable='off';end
            uilabel(left,'Text','Known temperature (K)');app.Temperature=uieditfield(left,'numeric','Value',300,'Limits',[1 Inf]);
            uilabel(left,'Text','Max evaluations');app.Budget=uispinner(left,'Value',500,'Limits',[1 2000]);
            uilabel(left,'Text','X columns (all = full map)');app.Columns=uieditfield(left,'text','Value','all');
            uilabel(left,'Text','Acquisition group ID');app.GroupEdit=uieditfield(left,'text','Placeholder','One ID per acquisition');
            b=uibutton(left,'Text','Import strategy','ButtonPushedFcn',@(~,~)app.importStrategy());b.Layout.Row=6;b.Layout.Column=1;
            b=uibutton(left,'Text','Restore defaults','ButtonPushedFcn',@(~,~)app.reset());b.Layout.Row=6;b.Layout.Column=2;
            app.ActiveLabel=uilabel(left,'Text','Current/default fitting settings active','FontColor',[.1 .3 .55]);app.ActiveLabel.Layout.Row=7;app.ActiveLabel.Layout.Column=[1 2];
            app.ParameterTable=uitable(left,'ColumnName',{'Enabled','Parameter','Initial value','Lower bound','Upper bound'}, ...
                'ColumnEditable',[true false true true true],'ColumnWidth',{40,290,85,85,85},'RowName',{});
            app.ParameterTable.Layout.Row=8;app.ParameterTable.Layout.Column=[1 2];app.refreshParameters();
            app.RunButton=uibutton(left,'Text','Start tuning','ButtonPushedFcn',@(~,~)app.run());app.RunButton.Layout.Row=9;app.RunButton.Layout.Column=1;
            b=uibutton(left,'Text','Stop and keep best','ButtonPushedFcn',@(~,~)app.stop());b.Layout.Row=9;b.Layout.Column=2;
            app.ApplyButton=uibutton(left,'Text','Apply to Stage1/2','Enable','off','ButtonPushedFcn',@(~,~)app.apply());app.ApplyButton.Layout.Row=10;app.ApplyButton.Layout.Column=1;
            app.ExportButton=uibutton(left,'Text','Save dataset as','Enable','off','ButtonPushedFcn',@(~,~)app.export());app.ExportButton.Layout.Row=10;app.ExportButton.Layout.Column=2;
            lab=uilabel(left,'Text','Auto-save every run, including failed trials');lab.Layout.Row=11;lab.Layout.Column=[1 2];
            app.OutputEdit=uieditfield(left,'text','Value',fullfile(fileparts(mfilename('fullpath')),'tuning_runs'));app.OutputEdit.Layout.Row=12;app.OutputEdit.Layout.Column=[1 2];
            b=uibutton(left,'Text','Sampling / heldout validation','ButtonPushedFcn',@(~,~)app.showJoint());b.Layout.Row=13;b.Layout.Column=[1 2];
            app.Status=uitextarea(left,'Editable','off','Value',{'Known temperature is used to score calibration candidates.','Run/import does not apply a strategy automatically. Validate before applying.'});
            app.Status.Layout.Row=14;app.Status.Layout.Column=[1 2];
            app.PlotTabs=uitabgroup(app.Container);
            overview=uitab(app.PlotTabs,'Title','Tuning progress');app.OverviewTab=overview;g=uigridlayout(overview,[2 2]);
            app.LossAxes=uiaxes(g);app.ZLPAxes=uiaxes(g);app.ValleyAxes=uiaxes(g);app.ResidualAxes=uiaxes(g);
            space=uitab(app.PlotTabs,'Title','Parameter space');g=uigridlayout(space,[2 1]);g.RowHeight={'1x','1.4x'};
            app.ProjectionAxes=uiaxes(g);app.HeatmapAxes=uiaxes(g);
            input=uitab(app.PlotTabs,'Title','Reference and features');app.ReferenceViewTab=input;g=uigridlayout(input,[2 1]);
            app.SpectrumAxes=uiaxes(g);app.FeatureTable=uitable(g,'ColumnName',{'Pre-fit feature','Value','Valid'},'ColumnWidth',{290,130,55},'RowName',{});
            trials=uitab(app.PlotTabs,'Title','Candidate history');g=uigridlayout(trials,[1 1]);
            app.HistoryTable=uitable(g,'ColumnName',{'Evaluation','Loss','Valley K','Error K','NRMSE','Time s','Status','Over-limit points','Failed points','Max residual ratio','Search objective'},'RowName',{});
            joint=uitab(app.PlotTabs,'Title','Training and heldout');app.JointTab=joint;g=uigridlayout(joint,[3 1]);g.RowHeight={110,'1x','1x'};
            samplePanel=uigridlayout(g,[3 4]);samplePanel.ColumnWidth={220,'1x',110,'1x'};samplePanel.RowHeight={28,28,30};
            samplePanel.Padding=[0 0 0 0];samplePanel.RowSpacing=6;
            b=uibutton(samplePanel,'Text','Sample training / heldout points from Data','ButtonPushedFcn',@(~,~)app.sampleCurrent());b.Layout.Row=3;b.Layout.Column=[1 4];
            if isempty(reference_callback),b.Enable='off';end
            lab=uilabel(samplePanel,'Text','Training points (same known T)');lab.Layout.Row=1;lab.Layout.Column=1;
            app.TrainCount=uispinner(samplePanel,'Value',32,'Limits',[2 100],'RoundFractionalValues','on');app.TrainCount.Layout.Row=1;app.TrainCount.Layout.Column=2;
            lab=uilabel(samplePanel,'Text','Heldout points');lab.Layout.Row=1;lab.Layout.Column=3;
            app.ValidationCount=uispinner(samplePanel,'Value',8,'Limits',[0 100],'RoundFractionalValues','on');app.ValidationCount.Layout.Row=1;app.ValidationCount.Layout.Column=4;
            lab=uilabel(samplePanel,'Text','Random seed');lab.Layout.Row=2;lab.Layout.Column=1;
            app.SampleSeed=uispinner(samplePanel,'Value',2206,'Limits',[0 2^32-1],'RoundFractionalValues','on');app.SampleSeed.Layout.Row=2;app.SampleSeed.Layout.Column=2;
            app.ParallelEnabled=uicheckbox(samplePanel,'Text','C++ parallel','Value',true);app.ParallelEnabled.Layout.Row=2;app.ParallelEnabled.Layout.Column=3;
            app.Workers=uispinner(samplePanel,'Value',EELS_Parallel.default_workers(),'Limits',[1 128],'RoundFractionalValues','on','Tooltip','Process workers; shared with Stage1/2');app.Workers.Layout.Row=2;app.Workers.Layout.Column=4;


            app.JointAxes=uiaxes(g);app.JointTable=uitable(g,'ColumnName',{'Role','Y','X','Baseline valley K','Best valley K','Absolute error K'},'RowName',{});
            weightsTab=uitab(app.PlotTabs,'Title','Fit weights');wg=uigridlayout(weightsTab,[2 1]);
            app.Stage1WeightAxes=uiaxes(wg);app.Stage2WeightAxes=uiaxes(wg);
            gateTab=uitab(app.PlotTabs,'Title','Residual gate');gg=uigridlayout(gateTab,[1 1]);app.ConstraintAxes=uiaxes(gg);
            costTab=uitab(app.PlotTabs,'Title','Cost components');cg=uigridlayout(costTab,[1 1]);app.CostComponentAxes=uiaxes(cg);
            paramTab=uitab(app.PlotTabs,'Title','Parameters');g=uigridlayout(paramTab,[2 1]);g.RowHeight={68,'1x'};
            presets=uigridlayout(g,[2 3]);presets.Padding=[0 0 0 0];
            uibutton(presets,'Text','Weights (10)','ButtonPushedFcn',@(~,~)app.setParameterPreset('weights'));
            uibutton(presets,'Text',sprintf('All available (up to %d)',numel(EELS_ReferenceTuner.parameters())),'ButtonPushedFcn',@(~,~)app.setParameterPreset('all'));
            uibutton(presets,'Text','Original (22)','ButtonPushedFcn',@(~,~)app.setParameterPreset('legacy'));
            uibutton(presets,'Text','Weights + Stage2 (22)','ButtonPushedFcn',@(~,~)app.setParameterPreset('joint'));
            uibutton(presets,'Text','Stage2 ZLP (12)','ButtonPushedFcn',@(~,~)app.setParameterPreset('stage2'));
            app.ParameterTable.Parent=g;app.ParameterTable.Layout.Row=2;app.ParameterTable.Layout.Column=1;
            app.ParameterTable.CellEditCallback=@(~,~)app.markCustomParameters();
            directionTab=uitab(app.PlotTabs,'Title','Direction probes');dg=uigridlayout(directionTab,[1 1]);
            app.DirectionTable=uitable(dg,'ColumnName',{'Parameter','Direction','Search improvement','Mean error improvement K'},'RowName',{});
            left.RowHeight{8}=152;
            searchPanel=uigridlayout(left,[5 2]);searchPanel.Layout.Row=8;searchPanel.Layout.Column=[1 2];
            searchPanel.Padding=[0 0 0 0];searchPanel.RowHeight={26,26,26,26,26};searchPanel.RowSpacing=4;
            uilabel(searchPanel,'Text','Search algorithm');
            app.Algorithm=uidropdown(searchPanel,'Items',{'Grouped + joint SA','Grouped + joint PSO','Surrogate + PSO + polish','All-selected PSO','Local GP'},'ItemsData',{'grouped_sa','grouped_pso','hybrid','pso','gp'},'Value','grouped_sa', ...
                'ValueChangedFcn',@(~,~)app.refreshSearchControls());
            app.Algorithm.Tooltip='SA uses the group order and budget of grouped PSO. It can accept uphill steps; only the best strategy passing the residual gate can be applied.';
            uilabel(searchPanel,'Text','PSO swarm size');
            app.SwarmSize=uispinner(searchPanel,'Value',12,'Limits',[2 200]);
            uilabel(searchPanel,'Text','SA temperature (unitless)');
            app.AnnealingTemperature=uieditfield(searchPanel,'numeric','Value',.1,'Limits',[.001 2]);
            app.AnnealingTemperature.Tooltip='Controls acceptance of worse trials, not sample temperature K. Default 0.1; the search radius controls step size.';
            uilabel(searchPanel,'Text','Max residual / baseline');app.FitRatioLimit=uieditfield(searchPanel,'numeric','Value',2,'Limits',[1 5]);
            uilabel(searchPanel,'Text','Initial search radius');app.SearchRadius=uieditfield(searchPanel,'numeric','Value',.12,'Limits',[.025 .5]);
            app.refreshSearchControls();
            invTab=uitab(app.PlotTabs,'Title','Full settings catalog');ig=uigridlayout(invTab,[2 1]);ig.RowHeight={30,'1x'};
            uilabel(ig,'Text','All settings are listed. Inactive and unregistered settings are excluded from search.');
            app.InventoryTable=uitable(ig,'Data',EELS_ParameterCatalog.inventory(app.Settings),'ColumnEditable',false);
            groupTab=uitab(app.PlotTabs,'Title','Group results');gg=uigridlayout(groupTab,[1 1]);app.GroupTable=uitable(gg,'ColumnName',{'Group','Evaluations','Before loss','After loss','Algorithm'},'RowName',{});
            modelTab=uitab(app.PlotTabs,'Title','Model configuration');mg=uigridlayout(modelTab,[3 1]);mg.RowHeight={30,30,'1x'};
            app.ModelChoice=uidropdown(mg,'Items',{'asym_pearson','asym_pvoigt'},'Value',app.Settings.zlp_model_type, ...
                'ValueChangedFcn',@(~,~)app.setModelConfiguration());
            uilabel(mg,'Text','Model choices and flags stay fixed during each run. Changes update available search parameters.');
            keys={'zlp_lowfreq_guard_enabled','zlp_signal_guard_enabled','stage1_zlp_core_fartail_enabled', ...
                'stage1_zlp_gain_tail_prior_enabled','stage1_zlp_gain_asym_shape_prior_enabled','stage1_zlp_gain_heavy_init_enabled', ...
                'em_lowfreq_dense_enabled','stage2_attention_enabled','stage2_attention_hardbe_only_enabled', ...
                'stage2_soft_dense_em_enabled','stage2_far_tail_downweight_enabled','stage2_far_tail_loss_budget_downweight_enabled', ...
                'stage2_zlp_refine_enabled','stage2_zlp_refine_scale_enabled','stage2_zlp_refine_tilt_enabled', ...
                'stage2_zlp_refine_shift_enabled','stage2_zlp_refine_width_enabled','stage2_zlp_refine_tail_enabled','stage2_zlp_refine_tail_asym_enabled'};
            values=cellfun(@(key)logical(app.Settings.(key)),keys);
            app.ModelFlags=uitable(mg,'Data',[keys',num2cell(values')],'ColumnName',{'Feature','Enabled'},'ColumnEditable',[false true], ...
                'ColumnWidth',{390,80},'CellEditCallback',@(~,~)app.setFeatureConfiguration());
            if isstruct(result)&&isfield(result,'ene_raw')&&isfield(result,'raw_data')
                app.setReference(result.ene_raw,result.raw_data,'Selected spectrum');
            end
        end
        function refreshSearchControls(app)
            app.SwarmSize.Enable='off';app.AnnealingTemperature.Enable='off';
            if ~app.Running
                if strcmp(app.Algorithm.Value,'grouped_sa'),app.AnnealingTemperature.Enable='on';end
                if ismember(app.Algorithm.Value,{'pso','grouped_pso','hybrid'}),app.SwarmSize.Enable='on';end
            end
        end
        function refreshParameters(app)
            p=EELS_ReferenceTuner.parameters();rows=cell(numel(p),5);
            old=app.ParameterTable.Data;
            if ~isempty(app.ModelChoice),app.ModelChoice.Value=app.Settings.zlp_model_type;end
            if ~isempty(app.ModelFlags)
                flags=app.ModelFlags.Data;for j=1:size(flags,1),flags{j,2}=logical(app.Settings.(flags{j,1}));end;app.ModelFlags.Data=flags;
            end
            for k=1:numel(p)
                v=app.Settings.(p(k).name);rows(k,:)={p(k).enabled,p(k).name,v,min(v,p(k).lower),max(v,p(k).upper)};
                if strcmp(app.ParameterPreset,'all'),rows{k,1}=true;
                elseif strcmp(app.ParameterPreset,'legacy'),rows{k,1}=k<=22;
                elseif strcmp(app.ParameterPreset,'weights'),rows{k,1}=k>=23&&k<=32;
                elseif strcmp(app.ParameterPreset,'stage2'),rows{k,1}=k>=33&&k<=44;
                elseif strcmp(app.ParameterPreset,'joint'),rows{k,1}=k>=23&&k<=44;
                elseif strcmp(app.ParameterPreset,'custom')&&size(old,1)==numel(p),rows{k,1}=old{k,1};end
            end
            available=EELS_ParameterCatalog.available(p,app.Settings);
            for j=1:numel(p),rows{j,1}=rows{j,1}&&available(j);end
            app.ParameterTable.Data=rows;
            if ~isempty(app.InventoryTable),app.InventoryTable.Data=EELS_ParameterCatalog.inventory(app.Settings);end
        end
        function setModelConfiguration(app)
            if app.Running,app.ModelChoice.Value=app.Settings.zlp_model_type;return;end
            app.Settings.zlp_model_type=app.ModelChoice.Value;app.ConfigurationOverrides.zlp_model_type=app.ModelChoice.Value;
        end
        function setFeatureConfiguration(app)
            if app.Running,return;end
            rows=app.ModelFlags.Data;
            for j=1:size(rows,1)
                app.Settings.(rows{j,1})=logical(rows{j,2});app.ConfigurationOverrides.(rows{j,1})=logical(rows{j,2});
            end
            app.refreshParameters();
        end
        function markCustomParameters(app),app.ParameterPreset='custom';end
        function setParameterPreset(app,mode)
            if app.Running,return;end
            app.ParameterPreset=mode;p=EELS_ReferenceTuner.parameters();rows=app.ParameterTable.Data;
            for k=1:numel(p)
                if strcmp(mode,'all'),rows{k,1}=true;
                elseif strcmp(mode,'legacy'),rows{k,1}=k<=22;
                elseif strcmp(mode,'weights'),rows{k,1}=k>=23&&k<=32;
                elseif strcmp(mode,'stage2'),rows{k,1}=k>=33&&k<=44;
                elseif strcmp(mode,'joint'),rows{k,1}=k>=23&&k<=44;
                else,rows{k,1}=p(k).enabled;end
            end
            available=EELS_ParameterCatalog.available(p,app.Settings);
            for j=1:numel(p),rows{j,1}=rows{j,1}&&available(j);end
            app.ParameterTable.Data=rows;
            if ~isempty(app.InventoryTable),app.InventoryTable.Data=EELS_ParameterCatalog.inventory(app.Settings);end
        end
        function setReference(app,e,y,source)
            if app.Running,error('EELS:TunerBusy','Wait for the current evaluation to stop.');end
            e=double(e(:));y=double(y);if isvector(y),y=y(:);end
            if size(y,1)~=numel(e)||any(~isfinite(y(:))),error('EELS:Reference','Energy rows must match finite spectrum columns.');end
            f=EELS_TuningDataset.features(e,y(:,1),app.Settings);
            app.ValidationSpectra=[];app.SampleInfo=struct();cla(app.JointAxes);app.JointTable.Data={};
            app.Energy=e;app.Spectrum=y;app.Source=source;app.Report=[];app.PendingSettings=[];
            app.ApplyButton.Enable='off';app.ExportButton.Enable='off';app.ProgressHistory=[];
            for ax=[app.LossAxes,app.ZLPAxes,app.ValleyAxes,app.ResidualAxes,app.ProjectionAxes,app.HeatmapAxes],cla(ax);end
            app.HistoryTable.Data={};app.DirectionTable.Data={};app.GroupTable.Data={};
            cla(app.Stage1WeightAxes);cla(app.Stage2WeightAxes);cla(app.ConstraintAxes);cla(app.CostComponentAxes);
            plot(app.SpectrumAxes,e,y);xlabel(app.SpectrumAxes,'Energy (meV)');title(app.SpectrumAxes,source,'Interpreter','none');
            app.FeatureTable.ColumnName={'Pre-fit features of reference 1','Value','Valid'};
            app.FeatureTable.Data=[f.names(:),num2cell(f.values(:)),num2cell(f.valid(:))];
            app.PlotTabs.SelectedTab=app.ReferenceViewTab;
            app.Status.Value={['Reference: ',source],'Reference updated. The strategy already applied to Stage1/2 is unchanged.'};
        end
        function useCurrent(app)
            if app.Running||isempty(app.ReferenceCallback),return;end
            try
                r=app.ReferenceCallback();
                if isfield(r,'cube'),app.setCubeReference(r.ene,r.cube,r.group);return;end
                app.setReference(r.ene,r.spectrum,r.source);
                if isfield(r,'group'),app.GroupEdit.Value=r.group;end
            catch ME,app.Status.Value={ME.message};end
        end
        function showJoint(app),app.PlotTabs.SelectedTab=app.JointTab;end
        function sampleCurrent(app)
            if app.Running||isempty(app.ReferenceCallback),return;end
            try
                r=app.ReferenceCallback();
                if ~isfield(r,'cube'),error('Import a full map into Data before sampling multiple points.');end
                app.setCubeReference(r.ene,r.cube,r.group);
            catch ME,app.LastError=ME.message;app.Status.Value={ME.message};end
        end
        function setCubeReference(app,e,cube,source)
            if app.Running,error('EELS:TunerBusy','Wait until tuning stops.');end
            [Y,V,info]=EELS_ReferenceTuner.sampleCube(cube,app.Columns.Value,app.TrainCount.Value,app.ValidationCount.Value,app.SampleSeed.Value);
            app.setReference(e,Y,sprintf('%s | %d train + %d heldout',source,size(Y,2),size(V,2)));
            app.ValidationSpectra=V;app.SampleInfo=info;app.GroupEdit.Value=source;app.LastError='';
            app.updateJointView();app.showJoint();
            app.Status.Value={sprintf('Fixed sample: %d training and %d heldout points.',size(Y,2),size(V,2)), ...
                'All points must share the known temperature. Sampling uses available points without replacement.', ...
                'C++ workers process spectra; heldout points are evaluated only after search.'};
        end
        function updateJointView(app)
            if ~isfield(app.SampleInfo,'train_yx'),return;end
            tr=app.SampleInfo.train_yx;va=app.SampleInfo.validation_yx;n=size(tr,1);m=size(va,1);
            before=NaN(n+m,1);after=before;
            if ~isempty(app.Report)
                before(1:n)=[app.Report.baseline.point_outcomes.valley_T_K];after(1:n)=[app.Report.best.point_outcomes.valley_T_K];
                if isfield(app.Report.validation,'before_K')
                    before(n+1:end)=app.Report.validation.before_K;after(n+1:end)=app.Report.validation.after_K;
                end
            end
            app.JointTable.Data=[repmat({'Train'},n,1),num2cell(tr);repmat({'Heldout'},m,1),num2cell(va)];
            app.JointTable.Data=[app.JointTable.Data,num2cell(before),num2cell(after),num2cell(abs(after-app.Temperature.Value))];
            cla(app.JointAxes);
            if isempty(app.Report)
                scatter(app.JointAxes,tr(:,2),tr(:,1),45,'b','filled');hold(app.JointAxes,'on');
                scatter(app.JointAxes,va(:,2),va(:,1),60,'r','x');hold(app.JointAxes,'off');
                title(app.JointAxes,'Fixed random sample positions');xlabel(app.JointAxes,'X');ylabel(app.JointAxes,'Y');
                legend(app.JointAxes,{'Train','Heldout'});app.JointAxes.YDir='reverse';
            else
                plot(app.JointAxes,1:n+m,before,'bo--',1:n+m,after,'ro-');hold(app.JointAxes,'on');
                yline(app.JointAxes,app.Temperature.Value,'k:');xline(app.JointAxes,n+.5,'k--');hold(app.JointAxes,'off');
                app.JointAxes.YDir='normal';title(app.JointAxes,'Per-point temperature: train | heldout');
                xlabel(app.JointAxes,'Point index');ylabel(app.JointAxes,'Raw valley (K)');legend(app.JointAxes,{'Before','Best','Known T','Holdout boundary'});
            end
        end
        function importSpectrum(app)
            if app.Running,return;end
            [name,path]=uigetfile({'*.mat;*.csv','Reference spectra MAT / CSV'});if isequal(name,0),return;end
            try,app.loadReferenceFile(fullfile(path,name));catch ME,app.LastError=ME.message;app.Status.Value={ME.message};end
        end
        function loadReferenceFile(app,file)
            if app.Running,error('EELS:TunerBusy','Wait until tuning stops.');end
            if endsWith(lower(file),'.csv')
                a=readmatrix(file);e=a(:,1);y=a(:,2:end);
            else
                a=load(file);
                if isfield(a,'data')&&(isobject(a.data)||isstruct(a.data))
                    [cube,e]=EELS_PearsonSparse_FullFit.load_datacube_mat(file);
                    app.setCubeReference(e,cube,file);return;
                elseif isfield(a,'ene')&&(isfield(a,'spectrum')||isfield(a,'reference'))
                    e=a.ene;if isfield(a,'spectrum'),y=a.spectrum;else,y=a.reference;end
                else
                    [cube,e]=EELS_PearsonSparse_FullFit.load_datacube_mat(file);
                    app.setCubeReference(e,cube,file);return;
                end
            end
            if ndims(y)==3
                app.setCubeReference(e,y,file);return;
            elseif ~isvector(y)&&size(y,1)==numel(e)
                cube=reshape(y',[1,size(y,2),size(y,1)]);app.setCubeReference(e,cube,file);return;
            end
            app.setReference(e,y,file);app.GroupEdit.Value=file;app.LastError='';
            app.Status.Value={'File contains one spectrum: using single-spectrum mode.', ...
                'Import a reference map or multiple spectra for joint tuning.'};
        end
        function pointProgress(app,done,total,label)
            app.Status.Value={sprintf('%s | %d/%d spectra completed',label,done,total), ...
                'Running Stage1/2 kernels; Stop cancels unfinished parallel spectra.'};
        end
        function run(app)
            if app.Running,return;end
            if isempty(app.Spectrum),app.Status.Value={'Import a reference spectrum first.'};return;end
            app.Running=true;app.StopRequested=false;app.LastError='';app.RunButton.Enable='off';app.ApplyButton.Enable='off';app.ExportButton.Enable='off';
            app.ParameterTable.Enable='off';app.Temperature.Enable='off';app.Budget.Enable='off';app.OutputEdit.Enable='off';app.GroupEdit.Enable='off';app.Algorithm.Enable='off';app.SwarmSize.Enable='off';
            app.FitRatioLimit.Enable='off';app.SearchRadius.Enable='off';
            app.ModelChoice.Enable='off';app.ModelFlags.Enable='off';
            app.AnnealingTemperature.Enable='off';
            app.ParallelEnabled.Enable='off';app.Workers.Enable='off';
            app.Report=[];app.PendingSettings=[];app.ProgressHistory=[];app.ProgressBaseline=[];app.LastDraw=0;
            app.LastExportDirectory='';cleanup=onCleanup(@()app.finish()); %#ok<NASGU>
            app.PlotTabs.SelectedTab=app.OverviewTab;
            try
                rows=app.ParameterTable.Data;p=struct('name',{},'lower',{},'upper',{});s=app.Settings;
                if ~isempty(app.SettingsCallback),s=app.SettingsCallback();end
                fields=fieldnames(app.ConfigurationOverrides);for j=1:numel(fields),s.(fields{j})=app.ConfigurationOverrides.(fields{j});end
                s.zlp_model_type=app.ModelChoice.Value;
                for k=1:size(rows,1)
                    s.(rows{k,2})=rows{k,3};
                    if rows{k,1}
                        p(end+1)=struct('name',rows{k,2},'lower',rows{k,4},'upper',rows{k,5}); %#ok<AGROW>
                    end
                end
                directory=EELS_TuningDataset.newDirectory(app.OutputEdit.Value);app.LastExportDirectory=directory;
                request=struct('energy_meV',app.Energy,'spectrum',app.Spectrum,'known_T_K',app.Temperature.Value, ...
                    'settings',s,'parameters',p,'source',app.Source,'acquisition_group',app.GroupEdit.Value,'budget',app.Budget.Value,'algorithm',app.Algorithm.Value,'swarm_size',app.SwarmSize.Value, ...
                    'max_fit_ratio',app.FitRatioLimit.Value,'search_radius',app.SearchRadius.Value,'annealing_temperature',app.AnnealingTemperature.Value, ...
                    'validation_spectra',app.ValidationSpectra,'sample_info',app.SampleInfo, ...
                    'parallel_enabled',app.ParallelEnabled.Value,'parallel_workers',app.Workers.Value); %#ok<NASGU>
                save(fullfile(directory,'request.mat'),'request');
                EELS_TuningDataset.writeJson(fullfile(directory,'status.json'),struct('status','running'));
                o=struct('parameters',p,'max_evaluations',app.Budget.Value, ...
                    'max_fit_ratio',app.FitRatioLimit.Value,'search_radius',app.SearchRadius.Value,'annealing_temperature',app.AnnealingTemperature.Value, ...
                    'algorithm',app.Algorithm.Value,'swarm_size',app.SwarmSize.Value, ...
                    'progress_callback',@(row,best)app.progress(row,best),'stop_callback',@()app.StopRequested, ...
                    'validation_spectra',app.ValidationSpectra,'sample_info',app.SampleInfo, ...
                    'parallel_enabled',app.ParallelEnabled.Value,'parallel_workers',app.Workers.Value, ...
                    'point_progress_callback',@(done,total,label)app.pointProgress(done,total,label));
                app.Status.Value={'Preparing tuning; starting the process pool can take a moment.'};drawnow;
                app.Report=EELS_ReferenceTuner.tune(app.Energy,app.Spectrum,app.Temperature.Value,s,o);
                app.Report.source=app.Source;app.Report.acquisition_group=app.GroupEdit.Value;
                app.PendingSettings=app.Report.best.settings;app.ApplyButton.Enable='on';app.ExportButton.Enable='on';
                app.renderReport();
                app.saveRun(directory);
                app.Status.Value={sprintf('Reference valley %.2f -> %.2f K; known %.2f K. Not yet applied.', ...
                    app.Report.baseline.valley.T_K,app.Report.best.valley.T_K,app.Report.known_T_K),['Dataset: ',directory]};
                if size(app.Spectrum,2)>1
                    app.Status.Value={sprintf('Training MAE (%d points) %.2f -> %.2f K. Not yet applied.',size(app.Spectrum,2), ...
                        mean([app.Report.baseline.point_outcomes.error_K]),mean([app.Report.best.point_outcomes.error_K])),['Dataset: ',directory]};
                    if isfield(app.Report.validation,'after_MAE_K')
                        app.Status.Value=[app.Status.Value;{sprintf('Heldout MAE %.2f -> %.2f K (excluded from selection)',app.Report.validation.before_MAE_K,app.Report.validation.after_MAE_K)}];
                    end
                end
                [accepted,regressed,failed]=EELS_TuningObjective.masks(app.Report.history);
                app.Status.Value=[app.Status.Value;{sprintf('%d candidates: %d eligible, %d over residual limit, %d with failures. See Residual gate.', ...
                    numel(accepted),nnz(accepted),nnz(regressed),nnz(failed))}];
            catch ME
                app.LastError=ME.message;app.LastErrorReport=getReport(ME,'extended','hyperlinks','off');app.Status.Value={['Run/export failed: ',ME.message]};
                if ~isempty(app.LastExportDirectory)
                    try
                        EELS_TuningDataset.writeJson(fullfile(app.LastExportDirectory,'status.json'), ...
                            struct('status','failed','identifier',ME.identifier,'message',ME.message));
                    catch
                    end
                end
            end
        end
        function progress(app,row,best)
            if ~isgraphics(app.Container),app.StopRequested=true;return;end
            if isempty(app.ProgressHistory),app.ProgressHistory=row;app.ProgressBaseline=best;
            else,app.ProgressHistory(end+1)=row;end
            app.Status.Value={sprintf('Evaluation %d: loss %.4g, valley %.2f K; best %.2f K',row.evaluation,row.score,row.valley_T_K,best.valley.T_K),row.status};
            if isfield(row,'regressed_count')
                phase=row.proposal;
                if strcmp(phase,'direction_probe'),phase='Bidirectional probes';
                elseif strcmp(phase,'measured_direction_combination'),phase='Combined improving directions';
                elseif strcmp(phase,'local_particleswarm'),phase='Local PSO';end
                app.Status.Value{2}=sprintf('%s: %d/%d points over residual limit, %d failed; search objective %.4g', ...
                    phase,row.regressed_count,numel(row.point_outcomes),row.failed_count,row.search_score);
            end
            if (now-app.LastDraw)*86400>.8
                p=app.selectedParameters();
                EELS_TuningPlots.render(app.axes(),app.ProgressHistory,app.ProgressBaseline,best,app.Temperature.Value,p);
                EELS_TuningPlots.weightProfiles(app.Stage1WeightAxes,app.Stage2WeightAxes,app.ProgressBaseline,best);
                EELS_TuningPlots.constraints(app.ConstraintAxes,app.ProgressHistory,app.FitRatioLimit.Value);
                app.updateHistory(app.ProgressHistory);app.LastDraw=now;drawnow limitrate;
            end
        end
        function p=selectedParameters(app)
            rows=app.ParameterTable.Data;p=struct('name',{},'lower',{},'upper',{});
            for k=1:size(rows,1)
                if rows{k,1},p(end+1)=struct('name',rows{k,2},'lower',rows{k,4},'upper',rows{k,5});end %#ok<AGROW>
            end
        end
        function ax=axes(app)
            ax=struct('loss',app.LossAxes,'zlp',app.ZLPAxes,'valley',app.ValleyAxes,'residual',app.ResidualAxes, ...
                'projection',app.ProjectionAxes,'parameters',app.HeatmapAxes);
        end
        function renderReport(app)
            r=app.Report;EELS_TuningPlots.render(app.axes(),r.history,r.baseline,r.best,r.known_T_K,r.parameters);app.updateHistory(r.history);
            EELS_TuningPlots.weightProfiles(app.Stage1WeightAxes,app.Stage2WeightAxes,r.baseline,r.best);
            EELS_TuningPlots.constraints(app.ConstraintAxes,r.history,r.objective_settings.max_fit_ratio);
            EELS_Stage2Tuning.plot(app.CostComponentAxes,r.best.result);
            app.updateJointView();
            app.GroupTable.Data={};
            if isfield(r.optimizer,'stages')
                a=r.optimizer.stages;app.GroupTable.Data=[{a.group}',num2cell([a.evaluations]'),num2cell([a.before]'),num2cell([a.after]')];
                if isfield(a,'method'),app.GroupTable.Data(:,5)={a.method}';else,app.GroupTable.Data(:,5)=repmat({'pso'},numel(a),1);end
            end
        end
        function updateHistory(app,h)
            rows=cell(numel(h),11);
            for k=1:numel(h)
                row=h(k);states={row.point_outcomes.status};
                regressed=nnz(strcmp(states,'fit_regression'));failed=nnz(~ismember(states,{'ok','fit_regression'}));
                ratio=max([row.point_outcomes.NRMSE]./max([h(1).point_outcomes.NRMSE],eps));search=NaN;
                if isfield(row,'search_score'),search=row.search_score;end
                rows(k,:)={row.evaluation,row.score,row.valley_T_K,row.error_K,row.NRMSE,row.seconds,row.status,regressed,failed,ratio,search};
            end
            app.HistoryTable.Data=rows;
            p=app.selectedParameters();if ~isempty(app.Report),p=app.Report.parameters;end
            directions={};
            for k=2:numel(h)
                if ~strcmp(h(k).proposal,'direction_probe'),continue;end
                j=find(abs(h(k).values-h(1).values)>1e-12);if numel(j)~=1,continue;end
                label='Increase';if h(k).values(j)<h(1).values(j),label='Decrease';end
                directions(end+1,:)={p(j).name,label,h(1).search_score-h(k).search_score,h(1).error_K-h(k).error_K}; %#ok<AGROW>
            end
            app.DirectionTable.Data=directions;
        end
        function saveRun(app,directory)
            EELS_TuningDataset.export(app.Report,directory);EELS_TuningPlots.export(app.Report,directory);
            state='complete';if app.Report.cancelled,state='cancelled_with_best';end
            EELS_TuningDataset.writeJson(fullfile(directory,'status.json'),struct('status',state,'applied_automatically',false));
        end
        function stop(app),app.StopRequested=true;end
        function finish(app)
            app.Running=false;
            if isgraphics(app.Container)
                app.ParallelEnabled.Enable='on';app.Workers.Enable='on';
                app.RunButton.Enable='on';app.ParameterTable.Enable='on';app.Temperature.Enable='on';
                app.Budget.Enable='on';app.OutputEdit.Enable='on';app.GroupEdit.Enable='on';app.Algorithm.Enable='on';app.SwarmSize.Enable='on';
                app.FitRatioLimit.Enable='on';app.SearchRadius.Enable='on';
                app.ModelChoice.Enable='on';app.ModelFlags.Enable='on';
                app.refreshSearchControls();
            end
        end
        function apply(app)
            if app.Running||isempty(app.PendingSettings),return;end
            s=app.PendingSettings;
            if ~isempty(app.ApplyCallback),app.ApplyCallback(s);end
            app.Settings=s;app.ConfigurationOverrides=struct();app.refreshParameters();
            app.ActiveLabel.Text='Tuned strategy active for Stage1/2';
            app.Status.Value={'Initialization rules and penalties applied. Each spectrum still optimizes its own shape.','Reference theta is exported as a descriptor, not a fixed shape or temperature for new spectra.'};
        end
        function reset(app)
            if app.Running,return;end
            if ~isempty(app.ResetCallback),app.ResetCallback();end
            app.Settings=EELS_PearsonSparse_FullFit.defaults(struct());app.ConfigurationOverrides=struct();app.refreshParameters();
            app.ActiveLabel.Text='Default settings restored';app.Status.Value={'Stage1/2 now use default initialization and penalties.'};
        end
        function importStrategy(app)
            if app.Running,return;end
            [name,path]=uigetfile('*.mat','Import fitting strategy or tuning report');if isequal(name,0),return;end
            try,app.loadStrategy(fullfile(path,name));catch ME,app.Status.Value={ME.message};end
        end
        function loadStrategy(app,file)
            if app.Running,error('EELS:TunerBusy','Cannot import while tuning.');end
            a=load(file);s=[];
            if isfield(a,'strategy'),s=a.strategy.settings;
            elseif isfield(a,'tuned_settings'),s=a.tuned_settings;
            elseif isfield(a,'report'),s=a.report.best.settings;end
            if ~isstruct(s)||isempty(fieldnames(s)),error('EELS:Strategy','File contains no fit strategy.');end
            s=EELS_PearsonSparse_FullFit.defaults(s);catalog=EELS_ReferenceTuner.parameters();
            for k=1:numel(catalog)
                v=s.(catalog(k).name);
                if ~isnumeric(v)||~isscalar(v)||~isfinite(v)||v<0||(v==0&&catalog(k).lower>0),error('EELS:Strategy','Invalid parameter %s.',catalog(k).name);end
            end
            s.initial_theta=[];s.stage1_hierarchy_context=struct();
            app.Report=[];app.ExportButton.Enable='off';app.HistoryTable.Data={};app.GroupTable.Data={};app.DirectionTable.Data={};
            cla(app.Stage1WeightAxes);cla(app.Stage2WeightAxes);cla(app.ConstraintAxes);cla(app.CostComponentAxes);
            for ax=[app.LossAxes,app.ZLPAxes,app.ValleyAxes,app.ResidualAxes,app.ProjectionAxes,app.HeatmapAxes],cla(ax);end
            app.PendingSettings=s;app.ApplyButton.Enable='on';
            app.Settings=s;app.ConfigurationOverrides=struct();app.refreshParameters();
            app.Status.Value={['Strategy loaded: ',file],'Not yet applied. Click Apply to Stage1/2 to activate.'};
        end
        function export(app)
            if isempty(app.Report)||app.Running,return;end
            root=uigetdir(app.OutputEdit.Value,'Save complete tuning dataset as');if isequal(root,0),return;end
            try
                directory=EELS_TuningDataset.newDirectory(root);app.saveRun(directory);app.Status.Value={['Saved: ',directory]};
            catch ME,app.Status.Value={ME.message};end
        end
        function close(app)
            if app.Running,app.StopRequested=true;app.Status.Value={'Stop requested. Save the dataset after the current candidate, then close.'};return;end
            delete(app);
        end
        function delete(app)
            if app.OwnsFigure
                if ~isempty(app.UIFigure)&&isgraphics(app.UIFigure),delete(app.UIFigure);end
            elseif ~isempty(app.Container)&&isgraphics(app.Container),delete(app.Container);
            end
        end
    end
end
