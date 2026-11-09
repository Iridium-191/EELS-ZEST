classdef EELS_ReferenceTuner
    methods (Static)
        function p=parameters()
            names={'zlp_fwhm_guard_factor','stage1_be_soft_scale','stage1_initial_width_scale', ...
                'stage1_zlp_gain_asym_width_ratio','stage1_zlp_gain_asym_shape_ratio', ...
                'stage1_zlp_gain_tail_prior_ratio','zlp_fwhm_core_factor', ...
                'stage1_zlp_core_data_weight','stage1_zlp_fartail_gain_boost', ...
                'stage1_zlp_gain_asym_shape_prior_weight','stage1_zlp_core_anchor_weight', ...
                'stage1_zlp_gain_width_boost','stage1_zlp_loss_width_scale', ...
                'lambda_pm_l1','lambda_pm_smooth','lambda_pm_l2','stage1_em_step_scale', ...
                'stage2_attention_start_scale','stage2_attention_span_scale', ...
                'stage2_attention_threshold','stage2_attention_fit_threshold','stage2_attention_peak_sigma_meV'};
            lo=[.3 .02 .3 .5 .5 .5 .25 .02 .3 .0001 .0005 .5 .25 0 1e-5 1e-7 .5 .25 .25 .03 .02 1];
            hi=[3.5 50 3 3 3 5 2.5 2 5 .15 .6 4 2.5 .01 .1 .01 3 2.5 1.6 .7 .6 20];
            p=struct('name',names,'lower',num2cell(lo),'upper',num2cell(hi), ...
                'enabled',num2cell(false(1,numel(names))));
            p=[p,EELS_FitWeights.parameters()];
            p=[p,EELS_Stage2Tuning.parameters()];
            p=[p,EELS_ParameterCatalog.extra()];
        end

        function report=tune(ene,spectrum,known_T,settings,options)
            if nargin<4,settings=struct();end
            if nargin<5,options=struct();end
            if ~isfield(options,'parameters')
                catalog=EELS_ReferenceTuner.parameters();active=EELS_ParameterCatalog.available(catalog,EELS_PearsonSparse_FullFit.defaults(settings));options.parameters=catalog(active);
            end
            if ~isfield(options,'max_evaluations'),options.max_evaluations=500;end
            if ~isfield(options,'algorithm'),options.algorithm='grouped_sa';end
            if ~isfield(options,'annealing_temperature'),options.annealing_temperature=.1;end
            if ~isfield(options,'swarm_size'),options.swarm_size=12;end
            if ~isfield(options,'max_fit_ratio'),options.max_fit_ratio=2;end
            if ~isfield(options,'search_radius'),options.search_radius=.12;end
            if ~isfield(options,'progress_callback'),options.progress_callback=[];end
            if ~isfield(options,'stop_callback'),options.stop_callback=[];end
            if ~isfield(options,'parallel_enabled'),options.parallel_enabled=true;end
            if ~isfield(options,'parallel_workers'),options.parallel_workers=EELS_Parallel.default_workers();end
            if ~isscalar(known_T)||~isfinite(known_T)||known_T<=0,error('EELS:TunerInput','Known temperature must be positive K.');end
            ene=double(ene(:));spectrum=double(spectrum);
            if isvector(spectrum),spectrum=spectrum(:);end
            spectra=spectrum;spectrum=spectra(:,1);ns=size(spectra,2);
            if ~isfield(options,'validation_spectra'),options.validation_spectra=[];end
            if ~isfield(options,'sample_info'),options.sample_info=struct();end
            validation=double(options.validation_spectra);
            if ~isempty(validation)&&(size(validation,1)~=numel(ene)||any(~isfinite(validation(:))))
                error('EELS:TunerInput','Invalid validation spectra.');
            end
            if numel(ene)~=numel(spectrum)||numel(ene)<20||any(~isfinite(ene))||any(~isfinite(spectra(:)))||any(diff(ene)<=0)
                error('EELS:TunerInput','Supply finite ascending energy (meV) and matching spectrum vectors.');
            end
            base=EELS_PearsonSparse_FullFit.defaults(settings);
            base=EELS_PearsonSparse_FullFit.strip_runtime_settings(base);
            base.initial_theta=[];base.stage1_hierarchy_context=struct();
            base.run_stage2=false;base.save_plots=false;base.save_diagnostic_plots=false;base.silent=true;
            if known_T<min(base.stage2_T_grid)||known_T>max(base.stage2_T_grid)
                error('EELS:TunerInput','Known temperature is outside the configured temperature scan.');
            end
            p=options.parameters;catalog=EELS_ReferenceTuner.parameters();
            if isempty(p)||~all(ismember({p.name},{catalog.name}))||numel(unique({p.name}))~=numel(p)
                error('EELS:TunerParameters','Choose unique parameters from the supported catalog.');
            end
            % A per-cube probe would otherwise replace the tuned EM grid.
            if any(strcmp({p.name},'stage1_em_step_scale'))
                base.stage1_fixed_EM=[];base.datacube_stage1_auto_sparse_pm_enabled=false;
            end
            if ns>1,base.datacube_stage1_auto_sparse_pm_enabled=false;end
            values=zeros(1,numel(p));
            EELS_TuningSearch.validateActive(p,base);
            for j=1:numel(p)
                if ~isfinite(p(j).lower)||~isfinite(p(j).upper)||p(j).lower<0||p(j).upper<=p(j).lower
                    error('EELS:TunerParameters','Parameter limits must be finite, nonnegative and ordered.');
                end
                values(j)=base.(p(j).name);
                if ~isfinite(values(j))||values(j)<p(j).lower||values(j)>p(j).upper
                    error('EELS:TunerParameters','Bounds for %s must contain its current value.',p(j).name);
                end
            end
            if ~isscalar(options.max_evaluations)||~isfinite(options.max_evaluations)||options.max_evaluations<1||options.max_evaluations>2000|| ...
                    ~isscalar(options.max_fit_ratio)||~isfinite(options.max_fit_ratio)||options.max_fit_ratio<1
                error('EELS:TunerInput','Invalid evaluation budget or acceptance limits.');
            end
            nmax=round(options.max_evaluations);
            if ~isscalar(options.search_radius)||~isfinite(options.search_radius)||options.search_radius<.025||options.search_radius>.5
                error('EELS:TunerInput','Search radius must be between 0.025 and 0.5 in normalized parameter space.');
            end
            algorithm=lower(char(options.algorithm));
            if ~ismember(algorithm,{'grouped_sa','hybrid','grouped_pso','pso','gp'}),error('EELS:TunerInput','Unknown search algorithm.');end
            if ismember(algorithm,{'grouped_sa','hybrid','grouped_pso'})
                groups=unique(cellfun(@EELS_ParameterCatalog.group,{p.name},'UniformOutput',false));
                if nmax<2*numel(groups)+3,error('EELS:GroupBudget','%d groups require at least %d evaluations; increase budget or select fewer groups.',numel(groups),2*numel(groups)+3);end
            end
            if ~isscalar(options.swarm_size)||~isfinite(options.swarm_size)||options.swarm_size<2||options.swarm_size>200
                error('EELS:TunerInput','Swarm size must be between 2 and 200.');
            end
            if ismember(algorithm,{'hybrid','grouped_pso','pso'})&&(exist('particleswarm','file')~=2||~license('test','GADS_Toolbox'))
                error('EELS:PSOUnavailable','MATLAB particleswarm requires a licensed Global Optimization Toolbox. Select GP to use the toolbox-independent search.');
            end
            if strcmp(algorithm,'hybrid')&&(exist('surrogateopt','file')~=2||exist('patternsearch','file')~=2)
                error('EELS:HybridUnavailable','Hybrid search requires surrogateopt and patternsearch; select grouped PSO or GP.');
            end
            if strcmp(algorithm,'grouped_sa')
                if exist('simulannealbnd','file')~=2||~license('test','GADS_Toolbox')
                    error('EELS:SAUnavailable','Simulated annealing requires Global Optimization Toolbox.');
                end
                if ~isscalar(options.annealing_temperature)||~isfinite(options.annealing_temperature)||options.annealing_temperature<.001||options.annealing_temperature>2
                    error('EELS:TunerInput','Annealing temperature must be between 0.001 and 2 (dimensionless).');
                end
            end
            if ~isscalar(options.parallel_enabled)||~isscalar(options.parallel_workers)|| ...
                    ~isfinite(options.parallel_workers)||options.parallel_workers<1||options.parallel_workers~=round(options.parallel_workers)
                error('EELS:TunerInput','Invalid parallel controls.');
            end
            options.parallel_workers=min(options.parallel_workers,max(ns,size(validation,2)));
            execution=struct('mode','serial_points','workers',1,'native_backend',base.native_backend);
            if options.parallel_enabled
                base.native_backend='cpp';
                execution=EELS_TuningParallel.prepare(options.parallel_workers);
            end
            cache=containers.Map('KeyType','char','ValueType','any');
            records=struct('evaluation',{},'values',{},'score',{},'valley_T_K',{},'reported_T_K',{}, ...
                'error_K',{},'NRMSE',{},'seconds',{},'status',{},'fit_snapshot',{}, ...
                'search_attempt',{},'proposal',{},'incumbent_values',{},'search_radius',{},'acquisition',{},'swarm_iteration',{},'point_outcomes',{}, ...
                'search_score',{},'feasible',{},'regressed_count',{},'failed_count',{},'max_fit_ratio',{},'constraint_violation',{},'annealing_temperature',{});
            best=[];baseline=[];baseline_nrmse=NaN(1,ns);best_values=values;cancelled=false;
            lookup=containers.Map('KeyType','char','ValueType','double');
            directions=ones(1,numel(p));
            optimizer=struct('name','zlp_trust_gp_v1','seed',2206, ...
                'description','Residual-ordered coordinate probes followed by local GP expected improvement.', ...
                'active_dimensions',numel(p));
            evaluate(values,struct('method','baseline','radius',NaN,'acquisition',NaN,'iteration',-1));
            if strcmp(algorithm,'gp')
                for attempt=2:nmax
                    drawnow limitrate;
                    if ~isempty(options.stop_callback)&&options.stop_callback(),cancelled=true;break;end
                    [candidate,search]=EELS_TuningSearch.propose(records,p,values,directions);search.iteration=-1;
                    try,evaluate(candidate,search);
                    catch ME,if cancelled,break;else,rethrow(ME);end;end
                end
                optimizer.initial_probe_directions=directions;
            else
                optimizer.name='measured_direction_local_matlab_pso_v1';
                optimizer.description='Measured directional probes then local low-inertia MATLAB particleswarm.';
                try
                    if strcmp(algorithm,'grouped_sa')
                        optimizer=EELS_GroupedSA.run(p,nmax,options.search_radius,options.annealing_temperature,@evaluate,@incumbent,@recordCount);
                    elseif strcmp(algorithm,'hybrid')
                        optimizer=EELS_HybridSearch.run(p,nmax,round(options.swarm_size),options.search_radius,@evaluate,@incumbent,@recordCount);
                    elseif strcmp(algorithm,'grouped_pso')
                        optimizer=EELS_GroupedPSO.run(p,nmax,round(options.swarm_size),options.search_radius,@evaluate,@incumbent,@recordCount);
                    else
                        optimizer=EELS_DirectedPSO.run(values,p,nmax,round(options.swarm_size),options.search_radius,@evaluate,records(1).search_score);
                    end
                catch ME
                    if cancelled,optimizer.stop_reason='user_stop';else,rethrow(ME);end
                end
            end
            if isempty(best),error('EELS:TunerCancelled','Stopped before baseline evaluation.');end
            report=struct('version','2.2','known_T_K',known_T,'energy_meV',ene,'spectrum',spectrum, ...
                'parameters',p,'baseline',baseline,'best',best,'history',records,'cancelled',cancelled, ...
                'validation_status','not_validated', ...
                'objective','mean absolute raw-valley temperature error in K');
            report.objective_settings=struct('loss_type','temperature_mae','loss_units','K', ...
                'max_fit_ratio',options.max_fit_ratio,'search_radius',options.search_radius, ...
                'random_seed',2206,'max_evaluations',options.max_evaluations);
            report.optimizer=optimizer;report.execution=execution;
            report.training_spectra=spectra;report.sample_info=options.sample_info;
            report.validation=struct('status','not_requested');
            report.validation_spectra=validation;
            report.objective_settings.training_count=ns;
            report.objective_settings.constraint_policy='No residual penalty in loss or search_score; exported best strategy must pass every point NRMSE gate.';
            report.objective_settings.search_transform='search_score=mean(abs(valley_T_K-known_T_K)); failed: 1e6+failed_count, raw loss=Inf';
            report.objective_settings.aggregate='arithmetic mean of absolute temperature errors, in K; no worst-point weighting';
            if ~isempty(validation)&&~cancelled
                report.validation=EELS_ReferenceTuner.validateHoldout(ene,validation,known_T,baseline.settings,best.settings,options);
                report.validation_status=report.validation.status;
            end


            function value=incumbent()
                value=struct('values',best_values,'score',best.score);
            end
            function n=recordCount(),n=numel(records);end
            function score=evaluate(candidate,search)
                if ~isfield(search,'annealing_temperature'),search.annealing_temperature=NaN;end
                integer=EELS_ParameterCatalog.isInteger({p.name});candidate(integer)=round(candidate(integer));
                key=sprintf('%.12g,',candidate);
                if isKey(lookup,key),score=lookup(key);return;end
                if numel(records)>=nmax
                    error('EELS:PSOBudget','Evaluation budget reached.');
                end
                k=numel(records)+1;incumbent_values=best_values;proposal=search.method;
                s=base;for j=1:numel(p),s.(p(j).name)=candidate(j);end
                s=EELS_PearsonSparse_FullFit.defaults(s);
                for j=1:numel(p)
                    if s.(p(j).name)~=candidate(j),error('EELS:TunerParameters','Solver overrides %s.',p(j).name);end
                end
                timer=tic;fit=[];v=[];score=Inf;status='ok';nr=NaN;reported=NaN;
                points=repmat(struct('valley_T_K',NaN,'error_K',NaN,'NRMSE',NaN, ...
                    'score',Inf,'status','not_evaluated','worker_pid',NaN,'native_version','', ...
                    'fit_ratio',NaN,'truth_excess',NaN),1,ns);
                batch=[];
                if options.parallel_enabled
                    try
                        batch=EELS_TuningParallel.fit(ene,spectra,s,known_T,options,sprintf('Candidate %d',k));
                    catch ME
                        if strcmp(ME.identifier,'EELS:PSOStop'),cancelled=true;end
                        rethrow(ME);
                    end
                end
                for point=1:ns
                    drawnow limitrate;
                    if ~isempty(options.stop_callback)&&options.stop_callback()
                        cancelled=true;error('EELS:PSOStop','User requested stop.');
                    end
                    try
                        if ~isempty(batch)
                            out=batch{point};
                            if ~strcmp(out.status,'ok'),error('EELS:TunerPoint','%s',out.status);end
                            current=out.result;valley=out.valley;nrmse=out.NRMSE;
                        else
                            stage1_indices=~startsWith({p.name},'stage2_');
                            cachekey=[sprintf('%d:',point),sprintf('%.12g,',candidate(stage1_indices))];
                            if isKey(cache,cachekey),first=cache(cachekey);
                            else
                                first=EELS_PearsonSparse_FullFit.solve_stage1(ene,spectra(:,point),s);
                                % Bound memory even with many pixels and evaluations.
                                if cache.Count>=2*ns,remove(cache,keys(cache));end
                                cache(cachekey)=first;
                            end
                            current=EELS_PearsonSparse_FullFit.promote_stage1_to_stage2(first,s);
                            valley=EELS_ReferenceTuner.valley(current,known_T);nrmse=current.metrics.NRMSE;
                        end
                        if ~isfinite(nrmse)||nrmse<0||~isfinite(valley.T_K)||~isfinite(valley.error_K)
                            error('EELS:TunerFit','Invalid residual or temperature metric.');
                        end
                        if k==1,baseline_nrmse(point)=nrmse;end
                        loss=abs(valley.T_K-known_T);
                        state='ok';
                        if nrmse>options.max_fit_ratio*max(baseline_nrmse(point),eps),state='fit_regression';end
                        points(point)=struct('valley_T_K',valley.T_K,'error_K',valley.error_K, ...
                            'NRMSE',nrmse,'score',loss,'status',state,'worker_pid',feature('getpid'),'native_version','', ...
                            'fit_ratio',nrmse/max(baseline_nrmse(point),eps),'truth_excess',valley.truth_excess);
                        if ~isempty(batch),points(point).worker_pid=out.worker_pid;points(point).native_version=out.native_version;end
                        if point==1,fit=current;v=valley;reported=current.stage2.T_fit;end
                    catch ME
                        if k==1,rethrow(ME);end
                        points(point).status=[ME.identifier,': ',ME.message];
                    end
                    if isempty(batch),EELS_TuningParallel.notify(options,point,ns,sprintf('Candidate %d',k));end
                end
                diagnosis=EELS_TuningObjective.aggregate(points,baseline_nrmse,options.max_fit_ratio);
                score=diagnosis.score;status=diagnosis.status;
                nr=mean([points.NRMSE]);
                tv=NaN;err=NaN;if ~isempty(v),tv=mean([points.valley_T_K]);err=mean([points.error_K]);end
                records(end+1)=struct('evaluation',numel(records)+1,'values',candidate,'score',score, ...
                    'valley_T_K',tv,'reported_T_K',reported,'error_K',err,'NRMSE',nr,'seconds',toc(timer),'status',status, ...
                    'fit_snapshot',EELS_TuningDataset.fitSnapshot(fit), ...
                    'search_attempt',k,'proposal',proposal,'incumbent_values',incumbent_values, ...
                    'search_radius',search.radius,'acquisition',search.acquisition,'swarm_iteration',search.iteration,'point_outcomes',points, ...
                    'search_score',diagnosis.search_score,'feasible',diagnosis.feasible, ...
                    'regressed_count',diagnosis.regressed_count,'failed_count',diagnosis.failed_count, ...
                    'max_fit_ratio',diagnosis.max_fit_ratio,'constraint_violation',diagnosis.constraint_violation, ...
                    'annealing_temperature',search.annealing_temperature); %#ok<AGROW>
                exported=s;exported.run_stage2=true;
                entry=struct('settings',exported,'result',fit,'valley',v,'score',score,'point_outcomes',points,'loss_units','K');
                if k==1,baseline=entry;directions=EELS_TuningSearch.residualDirections(fit,p);end
                if diagnosis.feasible&&(isempty(best)||score<best.score),best=entry;best_values=candidate;end
                if k==1&&isempty(best),error('EELS:TunerBaseline','Baseline does not have a finite feasible calibration loss.');end
                lookup(key)=diagnosis.search_score;
                if ~isempty(options.progress_callback),options.progress_callback(records(end),best);end
                score=diagnosis.search_score;
            end

        end

        function result=validateHoldout(e,Y,T,before,after,options)
            result=struct('status','heldout_evaluated','before_K',NaN(1,size(Y,2)), ...
                'after_K',NaN(1,size(Y,2)),'point_status',{repmat({'not_evaluated'},1,size(Y,2))});
            try
                a=EELS_TuningParallel.fit(e,Y,before,T,options,'Heldout baseline');
                b=EELS_TuningParallel.fit(e,Y,after,T,options,'Heldout best');
                for j=1:size(Y,2)
                    if strcmp(a{j}.status,'ok')&&strcmp(b{j}.status,'ok')
                        result.before_K(j)=a{j}.valley.T_K;result.after_K(j)=b{j}.valley.T_K;result.point_status{j}='ok';
                    else
                        result.point_status{j}=[a{j}.status,' | ',b{j}.status];result.status='heldout_has_failures';
                    end
                end
            catch ME
                if strcmp(ME.identifier,'EELS:PSOStop'),result.status='heldout_cancelled';else,rethrow(ME);end
            end
            result.before_MAE_K=mean(abs(result.before_K-T));
            result.after_MAE_K=mean(abs(result.after_K-T));
        end

        function [Y,V,info]=sampleCube(cube,columns,ntrain,nvalidation,seed)
            if ischar(columns)||isstring(columns)
                if strcmpi(strtrim(columns),'all'),columns=[1 size(cube,2)];
                else,columns=sscanf(columns,'%f')';end
            end
            if numel(columns)~=2||any(columns~=round(columns))||columns(1)<1||columns(2)>size(cube,2)||columns(1)>columns(2)
                error('EELS:TunerSample','X range must be all or two valid column indices.');
            end
            if ntrain<1||nvalidation<0||any([ntrain nvalidation seed]~=round([ntrain nvalidation seed]))||seed<0
                error('EELS:TunerSample','Counts and seed must be valid integers.');
            end
            [yy,xx]=ndgrid(1:size(cube,1),columns(1):columns(2));
            indices=sub2ind([size(cube,1),size(cube,2)],yy(:),xx(:));
            data=reshape(cube,[],size(cube,3));
            valid=all(isfinite(data(indices,:)),2)&sum(abs(data(indices,:)),2)>0;
            indices=indices(valid);
            requested=[ntrain nvalidation];
            if numel(indices)<ntrain+nvalidation
                nvalidation=min(nvalidation,floor(numel(indices)/5));ntrain=min(ntrain,numel(indices)-nvalidation);
            end
            if ntrain<1,error('EELS:TunerSample','No finite, nonempty spectra in the selected region.');end
            stream=RandStream('mt19937ar','Seed',seed);order=randperm(stream,numel(indices),ntrain+nvalidation);
            selected=indices(order);train=selected(1:ntrain);test=selected(ntrain+1:end);
            Y=double(data(train,:)');V=double(data(test,:)');
            [ty,tx]=ind2sub([size(cube,1),size(cube,2)],train);
            [vy,vx]=ind2sub([size(cube,1),size(cube,2)],test);
            info=struct('requested_counts',requested,'seed',seed,'map_size',[size(cube,1),size(cube,2)],'x_range',columns, ...
                'train_linear_indices',train,'validation_linear_indices',test, ...
                'train_yx',[ty(:),tx(:)],'validation_yx',[vy(:),vx(:)], ...
                'split_policy','Fixed random disjoint pixels; same known temperature; heldout never selects parameters.');
        end

        function v=valley(result,truth)
            [T,C]=EELS_PearsonSparse_FullFit.stage2_raw_profile(result.stage2);
            valid=isfinite(T)&isfinite(C);T=T(valid);C=C(valid);[T,order]=sort(T);C=C(order);
            if numel(T)<3||truth<min(T)||truth>max(T),error('EELS:TunerValley','Incomplete temperature profile.');end
            [minimum,j]=min(C);span=max(C)-minimum;
            if span<=1e-10*max(max(abs(C)),realmin),error('EELS:TunerValley','Flat cost profile cannot calibrate temperature.');end
            tied=T(C<=minimum+span*1e-9);
            % Score the reported raw minimum; retain tie ambiguity as a diagnostic.
            err=abs(T(j)-truth);
            v=struct('T_K',T(j),'error_K',err,'T_grid',T,'cost',C, ...
                'truth_excess',(interp1(T,C,truth)-minimum)/span, ...
                'minimum_width_K',max(tied)-min(tied),'tied_worst_error_K',max(abs(tied-truth)), ...
                'boundary_minimum',j==1||j==numel(T));
        end
    end
end
