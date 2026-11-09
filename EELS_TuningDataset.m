classdef EELS_TuningDataset
    % Versioned, label-separated contextual-bandit records for reference tuning.
    methods (Static)
        function f=features(e,y,s)
            e=double(e(:));y=double(y(:));
            valid=isfinite(e)&isfinite(y);e=e(valid);y=y(valid);
            [e,ix]=sort(e);y=y(ix);
            if numel(e)<7||any(diff(e)<=0),error('EELS:Features','Need at least 7 distinct finite energy samples.');end
            [~,policy]=EELS_ZLPWindow.apply(e,y,s);
            local=find(abs(e)<=s.adaptive_zlp_region_search_meV);
            if isempty(local),local=(1:numel(e))';end
            [peak,k]=max(y(local));center=e(local(k));x=e-center;
            scale=max(max(abs(y)),realmin);z=y/scale;positive=max(z,0);
            area=max(trapz(e,positive),eps);edge=abs(x)>=prctile(abs(x),80);
            d=z-movmedian(z,5);noise=1.4826*median(abs(d(edge)-median(d(edge))));
            h=policy.fwhm_meV;asym=NaN;
            if isfield(policy,'half_height_crossings_meV')
                a=policy.half_height_crossings_meV-center;asym=(a(2)+a(1))/max(a(2)-a(1),eps);
            end
            names={'sample_count','energy_step_meV','energy_min_meV','energy_max_meV', ...
                'peak_energy_meV','peak_intensity','intensity_abs_scale','negative_fraction', ...
                'positive_area_per_peak_meV','centroid_meV','zlp_fwhm_meV','half_width_asymmetry', ...
                'tail_noise_per_peak','peak_to_noise_proxy','roughness_per_peak','peak_plateau_fraction'};
            values=[numel(e),median(diff(e)),min(e),max(e),center,peak,scale,mean(y<0), ...
                area,trapz(e,x.*positive)/area,h,asym,noise,1/max(noise,eps), ...
                median(abs(diff(z))),mean(y(local)>=.999*peak)];
            bands=[0 10;10 30;30 80;80 150;150 250];
            for sign=[-1 1]
                for b=1:size(bands,1)
                    mask=sign*x>=bands(b,1)&sign*x<=bands(b,2);
                    side='gain';if sign==1,side='loss';end
                    names{end+1}=sprintf('%s_area_fraction_%d_%d',side,bands(b,1),bands(b,2)); %#ok<AGROW>
                    if nnz(mask)>1,v=trapz(e(mask),positive(mask))/area;else,v=NaN;end
                    values(end+1)=v; %#ok<AGROW>
                end
            end
            grid=(-250:2:250)';sampled=interp1(x,z,grid,'linear',NaN);mask=isfinite(sampled);
            sampled(~mask)=0;
            f=struct('schema','eels-observation-v1','names',{names},'values',values, ...
                'valid',isfinite(values),'spectrum_grid_meV',grid, ...
                'spectrum_normalized',sampled,'spectrum_valid',mask, ...
                'normalization','max(abs(raw intensity)); energy relative to local raw peak', ...
                'window_policy',policy,'fit_crop_meV',s.crop_range, ...
                'label_free',true);
        end

        function snap=fitSnapshot(fit)
            snap=struct();if isempty(fit),return;end
            names={'ene_raw','ene_corrected','raw_data','zlp_only','fitted_spectrum','residual', ...
                'norm_val','theta','zlp_params','metrics','diagnostic_flags','error_budget','iter_history', ...
                'T_uncertainty','reliability_score','zlp_region_policy','numerical_backend'};
            for k=1:numel(names)
                if isfield(fit,names{k}),snap.(names{k})=fit.(names{k});end
            end
            snap.fit_weights=EELS_FitWeights.curves(fit);
            if isfield(fit,'stage2')
                [t,c]=EELS_PearsonSparse_FullFit.stage2_raw_profile(fit.stage2);
                snap.raw_valley=struct('temperature_K',t,'cost',c);
                snap.temperature_profile=fit.stage2.temperature_profile;
                snap.zlp_refine=fit.stage2.zlp_refine;
                snap.zlp_delta_norm=fit.stage2.zlp_delta_norm;
            end
        end

        function [xy,model]=embedding(history,p)
            n=numel(history);xy=zeros(n,2);model=struct('method','PCA of bounds-normalized log/linear parameters', ...
                'center',[],'components',[],'explained_percent',[0 0]);
            if n==0,return;end
            X=EELS_TuningSearch.encode(vertcat(history.values),p);mu=mean(X,1);X=X-mu;
            [~,D,V]=svd(X,'econ');q=min(2,size(V,2));
            xy(:,1:q)=X*V(:,1:q);energy=diag(D).^2;
            explained=zeros(1,2);if sum(energy)>0,explained(1:q)=100*energy(1:q)/sum(energy);end
            model.center=mu;model.components=V(:,1:q);model.explained_percent=explained;
        end

        function strategy=strategy(report)
            settings=eels_normalize_settings(report.best.settings);
            settings=EELS_PearsonSparse_FullFit.strip_runtime_settings(settings);
            settings.initial_theta=[];settings.stage1_hierarchy_context=struct();
            strategy=struct('schema','eels-fit-strategy-v1','settings',settings, ...
                'parameter_names',{{report.parameters.name}}, ...
                'reference_features',EELS_TuningDataset.features(report.energy_meV,report.spectrum,settings), ...
                'reference_fit',EELS_TuningDataset.fitSnapshot(report.best.result), ...
                'reference_temperature_K',report.known_T_K, ...
                'validation_status',report.validation_status, ...
                'initialization_policy','Transfer initial-width rules and penalties; re-estimate each spectrum. Reference theta is diagnostic, not a fixed fit or temperature seed.');
        end

        function directory=newDirectory(root)
            if ~isfolder(root),mkdir(root);end
            id=char(java.util.UUID.randomUUID());
            directory=fullfile(root,['run_',char(datetime('now','Format','yyyyMMdd_HHmmss_SSS')),'_',id(1:8)]);
            mkdir(directory);
        end

        function dataset=export(report,directory)
            if ~isfolder(directory),mkdir(directory);end
            strategy=EELS_TuningDataset.strategy(report);tuned_settings=strategy.settings; %#ok<NASGU>
            observation=EELS_TuningDataset.features(report.energy_meV,report.spectrum,report.baseline.settings);
            hist=report.history;n=numel(hist);actions=vertcat(hist.values);
            valid=isfinite([hist.score])&strcmp({hist.status},'ok');
            rewards=-[hist.score];rewards(~valid)=-20;
            if isfield(hist,'search_score'),rewards=-[hist.search_score];end
            [xy,pca_model]=EELS_TuningDataset.embedding(hist,report.parameters);
            hashSpectra=report.spectrum;
            if isfield(report,'training_spectra'),hashSpectra=report.training_spectra;end
            spectrum_hash=EELS_TuningDataset.hashBytes(typecast([double(report.energy_meV(:));double(hashSpectra(:))],'uint8'));
            code=struct();root=fileparts(mfilename('fullpath'));
            for name={'EELS_ReferenceTuner.m','EELS_DirectedPSO.m','EELS_GroupedPSO.m','EELS_GroupedSA.m','EELS_HybridSearch.m','EELS_ParameterCatalog.m','EELS_Stage2Tuning.m','EELS_TuningObjective.m','EELS_FitWeights.m','EELS_TuningSearch.m','EELS_TuningParallel.m','EELS_Parallel.m','EELS_TuningDataset.m','EELS_PearsonSparse_FullFit.m','EELS_ZLPWindow.m','native/eels_native_mex.cpp','eels_native_v22_mex.mexw64'}
                file=fullfile(root,name{1});if ~isfile(file),continue;end
                fid=fopen(file,'rb');bytes=fread(fid,Inf,'*uint8');fclose(fid);
                code.(matlab.lang.makeValidName(name{1}))=EELS_TuningDataset.hashBytes(bytes);
            end
            source='';if isfield(report,'source'),source=report.source;end
            dataset=struct('schema','eels-tuning-bandit-v1','created_utc',char(datetime('now','TimeZone','UTC','Format','yyyy-MM-dd HH:mm:ss')), ...
                'spectrum_sha256',spectrum_hash,'source',source,'split_group',spectrum_hash, ...
                'observation',observation,'action_names',{{report.parameters.name}}, ...
                'action_bounds',[report.parameters.lower;report.parameters.upper], ...
                'actions',actions,'actions_log_normalized',EELS_TuningSearch.encode(actions,report.parameters), ...
                'reward',rewards(:),'valid',valid(:),'terminal',true(n,1), ...
                'labels',struct('known_temperature_K',report.known_T_K, ...
                    'valley_temperature_K',[hist.valley_T_K]','temperature_error_K',[hist.error_K]', ...
                    'NRMSE',[hist.NRMSE]','score',[hist.score]','seconds',[hist.seconds]'), ...
                'status',{{hist.status}},'embedding',xy,'embedding_model',pca_model, ...
                'code_sha256',code,'matlab_version',version, ...
                'task_type','offline contextual bandit; each candidate is a one-step terminal action on the same spectrum', ...
                'behavior_policy','Seed 2206; baseline then coordinate log steps, every fifth attempt log-uniform global. Adaptive search; propensity unavailable.', ...
                'reward_definition','-score for valid candidates; -20 for failed/rejected candidates, use valid mask explicitly', ...
                'split_warning','Keep all actions for a spectrum AND the same physical acquisition in one train/test group. Source-group assignment must be completed before benchmark splitting.');
            if isfield(report,'execution'),dataset.execution=report.execution;end
            if isfield(hist,'annealing_temperature'),dataset.annealing_temperature=[hist.annealing_temperature]';end
            if isfield(report,'objective_settings'),dataset.objective_settings=report.objective_settings;end
            if isfield(report,'optimizer')
                dataset.optimizer=report.optimizer;
                dataset.behavior_policy=[report.optimizer.name,': ',report.optimizer.description,' Propensity unavailable.'];
            end
            dataset.parameter_inventory=table2struct(EELS_ParameterCatalog.inventory(report.baseline.settings));
            dataset.action_groups=cellfun(@EELS_ParameterCatalog.group,{report.parameters.name},'UniformOutput',false);
            dataset.integer_actions=EELS_ParameterCatalog.isInteger({report.parameters.name});
            if isfield(report,'training_spectra')
                observations=cell(1,size(report.training_spectra,2));
                for j=1:numel(observations)
                    observations{j}=EELS_TuningDataset.features(report.energy_meV,report.training_spectra(:,j),report.baseline.settings);
                end
                dataset.group_observations=observations;dataset.sample_info=report.sample_info;
                dataset.schema='eels-tuning-bandit-v2';
                dataset.validation=report.validation;
                dataset.point_outcomes={hist.point_outcomes};
                dataset.action_transform='log(x+offset); offset=upper*1e-4 only for zero lower bounds; attention weight strength uses linear coordinates';
                dataset.task_type='offline contextual bandit; each action fits a fixed group of spectra; per-point outcomes retained';
                ntrain=size(report.training_spectra,2);nvalidation=size(report.validation_spectra,2);
                before=[report.baseline.point_outcomes.valley_T_K]';after=[report.best.point_outcomes.valley_T_K]';
                if isfield(report.validation,'before_K')
                    before=[before;report.validation.before_K(:)];after=[after;report.validation.after_K(:)];
                else
                    before=[before;NaN(nvalidation,1)];after=[after;NaN(nvalidation,1)];
                end
                coordinates=NaN(ntrain+nvalidation,2);
                if isfield(report.sample_info,'train_yx')
                    coordinates=[report.sample_info.train_yx;report.sample_info.validation_yx];
                end
                joint=table([repmat("train",ntrain,1);repmat("heldout",nvalidation,1)], ...
                    coordinates(:,1),coordinates(:,2),before,after,abs(after-report.known_T_K), ...
                    'VariableNames',{'role','Y','X','before_K','best_K','best_error_K'});
                writetable(joint,fullfile(directory,'joint_points.csv'));
            end
            dataset.acquisition_group='';
            if isfield(hist,'search_score')
                dataset.schema='eels-tuning-bandit-v3';
                dataset.labels.search_score=[hist.search_score]';
                dataset.labels.constraint_violation=[hist.constraint_violation]';
                dataset.labels.regressed_count=[hist.regressed_count]';
                dataset.labels.failed_count=[hist.failed_count]';
                dataset.labels.max_fit_ratio=[hist.max_fit_ratio]';
                dataset.reward_definition='-search_score; exact transform in objective_settings.search_transform; valid means eligible for strategy selection; raw calibration losses retained in labels.score';
            end
            if isfield(report,'acquisition_group'),dataset.acquisition_group=report.acquisition_group;end
            save(fullfile(directory,'tuning_run.mat'),'report','strategy','tuned_settings','dataset','-v7.3');
            save(fullfile(directory,'fit_strategy.mat'),'strategy','tuned_settings');
            save(fullfile(directory,'bandit_dataset.mat'),'dataset','-v7.3');
            EELS_TuningDataset.writeJson(fullfile(directory,'manifest.json'),rmfield(dataset, ...
                {'actions','actions_log_normalized','embedding','observation','reward','valid','terminal','labels','status'}));
            EELS_TuningDataset.writeJson(fullfile(directory,'observation.json'),observation);
            EELS_TuningDataset.writeJson(fullfile(directory,'strategy.json'),strategy);
            feature_table=table(string(observation.names(:)),observation.values(:),observation.valid(:), ...
                'VariableNames',{'feature','value','valid'});writetable(feature_table,fullfile(directory,'features.csv'));
            t=array2table(actions,'VariableNames',{report.parameters.name});
            t.evaluation=(1:n)';t.score=[hist.score]';t.reward=rewards(:);t.valid=valid(:);
            t.valley_K=[hist.valley_T_K]';t.error_K=[hist.error_K]';t.NRMSE=[hist.NRMSE]';
            t.seconds=[hist.seconds]';t.status=string({hist.status})';
            if isfield(hist,'swarm_iteration'),t.swarm_iteration=[hist.swarm_iteration]';end
            if isfield(hist,'annealing_temperature'),t.annealing_temperature=[hist.annealing_temperature]';end
            if isfield(hist,'search_score')
                for key={'search_score','constraint_violation','regressed_count','failed_count','max_fit_ratio'}
                    t.(key{1})=reshape([hist.(key{1})],[],1);
                end
            end
            writetable(t,fullfile(directory,'trials.csv'));
            writetable(EELS_ParameterCatalog.inventory(report.baseline.settings),fullfile(directory,'parameter_inventory.csv'));
            if isfield(report.optimizer,'stages')
                stage=report.optimizer.stages;
                st=table(string({stage.group})',[stage.evaluations]',[stage.before]',[stage.after]', ...
                    'VariableNames',{'group','evaluations','before_loss','after_loss'});
                if isfield(stage,'method'),st.method=string({stage.method})';end
                if isfield(stage,'accepted_uphill'),st.accepted_uphill=[stage.accepted_uphill]';end
                names={report.parameters.name};
                st.parameters=string(arrayfun(@(a)strjoin(names(a.indices),';'),stage,'UniformOutput',false))';
                writetable(st,fullfile(directory,'search_stages.csv'));
            end
            writetable(EELS_Stage2Tuning.profile(report.baseline.result),fullfile(directory,'baseline_cost_components.csv'));
            writetable(EELS_Stage2Tuning.profile(report.best.result),fullfile(directory,'best_cost_components.csv'));
            if isfield(report.optimizer,'probes')&&~isempty(report.optimizer.probes)
                writetable(struct2table(report.optimizer.probes),fullfile(directory,'search_directions.csv'));
            end
            pointRows=table();baselineNR=max([hist(1).point_outcomes.NRMSE]',eps);
            for k=1:n
                points=hist(k).point_outcomes;m=numel(points);
                rows=table(repmat(k,m,1),(1:m)',[points.valley_T_K]',[points.error_K]', ...
                    [points.NRMSE]',[points.NRMSE]'./baselineNR,[points.score]',string({points.status})', ...
                    'VariableNames',{'evaluation','training_point','valley_K','error_K','NRMSE','fit_ratio','loss','status'});
                pointRows=[pointRows;rows]; %#ok<AGROW>
            end
            writetable(pointRows,fullfile(directory,'trial_points.csv'));
            fid=fopen(fullfile(directory,'trials.jsonl'),'w','n','UTF-8');
            if fid<0,error('EELS:DatasetExport','Cannot create trial log.');end
            closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
            for k=1:n
                row=struct('evaluation',k,'observation_id',spectrum_hash,'action',actions(k,:), ...
                    'action_log_normalized',dataset.actions_log_normalized(k,:),'reward',rewards(k), ...
                    'valid',valid(k),'terminal',true,'label_temperature_K',report.known_T_K, ...
                    'outcome',hist(k));
                fprintf(fid,'%s\n',jsonencode(row));
            end
        end

        function writeJson(file,value)
            fid=fopen(file,'w','n','UTF-8');if fid<0,error('EELS:DatasetExport','Cannot write %s',file);end
            closer=onCleanup(@()fclose(fid)); %#ok<NASGU>
            fprintf(fid,'%s',jsonencode(value,'PrettyPrint',true));
        end

        function hash=hashBytes(bytes)
            md=java.security.MessageDigest.getInstance('SHA-256');md.update(uint8(bytes));
            hash=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
        end
    end
end
