classdef EELS_HybridSearch
    % Small conditional surrogates, full interaction PSO, then local polishing.
    methods (Static)
        function meta=run(p,budget,swarm,radius,evaluator,getBest,getCount)
            group=cellfun(@EELS_ParameterCatalog.group,{p.name},'UniformOutput',false);
            groups=unique(group,'stable');n=numel(groups);
            if budget<2*n+3,error('EELS:GroupBudget','%d groups need at least %d evaluations.',n,2*n+3);end
            saved=rng;cleanup=onCleanup(@()rng(saved));rng(2206,'twister'); %#ok<NASGU>
            startBest=getBest();initial=startBest.values;
            cache=containers.Map('KeyType','char','ValueType','double');cache(key(initial))=startBest.score;
            stages=struct('group',{},'evaluations',{},'before',{},'after',{},'indices',{},'method',{});
            localBudget=max(1,floor(.15*(budget-1)));
            jointBudget=max(1,floor(.25*(budget-1)));
            groupBudget=max(2*n,budget-1-localBudget-jointBudget);
            if groupBudget+localBudget+jointBudget>budget-1
                localBudget=1;jointBudget=budget-1-groupBudget-localBudget;
            end
            weight=ones(1,n);weight(strcmp(groups,'Attention and BE selection'))=2;
            allocation=2*ones(1,n);extra=groupBudget-2*n;
            shares=extra*weight/sum(weight);allocation=allocation+floor(shares);
            [~,order]=sort(shares-floor(shares),'descend');
            for k=1:extra-sum(floor(shares)),allocation(order(k))=allocation(order(k))+1;end
            priority={'stage2_attention_start_scale','stage2_attention_fit_threshold', ...
                'stage2_attention_threshold','stage2_attention_peak_sigma_meV', ...
                'stage2_hard_be_pair_threshold','stage2_attention_span_scale', ...
                'stage1_zlp_gain_asym_width_ratio','stage1_zlp_gain_asym_shape_ratio', ...
                'stage1_zlp_core_data_weight','zlp_fwhm_guard_factor','zlp_wide_component_l2', ...
                'lambda_pm_smooth','lambda_pm_l2','stage1_be_soft_scale','stage2_pm_smooth_weight', ...
                'stage2_gain_weight_ratio','stage2_lowfreq_weight_ratio','stage2_highfreq_weight_ratio', ...
                'stage2_zlp_refine_tail_asym_l2','stage2_zlp_refine_tilt_l2','stage2_zlp_refine_width_l2'};
            priorityIdx=[];
            for k=1:numel(priority),priorityIdx=[priorityIdx find(strcmp({p.name},priority{k}))];end %#ok<AGROW>
            priorityIdx=[priorityIdx,setdiff(1:numel(p),priorityIdx,'stable')];
            integer=EELS_ParameterCatalog.isInteger({p.name});
            for k=1:n
                idx=priorityIdx(strcmp(group(priorityIdx),groups{k})&~integer(priorityIdx));
                % Reserve roughly half each block for adaptive proposals after initialization.
                idx=idx(1:min(numel(idx),max(1,min(6,floor(allocation(k)/2)))));
                if isempty(idx)
                    runBlock(find(strcmp(group,groups{k})),allocation(k),groups{k},'pso');
                else
                    runBlock(idx,allocation(k),groups{k},'surrogate');
                end
            end
            % All selected dimensions participate here, including integers and unprobed controls.
            runBlock(1:numel(p),max(1,budget-getCount()-localBudget),'Joint interaction pass','pso');
            afterJoint=getBest();movement=abs(EELS_TuningSearch.encode(afterJoint.values,p)-EELS_TuningSearch.encode(initial,p));
            [~,rank]=sort(movement(priorityIdx),'descend');idx=priorityIdx(rank);idx=idx(~integer(idx));
            idx=idx(1:min(3,numel(idx)));
            if ~isempty(idx)&&getCount()<budget
                runBlock(idx,budget-getCount(),'Local refinement','patternsearch');
            end
            meta=struct('name','conditional_surrogate_joint_pso_pattern_v1','seed',2206, ...
                'description','Budgeted small conditional surrogates, all-selected PSO, and at most three-dimensional pattern search.', ...
                'active_dimensions',numel(p),'groups',{groups},'stages',stages, ...
                'evaluations',getCount(),'stop_reason','schedule_complete', ...
                'group_allocation',allocation,'attention_budget_multiplier',2, ...
                'history_reuse','Only exact complete candidates from this run; no cross-context projected losses.');
            function runBlock(idx,allowance,label,method)
                allowance=min(allowance,budget-getCount());if allowance<1,return;end
                before=getBest();anchor=before.values;start=getCount();limit=start+allowance;
                z0=EELS_TuningSearch.encode(anchor(idx),p(idx));lo=max(0,z0-radius);hi=min(1,z0+radius);
                try
                    switch method
                        case 'surrogate'
                            opts=optimoptions('surrogateopt','Display','off','PlotFcn',[], ...
                                'UseParallel',false,'MaxFunctionEvaluations',allowance+1, ...
                                'MinSurrogatePoints',numel(idx)+1, ...
                                'InitialPoints',struct('X',z0,'Fval',before.score));
                            surrogateopt(@sample,lo,hi,opts);
                        case 'patternsearch'
                            opts=optimoptions('patternsearch','Display','off','UseParallel',false, ...
                                'MaxFunctionEvaluations',allowance+1,'MaxIterations',allowance+1, ...
                                'InitialMeshSize',min(.04,radius/2),'MeshTolerance',1e-4);
                            patternsearch(@sample,z0,[],[],[],[],lo,hi,[],opts);
                        otherwise
                            EELS_DirectedPSO.run(anchor(idx),p(idx),allowance+1,swarm,radius,@sampleValues,before.score);
                    end
                catch ME
                    if ~strcmp(ME.identifier,'EELS:HybridBudget'),rethrow(ME);end
                end
                after=getBest();stages(end+1)=struct('group',label,'evaluations',getCount()-start, ...
                    'before',before.score,'after',after.score,'indices',idx,'method',method);
                function score=sample(z)
                    info=struct('method',method,'radius',radius,'acquisition',NaN,'iteration',-1);
                    score=sampleValues(EELS_TuningSearch.decode(z,p(idx)),info);
                end
                function score=sampleValues(x,info)
                    full=anchor;full(idx)=x;q=key(full);
                    if isKey(cache,q),score=cache(q);return;end
                    if getCount()>=limit,error('EELS:HybridBudget','Stage budget reached.');end
                    info.method=['group: ',label,' / ',info.method];score=evaluator(full,info);cache(q)=score;
                end
            end
            function q=key(x),q=sprintf('%.12g,',x);end
        end
    end
end
