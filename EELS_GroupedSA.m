classdef EELS_GroupedSA
    methods (Static)
        function meta=run(p,budget,radius,temperature,evaluator,getBest,getCount)
            group=cellfun(@EELS_ParameterCatalog.group,{p.name},'UniformOutput',false);
            groups=unique(group,'stable');n=numel(groups);
            if budget<2*n+3,error('EELS:GroupBudget','%d groups need at least %d evaluations.',n,2*n+3);end
            % Identical group scheduling and joint reservation to grouped PSO.
            reserve=max(2,min(max(3,round(.3*(budget-1))),budget-1-2*n));limit=budget-reserve;
            passes=1;if limit>=8*n,passes=2;end
            plans=repmat(groups,1,passes);
            stages=struct('group',{},'evaluations',{},'before',{},'after',{},'indices',{},'method',{},'accepted_uphill',{});
            saved=rng;cleanup=onCleanup(@()rng(saved));rng(2206,'twister'); %#ok<NASGU>
            cache=containers.Map('KeyType','char','ValueType','double');
            initial=getBest();cache(key(initial.values))=initial.score;
            for j=1:numel(plans)
                remaining=limit-getCount();if remaining<2,break;end
                allowance=max(2,floor(remaining/(numel(plans)-j+1)));
                runBlock(find(strcmp(group,plans{j})),allowance,plans{j});
            end
            if getCount()<budget,runBlock(1:numel(p),budget-getCount(),'Joint interaction pass');end
            meta=struct('name','grouped_matlab_simulated_annealing_v1','seed',2206, ...
                'description','Groupwise MATLAB simulated annealing followed by all-selected annealing; feasible incumbent exported independently of annealing chain.', ...
                'active_dimensions',numel(p),'groups',{groups},'stages',stages,'evaluations',getCount(), ...
                'stop_reason','schedule_complete','initial_temperature',temperature,'cooling_rate',.95, ...
                'score_scaling','search_score / max(abs(block incumbent loss),0.1)', ...
                'proposal','Gaussian per-coordinate step; reflection at global bounds; canonical integer moves');
            function runBlock(idx,allowance,label)
                before=getBest();anchor=before.values;start=getCount();stopAt=min(budget,start+allowance);
                z0=EELS_TuningSearch.encode(anchor(idx),p(idx));scale=max(abs(before.score),.1);
                currentTemperature=temperature;iteration=0;uphill=0;previous=before.score/scale;
                opts=optimoptions('simulannealbnd','Display','off','PlotFcn',[], ...
                    'AnnealingFcn',@proposal,'OutputFcn',@observe,'InitialTemperature',temperature, ...
                    'TemperatureFcn',@temperatureexp,'AcceptanceFcn',@acceptancesa, ...
                    'ReannealInterval',max(100,10*allowance), ...
                    'MaxFunctionEvaluations',10*allowance+10,'MaxIterations',10*allowance+10, ...
                    'MaxStallIterations',10*allowance+10,'FunctionTolerance',0,'HybridFcn',[]);
                try
                    simulannealbnd(@sample,z0,zeros(size(z0)),ones(size(z0)),opts);
                catch ME
                    if ~strcmp(ME.identifier,'EELS:SABudget'),rethrow(ME);end
                end
                after=getBest();stages(end+1)=struct('group',label,'evaluations',getCount()-start, ...
                    'before',before.score,'after',after.score,'indices',idx,'method','simulated_annealing','accepted_uphill',uphill);
                function z=proposal(v,~)
                    iteration=v.iteration;currentTemperature=max(v.temperature);
                    step=radius*sqrt(currentTemperature/temperature);
                    z=reshape(v.x,1,[])+step*randn(size(z0));
                    % Reflection avoids piling proposals onto bounds.
                    z=1-abs(mod(z,2)-1);
                    x=EELS_TuningSearch.decode(z,p(idx));int=EELS_ParameterCatalog.isInteger({p(idx).name});
                    % Adjacent integer moves remain possible even at a small radius.
                    if any(int)&&rand<.25
                        old=EELS_TuningSearch.decode(reshape(v.x,1,[]),p(idx));
                        x(int)=min(floor([p(idx(int)).upper]),max(ceil([p(idx(int)).lower]),old(int)+2*(rand> .5)-1));
                    end
                    z=EELS_TuningSearch.encode(x,p(idx));
                end
                function value=sample(z)
                    full=anchor;full(idx)=EELS_TuningSearch.decode(reshape(z,1,[]),p(idx));q=key(full);
                    if isKey(cache,q),value=cache(q)/scale;return;end
                    if getCount()>=stopAt,error('EELS:SABudget','Stage evaluation budget reached.');end
                    info=struct('method',['group: ',label,' / simulated_annealing'], ...
                        'radius',radius*sqrt(currentTemperature/temperature),'acquisition',NaN, ...
                        'iteration',iteration,'annealing_temperature',currentTemperature);
                    value=evaluator(full,info);cache(q)=value;value=value/scale;
                end
                function [stop,options,changed]=observe(options,v,flag)
                    changed=false;stop=getCount()>=stopAt;
                    if strcmp(flag,'iter')
                        if v.fval>previous+1e-12,uphill=uphill+1;end
                        previous=v.fval;
                    end
                end
            end
            function q=key(x),q=sprintf('%.12g,',x);end
        end
    end
end
