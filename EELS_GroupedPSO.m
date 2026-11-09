classdef EELS_GroupedPSO
    methods (Static)
        function meta=run(p,budget,swarm,radius,evaluator,getBest,getCount)
            group=cellfun(@EELS_ParameterCatalog.group,{p.name},'UniformOutput',false);
            groups=unique(group,'stable');n=numel(groups);
            if budget<2*n+3,error('EELS:GroupBudget','Selected %d groups need at least %d evaluations; increase budget or select fewer groups.',n,2*n+3);end
            % Baseline consumes one slot; reserve at least two new trials per group.
            reserve=max(2,min(max(3,round(.3*(budget-1))),budget-1-2*n));limit=budget-reserve;
            passes=1;if limit>=8*n,passes=2;end
            plans=repmat(groups,1,passes);stages=struct('group',{},'evaluations',{},'before',{},'after',{},'indices',{});
            for j=1:numel(plans)
                remaining=limit-getCount();if remaining<2,break;end
                allowance=max(2,floor(remaining/(numel(plans)-j+1)));
                idx=find(strcmp(group,plans{j}));runBlock(idx,allowance,plans{j});
            end
            % All selected dimensions return, including individually insensitive ones.
            if getCount()<budget,runBlock(1:numel(p),budget-getCount(),'Joint interaction pass');end
            meta=struct('name','grouped_matlab_pso_v1','seed',2206, ...
                'description','Groupwise low-inertia PSO, carrying the accepted incumbent, followed by an all-selected-parameter interaction pass.', ...
                'active_dimensions',numel(p),'groups',{groups},'stages',stages,'evaluations',getCount(),'stop_reason','schedule_complete');
            function runBlock(idx,allowance,label)
                before=getBest();anchor=before.values;start=getCount();
                EELS_DirectedPSO.run(anchor(idx),p(idx),allowance+1,swarm,radius,@evaluateBlock,before.score);
                after=getBest();stages(end+1)=struct('group',label,'evaluations',getCount()-start, ...
                    'before',before.score,'after',after.score,'indices',idx);
                function value=evaluateBlock(x,info)
                    full=anchor;full(idx)=x;info.method=['group: ',label,' / ',info.method];
                    value=evaluator(full,info);
                end
            end
        end
    end
end
