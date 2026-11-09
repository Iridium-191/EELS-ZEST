classdef EELS_DirectedPSO
    methods (Static)
        function meta=run(initial,p,budget,swarmSize,radius,evaluator,baselineScore)
            z0=EELS_TuningSearch.encode(initial,p);X=z0;F=baselineScore;
            probes=struct('parameter',{},'direction',{},'search_gain',{},'step',{});
            epochs=struct('radius',{},'before',{},'after',{},'evaluations',{});
            cache=containers.Map('KeyType','char','ValueType','double');cache(key(z0))=baselineScore;
            saved=rng;cleanup=onCleanup(@()rng(saved));rng(2206,'twister'); %#ok<NASGU>
            priority={'stage2_zlp_refine_tail_asym_l2','stage2_zlp_refine_tilt_l2','stage2_zlp_refine_width_l2','stage2_zlp_refine_shift_l2', ...
                'stage2_gain_weight_ratio','stage2_lowfreq_weight_ratio','stage2_highfreq_weight_ratio', ...
                'stage2_attention_weight_strength','stage1_shoulder_weight_ratio','stage1_dynamic_weight_scale', ...
                'stage1_vibration_weight_ratio','stage2_dynamic_weight_scale'};
            order=[];for j=1:numel(priority),order=[order find(strcmp({p.name},priority{j}))];end %#ok<AGROW>
            order=[order,setdiff(1:numel(p),order,'stable')];
            % Reserve most of a small budget for the swarm; probes are paired.
            count=min(numel(p),floor(.4*(budget-1)/2));step=min(.06,radius/2);
            epochEnd=budget;phase='direction_probe';epoch=0;
            for j=order(1:count)
                for sign=[-1 1]
                    z=z0;z(j)=min(1,max(0,z(j)+sign*step));
                    if abs(z(j)-z0(j))<1e-12,continue;end
                    value=sample(z);
                    probes(end+1)=struct('parameter',p(j).name,'direction',sign, ...
                        'search_gain',baselineScore-value,'step',abs(z(j)-z0(j))); %#ok<AGROW>
                end
            end
            % Combine only individually measured improving moves, at half step.
            if numel(F)<budget&&~isempty(probes)
                z=z0;
                for j=order(1:count)
                    rows=find(strcmp({probes.parameter},p(j).name));
                    if isempty(rows),continue;end
                    [gain,k]=max([probes(rows).search_gain]);
                    if gain>0,q=probes(rows(k));z(j)=z(j)+.5*q.direction*q.step;end
                end
                phase='measured_direction_combination';sample(min(1,max(0,z)));
            end
            while numel(F)<budget
                epoch=epoch+1;[previous,ib]=min(F);center=X(ib,:);start=numel(F);
                localLo=max(0,center-radius);localHi=min(1,center+radius);
                sizeSwarm=min(swarmSize,max(2,budget-start));
                epochEnd=min(budget,start+2*sizeSwarm);phase='local_particleswarm';
                swarm=repmat(center,sizeSwarm,1)+(.35*radius)*randn(sizeSwarm,numel(p));
                swarm=min(localHi,max(localLo,swarm));swarm(1,:)=center;
                [~,rank]=sort(F);slot=2;
                for k=rank
                    if k==ib||any(X(k,:)<localLo)||any(X(k,:)>localHi),continue;end
                    swarm(slot,:)=X(k,:);slot=slot+1;if slot>ceil(sizeSwarm/2),break;end
                end
                opts=optimoptions('particleswarm','Display','off','SwarmSize',sizeSwarm, ...
                    'InitialSwarmMatrix',swarm,'InitialSwarmSpan',radius, ...
                    'InertiaRange',[.2 .5],'SelfAdjustmentWeight',.8,'SocialAdjustmentWeight',1.2, ...
                    'MinNeighborsFraction',1,'MaxIterations',4,'MaxStallIterations',4, ...
                    'FunctionTolerance',1e-8,'UseParallel',false);
                try
                    particleswarm(@sample,numel(p),localLo,localHi,opts);
                catch ME
                    if ~strcmp(ME.identifier,'EELS:DirectedEpoch'),rethrow(ME);end
                end
                epochs(end+1)=struct('radius',radius,'before',previous,'after',min(F),'evaluations',numel(F)-start); %#ok<AGROW>
                if numel(F)==start,break;end
                if min(F)<previous-1e-5,radius=min(.35,1.3*radius);else,radius=max(.025,.6*radius);end
            end
            meta=struct('name','measured_direction_local_matlab_pso_v1','seed',2206, ...
                'description','Paired local sensitivity probes, measured improving combination, then bounded low-inertia MATLAB PSO with moving neighborhoods.', ...
                'active_dimensions',numel(p),'swarm_size',swarmSize,'probes',probes,'epochs',epochs, ...
                'stop_reason','evaluation_budget','evaluations',numel(F));
            if numel(F)<budget,meta.stop_reason='no_distinct_progress';end
            function s=key(z),s=sprintf('%.12g,',z);end
            function value=sample(z)
                z=min(1,max(0,z));z=EELS_TuningSearch.encode(EELS_TuningSearch.decode(z,p),p);q=key(z);
                if isKey(cache,q),value=cache(q);return;end
                if numel(F)>=epochEnd,error('EELS:DirectedEpoch','Local evaluation budget reached.');end
                value=evaluator(EELS_TuningSearch.decode(z,p),struct('method',phase,'radius',radius,'acquisition',NaN,'iteration',epoch));
                cache(q)=value;X(end+1,:)=z;F(end+1)=value;
            end
        end
    end
end
