classdef EELS_FitWeights
    % Stage-specific, bounded relative data weights. Neutral controls preserve
    % the previous objective; the weighted mean multiplier always equals one.
    methods (Static)
        function p=parameters()
            names={'stage1_shoulder_weight_ratio','stage1_vibration_weight_ratio', ...
                'stage1_tail_weight_ratio','stage1_dynamic_weight_scale', ...
                'stage2_gain_weight_ratio','stage2_lowfreq_weight_ratio', ...
                'stage2_highfreq_weight_ratio','stage2_tail_weight_ratio', ...
                'stage2_dynamic_weight_scale','stage2_attention_weight_strength'};
            lo=[.25 .25 .25 .25 .5 .25 .25 .25 .25 0];
            hi=[4 4 3 2.5 2.5 4 4 3 2.5 1];
            p=struct('name',names,'lower',num2cell(lo),'upper',num2cell(hi), ...
                'enabled',num2cell(true(1,numel(names))));
        end
        function s=defaults(s)
            p=EELS_FitWeights.parameters();
            for k=1:numel(p)
                if ~isfield(s,p(k).name),s.(p(k).name)=double(k<numel(p));end
            end
        end
        function multiplier=balance(base,profile)
            base=base(:);profile=profile(:);
            if any(~isfinite(profile)|profile<=0)||any(~isfinite(base)|base<=0)
                error('EELS:FitWeights','Data weights must be positive and finite.');
            end
            if all(profile==1),multiplier=ones(size(profile));return;end
            r=profile/(sum(base.*profile)/sum(base));
            % Contract toward one until both limits hold, preserving the mean.
            alpha=1;
            if min(r)<.25,alpha=min(alpha,.75/(1-min(r)));end
            if max(r)>4,alpha=min(alpha,3/(max(r)-1));end
            multiplier=1+alpha*(r-1);
        end
        function profile=regional(e,core,high,tail,s,stage)
            shoulder=(1-core).*(1-high).*(1-tail);
            vibration=(1-core).*high.*(1-tail);far=(1-core).*tail;
            if stage==1
                profile=exp(shoulder*log(s.stage1_shoulder_weight_ratio)+ ...
                    vibration*log(s.stage1_vibration_weight_ratio)+far*log(s.stage1_tail_weight_ratio));
            else
                gain=(1-core)./(1+exp(e(:)/4));
                profile=exp(shoulder*log(s.stage2_lowfreq_weight_ratio)+ ...
                    vibration*log(s.stage2_highfreq_weight_ratio)+far*log(s.stage2_tail_weight_ratio)+ ...
                    gain*log(s.stage2_gain_weight_ratio));
            end
        end
        function tf=stage2Active(s)
            p=EELS_FitWeights.parameters();tf=false;
            for k=5:9,tf=tf||s.(p(k).name)~=1;end
            tf=tf||s.stage2_attention_weight_strength~=0;
        end
        function [weights,ratio]=stage2(reference,regional,fit_alpha,feature,s)
            floor=max(.05,min(.95,s.stage2_fit_sigmoid_min_weight));
            feature_floor=max(.05,min(.95,s.stage2_attention_floor));
            target=(floor+(1-floor)*fit_alpha(:)).*(feature_floor+(1-feature_floor)*feature(:));
            % This smooth branch bypasses the legacy joint_floor=1 clamp.
            relative=regional(:).*max(target,.0025).^s.stage2_attention_weight_strength;
            ratio=EELS_FitWeights.balance(reference.^2,relative);
            weights=reference(:).*sqrt(ratio);
        end
        function w=curves(result)
            w=struct('energy_meV',result.ene_corrected(:),'stage1',[],'stage2',[]);
            if isfield(result,'stage1')&&isfield(result.stage1,'data_weights')
                w.stage1=result.stage1.data_weights(:);
            end
            if isfield(result,'stage2')&&isfield(result.stage2,'attention')&& ...
                    isfield(result.stage2.attention,'data_weight_policy')
                w.stage2=result.stage2.attention.data_weight_policy.effective_weights(:);
            end
        end
    end
end
