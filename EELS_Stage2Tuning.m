classdef EELS_Stage2Tuning
    methods (Static)
        function p=parameters()
            modes={'scale','tilt','shift','width','tail','tail_asym'};
            p=struct('name',{},'lower',{},'upper',{},'enabled',{});
            low=[.01 .01 .1 .01 .01 .01];high=[.3 .4 4 .4 .5 .5];
            for j=1:6
                suffix='_bound';if j==3,suffix='_bound_meV';end
                p(end+1)=struct('name',['stage2_zlp_refine_',modes{j},suffix], ...
                    'lower',low(j),'upper',high(j),'enabled',true); %#ok<AGROW>
                p(end+1)=struct('name',['stage2_zlp_refine_',modes{j},'_l2'], ...
                    'lower',.005,'upper',5,'enabled',true); %#ok<AGROW>
            end
        end
        function t=profile(result)
            q=result.stage2.temperature_profile;
            keys={'cost_raw','joint_data_mse','attention_data_mse','regularization_cost','inner_objective', ...
                'inner_data_sum','inner_regularization_sum','inner_l1_cost','ZLP_delta_rms'};
            t=table(q.T_grid(:),'VariableNames',{'temperature_K'});
            for j=1:numel(keys),if isfield(q,keys{j}),t.(keys{j})=q.(keys{j})(:);end;end
        end
        function plot(ax,fit)
            t=EELS_Stage2Tuning.profile(fit);cla(ax);hold(ax,'on');
            keys={'cost_raw','joint_data_mse','regularization_cost','inner_objective'};
            for j=1:numel(keys)
                if ~ismember(keys{j},t.Properties.VariableNames),continue;end
                c=t.(keys{j});[v,k]=min(c);span=max(c)-v;
                label=sprintf('%s: min %.2f K',keys{j},t.temperature_K(k));
                if span<=eps*max(1,max(abs(c))),label=[keys{j},': flat'];end
                plot(ax,t.temperature_K,(c-v)/max(span,realmin),'DisplayName',label);
            end
            hold(ax,'off');legend(ax,'show','Interpreter','none','Location','best');grid(ax,'on');
            title(ax,'Training point 1: minima of separate cost components');
            xlabel(ax,'Temperature (K)');ylabel(ax,'Each curve normalized separately (not relative strength)');
        end
    end
end
