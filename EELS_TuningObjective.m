classdef EELS_TuningObjective
    % Keep calibration loss, search feedback and strategy eligibility distinct.
    methods (Static)
        function d=aggregate(points,baselineNRMSE,maxRatio)
            ratio=[points.NRMSE]./max(baselineNRMSE,eps);
            failed=~isfinite(ratio)|~isfinite([points.valley_T_K])| ...
                ~ismember({points.status},{'ok','fit_regression'});
            excess=max(0,ratio/maxRatio-1);excess(failed)=0;
            regressed=~failed & ratio>maxRatio;
            loss=[points.score];score=mean(loss);
            violation=mean(excess)+.25*max(excess);
            feasible=~any(failed|regressed)&&isfinite(score);
            if any(failed)||~isfinite(score)
                score=Inf;search=1e6+nnz(failed);
                status='fit_failed';
            else
                % Pure temperature MAE; residual eligibility is a separate gate.
                search=score;
                if feasible,status='ok';else,status='fit_regression';end
            end
            d=struct('score',score,'search_score',search,'feasible',feasible, ...
                'regressed_count',nnz(regressed),'failed_count',nnz(failed), ...
                'max_fit_ratio',max(ratio),'constraint_violation',violation,'status',status);
        end

        function [accepted,regressed,failed]=masks(history)
            accepted=isfinite([history.score])&strcmp({history.status},'ok');
            regressed=false(size(accepted));failed=regressed;
            for k=1:numel(history)
                if isfield(history,'point_outcomes')&&~isempty(history(k).point_outcomes)
                    states={history(k).point_outcomes.status};
                    regressed(k)=any(strcmp(states,'fit_regression'));
                    failed(k)=any(~ismember(states,{'ok','fit_regression'}));
                else
                    regressed(k)=strcmp(history(k).status,'fit_regression');
                    failed(k)=~accepted(k)&&~regressed(k);
                end
                if ~accepted(k)&&~regressed(k)&&~failed(k),failed(k)=true;end
            end
        end
    end
end
