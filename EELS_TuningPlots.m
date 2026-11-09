classdef EELS_TuningPlots
    methods (Static)
        function render(ax,history,baseline,best,truth,p)
            if isempty(history)||isempty(best),return;end
            n=numel(history);loss=[history.score];ok=isfinite(loss);
            [accepted,regressed,failed]=EELS_TuningObjective.masks(history);
            inc=loss;inc(~accepted)=Inf;
            for k=2:n,inc(k)=min(inc(k-1),inc(k));end
            cla(ax.loss);trace=loss;trace(~ok)=NaN;
            plot(ax.loss,1:n,trace,'-','Color',[.65 .7 .8],'HandleVisibility','off');hold(ax.loss,'on');
            plot(ax.loss,find(accepted),loss(accepted),'bo','DisplayName','Accepted');
            plot(ax.loss,find(ok&~accepted),loss(ok&~accepted),'rx','MarkerSize',8,'DisplayName','NRMSE gate exceeded');
            missing=~ok;top=1;
            if any(ok),top=max(loss(ok))+max(.1,.12*(max(loss(ok))-min(loss(ok))));end
            if any(missing)
                plot(ax.loss,find(missing),top*ones(1,nnz(missing)),'kv','DisplayName','No finite loss (marker only)');
            end
            plot(ax.loss,1:n,inc,'-','Color',[.85 .2 .1],'LineWidth',1.7,'DisplayName','Best accepted');hold(ax.loss,'off');
            title(ax.loss,sprintf('%d candidates: %d accepted / %d regressed / %d failed',n,nnz(accepted),nnz(regressed),nnz(failed)));
            xlabel(ax.loss,'Candidate evaluation');ylabel(ax.loss,'Calibration loss (before acceptance gate)');
            if isfield(best,'loss_units')&&strcmp(best.loss_units,'K'),ylabel(ax.loss,'Temperature MAE (K)');end
            legend(ax.loss,'show','Location','best');grid(ax.loss,'on');
            cla(ax.zlp);a=baseline.result;b=best.result;
            plot(ax.zlp,a.ene_corrected,a.raw_data/max(a.norm_val,eps),'Color',[.65 .65 .65]);hold(ax.zlp,'on');
            plot(ax.zlp,a.ene_corrected,a.zlp_only/max(a.norm_val,eps),'b--','LineWidth',1.2);
            plot(ax.zlp,b.ene_corrected,b.zlp_only/max(b.norm_val,eps),'r-','LineWidth',1.2);hold(ax.zlp,'off');
            xlim(ax.zlp,[-65 65]);title(ax.zlp,'ZLP: before / best');xlabel(ax.zlp,'Corrected energy (meV)');ylabel(ax.zlp,'Normalized intensity');
            legend(ax.zlp,{'Reference','Before','Best'},'Location','best');grid(ax.zlp,'on');
            cla(ax.valley);
            a=baseline.valley;b=best.valley;
            plot(ax.valley,a.T_grid,EELS_TuningPlots.scaleCost(a.cost),'b--','LineWidth',1.2);hold(ax.valley,'on');
            plot(ax.valley,b.T_grid,EELS_TuningPlots.scaleCost(b.cost),'r-','LineWidth',1.2);
            xline(ax.valley,truth,'k:');hold(ax.valley,'off');
            tmin=max(min([a.T_grid(:);b.T_grid(:)]),min([truth,a.T_K,b.T_K])-120);
            tmax=min(max([a.T_grid(:);b.T_grid(:)]),max([truth,a.T_K,b.T_K])+160);
            if tmax>tmin,xlim(ax.valley,[tmin tmax]);end
            title(ax.valley,sprintf('Raw valley %.2f -> %.2f K',a.T_K,b.T_K));
            xlabel(ax.valley,'Temperature (K)');ylabel(ax.valley,'Cost above minimum / span');
            legend(ax.valley,{'Before','Best','Known T'},'Location','best');grid(ax.valley,'on');
            cla(ax.residual);a=baseline.result;b=best.result;
            plot(ax.residual,a.ene_corrected,a.residual/max(a.norm_val,eps),'b--');hold(ax.residual,'on');
            plot(ax.residual,b.ene_corrected,b.residual/max(b.norm_val,eps),'r-');yline(ax.residual,0,'k:');hold(ax.residual,'off');
            title(ax.residual,'Residual: before / best');xlabel(ax.residual,'Corrected energy (meV)');ylabel(ax.residual,'Residual / peak');grid(ax.residual,'on');
            legend(ax.residual,{'Before','Best'},'Location','best');
            if isfield(best,'point_outcomes')&&numel(best.point_outcomes)>1
                title(ax.zlp,'Training point 1: ZLP before / best');
                title(ax.valley,sprintf('Training point 1 valley %.2f -> %.2f K',baseline.valley.T_K,best.valley.T_K));
                title(ax.residual,'Training point 1: residual before / best');
            end
            if isfield(ax,'projection')
                [xy,m]=EELS_TuningDataset.embedding(history,p);cla(ax.projection);
                colors=loss(ok);scatter(ax.projection,xy(ok,1),xy(ok,2),32,colors(:),'filled');hold(ax.projection,'on');
                if any(~ok),scatter(ax.projection,xy(~ok,1),xy(~ok,2),30,'kx');end
                if any(regressed&ok),scatter(ax.projection,xy(regressed&ok,1),xy(regressed&ok,2),55,'rx');end
                eligible=loss;eligible(~accepted)=Inf;[~,ib]=min(eligible);
                if any(accepted),plot(ax.projection,xy(ib,1),xy(ib,2),'rp','MarkerSize',13,'LineWidth',1.5);end
                hold(ax.projection,'off');
                title(ax.projection,'PCA: color=loss; x=rejected; star=best accepted');
                xlabel(ax.projection,sprintf('PC1 %.1f%%',m.explained_percent(1)));ylabel(ax.projection,sprintf('PC2 %.1f%%',m.explained_percent(2)));
                grid(ax.projection,'on');colorbar(ax.projection);
            end
            if isfield(ax,'parameters')
                z=EELS_TuningSearch.encode(vertcat(history.values),p);
                imagesc(ax.parameters,1:n,1:numel(p),z');clim(ax.parameters,[0 1]);
                ax.parameters.YTick=1:numel(p);ax.parameters.YTickLabel=strrep({p.name},'_',' ');ax.parameters.TickLabelInterpreter='none';
                title(ax.parameters,'Parameter path (bounds normalized)');xlabel(ax.parameters,'Candidate evaluation');colorbar(ax.parameters);
            end
        end
        function c=scaleCost(c)
            c=(c-min(c))/max(max(c)-min(c),eps);
        end
        function constraints(ax,history,limit)
            n=numel(history);m=numel(history(1).point_outcomes);ratios=NaN(m,n);
            base=max([history(1).point_outcomes.NRMSE]',eps);
            for k=1:n,ratios(:,k)=[history(k).point_outcomes.NRMSE]'./base;end
            cla(ax);imagesc(ax,1:n,1:m,ratios);hold(ax,'on');
            [point,candidate]=find(ratios>limit);plot(ax,candidate,point,'kx','MarkerSize',5);
            [point,candidate]=find(~isfinite(ratios));plot(ax,candidate,point,'rs','MarkerSize',6);hold(ax,'off');
            title(ax,sprintf('NRMSE / baseline: x > %.2f; red square = missing',limit));
            xlabel(ax,'Candidate evaluation');ylabel(ax,'Training point index');
            cb=colorbar(ax);cb.Label.String='Unweighted NRMSE ratio';
        end
        function weightProfiles(ax1,ax2,before,after)
            a=EELS_FitWeights.curves(before.result);b=EELS_FitWeights.curves(after.result);
            for stage=1:2
                if stage==1,ax=ax1;else,ax=ax2;end
                key=sprintf('stage%d',stage);cla(ax);hold(ax,'on');
                if ~isempty(a.(key)),plot(ax,a.energy_meV,a.(key),'b--','LineWidth',1.2);end
                if ~isempty(b.(key)),plot(ax,b.energy_meV,b.(key),'r-','LineWidth',1.2);end
                hold(ax,'off');ax.YScale='log';grid(ax,'on');
                title(ax,sprintf('Training point 1: Stage%d effective weights',stage));
                xlabel(ax,'Corrected energy (meV)');ylabel(ax,'Squared-residual weight w(E)');
                if ~isempty(a.(key))&&~isempty(b.(key)),legend(ax,{'Before','Best'},'Location','best');end
            end
        end
        function export(report,directory)
            f=figure('Visible','off','Position',[30 30 1500 1050]);cleaner=onCleanup(@()close(f)); %#ok<NASGU>
            tiledlayout(f,3,2,'TileSpacing','compact');
            ax=struct('loss',nexttile,'zlp',nexttile,'valley',nexttile,'residual',nexttile,'projection',nexttile,'parameters',nexttile);
            EELS_TuningPlots.render(ax,report.history,report.baseline,report.best,report.known_T_K,report.parameters);
            exportgraphics(f,fullfile(directory,'tuning_dashboard.png'),'Resolution',140);
            clf(f);a=axes(f);limit=report.objective_settings.max_fit_ratio;
            EELS_TuningPlots.constraints(a,report.history,limit);
            exportgraphics(f,fullfile(directory,'fit_constraints.png'),'Resolution',140);
            if isfield(report,'training_spectra')&&size(report.training_spectra,2)>1
                clf(f);a=axes(f);n=numel(report.best.point_outcomes);
                before=[report.baseline.point_outcomes.valley_T_K];after=[report.best.point_outcomes.valley_T_K];
                if isfield(report.validation,'before_K')
                    before=[before report.validation.before_K];after=[after report.validation.after_K];
                end
                plot(a,1:numel(before),before,'bo--',1:numel(after),after,'ro-');hold(a,'on');
                yline(a,report.known_T_K,'k:');xline(a,n+.5,'k--');
                title(a,'Fixed training points | independent heldout points');xlabel(a,'Point index');ylabel(a,'Raw valley temperature (K)');
                legend(a,{'Before','Best','Known T','Heldout boundary'});grid(a,'on');
                exportgraphics(f,fullfile(directory,'joint_validation.png'),'Resolution',140);
            end
            clf(f);ax1=subplot(2,1,1,'Parent',f);ax2=subplot(2,1,2,'Parent',f);
            EELS_TuningPlots.weightProfiles(ax1,ax2,report.baseline,report.best);
            exportgraphics(f,fullfile(directory,'fit_weights.png'),'Resolution',140);
            rows=table();entries={report.baseline,report.best};labels={'before','best'};
            for k=1:2
                w=EELS_FitWeights.curves(entries{k}.result);
                for stage=1:2
                    v=w.(sprintf('stage%d',stage));if isempty(v),continue;end
                    t=table(repmat(string(labels{k}),numel(v),1),repmat(stage,numel(v),1),w.energy_meV,v, ...
                        'VariableNames',{'strategy','stage','energy_meV','effective_weight'});rows=[rows;t]; %#ok<AGROW>
                end
            end
            if ~isempty(rows),writetable(rows,fullfile(directory,'fit_weights.csv'));end
        end
    end
end
