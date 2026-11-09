classdef EELS_TuningSearch
    % Small-budget, toolbox-independent GP search in normalized log space.
    methods (Static)
        function z=encode(x,p)
            offset=([p.lower]==0).*[p.upper]*1e-4;
            lo=log([p.lower]+offset);span=log([p.upper]+offset)-lo;
            z=(log(x+offset)-lo)./span;
            linear=strcmp({p.name},'stage2_attention_weight_strength');
            if any(linear),z(:,linear)=(x(:,linear)-[p(linear).lower])./([p(linear).upper]-[p(linear).lower]);end
        end
        function x=decode(z,p)
            offset=([p.lower]==0).*[p.upper]*1e-4;
            lo=log([p.lower]+offset);span=log([p.upper]+offset)-lo;
            x=min([p.upper],max([p.lower],exp(lo+z.*span)-offset));
            linear=strcmp({p.name},'stage2_attention_weight_strength');
            if any(linear),x(:,linear)=[p(linear).lower]+z(:,linear).*([p(linear).upper]-[p(linear).lower]);end
            lower=repmat([p.lower],size(z,1),1);upper=repmat([p.upper],size(z,1),1);
            x(z<=0)=lower(z<=0);x(z>=1)=upper(z>=1);
            integer=EELS_ParameterCatalog.isInteger({p.name});x(:,integer)=round(x(:,integer));
        end
        function validateActive(p,s)
            names={p.name};
            for k=1:numel(names)
                reason=EELS_ParameterCatalog.inactive(names{k},s);
                if ~isempty(reason),error('EELS:TunerInactive','%s: %s',names{k},reason);end
            end
            modes={'scale','tilt','shift','width','tail','tail_asym'};
            for j=1:numel(modes)
                prefix=['stage2_zlp_refine_',modes{j}];
                if any(strcmp(names,[prefix,'_l2'])|strcmp(names,[prefix,'_bound'])|strcmp(names,[prefix,'_bound_meV']))&& ...
                        (~s.stage2_zlp_refine_enabled||~s.([prefix,'_enabled']))
                    error('EELS:TunerInactive','Enable %s before tuning its controls.',prefix);
                end
            end
            if any(startsWith(names,'stage2_attention_'))&&~s.stage2_attention_enabled
                error('EELS:TunerInactive','Attention parameters require attention enabled.');
            end
            if any(strcmp(names,'lambda_pm_l1'))&&exist('quadprog','file')~=2
                error('EELS:TunerInactive','PM L1 requires quadprog.');
            end
            groups={ ...
                {'zlp_fwhm_guard_factor','zlp_fwhm_core_factor'},'adaptive_zlp_region_enabled'; ...
                {'stage1_zlp_gain_asym_width_ratio','stage1_zlp_gain_asym_shape_ratio','stage1_zlp_gain_asym_shape_prior_weight'},'stage1_zlp_gain_asym_shape_prior_enabled'; ...
                {'stage1_zlp_gain_tail_prior_ratio'},'stage1_zlp_gain_tail_prior_enabled'; ...
                {'stage1_zlp_core_data_weight','stage1_zlp_core_anchor_weight','stage1_zlp_fartail_gain_boost'},'stage1_zlp_core_fartail_enabled'; ...
                {'stage1_zlp_gain_width_boost','stage1_zlp_loss_width_scale'},'stage1_zlp_gain_heavy_init_enabled'};
            for k=1:size(groups,1)
                if any(ismember(names,groups{k,1}))&&~s.(groups{k,2})
                    error('EELS:TunerInactive','Selected ZLP controls require %s=true.',groups{k,2});
                end
            end
            if any(ismember(names,groups{1,1}))&&~strcmpi(s.adaptive_zlp_region_method,'fwhm')
                error('EELS:TunerInactive','FWHM controls require adaptive_zlp_region_method=fwhm.');
            end
            if any(ismember(names,{'stage1_zlp_gain_asym_width_ratio','stage1_zlp_gain_asym_shape_ratio'}))&&s.stage1_zlp_gain_asym_shape_prior_weight<=0
                error('EELS:TunerInactive','Asymmetry target ratios need a nonzero shape-prior weight.');
            end
            if any(strcmp(names,'stage1_zlp_gain_tail_prior_ratio'))&&s.stage1_zlp_gain_tail_prior_weight<=0
                error('EELS:TunerInactive','Tail ratio needs a nonzero tail-prior weight.');
            end
        end

        function direction=residualDirections(fit,p)
            direction=ones(1,numel(p));
            e=fit.ene_corrected(:);r=fit.residual(:);low=12;
            if isfield(fit.zlp_region_policy,'fwhm_meV')&&isfinite(fit.zlp_region_policy.fwhm_meV)
                low=max(low,fit.zlp_region_policy.fwhm_meV);
            end
            high=min(100,3*low);left=e< -low&e> -high;right=e>low&e<high;
            if nnz(left)<3||nnz(right)<3,return;end
            % Only order probes; always test the opposite direction too.
            sign_left=sign(median(r(left))-median(r(right)));if sign_left==0,sign_left=1;end
            j=ismember({p.name},{'stage1_zlp_gain_asym_width_ratio','stage1_zlp_gain_asym_shape_ratio', ...
                'stage1_zlp_gain_tail_prior_ratio','stage1_zlp_gain_width_boost','stage1_zlp_fartail_gain_boost'});
            direction(j)=sign_left;
        end

        function [candidate,meta]=propose(h,p,initial,directions)
            d=numel(p);n=numel(h);
            meta=struct('method','baseline','radius',.12,'acquisition',NaN);
            if n==0,candidate=initial;return;end
            X=EELS_TuningSearch.encode(vertcat(h.values),p);y=[h.score]';
            if isfield(h,'search_score'),y=[h.search_score]';end
            [~,ib]=min(y);center=X(ib,:);
            if n<=2*d
                j=1+floor((n-1)/2);sgn=directions(j)*(1-2*mod(n-1,2));
                q=center;
                if mod(n,2)==0,q=EELS_TuningSearch.encode(h(end).incumbent_values,p);end
                q(j)=min(1,max(0,q(j)+sgn*.12));
                if min(sum((X-q).^2,2))>1e-16
                    candidate=EELS_TuningSearch.decode(q,p);meta.method='residual_ordered_probe';return;
                end
            end
            incumbent=cummin(y);gain=[Inf;incumbent(1:end-1)-incumbent(2:end)];
            last=find(gain>1e-4,1,'last');stalled=n-last;
            radius=max(.025,.18*.65^floor(stalled/4));meta.radius=radius;
            stream=RandStream('mt19937ar','Seed',2206+n);
            Q=repmat(center,400+2*d,1);
            for k=1:400
                [~,order]=sort(rand(stream,1,d));j=order(1:min(d,1+mod(k,3)));
                Q(k,j)=center(j)+(2*rand(stream,1,numel(j))-1)*radius;
            end
            for j=1:d,Q(400+2*j-1,j)=center(j)-radius;Q(400+2*j,j)=center(j)+radius;end
            Q=min(1,max(0,Q));
            if stalled>=8,Q=[Q;rand(stream,80,d)];end %#ok<AGROW>
            distance=EELS_TuningSearch.distance(Q,X);Q=Q(min(distance,[],2)>1e-14,:);
            if isempty(Q),error('EELS:TunerSearch','No distinct candidates in the search region.');end
            keep=unique([round(linspace(1,n,min(n,160))),ib]);X=X(keep,:);y=y(keep);
            finite=y(isfinite(y));penalty=max(finite)+max(1,std(finite));y(~isfinite(y))=penalty;
            mu=mean(y);scale=max(std(y),.1);z=(y-mu)/scale;
            K=exp(-EELS_TuningSearch.distance(X,X)/(2*.25^2));
            L=chol(K+1e-5*eye(size(K)),'lower');alpha=L'\(L\z);
            Ks=exp(-EELS_TuningSearch.distance(X,Q)/(2*.25^2));
            predicted=Ks'*alpha;v=L\Ks;sigma=sqrt(max(1-sum(v.^2,1)',1e-12));
            improvement=(min(finite)-mu)/scale-predicted-.01;
            u=improvement./sigma;ei=improvement.*(.5*(1+erf(u/sqrt(2))))+sigma.*exp(-u.^2/2)/sqrt(2*pi);
            [meta.acquisition,k]=max(ei);candidate=EELS_TuningSearch.decode(Q(k,:),p);meta.method='gp_expected_improvement';
        end

        function D=distance(X,Y)
            D=max(0,sum(X.^2,2)+sum(Y.^2,2)'-2*(X*Y'));
        end
    end
end
