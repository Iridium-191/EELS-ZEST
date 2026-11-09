classdef EELS_Native
    % Compiled numerical kernels; MATLAB owns policy, diagnostics and the API.
    methods (Static)
        function tf = enabled(settings)
            mode = 'auto';
            if isfield(settings,'native_backend'), mode = lower(char(settings.native_backend)); end
            if ~ismember(mode, {'auto','cpp','matlab'})
                error('EELS:NativeBackend','native_backend must be auto, cpp, or matlab.');
            end
            tf = ~strcmp(mode,'matlab') && exist('eels_native_v22_mex','file') == 3;
            if strcmp(mode,'cpp') && ~tf
                error('EELS:NativeMissing','Run build_eels_native before selecting native_backend=''cpp''.');
            end
        end

        function tf = supports_stage1(settings)
            % auto uses the compiled objective when available; cpp additionally
            % requires an installed MEX. matlab keeps the reference trajectory.
            tf = EELS_Native.enabled(settings) && ~settings.tail_overshoot_constraint && ...
                ismember(lower(char(settings.inner_solver)),{'auto','lsqlin','lsqnonneg'});
        end

        function c = stage1_context(ene, y, EM, settings)
            % Invariants are prepared once per nonlinear fit, not per theta.
            [s,~] = EELS_PearsonSparse_FullFit.adaptive_zlp_region_settings(ene,y,settings);
            n = numel(ene); nz = s.num_zlp_components; ns = numel(EM); p = nz+2*ns;
            layout = struct('zlp',1:nz,'pm_s',nz+(1:ns), ...
                'pm_a',nz+ns+(1:ns),'bg',[]);
            static = s; static.lambda_be_soft = 0;
            [R,rhs] = EELS_PearsonSparse_FullFit.regularization_rows(EM,layout,450,p,static);
            guards = static;
            guards.lambda_pm_smooth = 0; guards.lambda_pm_l2 = 0;
            guards.zlp_wide_component_l2 = 0;
            glayout = struct('zlp',1:n,'pm_s',[],'pm_a',[],'bg',[]);
            % These rows are linear in the ZLP shapes. An identity input
            % extracts the exact operator, including paired gain-tail rows.
            [P,prhs] = EELS_PearsonSparse_FullFit.regularization_rows( ...
                [],glayout,450,n,guards,speye(n),ene,y);
            kernel = 0;
            if strcmpi(s.kernel_mode,'full_zlp'), kernel=1; end
            if strcmpi(s.kernel_mode,'fitted_zlp'), kernel=2; end
            shape_weight = s.stage1_zlp_gain_asym_shape_prior_weight;
            if ~s.stage1_zlp_gain_asym_shape_prior_enabled, shape_weight=0; end
            [lb,ub] = EELS_PearsonSparse_FullFit.theta_bounds(s);
            [prior,center_scale] = EELS_Native.hierarchy_prior(s);
            c = struct('energy',double(ene(:)),'target',double(y(:)), ...
                'sqrtW',EELS_PearsonSparse_FullFit.build_weights(ene,y,s), ...
                'EM',double(EM(:)),'regular',full(R),'rhs',rhs, ...
                'guard',sparse(P),'guard_rhs',prhs,'lower',lb,'upper',ub, ...
                'prior',prior,'options',double([nz, ...
                EELS_PearsonSparse_FullFit.is_pvoigt_model(s.zlp_model_type), ...
                kernel,s.lambda_be_soft,shape_weight,s.stage1_zlp_gain_asym_prior_start_index, ...
                s.stage1_zlp_gain_asym_width_ratio,s.stage1_zlp_gain_asym_shape_ratio, ...
                center_scale,max(s.inner_solver_max_iter,10*p)]));
        end

        function [prior,center_scale] = hierarchy_prior(s)
            prior=zeros(0,4);center_scale=1;nphi=5*s.num_zlp_components;
            ctx=s.stage1_hierarchy_context;
            if ~isstruct(ctx)||~isfield(ctx,'enabled')||~ctx.enabled,return;end
            if isfield(ctx,'center_scale_meV'),center_scale=max(abs(double(ctx.center_scale_meV)),eps);end
            data_cost=1;
            if isfield(ctx,'data_cost')&&isfinite(ctx.data_cost)&&ctx.data_cost>0
                data_cost=max(double(ctx.data_cost),1e-10);
            end
            if isfield(ctx,'template_phi')&&numel(ctx.template_phi)==nphi
                scale=EELS_PearsonSparse_FullFit.hierarchy_scale_vector(ctx,'template_scale',nphi);
                membership=1;
                if isfield(ctx,'membership')&&isfinite(ctx.membership),membership=min(max(double(ctx.membership),0),1);end
                frac=s.stage1_cluster_template_cost_fraction*membership;
                if isfield(ctx,'template_cost_fraction')&&isfinite(ctx.template_cost_fraction)
                    frac=max(0,double(ctx.template_cost_fraction))*membership;
                end
                if frac>0
                    prior=[double(ctx.template_phi(:)),scale(:), ...
                        repmat(sqrt(frac*data_cost/nphi),nphi,1),zeros(nphi,1)];
                end
            end
            if s.stage1_spatial_zlp_enabled&&isfield(ctx,'neighbor_phi')&&size(ctx.neighbor_phi,2)==nphi
                phi=double(ctx.neighbor_phi);weights=ones(size(phi,1),1);
                if isfield(ctx,'neighbor_weight')&&numel(ctx.neighbor_weight)==size(phi,1)
                    weights=max(double(ctx.neighbor_weight(:)),0);
                end
                keep=all(isfinite(phi),2)&isfinite(weights)&weights>0;phi=phi(keep,:);weights=weights(keep);
                scale=EELS_PearsonSparse_FullFit.hierarchy_scale_vector(ctx,'spatial_scale',nphi);
                denom=max(nphi*sum(weights),eps);
                for j=1:size(phi,1)
                    factor=sqrt(max(s.stage1_spatial_zlp_weight,0)*data_cost*weights(j)/denom);
                    prior=[prior;phi(j,:)',scale(:),repmat(factor,nphi,1), ...
                        repmat(max(s.stage1_spatial_zlp_huber_delta,eps),nphi,1)]; %#ok<AGROW>
                end
            end
        end

        function [r,J] = objective(theta,c)
            if nargout>1
                [r,J]=eels_native_v22_mex('stage1',double(theta),c);
            else
                r=eels_native_v22_mex('stage1',double(theta),c);
            end
        end

        function [coef,status] = stage2_coefficients(T,basis,hs,ha,ss,sa,zlp,y,w,settings)
            s=settings;s.stage2_soft_be_weight=0;
            [R,rhs]=EELS_PearsonSparse_FullFit.stage2_regularization_rows(basis,0,s,[],[]);
            nh=size(hs,2);ns=size(ss,2);base=nh+2*ns;nz=size(zlp.B,2);
            lo=zeros(base,1);hi=Inf(base,1);
            if nz>0
                R=[R,zeros(size(R,1),nz)];
                [rz,dz]=EELS_PearsonSparse_FullFit.stage2_zlp_refine_regularization(zlp,base);
                R=[R;rz];rhs=[rhs;dz];lo=[lo;zlp.lb];hi=[hi;zlp.ub];
            end
            c=struct('basis',[hs,ha,ss,sa,zlp.B],'regular',full(R), ...
                'target',y(:),'sqrtW',w(:),'rhs',rhs,'EM_hard',double(basis.EM_hard(:)), ...
                'EM_soft',double(basis.EM_soft(:)),'temperatures',double(T(:)), ...
                'lower',lo,'upper',hi,'soft_be',double(settings.stage2_soft_be_weight), ...
                'max_iter',double(max(settings.inner_solver_max_iter,10*numel(lo))));
            [coef,status]=eels_native_v22_mex('stage2',c);
        end
    end
end
