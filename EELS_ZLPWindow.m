classdef EELS_ZLPWindow
    methods (Static)
        function [s,p]=apply(e,y,s)
            p=struct('enabled',true,'method','fwhm_fallback_static', ...
                'peak_center_meV',NaN,'fwhm_meV',NaN,'core_half_width_meV',NaN, ...
                'guard_start_meV',s.zlp_signal_guard_start_meV, ...
                'guard_end_meV',s.zlp_signal_guard_end_meV, ...
                'core_downweight_width_meV',s.core_downweight_width,'reason','unresolved_half_height');
            e=double(e(:));y=double(y(:));
            if numel(e)~=numel(y)||numel(e)<7||any(~isfinite(e))||any(diff(e)<=0),return;end
            good=isfinite(y);search=find(good&abs(e)<=s.adaptive_zlp_region_search_meV);
            if numel(search)<5,return;end
            % Three-point median suppresses isolated spikes without a broad
            % sample-count smoothing kernel that inflates narrow ZLPs.
            ys=movmedian(y,3,'omitnan');[peak,k]=max(ys(search));ip=search(k);
            floor_y=min(ys(search));level=floor_y+0.5*(peak-floor_y);
            if ~isfinite(level)||peak<=floor_y,return;end
            left=ip;right=ip;
            while left>search(1)&&isfinite(ys(left-1))&&ys(left)>level,left=left-1;end
            while right<search(end)&&isfinite(ys(right+1))&&ys(right)>level,right=right+1;end
            if left==ip||right==ip||ys(left)>level||ys(right)>level,return;end
            xl=e(left)+(level-ys(left))*(e(left+1)-e(left))/(ys(left+1)-ys(left));
            xr=e(right-1)+(level-ys(right-1))*(e(right)-e(right-1))/(ys(right)-ys(right-1));
            width=xr-xl;
            if ~isfinite(width)||width<=0,return;end
            half=max(e(ip)-xl,xr-e(ip));
            lower=max(s.zlp_fwhm_min_halfwidth_meV,median(diff(e)));
            upper=s.zlp_fwhm_max_halfwidth_meV;
            if ~isfinite(lower)||~isfinite(upper)||upper<lower|| ...
                    s.zlp_fwhm_guard_factor<=0||s.zlp_fwhm_core_factor<=0
                error('EELS:ZLPBounds','Invalid FWHM window factors or bounds.');
            end
            guard=min(max(s.zlp_fwhm_guard_factor*half,lower),upper);
            core=min(max(s.zlp_fwhm_core_factor*half,lower),upper);
            s.zlp_signal_guard_start_meV=guard;
            s.zlp_lowfreq_guard_start_meV=guard;
            s.stage1_zlp_core_anchor_width_meV=core;
            if s.adaptive_zlp_region_adjust_weights,s.core_downweight_width=core;end
            p.method='fwhm';p.reason='';p.peak_center_meV=e(ip);p.fwhm_meV=width;
            p.half_height_crossings_meV=[xl xr];p.core_half_width_meV=half;
            p.guard_start_meV=guard;p.core_downweight_width_meV=s.core_downweight_width;
        end
    end
end
