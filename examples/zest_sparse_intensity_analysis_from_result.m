% Sparse phonon-mode source-intensity analysis from a fitted datacube result.
%
% Required variable:
%   cuberes : output from eels_fit_datacube or promoted Stage 2 result
%
% This script creates a signed source-energy axis where gain/anti-Stokes is
% negative and loss/Stokes is positive, then exports a DOS-like intensity cube.
%
axis_meV = -250:1:250;
[ny, nx] = size(cuberes.R2_map);
pm_cube = NaN(ny, nx, numel(axis_meV));
valid_mask = false(ny, nx);

for ii = 1:numel(cuberes.results)
    res = cuberes.results{ii};
    if isempty(res) || ~isstruct(res) || isfield(res, 'error')
        continue;
    end
    if isfield(res, 'stage2') && isfield(res.stage2, 'EM_final')
        em = res.stage2.EM_final(:);
        ps = res.stage2.PM_stokes(:);
        pa = res.stage2.PM_antistokes(:);
    else
        em = res.EM_final(:);
        ps = res.PM_stokes(:);
        pa = res.PM_antistokes(:);
    end
    signed_em = [-flipud(em); em];
    signed_pm = [flipud(pa); ps];
    pix = cuberes.sample_pixels(ii);
    [yy, xx] = ind2sub([ny, nx], pix);
    pm_cube(yy, xx, :) = interp1(signed_em, signed_pm, axis_meV, 'linear', 0);
    valid_mask(yy, xx) = true;
end

mean_PM = squeeze(mean(mean(pm_cube, 1, 'omitnan'), 2, 'omitnan'));
eels_zest_pm_cube = struct('axis_meV', axis_meV, 'pm_cube', pm_cube, ...
    'mean_PM', mean_PM, 'valid_mask', valid_mask);

figure;
plot(axis_meV, mean_PM, 'k-', 'LineWidth', 1.5);
grid on;
xlabel('Signed EM (meV)');
ylabel('Mean sparse phonon source intensity');
title('EELS-ZEST mean signed sparse phonon source intensity');
