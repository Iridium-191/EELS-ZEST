function plot_example_figures()
%PLOT_EXAMPLE_FIGURES Publication-style SVG figures for the example data.

    script_dir = fileparts(mfilename('fullpath'));
    data_file = fullfile(script_dir, 'example_figure_data.mat');
    if ~exist(data_file, 'file')
        error('Missing %s. Run export_example_figure_data first.', data_file);
    end
    S = load(data_file);

    out_dir = script_dir;
    set(groot, 'defaultAxesFontName', 'Arial');
    set(groot, 'defaultTextFontName', 'Arial');
    set(groot, 'defaultAxesFontSize', 8);

    fig1 = local_figure1(S);
    local_export_figure(fig1, out_dir, 'Figure1_ZLP_sparse_reconstruction');
    close(fig1);

    fig2 = local_figure2(S);
    local_export_figure(fig2, out_dir, 'Figure2_attention_BE_thermometry');
    close(fig2);

    local_export_panel_csv(S, out_dir);
    fprintf('Wrote SVG/PNG figures and CSV panel data to %s\n', out_dir);
end

function fig = local_figure1(S)
    ene = S.ene_fit(:);
    raw = S.raw_original(:);
    den = S.denoised_used(:);
    model = S.model(:);
    zlp = S.zlp(:);
    vib = S.vib(:);
    bg = S.bg(:);
    EM = S.EM(:);
    PMs = S.PM_stokes(:);
    PMa = S.PM_antistokes(:);

    c_raw = [207 224 240] ./ 255;   % #CEDFEF
    c_den = [188 206 222] ./ 255;   % softened cool blue-gray
    c_fit = [132 129 186] ./ 255;   % #8481BA
    c_zlp = [116 178 103] ./ 255;   % softened green
    c_vib = [245 189 129] ./ 255;   % softened orange
    c_gain = [146 194 221] ./ 255;  % #92C2DD
    c_dos = [145 184 133] ./ 255;   % muted sage green
    c_grid = [239 243 247] ./ 255;

    fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 18.0 7.8], 'Renderer', 'painters');
    ax1 = axes(fig, 'Position', [0.07 0.16 0.39 0.73]);
    mask = abs(ene) <= 145;
    h_raw = plot(ax1, ene(mask), raw(mask), 'o', 'Color', c_raw, 'MarkerSize', 3.6, ...
        'LineWidth', 0.45, 'MarkerFaceColor', 'none'); hold(ax1, 'on');
    h_den = plot(ax1, ene(mask), den(mask), '-', 'Color', c_den, 'LineWidth', 0.55);
    h_fit = plot(ax1, ene(mask), model(mask), '-', 'Color', c_fit, 'LineWidth', 1.25);
    h_zlp = plot(ax1, ene(mask), zlp(mask), '-', 'Color', c_zlp, 'LineWidth', 1.10);
    h_vib = plot(ax1, ene(mask), vib(mask), '-', 'Color', c_vib, 'LineWidth', 1.10);
    yline(ax1, 0, '-', 'Color', [0.50 0.50 0.48], 'LineWidth', 0.45);
    signal_mask = abs(ene) >= 8 & abs(ene) <= 145;
    signal_ymax = max([vib(signal_mask); den(signal_mask) - zlp(signal_mask) - bg(signal_mask); eps], [], 'omitnan');
    y_upper = max(1.1 .* signal_ymax, 1);
    local_plot_pm_sticks(ax1, EM, PMs, PMa, c_vib, c_gain, 0, ...
        0.26 .* y_upper, 145, ene, max(vib, 0));
    grid(ax1, 'on'); box(ax1, 'on');
    ax1.GridColor = c_grid;
    ax1.GridAlpha = 0.45;
    ax1.LineWidth = 0.7;
    ax1.TickDir = 'out';
    xlabel(ax1, 'Energy loss (meV)');
    ylabel(ax1, 'Intensity (counts)');
    xlim(ax1, [-145 145]);
    ylim(ax1, [-0.05 .* y_upper, y_upper]);
    legend(ax1, [h_raw, h_den, h_fit, h_zlp, h_vib], ...
        {'raw', 'SVD', 'fit', 'ZLP', 'phonon'}, ...
        'Location', 'northwest', 'Box', 'off', 'FontSize', 6.3);
    text(ax1, 0.985, 0.94, sprintf('T = %.0f K', S.T_fit), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'FontSize', 6.8, 'Color', [0.22 0.22 0.22], ...
        'BackgroundColor', 'w', 'Margin', 1);
    text(ax1, 0.02, 0.98, 'a', 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', 10, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'w', 'Margin', 0.2);

    pos1 = ax1.Position;
    axin = axes(fig, 'Position', [pos1(1) + 0.72 * pos1(3), pos1(2) + 0.52 * pos1(4), ...
        0.23 * pos1(3), 0.30 * pos1(4)]);
    inset_mask = abs(ene) <= 200;
    semilogy(axin, ene(inset_mask), local_positive_floor(zlp(inset_mask)), ...
        '-', 'Color', c_zlp, 'LineWidth', 0.85); hold(axin, 'on');
    semilogy(axin, ene(inset_mask), local_positive_floor(abs(vib(inset_mask))), ...
        '-', 'Color', c_vib, 'LineWidth', 0.85);
    sig = den(:) - zlp(:) - bg(:);
    semilogy(axin, ene(inset_mask), local_positive_floor(abs(sig(inset_mask))), ...
        '-', 'Color', [0.15 0.15 0.15], 'LineWidth', 0.60);
    xlim(axin, [-160 160]); grid(axin, 'on'); box(axin, 'on');
    axin.LineWidth = 0.7;
    set(axin, 'FontSize', 5.8, 'TickDir', 'out');

    ax2 = axes(fig, 'Position', [0.54 0.16 0.31 0.73]);
    xgrid = linspace(0, 118, 900).';
    lr_curve = local_lr_deconvolution_curve(S, xgrid);
    lr_sta_curve = local_sta_envelope_removed(lr_curve, xgrid);
    gan_wz = local_interp_and_norm(S.energy_gan_wurtzite(:), S.dos_gan_wurtzite(:), xgrid);
    lr_sta_curve = max(local_norm01(lr_sta_curve, xgrid >= 10 & xgrid <= 112), 0);

    plot(ax2, xgrid, lr_sta_curve, '-', 'Color', c_fit, 'LineWidth', 1.35); hold(ax2, 'on');
    plot(ax2, xgrid, gan_wz, '-', 'Color', c_dos, 'LineWidth', 1.65);
    grid(ax2, 'on'); box(ax2, 'on');
    ax2.GridColor = c_grid;
    ax2.GridAlpha = 0.45;
    ax2.LineWidth = 0.7;
    ax2.TickDir = 'out';
    xlabel(ax2, 'Phonon energy (meV)');
    ylabel(ax2, 'Normalized intensity');
    xlim(ax2, [0 118]); ylim(ax2, [-0.04 1.18]);
    legend(ax2, {'LR \rightarrow STA', 'w-GaN DOS'}, ...
        'Location', 'northwest', 'Box', 'off', 'FontSize', 6.8);
    text(ax2, 0.02, 0.98, 'b', 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', 10, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'w', 'Margin', 0.2);
end

function fig = local_figure2(S)
    ene = S.ene_fit(:);
    de_ene = S.deconv_ene(:);
    target = S.deconv_target(:);
    att = S.att_feature(:);
    fit_mask = logical(S.fit_mask(:));
    Tfit = S.T_fit;

    c_loss = [252 128 2] ./ 255;    % #FC8002
    c_gain = [73 149 198] ./ 255;   % #4995C6
    c_att = [250 199 179] ./ 255;   % #FAC7B3
    c_fit = [97 64 153] ./ 255;     % #614099
    c_grid = [239 243 247] ./ 255;

    be_x = S.be_EM(:);
    be_sel = logical(S.be_selected(:));
    good_be = isfinite(be_x) & be_x > 0 & isfinite(double(be_sel));
    fit_lo = 40;
    fit_hi = 100;
    if any(good_be & be_sel)
        fit_lo = min(be_x(good_be & be_sel));
        fit_hi = max(be_x(good_be & be_sel));
    end

    fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'centimeters', ...
        'Position', [2 2 18.2 6.6], 'Renderer', 'painters');

    ax1 = axes(fig, 'Position', [0.06 0.20 0.25 0.68]);
    x = linspace(0, 150, 650).';
    loss = local_interp_positive(de_ene, target, x);
    gain = interp1(de_ene, target, -x, 'linear', NaN);
    scale = max([prctile(abs(loss(isfinite(loss))), 98), prctile(abs(gain(isfinite(gain))), 98), eps]);
    loss_n = smoothdata(loss ./ scale, 'movmean', 7);
    gain_n = smoothdata(gain ./ scale, 'movmean', 7);
    att_pos = interp1(ene, att, x, 'linear', 0);
    att_neg = interp1(ene, att, -x, 'linear', 0);
    att_pair = max(att_pos, att_neg);
    fit_window_mask = x >= fit_lo & x <= fit_hi;
    local_shade_mask_bands(ax1, x, fit_window_mask, [-0.10 1.08], c_att, 0.18);
    hold(ax1, 'on');
    plot(ax1, x, loss_n, '-', 'Color', c_loss, 'LineWidth', 1.20);
    plot(ax1, x, gain_n, '-', 'Color', c_gain, 'LineWidth', 1.20);
    plot(ax1, x, 0.98 .* local_norm01(att_pair, x > 0), '-', 'Color', c_att .* 0.85, 'LineWidth', 0.85);
    yline(ax1, 0, ':', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.6);
    grid(ax1, 'on'); box(ax1, 'on');
    ax1.GridColor = c_grid;
    ax1.GridAlpha = 0.45;
    ax1.LineWidth = 0.7;
    ax1.TickDir = 'out';
    xlabel(ax1, '|Energy| (meV)');
    ylabel(ax1, 'Normalized signal');
    xlim(ax1, [0 150]); ylim(ax1, [-0.10 1.08]);
    legend(ax1, {'fit window', 'loss', 'gain', 'attention'}, ...
        'Location', 'northeast', 'Box', 'off', 'FontSize', 6.5);
    text(ax1, 0.02, 0.98, 'a', 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', 10, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'w', 'Margin', 0.2);

    ax2 = axes(fig, 'Position', [0.39 0.20 0.25 0.68]);
    be_x = S.be_EM(:);
    be_y = S.be_ln_ratio(:);
    be_w = S.be_weight(:);
    be_sel = logical(S.be_selected(:));
    good = isfinite(be_x) & isfinite(be_y) & isfinite(be_w) & be_x > 0;
    be_x = be_x(good); be_y = be_y(good); be_w = be_w(good); be_sel = be_sel(good);
    hold(ax2, 'on');
    if any(~be_sel)
        scatter(ax2, be_x(~be_sel), be_y(~be_sel), 18 + 42 .* be_w(~be_sel), ...
            [0.70 0.70 0.76], 'filled', 'MarkerFaceAlpha', 0.42, 'MarkerEdgeColor', 'none');
    end
    scatter(ax2, be_x(be_sel), be_y(be_sel), 20 + 55 .* be_w(be_sel), ...
        c_fit, 'filled', 'MarkerFaceAlpha', 0.75, 'MarkerEdgeColor', 'w', 'LineWidth', 0.35);
    xx = linspace(0, 100, 200);
    plot(ax2, xx, S.be_slope .* xx, '-', 'Color', [0.05 0.05 0.05], 'LineWidth', 1.1);
    grid(ax2, 'on'); box(ax2, 'on');
    ax2.GridColor = c_grid;
    ax2.GridAlpha = 0.45;
    ax2.LineWidth = 0.7;
    ax2.TickDir = 'out';
    xlabel(ax2, 'Sparse source energy (meV)');
    ylabel(ax2, 'ln(P_{loss}/P_{gain})');
    xlim(ax2, [0 100]);
    ymax = max([be_y(:); S.be_slope .* xx(:)], [], 'omitnan');
    ylim(ax2, [0 max(1.2, 1.12 .* ymax)]);
    text(ax2, 0.18, 0.92, sprintf('T = %.0f K', Tfit), ...
        'Units', 'normalized', 'FontSize', 7.2, 'VerticalAlignment', 'top', ...
        'BackgroundColor', 'w', 'Margin', 2, 'Interpreter', 'none');
    text(ax2, 0.02, 0.98, 'b', 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', 10, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'w', 'Margin', 0.2);

    ax3 = axes(fig, 'Position', [0.72 0.20 0.25 0.68]);
    T = S.T_grid(:);
    cost_raw = S.cost_raw(:);
    finite = isfinite(T) & isfinite(cost_raw);
    T = T(finite); cost_raw = cost_raw(finite);
    [T, ord] = sort(T);
    cost_raw = cost_raw(ord);
    valley_mask = T >= 600 & T <= 1000;
    T = T(valley_mask);
    cost_raw = cost_raw(valley_mask);
    cden = max(max(cost_raw, [], 'omitnan') - min(cost_raw, [], 'omitnan'), eps);
    raw_delta = (cost_raw - min(cost_raw, [], 'omitnan')) ./ cden;
    Td = linspace(min(T, [], 'omitnan'), max(T, [], 'omitnan'), 480).';
    if numel(raw_delta) >= 3
        smooth_delta = interp1(T, raw_delta, Td, 'pchip');
    else
        smooth_delta = interp1(T, raw_delta, Td, 'linear');
    end
    plot(ax3, T, raw_delta, 'o', 'Color', [0.60 0.61 0.66], ...
        'MarkerSize', 3.5, 'LineWidth', 0.65); hold(ax3, 'on');
    plot(ax3, Td, smooth_delta, '-', 'Color', c_fit, 'LineWidth', 1.55);
    xline(ax3, Tfit, '-', sprintf(' T = %.0f K', Tfit), 'Color', [0.16 0.16 0.16], ...
        'LineWidth', 0.9, 'LabelVerticalAlignment', 'bottom', 'FontSize', 7);
    ylabel(ax3, 'Normalized cost');
    ylim(ax3, [-0.03, max(1.02, max(smooth_delta, [], 'omitnan') .* 1.08)]);
    grid(ax3, 'on'); box(ax3, 'on');
    ax3.GridColor = c_grid;
    ax3.GridAlpha = 0.45;
    ax3.LineWidth = 0.7;
    ax3.TickDir = 'out';
    xlabel(ax3, 'Temperature candidate (K)');
    legend(ax3, {'grid', 'smooth', 'T'}, ...
        'Location', 'northwest', 'Box', 'off', 'FontSize', 6.4);
    xlim(ax3, [670 1000]);
    text(ax3, 0.02, 0.98, 'c', 'Units', 'normalized', 'FontWeight', 'bold', ...
        'FontSize', 10, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'BackgroundColor', 'w', 'Margin', 0.2);
end

function local_export_figure(fig, out_dir, basename)
    svg_file = fullfile(out_dir, [basename '.svg']);
    png_file = fullfile(out_dir, [basename '.png']);
    tif_file = fullfile(out_dir, [basename '.tif']);
    set(fig, 'InvertHardcopy', 'off');
    print(fig, svg_file, '-dsvg', '-painters');
    print(fig, png_file, '-dpng', '-r360');
    print(fig, tif_file, '-dtiff', '-r360');
    fprintf('Wrote %s\n', svg_file);
    fprintf('Wrote %s\n', png_file);
    fprintf('Wrote %s\n', tif_file);
end

function local_export_panel_csv(S, out_dir)
    ene = S.ene_fit(:);
    T1 = table(ene, S.raw_original(:), S.denoised_used(:), S.model(:), S.zlp(:), S.vib(:), ...
        'VariableNames', {'energy_meV', 'raw_counts', 'svd_denoised_counts', ...
        'stage2_fit_counts', 'zlp_counts', 'phonon_reconstruction_counts'});
    writetable(T1, fullfile(out_dir, 'Figure1a_spectrum_components.csv'));

    xgrid = linspace(0, 118, 900).';
    lr_curve = local_lr_deconvolution_curve(S, xgrid);
    lr_sta_curve = max(local_norm01(local_sta_envelope_removed(lr_curve, xgrid), xgrid >= 10), 0);
    gan_wz = local_interp_and_norm(S.energy_gan_wurtzite(:), S.dos_gan_wurtzite(:), xgrid);
    T2 = table(xgrid, lr_sta_curve, gan_wz, ...
        'VariableNames', {'energy_meV', 'LR_then_STA', 'wurtzite_GaN_DOS_mp804'});
    writetable(T2, fullfile(out_dir, 'Figure1b_GaN_DOS_comparison.csv'));

    x = linspace(0, 210, 850).';
    att_pos = interp1(ene, S.att_feature(:), x, 'linear', 0);
    att_neg = interp1(ene, S.att_feature(:), -x, 'linear', 0);
    T3 = table(x, local_interp_positive(S.deconv_ene(:), S.deconv_target(:), x), ...
        interp1(S.deconv_ene(:), S.deconv_target(:), -x, 'linear', NaN), max(att_pos, att_neg), ...
        'VariableNames', {'abs_energy_meV', 'loss_EEL_zlp_subtracted', ...
        'gain_EEG_zlp_subtracted', 'attention_pair'});
    writetable(T3, fullfile(out_dir, 'Figure2a_attention_gain_loss.csv'));

    T4 = table(S.be_EM(:), S.be_ln_ratio(:), S.be_weight(:), logical(S.be_selected(:)), ...
        'VariableNames', {'EM_meV', 'ln_PM_loss_gain', 'weight', 'selected'});
    writetable(T4, fullfile(out_dir, 'Figure2b_BE_linearity.csv'));

    T5 = table(S.T_grid(:), S.cost_raw(:), S.cost_smooth(:), S.linearity_R2(:), ...
        S.pm_smooth_penalty(:), S.zlp_delta_rms(:), S.effective_pair_count(:), ...
        'VariableNames', {'T_grid_K', 'cost_raw', 'cost_smooth', 'linearity_R2', ...
        'PM_smooth_penalty', 'ZLP_delta_rms', 'effective_pair_count'});
    writetable(T5, fullfile(out_dir, 'Figure2c_cost_valley.csv'));
end

function y = local_positive_floor(y)
    y = abs(y(:));
    good = isfinite(y) & y > 0;
    if any(good)
        floor_val = max(min(y(good)) .* 0.5, max(y(good)) .* 1e-7);
    else
        floor_val = 1e-12;
    end
    y(~isfinite(y) | y <= 0) = floor_val;
end

function local_plot_pm_sticks(ax, EM, PMs, PMa, c_stokes, c_gain, base, height, xmax, ene_ref, signal_ref)
    EM = EM(:); PMs = max(PMs(:), 0); PMa = max(PMa(:), 0);
    if isempty(EM), return; end
    if nargin < 7 || isempty(base), base = 0; end
    if nargin < 8 || isempty(height), height = 1; end
    if nargin < 9 || isempty(xmax), xmax = max(EM, [], 'omitnan'); end
    if nargin < 10, ene_ref = []; end
    if nargin < 11, signal_ref = []; end
    scale = height ./ max([PMs; PMa; eps], [], 'omitnan');
    keep = EM <= xmax & isfinite(EM);
    for ii = find(keep(:).')
        hs = PMs(ii) .* scale;
        ha = PMa(ii) .* scale;
        if ~isempty(ene_ref) && ~isempty(signal_ref)
            cap_s = 0.72 .* interp1(ene_ref(:), signal_ref(:), EM(ii), 'linear', NaN);
            cap_a = 0.72 .* interp1(ene_ref(:), signal_ref(:), -EM(ii), 'linear', NaN);
            if isfinite(cap_s) && cap_s > 0, hs = min(hs, cap_s); end
            if isfinite(cap_a) && cap_a > 0, ha = min(ha, cap_a); end
        end
        if hs > 0
            line(ax, [EM(ii) EM(ii)], [base base + hs], 'Color', c_stokes, 'LineWidth', 0.45);
        end
        if ha > 0
            line(ax, [-EM(ii) -EM(ii)], [base base + ha], 'Color', c_gain, 'LineWidth', 0.45);
        end
    end
end

function y = local_interp_positive(x, y0, xgrid)
    x = x(:); y0 = y0(:);
    good = isfinite(x) & isfinite(y0) & x >= 0;
    if nnz(good) < 2
        y = NaN(size(xgrid));
    else
        [xu, ia] = unique(x(good), 'stable');
        yu = y0(good);
        yu = yu(ia);
        y = interp1(xu, yu, xgrid, 'linear', NaN);
    end
end

function y = local_gaussian_broaden(centers, amp, xgrid, sigma)
    centers = centers(:);
    amp = max(amp(:), 0);
    xgrid = xgrid(:);
    y = zeros(size(xgrid));
    good = isfinite(centers) & isfinite(amp) & centers >= 0 & amp > 0;
    centers = centers(good);
    amp = amp(good);
    if isempty(centers), return; end
    sigma = max(sigma, eps);
    for ii = 1:numel(centers)
        y = y + amp(ii) .* exp(-0.5 .* ((xgrid - centers(ii)) ./ sigma) .^ 2);
    end
end

function y = local_interp_and_norm(x, y0, xgrid)
    x = x(:); y0 = max(y0(:), 0);
    good = isfinite(x) & isfinite(y0) & x >= 0;
    if nnz(good) < 2
        y = NaN(size(xgrid));
        return;
    end
    [xu, ia] = unique(x(good), 'stable');
    yu = y0(good);
    yu = yu(ia);
    y = interp1(xu, yu, xgrid(:), 'linear', 0);
    y = smoothdata(y, 'movmean', 7);
    y = local_norm01(y, xgrid(:) >= 0);
end

function y = local_sta_envelope_removed(signal_curve, energy_meV)
    signal_curve = max(signal_curve(:), 0);
    energy_meV = energy_meV(:);
    % Vibrational EELS carries a strong low-energy transfer envelope. For
    % a DOS-like comparison, remove that smooth envelope while preserving
    % peak positions. The exponent is intentionally mild and only affects
    % Fig.1b visualization.
    envelope_inverse = ((max(energy_meV, 0) + 8) ./ 88) .^ 1.8;
    y = signal_curve .* envelope_inverse;
end

function y = local_lr_deconvolution_curve(S, xgrid)
    ene = S.ene_fit(:);
    target = max(S.denoised_used(:) - S.zlp(:) - S.bg(:), 0);
    zlp = max(S.zlp(:), 0);
    n = numel(ene);
    if n < 16 || numel(target) ~= n || numel(zlp) ~= n
        y = zeros(size(xgrid(:)));
        return;
    end
    eabs = abs(ene);
    lo = 6;
    hi = 112;
    edge = 7;
    taper = min(1, max(0, (eabs - lo) ./ max(edge, eps))) .* ...
        min(1, max(0, (hi - eabs) ./ max(edge, eps)));
    positive = target(:) .* taper(:);
    peak = max(positive, [], 'omitnan');
    if ~isfinite(peak) || peak <= 0
        y = zeros(size(xgrid(:)));
        return;
    end
    psf = zlp(:);
    psf(~isfinite(psf) | psf < 0) = 0;
    psf = psf ./ max(sum(psf, 'omitnan'), eps);
    psf_flip = flipud(psf);
    floor_val = max(peak .* 1e-6, eps);
    latent = max(positive, floor_val .* double(taper > 0.01));
    latent = latent .* taper;
    for ii = 1:36
        estimate = conv(latent, psf, 'same');
        estimate = max(estimate, floor_val);
        ratio = positive ./ estimate;
        ratio(~isfinite(ratio)) = 0;
        corr = conv(ratio, psf_flip, 'same');
        corr(~isfinite(corr) | corr < 0) = 0;
        latent = latent .* max(corr, 0) .^ 0.55;
        latent(~isfinite(latent) | latent < 0) = 0;
        latent = smoothdata(latent, 'gaussian', 11) .* taper;
    end
    y = interp1(ene, latent, xgrid(:), 'linear', 0);
    y = max(y, 0);
    y = y .* local_soft_high_energy_gate(xgrid(:), 86, 98, 0.005);
end

function gate = local_soft_high_energy_gate(x, start_meV, stop_meV, floor_value)
    x = x(:);
    gate = ones(size(x));
    mid = x > start_meV & x < stop_meV;
    gate(mid) = floor_value + (1 - floor_value) .* ...
        0.5 .* (1 + cos(pi .* (x(mid) - start_meV) ./ max(stop_meV - start_meV, eps)));
    gate(x >= stop_meV) = floor_value;
end

function y = local_norm01(y, mask)
    y = y(:);
    if nargin < 2 || isempty(mask)
        mask = isfinite(y);
    else
        mask = logical(mask(:)) & isfinite(y);
    end
    if nnz(mask) < 2
        y(:) = 0;
        return;
    end
    y = y - min(y(mask), [], 'omitnan');
    den = max(y(mask), [], 'omitnan');
    if ~isfinite(den) || den <= 0
        den = max(abs(y(mask)), [], 'omitnan');
    end
    if ~isfinite(den) || den <= 0
        y(:) = 0;
    else
        y = y ./ den;
    end
end

function local_shade_mask_bands(ax, x, mask, ylim_vals, color, alpha)
    bands = local_mask_to_bands(x(:), logical(mask(:)));
    hold(ax, 'on');
    for ii = 1:size(bands, 1)
        patch(ax, [bands(ii, 1) bands(ii, 2) bands(ii, 2) bands(ii, 1)], ...
            [ylim_vals(1) ylim_vals(1) ylim_vals(2) ylim_vals(2)], color, ...
            'FaceAlpha', alpha, 'EdgeColor', 'none', 'HandleVisibility', 'on');
    end
end

function bands = local_mask_to_bands(x, mask)
    x = x(:);
    mask = logical(mask(:));
    if isempty(x) || ~any(mask)
        bands = zeros(0, 2);
        return;
    end
    d = diff([false; mask; false]);
    starts = find(d == 1);
    stops = find(d == -1) - 1;
    bands = [x(starts), x(stops)];
    if size(bands, 1) > 1
        merged = bands(1, :);
        for ii = 2:size(bands, 1)
            if bands(ii, 1) - merged(end, 2) <= 3
                merged(end, 2) = bands(ii, 2);
            else
                merged(end + 1, :) = bands(ii, :); %#ok<AGROW>
            end
        end
        bands = merged;
    end
end
