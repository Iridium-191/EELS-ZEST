function export_example_figure_data()
%EXPORT_EXAMPLE_FIGURE_DATA Build figure-ready data from the bundled example.
%
% Workflow:
%   example_data.mat -> crop [-250,250] meV -> SVD/PCA denoise ->
%   pick a high-signal single spectrum -> EELS-ZEST 2.0 Stage1/Stage2 ->
%   PM-guided deconvolution and BE-linearity diagnostics.

    script_dir = fileparts(mfilename('fullpath'));
    zest_dir = fileparts(fileparts(script_dir));
    workspace_dir = getenv('EELS_ZEST_WORKSPACE_ROOT');
    if isempty(workspace_dir), workspace_dir = fileparts(zest_dir); end
    addpath(zest_dir);

    input_file = fullfile(zest_dir, 'examples', 'data', 'example_data.mat');
    phonon_dir = getenv('EELS_ZEST_PHONON_DIR');
    if isempty(phonon_dir), phonon_dir = fullfile(workspace_dir, 'phonon'); end
    gan_wurtzite_file = fullfile(phonon_dir, 'mp-804.json');
    gan_zincblende_file = fullfile(phonon_dir, 'mp-830.json');
    gan_rocksalt_file = fullfile(phonon_dir, 'mp-2853.json');

    out_dir = script_dir;
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    fprintf('Loading %s\n', input_file);
    [cube, ene, info] = EELS_PearsonSparse_FullFit.load_datacube_mat(input_file);
    ene = double(ene(:));

    crop_range = [-250, 250];
    crop_mask = ene >= crop_range(1) & ene <= crop_range(2);
    if nnz(crop_mask) < 50
        error('Energy crop [%g,%g] meV leaves too few channels.', crop_range(1), crop_range(2));
    end
    cube_crop = double(cube(:, :, crop_mask));
    ene_crop = ene(crop_mask);

    svd_rank = 12;
    svd_smooth_points = 1;
    [cube_svd, svd_info] = local_svd_denoise(cube_crop, svd_rank, svd_smooth_points);
    fprintf('SVD denoise complete: rank %d, smoothing %d\n', ...
        svd_info.rank, svd_info.energy_smooth_points);

    candidates = local_rank_candidate_pixels(cube_svd, ene_crop, 18);
    if isempty(candidates)
        error('No candidate spectra found after SVD denoise.');
    end

    settings = struct();
    settings.run_stage2 = true;
    settings.save_plots = false;
    settings.stage1_speed_profile = 'v20_fast';
    settings.stage1_zlp_align_enabled = false;
    settings.pm_deconv_max_centers = 1000;
    settings.pm_deconv_min_pm_fraction = 0.018;
    settings.pm_deconv_min_spacing_meV = 2.5;
    settings.pm_deconv_kernel_width_scale = 0.75;
    settings.pm_deconv_pm_prior_weight = 0.06;
    settings.pm_deconv_l2_weight = 0.004;
    settings.pm_deconv_be_prior_weight = 0.045;

    best = struct('score', -Inf, 'idx', NaN, 'yy', NaN, 'xx', NaN, ...
        'res', [], 'deconv', [], 'befit', []);

    fprintf('Trying %d high-signal candidate spectra...\n', numel(candidates));
    for kk = 1:numel(candidates)
        pix = candidates(kk);
        [yy, xx] = ind2sub(size(cube_svd, [1 2]), pix);
        spec_svd = squeeze(cube_svd(yy, xx, :));
        try
            evalc('res = EELS_PearsonSparse_FullFit.solve_spectrum(ene_crop, spec_svd, settings);');
        catch ME
            fprintf('  candidate %02d pixel (%d,%d) failed: %s\n', kk, yy, xx, ME.message);
            continue;
        end
        try
            deconv = EELS_PearsonSparse_FullFit.pm_guided_deconvolution(res, settings);
        catch ME
            deconv = struct('available', false, 'reason', ME.message);
        end
        try
            befit = EELS_PearsonSparse_FullFit.hard_be_ln_ratio_fit(res);
        catch ME
            befit = struct('available', false, 'reason', ME.message);
        end

        if ~(isfield(res, 'stage2') && isfield(res.stage2, 'status') && strcmp(res.stage2.status, 'ok'))
            fprintf('  candidate %02d pixel (%d,%d): Stage2 not ok\n', kk, yy, xx);
            continue;
        end
        if ~isfield(res.stage2, 'T_fit') || ~isfinite(res.stage2.T_fit)
            fprintf('  candidate %02d pixel (%d,%d): no finite Stage2 T\n', kk, yy, xx);
            continue;
        end

        r2 = local_get_numeric(res, 'R2', NaN);
        nrmse = local_get_numeric(res, 'NRMSE', NaN);
        rel = local_get_numeric(res.stage2, 'reliability_score', 0);
        att_frac = 0;
        if isfield(res.stage2, 'attention') && isfield(res.stage2.attention, 'selected_fraction')
            att_frac = res.stage2.attention.selected_fraction;
        end
        deconv_ok = isfield(deconv, 'available') && deconv.available;
        be_ok = isfield(befit, 'available') && befit.available;
        deconv_r2 = local_get_numeric(deconv, 'R2_weighted', -1);
        be_r2 = local_get_numeric(befit, 'R2_origin', -1);
        score = 1.0 .* max(r2, 0) - 0.25 .* max(nrmse, 0) + ...
            0.35 .* max(rel, 0) + 0.18 .* max(att_frac, 0) + ...
            0.28 .* max(deconv_r2, 0) + 0.20 .* max(be_r2, 0) + ...
            0.35 .* double(deconv_ok) + 0.20 .* double(be_ok);

        fprintf('  candidate %02d pixel (%d,%d): T %.1f K, R2 %.5f, NRMSE %.4g, deconv %d, BE %d, score %.3f\n', ...
            kk, yy, xx, res.stage2.T_fit, r2, nrmse, deconv_ok, be_ok, score);

        if score > best.score
            best.score = score;
            best.idx = pix;
            best.yy = yy;
            best.xx = xx;
            best.res = res;
            best.deconv = deconv;
            best.befit = befit;
        end
        if kk >= 8 && best.score > 2.4 && isfield(best.deconv, 'available') && best.deconv.available && ...
                isfield(best.befit, 'available') && best.befit.available
            break;
        end
    end

    if ~isfinite(best.score)
        error('No usable Stage2 spectrum found among candidate spectra.');
    end

    yy = best.yy;
    xx = best.xx;
    res = best.res;
    deconv = best.deconv;
    befit = best.befit;
    raw_original = squeeze(cube_crop(yy, xx, :));
    denoised_used = squeeze(cube_svd(yy, xx, :));

    n = numel(res.ene_corrected);
    ene_fit = res.ene_corrected(:);
    norm_val = local_get_numeric(res, 'norm_val', 1);
    if ~isfinite(norm_val) || norm_val <= 0, norm_val = 1; end

    model = local_stage2_or_result_vector(res, 'model_norm', 'fitted_spectrum', n) .* norm_val;
    zlp = local_stage2_or_result_vector(res, 'zlp_norm', 'zlp_only', n, norm_val) .* norm_val;
    vib = local_stage2_or_result_vector(res, 'vib_norm', 'vib_only', n, norm_val) .* norm_val;
    bg = local_stage2_or_result_vector(res, 'bg_norm', 'bg_only', n, norm_val) .* norm_val;
    residual = denoised_used(:) - model(:);

    EM = local_stage2_vector(res, 'EM_final');
    PM_stokes = local_stage2_vector(res, 'PM_stokes');
    PM_antistokes = local_stage2_vector(res, 'PM_antistokes');

    att_feature = zeros(n, 1);
    att_fit_alpha = zeros(n, 1);
    att_active = false(n, 1);
    fit_mask = false(n, 1);
    if isfield(res.stage2, 'attention')
        A = res.stage2.attention;
        if isfield(A, 'feature_score') && numel(A.feature_score) >= n
            att_feature = A.feature_score(1:n);
        end
        if isfield(A, 'fit_alpha') && numel(A.fit_alpha) >= n
            att_fit_alpha = A.fit_alpha(1:n);
        end
        if isfield(A, 'active_mask') && numel(A.active_mask) >= n
            att_active = logical(A.active_mask(1:n));
        end
    end
    if isfield(res.stage2, 'fit_mask') && numel(res.stage2.fit_mask) >= n
        fit_mask = logical(res.stage2.fit_mask(1:n));
    end

    profile = res.stage2.temperature_profile;
    T_grid = profile.T_grid(:);
    cost_raw = profile.cost_raw(:);
    cost_smooth = profile.cost_smooth(:);
    linearity_R2 = profile.linearity_R2(:);
    pm_smooth_penalty = profile.PM_smooth_penalty(:);
    zlp_delta_rms = profile.ZLP_delta_rms(:);
    effective_pair_count = profile.effective_pair_count(:);

    deconv_available = isfield(deconv, 'available') && deconv.available;
    if deconv_available
        deconv_ene = deconv.ene(:);
        deconv_target = deconv.target_signal(:);
        deconv_reconvolved = deconv.reconvolved_signal(:);
        deconv_EM = deconv.EM_meV(:);
        deconv_coef_stokes = deconv.coef_stokes(:);
        deconv_coef_antistokes = deconv.coef_antistokes(:);
    else
        deconv_ene = ene_fit;
        deconv_target = (denoised_used(:) - zlp(:) - bg(:)) ./ norm_val;
        deconv_reconvolved = vib(:) ./ norm_val;
        deconv_EM = EM(:);
        deconv_coef_stokes = PM_stokes(:);
        deconv_coef_antistokes = PM_antistokes(:);
    end

    be_available = isfield(befit, 'available') && befit.available;
    if be_available
        be_EM = befit.EM_meV(:);
        be_ln_ratio = befit.ln_loss_gain(:);
        be_weight = befit.weight(:);
        be_selected = logical(befit.selected_mask(:));
        be_slope = befit.slope_per_meV;
        be_R2 = befit.R2_origin;
    else
        be_EM = EM(:);
        be_ln_ratio = log(max(PM_stokes(:), eps) ./ max(PM_antistokes(:), eps));
        be_weight = sqrt(max(PM_stokes(:), 0) .* max(PM_antistokes(:), 0));
        be_weight = be_weight ./ max(max(be_weight), eps);
        be_selected = true(size(be_EM));
        be_slope = 1 ./ (0.08617333262 .* res.stage2.T_fit);
        be_R2 = NaN;
    end

    [energy_gan_wurtzite, dos_gan_wurtzite, gan_wurtzite_meta] = local_load_materials_project_phonon_json(gan_wurtzite_file);
    [energy_gan_zincblende, dos_gan_zincblende, gan_zincblende_meta] = local_load_materials_project_phonon_json(gan_zincblende_file);
    [energy_gan_rocksalt, dos_gan_rocksalt, gan_rocksalt_meta] = local_load_materials_project_phonon_json(gan_rocksalt_file);

    T_fit = res.stage2.T_fit;
    R2 = local_get_numeric(res, 'R2', NaN);
    NRMSE = local_get_numeric(res, 'NRMSE', NaN);
    reliability = local_get_numeric(res.stage2, 'reliability_score', NaN);
    kb_mev_K = 0.08617333262;

    out_mat = fullfile(out_dir, 'example_figure_data.mat');
    save(out_mat, ...
        'ene_fit', 'raw_original', 'denoised_used', 'model', 'zlp', 'vib', 'bg', 'residual', ...
        'EM', 'PM_stokes', 'PM_antistokes', 'att_feature', 'att_fit_alpha', 'att_active', 'fit_mask', ...
        'T_grid', 'cost_raw', 'cost_smooth', 'linearity_R2', 'pm_smooth_penalty', ...
        'zlp_delta_rms', 'effective_pair_count', 'deconv_ene', 'deconv_target', ...
        'deconv_reconvolved', 'deconv_EM', 'deconv_coef_stokes', 'deconv_coef_antistokes', ...
        'be_EM', 'be_ln_ratio', 'be_weight', 'be_selected', 'be_slope', 'be_R2', ...
        'energy_gan_wurtzite', 'dos_gan_wurtzite', 'gan_wurtzite_meta', ...
        'energy_gan_zincblende', 'dos_gan_zincblende', 'gan_zincblende_meta', ...
        'energy_gan_rocksalt', 'dos_gan_rocksalt', 'gan_rocksalt_meta', ...
        'T_fit', 'R2', 'NRMSE', 'reliability', 'kb_mev_K', 'yy', 'xx', 'svd_info', 'info', ...
        'deconv_available', 'be_available', '-v7');

    out_res = fullfile(out_dir, 'example_selected_fit_result.mat');
    save(out_res, 'res', 'deconv', 'befit', 'settings', 'raw_original', 'denoised_used', ...
        'ene_crop', 'yy', 'xx', '-v7.3');

    summary_file = fullfile(out_dir, 'example_figure_summary.txt');
    fid = fopen(summary_file, 'w');
    if fid >= 0
        fprintf(fid, 'EELS-ZEST 2.0 example figure data\n');
        fprintf(fid, 'input_file: examples/data/example_data.mat\n');
        fprintf(fid, 'phonon_GaN_wurtzite_file: %s\n', gan_wurtzite_file);
        fprintf(fid, 'phonon_GaN_zincblende_file: %s\n', gan_zincblende_file);
        fprintf(fid, 'phonon_GaN_rocksalt_file: %s\n', gan_rocksalt_file);
        fprintf(fid, 'selected_pixel_yx: [%d, %d]\n', yy, xx);
        fprintf(fid, 'svd_rank: %d\n', svd_info.rank);
        fprintf(fid, 'energy_crop_meV: [%.3f, %.3f]\n', crop_range(1), crop_range(2));
        fprintf(fid, 'T_fit_K: %.8g\n', T_fit);
        fprintf(fid, 'R2: %.8g\n', R2);
        fprintf(fid, 'NRMSE: %.8g\n', NRMSE);
        fprintf(fid, 'reliability: %.8g\n', reliability);
        fprintf(fid, 'deconv_available: %d\n', deconv_available);
        if deconv_available
            fprintf(fid, 'deconv_R2_weighted: %.8g\n', deconv.R2_weighted);
            fprintf(fid, 'deconv_NRMSE_weighted: %.8g\n', deconv.NRMSE_weighted);
        end
        fprintf(fid, 'be_available: %d\n', be_available);
        fprintf(fid, 'be_R2_origin: %.8g\n', be_R2);
        fprintf(fid, 'figure_data: %s\n', out_mat);
        fclose(fid);
    end

    fprintf('Wrote %s\n', out_mat);
    fprintf('Wrote %s\n', out_res);
    fprintf('Wrote %s\n', summary_file);
end

function [cube_denoised, info_out] = local_svd_denoise(cube, requested_rank, smooth_points)
    [ny, nx, ne] = size(cube);
    npix = ny * nx;
    default_rank = min([12, max(1, npix - 1), max(1, ne - 1)]);
    if nargin < 2 || isempty(requested_rank) || ~isfinite(requested_rank)
        rk = default_rank;
    else
        rk = max(1, round(double(requested_rank)));
    end
    rk = min(rk, min(npix, ne));
    if nargin < 3 || isempty(smooth_points) || ~isfinite(smooth_points)
        smooth_points = 1;
    end
    smooth_points = max(1, round(double(smooth_points)));

    X = reshape(double(cube), npix, ne);
    finite = isfinite(X);
    col_ref = median(X, 1, 'omitnan');
    col_ref(~isfinite(col_ref)) = 0;
    for jj = 1:ne
        bad = ~finite(:, jj);
        if any(bad), X(bad, jj) = col_ref(jj); end
    end
    mu = mean(X, 1, 'omitnan');
    mu(~isfinite(mu)) = 0;
    Xc = X - mu;

    try
        if rk < min(size(Xc)) - 1
            [U, S, V] = svds(Xc, rk);
        else
            [U, S, V] = svd(Xc, 'econ');
            U = U(:, 1:rk); S = S(1:rk, 1:rk); V = V(:, 1:rk);
        end
    catch
        [U, S, V] = svd(Xc, 'econ');
        rk = min(rk, size(S, 1));
        U = U(:, 1:rk); S = S(1:rk, 1:rk); V = V(:, 1:rk);
    end
    Xden = U * S * V' + mu;
    if smooth_points > 1 && ne >= smooth_points
        if mod(smooth_points, 2) == 0, smooth_points = smooth_points + 1; end
        Xden = movmean(Xden, smooth_points, 2, 'omitnan');
    end
    cube_denoised = reshape(Xden, ny, nx, ne);
    info_out = struct('enabled', true, 'method', 'low_rank_svd_pca_gui_equivalent', ...
        'rank', rk, 'default_rank', default_rank, 'energy_smooth_points', smooth_points, ...
        'input_size', [ny, nx, ne]);
end

function candidates = local_rank_candidate_pixels(cube_svd, ene, max_count)
    [ny, nx, ne] = size(cube_svd);
    X = reshape(double(cube_svd), ny * nx, ne);
    ene = ene(:).';
    loss_mask = ene >= 35 & ene <= 190;
    gain_mask = ene <= -35 & ene >= -190;
    tail_mask = abs(ene) >= 210 & abs(ene) <= 245;
    core_mask = abs(ene) <= 8;
    if nnz(tail_mask) < 5
        tail_mask = abs(ene) >= 180;
    end
    base = median(X(:, tail_mask), 2, 'omitnan');
    loss_sig = sum(max(X(:, loss_mask) - base, 0), 2, 'omitnan') ./ max(nnz(loss_mask), 1);
    gain_sig = sum(max(X(:, gain_mask) - base, 0), 2, 'omitnan') ./ max(nnz(gain_mask), 1);
    core_amp = max(X(:, core_mask), [], 2, 'omitnan');
    tail_noise = 1.4826 .* median(abs(X(:, tail_mask) - median(X(:, tail_mask), 2, 'omitnan')), 2, 'omitnan');
    sym = min(loss_sig, gain_sig) ./ max(max(loss_sig, gain_sig), eps);
    score = loss_sig .* (0.35 + 0.65 .* max(sym, 0)) ./ max(tail_noise, eps);
    score = score .* double(isfinite(score) & isfinite(core_amp) & core_amp > 0);

    [~, order] = sort(score, 'descend', 'MissingPlacement', 'last');
    order = order(isfinite(score(order)) & score(order) > 0);
    if isempty(order)
        candidates = [];
        return;
    end
    keep = false(size(order));
    coords = zeros(numel(order), 2);
    for ii = 1:numel(order)
        [yy, xx] = ind2sub([ny, nx], order(ii));
        coords(ii, :) = [yy, xx];
    end
    selected = zeros(1, 0);
    min_dist = 4;
    for ii = 1:numel(order)
        if numel(selected) >= max_count, break; end
        yy = coords(ii, 1); xx = coords(ii, 2);
        if yy <= 1 || yy >= ny || xx <= 1 || xx >= nx
            continue;
        end
        if isempty(selected)
            selected(end + 1) = order(ii); %#ok<AGROW>
            keep(ii) = true; %#ok<NASGU>
            continue;
        end
        [sy, sx] = ind2sub([ny, nx], selected);
        if all(hypot(double(sy(:)) - yy, double(sx(:)) - xx) >= min_dist)
            selected(end + 1) = order(ii); %#ok<AGROW>
            keep(ii) = true; %#ok<NASGU>
        end
    end
    if isempty(selected)
        selected = order(1:min(max_count, numel(order))).';
    end
    candidates = selected(:).';
end

function val = local_get_numeric(s, field, default_value)
    val = default_value;
    if isstruct(s) && isfield(s, field)
        tmp = s.(field);
        if isnumeric(tmp) && isscalar(tmp) && isfinite(tmp)
            val = double(tmp);
        end
    end
end

function v = local_stage2_vector(res, name)
    v = zeros(0, 1);
    if isfield(res, 'stage2') && isstruct(res.stage2) && isfield(res.stage2, name)
        tmp = res.stage2.(name);
        if isnumeric(tmp)
            v = double(tmp(:));
        end
    end
end

function v = local_stage2_or_result_vector(res, stage2_name, fallback_name, n, fallback_norm)
    if nargin < 5, fallback_norm = 1; end
    v = zeros(n, 1);
    if isfield(res, 'stage2') && isstruct(res.stage2) && isfield(res.stage2, stage2_name)
        tmp = res.stage2.(stage2_name);
        if isnumeric(tmp) && numel(tmp) >= n
            v = double(tmp(1:n));
            return;
        end
    end
    if isfield(res, fallback_name)
        tmp = res.(fallback_name);
        if isnumeric(tmp) && numel(tmp) >= n
            v = double(tmp(1:n)) ./ max(fallback_norm, eps);
        end
    end
end

function [energy_meV, dos, meta] = local_load_materials_project_phonon_json(filename)
    raw = fileread(filename);
    data = jsondecode(raw);
    phonon = data.phonon;
    freq_cm1 = double(phonon.dos_frequencies(:));
    dos = double(phonon.ph_dos(:));
    cm1_per_meV = 8.065543937;
    energy_meV = freq_cm1 ./ cm1_per_meV;
    keep = isfinite(energy_meV) & isfinite(dos) & energy_meV >= 0;
    energy_meV = energy_meV(keep);
    dos = dos(keep);
    if ~isempty(dos)
        dos = max(dos, 0);
    end
    meta = struct();
    meta.source_file = filename;
    if isfield(data, 'metadata')
        md = data.metadata;
        if isfield(md, 'material_id'), meta.material_id = md.material_id; end
        if isfield(md, 'formula'), meta.formula = md.formula; end
        if isfield(md, 'point_group'), meta.point_group = md.point_group; end
        if isfield(md, 'nsites'), meta.nsites = md.nsites; end
        if isfield(md, 'qpoints_grid'), meta.qpoints_grid = md.qpoints_grid; end
    end
    meta.frequency_unit = 'cm^-1';
    meta.energy_unit = 'meV';
    meta.energy_conversion = 'meV = cm^-1 / 8.065543937';
end
