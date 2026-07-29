classdef EELS_PearsonSparse_FullFit
    % Full-spectrum Pearson-ZLP + sparse convolution prototype.
    %
    % Stage 1:
    %   I(E) = I_ZLP(E) + sum PM_s(i)*K(E-EM_i) + sum PM_a(i)*K(E+EM_i)
    %   PM_s and PM_a are independent reconstruction variables.
    %   BE detailed balance is only a soft penalty in the optimizer.
    %
    % Stage 2:
    %   Placeholder temperature-window hard-BE refinement using Stage-1 PM as
    %   the initial source. The first development target is Stage 1 fidelity.

    methods (Static)
        function settings = defaults(settings)
            if nargin < 1 || isempty(settings), settings = struct(); end
            if ~isfield(settings, 'crop_range'), settings.crop_range = [-Inf, Inf]; end
            if ~isfield(settings, 'num_zlp_components'), settings.num_zlp_components = 4; end
            if ~isfield(settings, 'zlp_model_type'), settings.zlp_model_type = 'asym_pearson'; end
            if ~isfield(settings, 'pearson_m_bounds'), settings.pearson_m_bounds = [1.60, 16]; end
            if ~isfield(settings, 'pvoigt_eta_bounds'), settings.pvoigt_eta_bounds = [0.02, 0.65]; end
            if ~isfield(settings, 'zlp_center_limit'), settings.zlp_center_limit = 3.0; end
            if ~isfield(settings, 'zlp_width_bounds'), settings.zlp_width_bounds = [0.18, 42]; end
            if ~isfield(settings, 'zlp_component_width_max_meV'), settings.zlp_component_width_max_meV = [3.3, 7.0, 18, 58]; end
            if ~isfield(settings, 'zlp_component_width_min_meV'), settings.zlp_component_width_min_meV = [0.18, 0.35, 0.8, 4.0]; end
            if ~isfield(settings, 'zlp_wide_component_l2'), settings.zlp_wide_component_l2 = 3e-3; end
            if ~isfield(settings, 'zlp_wide_component_start_index'), settings.zlp_wide_component_start_index = 3; end
            if ~isfield(settings, 'zlp_lowfreq_guard_enabled'), settings.zlp_lowfreq_guard_enabled = true; end
            if ~isfield(settings, 'zlp_lowfreq_guard_start_meV'), settings.zlp_lowfreq_guard_start_meV = 14; end
            if ~isfield(settings, 'zlp_lowfreq_guard_end_meV'), settings.zlp_lowfreq_guard_end_meV = 58; end
            if ~isfield(settings, 'zlp_lowfreq_guard_weight'), settings.zlp_lowfreq_guard_weight = 0.05; end
            if ~isfield(settings, 'zlp_lowfreq_guard_taper_meV'), settings.zlp_lowfreq_guard_taper_meV = 4; end
            if ~isfield(settings, 'zlp_signal_guard_enabled'), settings.zlp_signal_guard_enabled = true; end
            if ~isfield(settings, 'zlp_signal_guard_start_meV'), settings.zlp_signal_guard_start_meV = 16; end
            if ~isfield(settings, 'zlp_signal_guard_end_meV'), settings.zlp_signal_guard_end_meV = 120; end
            if ~isfield(settings, 'zlp_signal_guard_weight'), settings.zlp_signal_guard_weight = 0.20; end
            if ~isfield(settings, 'zlp_signal_guard_taper_meV'), settings.zlp_signal_guard_taper_meV = 5; end
            if ~isfield(settings, 'zlp_signal_guard_tail_release_start_meV'), settings.zlp_signal_guard_tail_release_start_meV = 95; end
            if ~isfield(settings, 'zlp_signal_guard_tail_release_width_meV'), settings.zlp_signal_guard_tail_release_width_meV = 10; end
            if ~isfield(settings, 'zlp_signal_guard_tail_floor'), settings.zlp_signal_guard_tail_floor = 0.05; end
            if ~isfield(settings, 'zlp_signal_guard_feature_weighted'), settings.zlp_signal_guard_feature_weighted = true; end
            if ~isfield(settings, 'zlp_signal_guard_min_profile'), settings.zlp_signal_guard_min_profile = 0.04; end
            if ~isfield(settings, 'zlp_signal_guard_feature_smooth_points'), settings.zlp_signal_guard_feature_smooth_points = 11; end
            if ~isfield(settings, 'adaptive_zlp_region_enabled'), settings.adaptive_zlp_region_enabled = true; end
            if ~isfield(settings, 'adaptive_zlp_region_adjust_weights'), settings.adaptive_zlp_region_adjust_weights = true; end
            if ~isfield(settings, 'adaptive_zlp_region_search_meV'), settings.adaptive_zlp_region_search_meV = 45; end
            if ~isfield(settings, 'adaptive_zlp_core_fraction'), settings.adaptive_zlp_core_fraction = 0.28; end
            if ~isfield(settings, 'adaptive_zlp_guard_margin_meV'), settings.adaptive_zlp_guard_margin_meV = 6; end
            if ~isfield(settings, 'adaptive_zlp_guard_min_start_meV'), settings.adaptive_zlp_guard_min_start_meV = 12; end
            if ~isfield(settings, 'adaptive_zlp_guard_max_start_meV'), settings.adaptive_zlp_guard_max_start_meV = 28; end
            if ~isfield(settings, 'adaptive_zlp_core_min_meV'), settings.adaptive_zlp_core_min_meV = 8; end
            if ~isfield(settings, 'adaptive_zlp_core_max_meV'), settings.adaptive_zlp_core_max_meV = 35; end
            if ~isfield(settings, 'adaptive_zlp_guard_end_fraction'), settings.adaptive_zlp_guard_end_fraction = 0.82; end
            if ~isfield(settings, 'eshift_range'), settings.eshift_range = [0, 0]; end
            if ~isfield(settings, 'stage1_zlp_align_enabled'), settings.stage1_zlp_align_enabled = false; end
            if ~isfield(settings, 'stage1_zlp_align_max_shift_meV'), settings.stage1_zlp_align_max_shift_meV = 2; end
            if ~isfield(settings, 'stage1_zlp_align_search_meV'), settings.stage1_zlp_align_search_meV = 12; end
            if ~isfield(settings, 'stage1_zlp_align_smooth_points'), settings.stage1_zlp_align_smooth_points = 5; end
            if ~isfield(settings, 'stage1_zlp_align_min_peak_fraction'), settings.stage1_zlp_align_min_peak_fraction = 0.25; end
            if ~isfield(settings, 'stage1_zlp_core_fartail_enabled'), settings.stage1_zlp_core_fartail_enabled = true; end
            if ~isfield(settings, 'stage1_zlp_core_anchor_width_meV'), settings.stage1_zlp_core_anchor_width_meV = 10; end
            if ~isfield(settings, 'stage1_zlp_core_data_weight'), settings.stage1_zlp_core_data_weight = 0.38; end
            if ~isfield(settings, 'stage1_zlp_core_anchor_weight'), settings.stage1_zlp_core_anchor_weight = 0.055; end
            if ~isfield(settings, 'stage1_zlp_fartail_start_meV'), settings.stage1_zlp_fartail_start_meV = 145; end
            if ~isfield(settings, 'stage1_zlp_fartail_end_meV'), settings.stage1_zlp_fartail_end_meV = 280; end
            if ~isfield(settings, 'stage1_zlp_fartail_taper_meV'), settings.stage1_zlp_fartail_taper_meV = 14; end
            if ~isfield(settings, 'stage1_zlp_fartail_data_weight'), settings.stage1_zlp_fartail_data_weight = 0.30; end
            if ~isfield(settings, 'stage1_zlp_fartail_anchor_weight'), settings.stage1_zlp_fartail_anchor_weight = 0.007; end
            if ~isfield(settings, 'stage1_zlp_fartail_anchor_side'), settings.stage1_zlp_fartail_anchor_side = 'gain'; end
            if ~isfield(settings, 'stage1_zlp_fartail_gain_boost'), settings.stage1_zlp_fartail_gain_boost = 1.25; end
            if ~isfield(settings, 'stage1_zlp_gain_tail_prior_enabled'), settings.stage1_zlp_gain_tail_prior_enabled = true; end
            if ~isfield(settings, 'stage1_zlp_gain_tail_prior_weight'), settings.stage1_zlp_gain_tail_prior_weight = 0.002; end
            if ~isfield(settings, 'stage1_zlp_gain_tail_prior_ratio'), settings.stage1_zlp_gain_tail_prior_ratio = 1.18; end
            if ~isfield(settings, 'stage1_zlp_gain_tail_pair_step_meV'), settings.stage1_zlp_gain_tail_pair_step_meV = 8; end
            if ~isfield(settings, 'stage1_zlp_gain_asym_shape_prior_enabled'), settings.stage1_zlp_gain_asym_shape_prior_enabled = true; end
            if ~isfield(settings, 'stage1_zlp_gain_asym_shape_prior_weight'), settings.stage1_zlp_gain_asym_shape_prior_weight = 0.008; end
            if ~isfield(settings, 'stage1_zlp_gain_asym_prior_start_index'), settings.stage1_zlp_gain_asym_prior_start_index = 2; end
            if ~isfield(settings, 'stage1_zlp_gain_asym_width_ratio'), settings.stage1_zlp_gain_asym_width_ratio = 1.06; end
            if ~isfield(settings, 'stage1_zlp_gain_asym_shape_ratio'), settings.stage1_zlp_gain_asym_shape_ratio = 1.04; end
            if ~isfield(settings, 'stage1_zlp_gain_heavy_init_enabled'), settings.stage1_zlp_gain_heavy_init_enabled = true; end
            if ~isfield(settings, 'stage1_zlp_gain_width_boost'), settings.stage1_zlp_gain_width_boost = 1.45; end
            if ~isfield(settings, 'stage1_zlp_loss_width_scale'), settings.stage1_zlp_loss_width_scale = 0.92; end
            if ~isfield(settings, 'stage1_zlp_gain_shape_heavy_factor'), settings.stage1_zlp_gain_shape_heavy_factor = 0.86; end
            if ~isfield(settings, 'stage1_zlp_loss_shape_light_factor'), settings.stage1_zlp_loss_shape_light_factor = 1.12; end
            settings.use_bg = false;
            settings.bg_order = 0;
            if ~isfield(settings, 'em_min_meV'), settings.em_min_meV = 8; end
            if ~isfield(settings, 'em_max_meV'), settings.em_max_meV = []; end
            if ~isfield(settings, 'em_step_meV'), settings.em_step_meV = 6; end
            if ~isfield(settings, 'em_lowfreq_dense_enabled'), settings.em_lowfreq_dense_enabled = true; end
            if ~isfield(settings, 'em_lowfreq_dense_range_meV'), settings.em_lowfreq_dense_range_meV = [8, 72]; end
            if ~isfield(settings, 'em_lowfreq_dense_step_factor'), settings.em_lowfreq_dense_step_factor = 0.55; end
            if ~isfield(settings, 'em_lowfreq_dense_max_step_meV'), settings.em_lowfreq_dense_max_step_meV = 4.0; end
            if ~isfield(settings, 'em_exclude_core_meV'), settings.em_exclude_core_meV = 6; end
            if ~isfield(settings, 'stage1_fixed_EM'), settings.stage1_fixed_EM = []; end
            if ~isfield(settings, 'datacube_stage1_auto_sparse_pm_enabled'), settings.datacube_stage1_auto_sparse_pm_enabled = true; end
            if ~isfield(settings, 'stage1_auto_sparse_pm_min_count'), settings.stage1_auto_sparse_pm_min_count = 18; end
            if ~isfield(settings, 'stage1_auto_sparse_pm_max_count'), settings.stage1_auto_sparse_pm_max_count = 56; end
            if ~isfield(settings, 'stage1_auto_sparse_pm_step'), settings.stage1_auto_sparse_pm_step = 6; end
            if ~isfield(settings, 'stage1_auto_sparse_pm_R2_tolerance'), settings.stage1_auto_sparse_pm_R2_tolerance = 6e-4; end
            if ~isfield(settings, 'stage1_auto_sparse_pm_min_R2'), settings.stage1_auto_sparse_pm_min_R2 = 0.992; end
            if ~isfield(settings, 'stage1_auto_sparse_pm_min_spacing_meV'), settings.stage1_auto_sparse_pm_min_spacing_meV = 3.2; end
            if ~isfield(settings, 'stage1_auto_sparse_pm_feature_threshold'), settings.stage1_auto_sparse_pm_feature_threshold = 0.10; end
            if ~isfield(settings, 'stage1_auto_sparse_pm_protect_max_centers'), settings.stage1_auto_sparse_pm_protect_max_centers = 20; end
            if ~isfield(settings, 'stage1_auto_sparse_pm_probe_pixel'), settings.stage1_auto_sparse_pm_probe_pixel = []; end
            if ~isfield(settings, 'kernel_mode'), settings.kernel_mode = 'fitted_zlp'; end
            if ~isfield(settings, 'lambda_be_soft'), settings.lambda_be_soft = 1e-3; end
            if ~isfield(settings, 'lambda_pm_smooth'), settings.lambda_pm_smooth = 3e-3; end
            if ~isfield(settings, 'lambda_pm_l2'), settings.lambda_pm_l2 = 2e-5; end
            if ~isfield(settings, 'lambda_pm_l1'), settings.lambda_pm_l1 = 3e-5; end
            settings.lambda_bg_l2 = 0;
            if ~isfield(settings, 'inner_solver'), settings.inner_solver = 'auto'; end
            if ~isfield(settings, 'nonlinear_inner_solver'), settings.nonlinear_inner_solver = 'lsqlin'; end
            if ~isfield(settings, 'nonlinear_tail_overshoot_constraint'), settings.nonlinear_tail_overshoot_constraint = false; end
            if ~isfield(settings, 'inner_solver_max_iter'), settings.inner_solver_max_iter = 300; end
            if ~isfield(settings, 'fast_nnls_refine_passes'), settings.fast_nnls_refine_passes = 2; end
            if ~isfield(settings, 'tail_overshoot_constraint'), settings.tail_overshoot_constraint = true; end
            if ~isfield(settings, 'tail_overshoot_start_meV'), settings.tail_overshoot_start_meV = 180; end
            if ~isfield(settings, 'tail_overshoot_noise_sigma'), settings.tail_overshoot_noise_sigma = 4; end
            if ~isfield(settings, 'tail_overshoot_rel_margin'), settings.tail_overshoot_rel_margin = 0.018; end
            if ~isfield(settings, 'tail_overshoot_abs_margin'), settings.tail_overshoot_abs_margin = 2e-4; end
            if ~isfield(settings, 'T_soft_init'), settings.T_soft_init = 450; end
            if ~isfield(settings, 'T_soft_bounds'), settings.T_soft_bounds = [80, 2500]; end
            if ~isfield(settings, 'core_downweight_width'), settings.core_downweight_width = 12; end
            if ~isfield(settings, 'core_weight'), settings.core_weight = 0.22; end
            if ~isfield(settings, 'shoulder_weight'), settings.shoulder_weight = 0.72; end
            if ~isfield(settings, 'signal_weight'), settings.signal_weight = 1.45; end
            if ~isfield(settings, 'highfreq_weight'), settings.highfreq_weight = 2.2; end
            if ~isfield(settings, 'highfreq_start_meV'), settings.highfreq_start_meV = 70; end
            if ~isfield(settings, 'tail_weight'), settings.tail_weight = 2.5; end
            if ~isfield(settings, 'tail_start_fraction'), settings.tail_start_fraction = 0.62; end
            if ~isfield(settings, 'dynamic_weight_strength'), settings.dynamic_weight_strength = 0.8; end
            if ~isfield(settings, 'dynamic_weight_max_gain'), settings.dynamic_weight_max_gain = 10; end
            if ~isfield(settings, 'dynamic_weight_noise_sigma'), settings.dynamic_weight_noise_sigma = 4; end
            if ~isfield(settings, 'max_nonlinear_iter'), settings.max_nonlinear_iter = 45; end
            if ~isfield(settings, 'stage1_early_stop_enabled'), settings.stage1_early_stop_enabled = true; end
            if ~isfield(settings, 'stage1_early_stop_min_rounds'), settings.stage1_early_stop_min_rounds = 2; end
            if ~isfield(settings, 'stage1_early_stop_R2'), settings.stage1_early_stop_R2 = 0.99995; end
            if ~isfield(settings, 'stage1_early_stop_weighted_NRMSE'), settings.stage1_early_stop_weighted_NRMSE = 0.0025; end
            if ~isfield(settings, 'stage1_early_stop_highfreq_NRMSE'), settings.stage1_early_stop_highfreq_NRMSE = 0.0085; end
            if ~isfield(settings, 'stage1_warm_start_skip_coarse'), settings.stage1_warm_start_skip_coarse = true; end
            if ~isfield(settings, 'stage1_strict_polish_enabled'), settings.stage1_strict_polish_enabled = true; end
            if ~isfield(settings, 'stage1_strict_polish_iter'), settings.stage1_strict_polish_iter = 7; end
            if ~isfield(settings, 'stage1_strict_polish_inner_solver'), settings.stage1_strict_polish_inner_solver = 'lsqnonneg'; end
            if ~isfield(settings, 'stage1_strict_polish_accept_factor'), settings.stage1_strict_polish_accept_factor = 1.03; end
            if ~isfield(settings, 'initial_theta'), settings.initial_theta = []; end
            if ~isfield(settings, 'optimization_max_points'), settings.optimization_max_points = Inf; end
            if ~isfield(settings, 'optimization_keep_core_meV'), settings.optimization_keep_core_meV = 35; end
            if ~isfield(settings, 'silent'), settings.silent = true; end
            if ~isfield(settings, 'save_plots'), settings.save_plots = true; end
            if ~isfield(settings, 'save_diagnostic_plots'), settings.save_diagnostic_plots = true; end
            if ~isfield(settings, 'max_saved_spectrum_plots'), settings.max_saved_spectrum_plots = 8; end
            if ~isfield(settings, 'force_hidden_figures'), settings.force_hidden_figures = false; end
            if ~isfield(settings, 'close_parallel_pool_after_run'), settings.close_parallel_pool_after_run = false; end
            if ~isfield(settings, 'ci_window_points'), settings.ci_window_points = 31; end
            if ~isfield(settings, 'ci_min_noise_fraction'), settings.ci_min_noise_fraction = 2e-4; end
            if ~isfield(settings, 'ci_level'), settings.ci_level = 0.95; end
            if ~isfield(settings, 'run_stage2'), settings.run_stage2 = false; end
            if ~isfield(settings, 'stage2_window_meV'), settings.stage2_window_meV = [40, 190]; end
            if ~isfield(settings, 'stage2_attention_window_meV'), settings.stage2_attention_window_meV = [40, 230]; end
            if ~isfield(settings, 'stage2_temperature_algorithm'), settings.stage2_temperature_algorithm = 'grid_scan'; end
            if ~isfield(settings, 'stage2_T_grid'), settings.stage2_T_grid = 80:20:1800; end
            if ~isfield(settings, 'stage2_scan_inner_solver'), settings.stage2_scan_inner_solver = 'lsqnonneg'; end
            if ~isfield(settings, 'stage2_scan_tail_overshoot_constraint'), settings.stage2_scan_tail_overshoot_constraint = false; end
            if ~isfield(settings, 'stage2_adaptive_grid_scan_enabled'), settings.stage2_adaptive_grid_scan_enabled = true; end
            if ~isfield(settings, 'stage2_adaptive_grid_min_full_points'), settings.stage2_adaptive_grid_min_full_points = 25; end
            if ~isfield(settings, 'stage2_adaptive_grid_coarse_step_K'), settings.stage2_adaptive_grid_coarse_step_K = 100; end
            if ~isfield(settings, 'stage2_adaptive_grid_fine_half_width_K'), settings.stage2_adaptive_grid_fine_half_width_K = 160; end
            if ~isfield(settings, 'stage2_temperature_refine_enabled'), settings.stage2_temperature_refine_enabled = true; end
            if ~isfield(settings, 'stage2_temperature_refine_half_width_K'), settings.stage2_temperature_refine_half_width_K = 20; end
            if ~isfield(settings, 'stage2_temperature_refine_initial_step_K'), settings.stage2_temperature_refine_initial_step_K = 10; end
            if ~isfield(settings, 'stage2_temperature_refine_tol_K'), settings.stage2_temperature_refine_tol_K = 1; end
            if ~isfield(settings, 'stage2_temperature_refine_max_evals'), settings.stage2_temperature_refine_max_evals = 8; end
            if ~isfield(settings, 'stage2_pm_smooth_weight'), settings.stage2_pm_smooth_weight = 1.2e-3; end
            if ~isfield(settings, 'stage2_pm_smooth_max_gap_meV'), settings.stage2_pm_smooth_max_gap_meV = 6; end
            if ~isfield(settings, 'stage2_attention_hardbe_only_enabled'), settings.stage2_attention_hardbe_only_enabled = true; end
            if ~isfield(settings, 'stage2_hard_be_pair_threshold'), settings.stage2_hard_be_pair_threshold = 0.35; end
            if ~isfield(settings, 'stage2_hard_be_pair_floor'), settings.stage2_hard_be_pair_floor = 0.10; end
            if ~isfield(settings, 'stage2_hard_dense_em_enabled'), settings.stage2_hard_dense_em_enabled = false; end
            if ~isfield(settings, 'stage2_hard_dense_min_meV'), settings.stage2_hard_dense_min_meV = []; end
            if ~isfield(settings, 'stage2_hard_dense_step_meV'), settings.stage2_hard_dense_step_meV = 4.0; end
            if ~isfield(settings, 'stage2_hard_dense_extra_meV'), settings.stage2_hard_dense_extra_meV = 0; end
            if ~isfield(settings, 'stage2_hard_dense_max_meV'), settings.stage2_hard_dense_max_meV = []; end
            if ~isfield(settings, 'stage2_soft_dense_em_enabled'), settings.stage2_soft_dense_em_enabled = true; end
            if ~isfield(settings, 'stage2_soft_dense_min_meV'), settings.stage2_soft_dense_min_meV = []; end
            if ~isfield(settings, 'stage2_soft_dense_step_meV'), settings.stage2_soft_dense_step_meV = 4.0; end
            if ~isfield(settings, 'stage2_soft_dense_extra_meV'), settings.stage2_soft_dense_extra_meV = 30; end
            if ~isfield(settings, 'stage2_soft_dense_max_meV'), settings.stage2_soft_dense_max_meV = []; end
            if ~isfield(settings, 'stage2_soft_dense_hard_exclusion_meV'), settings.stage2_soft_dense_hard_exclusion_meV = 1.25; end
            if ~isfield(settings, 'stage2_soft_dense_pair_score_max'), settings.stage2_soft_dense_pair_score_max = settings.stage2_hard_be_pair_threshold; end
            if ~isfield(settings, 'stage2_soft_pm_smooth_weight'), settings.stage2_soft_pm_smooth_weight = 0; end
            if ~isfield(settings, 'stage2_soft_pm_l2_weight'), settings.stage2_soft_pm_l2_weight = 1e-8; end
            if ~isfield(settings, 'stage2_soft_pm_stage1_prior_weight'), settings.stage2_soft_pm_stage1_prior_weight = 0; end
            if ~isfield(settings, 'stage2_soft_be_weight'), settings.stage2_soft_be_weight = 0; end
            if ~isfield(settings, 'stage2_pm_ratio_plot_hide_flyers'), settings.stage2_pm_ratio_plot_hide_flyers = true; end
            if ~isfield(settings, 'stage2_pm_ratio_plot_flyer_min_residual'), settings.stage2_pm_ratio_plot_flyer_min_residual = 1.25; end
            if ~isfield(settings, 'stage2_pm_ratio_plot_flyer_weight_fraction'), settings.stage2_pm_ratio_plot_flyer_weight_fraction = 0.30; end
            if ~isfield(settings, 'stage2_regularization_cost_weight'), settings.stage2_regularization_cost_weight = 1.0; end
            if ~isfield(settings, 'stage2_far_tail_downweight_enabled'), settings.stage2_far_tail_downweight_enabled = true; end
            if ~isfield(settings, 'stage2_far_tail_weight_scale'), settings.stage2_far_tail_weight_scale = 0.12; end
            if ~isfield(settings, 'stage2_far_tail_downweight_start_meV'), settings.stage2_far_tail_downweight_start_meV = 220; end
            if ~isfield(settings, 'stage2_far_tail_downweight_width_meV'), settings.stage2_far_tail_downweight_width_meV = 16; end
            if ~isfield(settings, 'stage2_far_tail_loss_budget_downweight_enabled'), settings.stage2_far_tail_loss_budget_downweight_enabled = true; end
            if ~isfield(settings, 'stage2_far_tail_loss_budget_scale'), settings.stage2_far_tail_loss_budget_scale = 0.25; end
            if ~isfield(settings, 'stage2_attention_enabled'), settings.stage2_attention_enabled = true; end
            if ~isfield(settings, 'stage2_attention_floor'), settings.stage2_attention_floor = 0.18; end
            if ~isfield(settings, 'stage2_attention_gain'), settings.stage2_attention_gain = 9; end
            if ~isfield(settings, 'stage2_attention_peak_sigma_meV'), settings.stage2_attention_peak_sigma_meV = 5; end
            if ~isfield(settings, 'stage2_attention_threshold'), settings.stage2_attention_threshold = 0.18; end
            if ~isfield(settings, 'stage2_attention_fit_enabled'), settings.stage2_attention_fit_enabled = true; end
            if ~isfield(settings, 'stage2_attention_fit_threshold'), settings.stage2_attention_fit_threshold = 0.10; end
            if ~isfield(settings, 'stage2_attention_fit_expand_meV'), settings.stage2_attention_fit_expand_meV = 8; end
            if ~isfield(settings, 'stage2_attention_tail_enabled'), settings.stage2_attention_tail_enabled = true; end
            if ~isfield(settings, 'stage2_attention_tail_min_abs_meV'), settings.stage2_attention_tail_min_abs_meV = 70; end
            if ~isfield(settings, 'stage2_attention_tail_signal_threshold'), settings.stage2_attention_tail_signal_threshold = 0.10; end
            if ~isfield(settings, 'stage2_attention_tail_soft_width'), settings.stage2_attention_tail_soft_width = 0.05; end
            if ~isfield(settings, 'stage2_attention_tail_pair_release_threshold'), settings.stage2_attention_tail_pair_release_threshold = 0.10; end
            if ~isfield(settings, 'stage2_allow_aggressive_antistokes'), settings.stage2_allow_aggressive_antistokes = false; end
            if ~isfield(settings, 'stage2_antistokes_attention_boost'), settings.stage2_antistokes_attention_boost = 1.35; end
            if ~isfield(settings, 'stage2_antistokes_fit_weight_boost'), settings.stage2_antistokes_fit_weight_boost = 1.45; end
            if ~isfield(settings, 'stage2_antistokes_noise_sigma'), settings.stage2_antistokes_noise_sigma = 1.8; end
            if ~logical(settings.stage2_allow_aggressive_antistokes)
                settings.stage2_antistokes_attention_boost = min(settings.stage2_antistokes_attention_boost, 1.35);
                settings.stage2_antistokes_fit_weight_boost = min(settings.stage2_antistokes_fit_weight_boost, 1.45);
                settings.stage2_antistokes_noise_sigma = max(settings.stage2_antistokes_noise_sigma, 1.8);
            end
            if ~isfield(settings, 'stage2_fit_sigmoid_taper_meV'), settings.stage2_fit_sigmoid_taper_meV = 10; end
            if ~isfield(settings, 'stage2_fit_sigmoid_min_weight'), settings.stage2_fit_sigmoid_min_weight = 0.08; end
            if ~isfield(settings, 'stage2_transition_fit_enabled'), settings.stage2_transition_fit_enabled = false; end
            if ~isfield(settings, 'stage2_transition_min_meV'), settings.stage2_transition_min_meV = 14; end
            if ~isfield(settings, 'stage2_zlp_refine_enabled'), settings.stage2_zlp_refine_enabled = true; end
            if ~isfield(settings, 'stage2_zlp_refine_scale_enabled'), settings.stage2_zlp_refine_scale_enabled = true; end
            if ~isfield(settings, 'stage2_zlp_refine_tilt_enabled'), settings.stage2_zlp_refine_tilt_enabled = true; end
            if ~isfield(settings, 'stage2_zlp_refine_shift_enabled'), settings.stage2_zlp_refine_shift_enabled = true; end
            if ~isfield(settings, 'stage2_zlp_refine_width_enabled'), settings.stage2_zlp_refine_width_enabled = true; end
            if ~isfield(settings, 'stage2_zlp_refine_tail_enabled'), settings.stage2_zlp_refine_tail_enabled = true; end
            if ~isfield(settings, 'stage2_zlp_refine_tail_asym_enabled'), settings.stage2_zlp_refine_tail_asym_enabled = true; end
            if ~isfield(settings, 'stage2_zlp_refine_scale_bound'), settings.stage2_zlp_refine_scale_bound = 0.10; end
            if ~isfield(settings, 'stage2_zlp_refine_tilt_bound'), settings.stage2_zlp_refine_tilt_bound = 0.15; end
            if ~isfield(settings, 'stage2_zlp_refine_shift_bound_meV'), settings.stage2_zlp_refine_shift_bound_meV = 2.00; end
            if ~isfield(settings, 'stage2_zlp_refine_width_bound'), settings.stage2_zlp_refine_width_bound = 0.15; end
            if ~isfield(settings, 'stage2_zlp_refine_tail_bound'), settings.stage2_zlp_refine_tail_bound = 0.20; end
            if ~isfield(settings, 'stage2_zlp_refine_tail_asym_bound'), settings.stage2_zlp_refine_tail_asym_bound = 0.20; end
            if ~isfield(settings, 'stage2_zlp_refine_scale_l2'), settings.stage2_zlp_refine_scale_l2 = 1.0; end
            if ~isfield(settings, 'stage2_zlp_refine_tilt_l2'), settings.stage2_zlp_refine_tilt_l2 = 1.0; end
            if ~isfield(settings, 'stage2_zlp_refine_shift_l2'), settings.stage2_zlp_refine_shift_l2 = 0.20; end
            if ~isfield(settings, 'stage2_zlp_refine_width_l2'), settings.stage2_zlp_refine_width_l2 = 0.60; end
            if ~isfield(settings, 'stage2_zlp_refine_tail_l2'), settings.stage2_zlp_refine_tail_l2 = 0.60; end
            if ~isfield(settings, 'stage2_zlp_refine_tail_asym_l2'), settings.stage2_zlp_refine_tail_asym_l2 = 0.60; end
            if ~isfield(settings, 'stage2_zlp_refine_tilt_scale_meV'), settings.stage2_zlp_refine_tilt_scale_meV = 65; end
            if ~isfield(settings, 'stage2_zlp_refine_tail_scale_meV'), settings.stage2_zlp_refine_tail_scale_meV = 75; end
            if ~isfield(settings, 'stage2_ratio_min_pairs'), settings.stage2_ratio_min_pairs = 5; end
            if ~isfield(settings, 'stage2_ratio_reliable_spread_K'), settings.stage2_ratio_reliable_spread_K = 420; end
            if ~isfield(settings, 'stage2_ratio_cost_weight_sum_ref'), settings.stage2_ratio_cost_weight_sum_ref = 0.35; end
            if ~isfield(settings, 'stage2_ratio_best_line_enabled'), settings.stage2_ratio_best_line_enabled = true; end
            if ~isfield(settings, 'stage2_ratio_best_line_min_points'), settings.stage2_ratio_best_line_min_points = 6; end
            if ~isfield(settings, 'stage2_ratio_best_line_min_span_meV'), settings.stage2_ratio_best_line_min_span_meV = 22; end
            if ~isfield(settings, 'stage2_ratio_best_line_ref_span_meV'), settings.stage2_ratio_best_line_ref_span_meV = 55; end
            if ~isfield(settings, 'stage2_ratio_best_line_min_R2'), settings.stage2_ratio_best_line_min_R2 = 0.955; end
            if ~isfield(settings, 'stage2_ratio_best_line_edge_points'), settings.stage2_ratio_best_line_edge_points = 3; end
            if ~isfield(settings, 'stage2_ratio_best_line_right_edge_median_residual_max'), settings.stage2_ratio_best_line_right_edge_median_residual_max = 0.55; end
            if ~isfield(settings, 'stage2_ratio_best_line_right_edge_max_residual_max'), settings.stage2_ratio_best_line_right_edge_max_residual_max = 1.10; end
            if ~isfield(settings, 'stage2_ratio_best_line_allow_all_fallback'), settings.stage2_ratio_best_line_allow_all_fallback = false; end
            if ~isfield(settings, 'stage2_report_ratio_source'), settings.stage2_report_ratio_source = 'stage1_pm_ratio'; end
            if ~isfield(settings, 'stage2_report_ratio_require_reliable'), settings.stage2_report_ratio_require_reliable = true; end
            if ~isfield(settings, 'stage2_boundary_margin_K'), settings.stage2_boundary_margin_K = 40; end
            if ~isfield(settings, 'stage2_outlier_refit_top_loss_count'), settings.stage2_outlier_refit_top_loss_count = 3; end
            if ~isfield(settings, 'stage2_outlier_refit_loss_weight_boost'), settings.stage2_outlier_refit_loss_weight_boost = 1.35; end
            if ~isfield(settings, 'stage2_outlier_refit_attention_boost'), settings.stage2_outlier_refit_attention_boost = 1.25; end
            if ~isfield(settings, 'stage2_outlier_refit_l1_boost'), settings.stage2_outlier_refit_l1_boost = 1.45; end
            if ~isfield(settings, 'stage2_outlier_refit_l2_boost'), settings.stage2_outlier_refit_l2_boost = 2.0; end
            if ~isfield(settings, 'stage2_outlier_refit_smooth_boost'), settings.stage2_outlier_refit_smooth_boost = 2.2; end
            if ~isfield(settings, 'stage2_outlier_refit_low_R2_extra_boost'), settings.stage2_outlier_refit_low_R2_extra_boost = 1.35; end
            if ~isfield(settings, 'stage2_outlier_refit_max_l1'), settings.stage2_outlier_refit_max_l1 = 2.5e-5; end
            if ~isfield(settings, 'stage2_outlier_refit_max_l2'), settings.stage2_outlier_refit_max_l2 = 5.0e-6; end
            if ~isfield(settings, 'stage2_outlier_refit_max_smooth'), settings.stage2_outlier_refit_max_smooth = 1.8e-3; end
            if ~isfield(settings, 'test_output_dir'), settings.test_output_dir = fullfile('pearson_sparse_fullfit', 'test_outputs'); end
            if ~isfield(settings, 'test_indices'), settings.test_indices = []; end
            if ~isfield(settings, 'test_max_spectra'), settings.test_max_spectra = 6; end
            if ~isfield(settings, 'datacube_sample_pixels'), settings.datacube_sample_pixels = 24; end
            if ~isfield(settings, 'datacube_sample_strategy'), settings.datacube_sample_strategy = 'uniform_grid'; end
            if ~isfield(settings, 'datacube_random_seed'), settings.datacube_random_seed = []; end
            if ~isfield(settings, 'datacube_use_warm_start'), settings.datacube_use_warm_start = true; end
            if ~isfield(settings, 'datacube_parallel_enabled'), settings.datacube_parallel_enabled = false; end
            if ~isfield(settings, 'datacube_parallel_workers'), settings.datacube_parallel_workers = 0; end
            if ~isfield(settings, 'datacube_parallel_min_pixels'), settings.datacube_parallel_min_pixels = 8; end
            if ~isfield(settings, 'datacube_parallel_chunks'), settings.datacube_parallel_chunks = 0; end
            if ~isfield(settings, 'datacube_parallel_warm_start_within_chunk'), settings.datacube_parallel_warm_start_within_chunk = false; end
            if ~isfield(settings, 'datacube_fit_retries'), settings.datacube_fit_retries = 1; end
            if ~isfield(settings, 'datacube_retry_with_strict_nnls'), settings.datacube_retry_with_strict_nnls = true; end
            if ~isfield(settings, 'datacube_checkpoint_enabled'), settings.datacube_checkpoint_enabled = false; end
            if ~isfield(settings, 'datacube_resume_checkpoint'), settings.datacube_resume_checkpoint = false; end
            if ~isfield(settings, 'datacube_checkpoint_every'), settings.datacube_checkpoint_every = 25; end
            if ~isfield(settings, 'datacube_benchmark_repeats'), settings.datacube_benchmark_repeats = 1; end
            if ~isfield(settings, 'datacube_progress_every'), settings.datacube_progress_every = 25; end
            if ~isfield(settings, 'datacube_progress_callback'), settings.datacube_progress_callback = []; end
            if ~isfield(settings, 'datacube_stop_callback'), settings.datacube_stop_callback = []; end
            if ~isfield(settings, 'pm_deconv_min_em_meV'), settings.pm_deconv_min_em_meV = 8; end
            if ~isfield(settings, 'pm_deconv_max_em_meV'), settings.pm_deconv_max_em_meV = 180; end
            if ~isfield(settings, 'pm_deconv_min_pm_fraction'), settings.pm_deconv_min_pm_fraction = 0.025; end
            if ~isfield(settings, 'pm_deconv_max_centers'), settings.pm_deconv_max_centers = 36; end
            if ~isfield(settings, 'pm_deconv_min_spacing_meV'), settings.pm_deconv_min_spacing_meV = 1.2; end
            if ~isfield(settings, 'pm_deconv_data_weight_floor'), settings.pm_deconv_data_weight_floor = 0.15; end
            if ~isfield(settings, 'pm_deconv_pm_prior_weight'), settings.pm_deconv_pm_prior_weight = 0.04; end
            if ~isfield(settings, 'pm_deconv_be_prior_weight'), settings.pm_deconv_be_prior_weight = 0.025; end
            if ~isfield(settings, 'pm_deconv_l2_weight'), settings.pm_deconv_l2_weight = 2e-5; end
            if ~isfield(settings, 'pm_deconv_kernel_width_scale'), settings.pm_deconv_kernel_width_scale = 1; end
            if ~isfield(settings, 'pm_interp_min_em_meV'), settings.pm_interp_min_em_meV = settings.pm_deconv_min_em_meV; end
            if ~isfield(settings, 'pm_interp_max_em_meV'), settings.pm_interp_max_em_meV = settings.pm_deconv_max_em_meV; end
            if ~isfield(settings, 'pm_interp_min_pm_fraction'), settings.pm_interp_min_pm_fraction = settings.pm_deconv_min_pm_fraction; end
            if ~isfield(settings, 'pm_interp_max_centers'), settings.pm_interp_max_centers = 60; end
            if ~isfield(settings, 'pm_interp_min_spacing_meV'), settings.pm_interp_min_spacing_meV = 0.8; end
            if ~isfield(settings, 'pm_interp_peak_sigma_meV'), settings.pm_interp_peak_sigma_meV = 1.5; end
            if ~isfield(settings, 'pm_interp_method'), settings.pm_interp_method = 'linear'; end
            if ~isfield(settings, 'pm_interp_data_weight_floor'), settings.pm_interp_data_weight_floor = 0.15; end
            if ~isfield(settings, 'pm_interp_scale_mode'), settings.pm_interp_scale_mode = 'global_nnls'; end
            if ~isfield(settings, 'rl_deconv_min_energy_meV'), settings.rl_deconv_min_energy_meV = settings.pm_deconv_min_em_meV; end
            if ~isfield(settings, 'rl_deconv_max_energy_meV'), settings.rl_deconv_max_energy_meV = settings.pm_deconv_max_em_meV; end
            if ~isfield(settings, 'rl_deconv_iterations'), settings.rl_deconv_iterations = 60; end
            if ~isfield(settings, 'rl_deconv_kernel_width_scale'), settings.rl_deconv_kernel_width_scale = 1; end
            if ~isfield(settings, 'rl_deconv_edge_taper_meV'), settings.rl_deconv_edge_taper_meV = 6; end
            if ~isfield(settings, 'rl_deconv_latent_smooth_sigma_meV'), settings.rl_deconv_latent_smooth_sigma_meV = 0.35; end
            if ~isfield(settings, 'rl_deconv_floor_fraction'), settings.rl_deconv_floor_fraction = 1e-6; end
            if ~isfield(settings, 'rl_deconv_data_weight_floor'), settings.rl_deconv_data_weight_floor = 0.15; end
            if ~isfield(settings, 'rl_deconv_final_scale'), settings.rl_deconv_final_scale = true; end
            if ~isfield(settings, 'datacube_y_bin'), settings.datacube_y_bin = 1; end
            if ~isfield(settings, 'datacube_x_bin'), settings.datacube_x_bin = 1; end
            if ~isfield(settings, 'datacube_energy_bin'), settings.datacube_energy_bin = 1; end
            if isfield(settings, 'spatial_bin')
                settings = EELS_PearsonSparse_FullFit.apply_spatial_bin_alias(settings, settings.spatial_bin);
            end
            if isfield(settings, 'datacube_spatial_bin')
                settings = EELS_PearsonSparse_FullFit.apply_spatial_bin_alias(settings, settings.datacube_spatial_bin);
            end
            if isfield(settings, 'energy_bin')
                settings.datacube_energy_bin = settings.energy_bin;
            end
            if ~isfield(settings, 'spatial_neighbor_enabled'), settings.spatial_neighbor_enabled = true; end
            if ~isfield(settings, 'spatial_neighbor_weight'), settings.spatial_neighbor_weight = 0.25; end
            if ~isfield(settings, 'spatial_interface_attention_weight'), settings.spatial_interface_attention_weight = 0.65; end
            if ~isfield(settings, 'spatial_edge_sigma'), settings.spatial_edge_sigma = 0.18; end
            if ~isfield(settings, 'spatial_spectrum_smooth_weight'), settings.spatial_spectrum_smooth_weight = 0.20; end
            if ~isfield(settings, 'spatial_attention_profile_bins'), settings.spatial_attention_profile_bins = 12; end
            if ~isfield(settings, 'spatial_peak_shift_weight'), settings.spatial_peak_shift_weight = 0.26; end
            if ~isfield(settings, 'spatial_signal_peak_shift_weight'), settings.spatial_signal_peak_shift_weight = 0.34; end
            if ~isfield(settings, 'spatial_bad_pixel_z'), settings.spatial_bad_pixel_z = 3.5; end
            if ~isfield(settings, 'spatial_bad_R2_enabled'), settings.spatial_bad_R2_enabled = true; end
            if ~isfield(settings, 'spatial_bad_R2_abs_threshold'), settings.spatial_bad_R2_abs_threshold = 0.9985; end
            if ~isfield(settings, 'spatial_bad_R2_local_drop'), settings.spatial_bad_R2_local_drop = 3e-4; end
            if ~isfield(settings, 'spatial_refit_low_R2_priority_bonus'), settings.spatial_refit_low_R2_priority_bonus = 0.45; end
            if ~isfield(settings, 'spatial_refit_accept_min_R2_improvement'), settings.spatial_refit_accept_min_R2_improvement = 1.5e-4; end
            if ~isfield(settings, 'spatial_refit_accept_min_NRMSE_improvement'), settings.spatial_refit_accept_min_NRMSE_improvement = 8e-4; end
            if ~isfield(settings, 'spatial_T_floor_K'), settings.spatial_T_floor_K = 35; end
            if ~isfield(settings, 'spatial_reliability_floor'), settings.spatial_reliability_floor = 0.25; end
            if ~isfield(settings, 'spatial_repair_isolated'), settings.spatial_repair_isolated = true; end
            if ~isfield(settings, 'spatial_post_smooth_enabled'), settings.spatial_post_smooth_enabled = true; end
            if ~isfield(settings, 'spatial_post_smooth_iterations'), settings.spatial_post_smooth_iterations = 2; end
            if ~isfield(settings, 'spatial_post_smooth_weight'), settings.spatial_post_smooth_weight = 0.28; end
            if ~isfield(settings, 'spatial_post_smooth_loss_quantile'), settings.spatial_post_smooth_loss_quantile = 0.76; end
            if ~isfield(settings, 'spatial_post_smooth_reliability_floor'), settings.spatial_post_smooth_reliability_floor = 0.45; end
            if ~isfield(settings, 'spatial_post_smooth_interface_release'), settings.spatial_post_smooth_interface_release = 0.75; end
            if ~isfield(settings, 'spatial_refit_enabled'), settings.spatial_refit_enabled = false; end
            if ~isfield(settings, 'spatial_refit_iterations'), settings.spatial_refit_iterations = 2; end
            if ~isfield(settings, 'spatial_refit_loss_quantile'), settings.spatial_refit_loss_quantile = 0.78; end
            if ~isfield(settings, 'spatial_refit_max_fraction'), settings.spatial_refit_max_fraction = 0.12; end
            if ~isfield(settings, 'spatial_refit_max_pixels_per_iter'), settings.spatial_refit_max_pixels_per_iter = 60; end
            if ~isfield(settings, 'spatial_refit_min_temp_delta_K'), settings.spatial_refit_min_temp_delta_K = 85; end
            if ~isfield(settings, 'spatial_refit_low_profile_delta_max'), settings.spatial_refit_low_profile_delta_max = 0.32; end
            if ~isfield(settings, 'spatial_refit_peak_shift_max'), settings.spatial_refit_peak_shift_max = 0.20; end
            if ~isfield(settings, 'spatial_refit_max_R2_drop'), settings.spatial_refit_max_R2_drop = 2.5e-4; end
            if ~isfield(settings, 'spatial_refit_max_NRMSE_increase'), settings.spatial_refit_max_NRMSE_increase = 0.004; end
            if ~isfield(settings, 'spatial_refit_accept_min_improvement_K'), settings.spatial_refit_accept_min_improvement_K = 12; end
            if ~isfield(settings, 'spatial_region_refit_enabled'), settings.spatial_region_refit_enabled = true; end
            if ~isfield(settings, 'spatial_region_refit_min_pixels'), settings.spatial_region_refit_min_pixels = 8; end
            if ~isfield(settings, 'spatial_region_refit_weight'), settings.spatial_region_refit_weight = 0.55; end
            if ~isfield(settings, 'spatial_region_refit_delta_K'), settings.spatial_region_refit_delta_K = 80; end
            if ~isfield(settings, 'spatial_region_refit_profile_delta_max'), settings.spatial_region_refit_profile_delta_max = 0.38; end
            if ~isfield(settings, 'spatial_region_post_smooth_enabled'), settings.spatial_region_post_smooth_enabled = true; end
            if ~isfield(settings, 'spatial_region_post_smooth_iterations'), settings.spatial_region_post_smooth_iterations = 2; end
            if ~isfield(settings, 'spatial_region_post_smooth_weight'), settings.spatial_region_post_smooth_weight = 0.0; end
            if ~isfield(settings, 'spatial_region_post_smooth_outlier_weight'), settings.spatial_region_post_smooth_outlier_weight = 0.55; end
            if ~isfield(settings, 'spatial_region_post_smooth_z'), settings.spatial_region_post_smooth_z = 1.25; end
            if ~isfield(settings, 'spatial_region_post_smooth_max_step_K'), settings.spatial_region_post_smooth_max_step_K = 75; end
            if ~isfield(settings, 'spatial_region_post_smooth_min_region_p90_K'), settings.spatial_region_post_smooth_min_region_p90_K = 90; end
            if ~isfield(settings, 'spatial_region_post_smooth_min_global_p90_K'), settings.spatial_region_post_smooth_min_global_p90_K = 90; end
            if ~isfield(settings, 'stage1_speed_profile'), settings.stage1_speed_profile = 'v20_fast'; end
            if ~isfield(settings, 'stage1_rounds')
                switch lower(char(string(settings.stage1_speed_profile)))
                    case {'full', 'quality', 'v19_full', 'legacy_full'}
                        settings.stage1_rounds = [ ...
                            struct('name', 'round1_coarse', 'em_step_meV', 10, 'lambda_be_soft', 2e-3, 'max_nonlinear_iter', 35); ...
                            struct('name', 'round2_dense', 'em_step_meV', 6, 'lambda_be_soft', 1e-3, 'max_nonlinear_iter', 40); ...
                            struct('name', 'round3_dense_weakBE', 'em_step_meV', 4, 'lambda_be_soft', 5e-4, 'max_nonlinear_iter', 45)];
                    otherwise
                        settings.stage1_rounds = [ ...
                            struct('name', 'round1_mid_v20fast', 'em_step_meV', 8, 'lambda_be_soft', 1.5e-3, 'max_nonlinear_iter', 24); ...
                            struct('name', 'round2_dense_v20fast', 'em_step_meV', 4, 'lambda_be_soft', 5e-4, 'max_nonlinear_iter', 28)];
                end
            end
            settings.stage2_version = '2.0-attention-hardBE-softPM-stage2';
            settings.stage2_temperature_algorithm = 'grid_scan';
            settings.stage2_scan_inner_solver = 'lsqlin';
            settings.stage2_pm_smooth_weight = max(settings.stage2_pm_smooth_weight, 1.2e-3);
            if ~isfield(settings, 'stage2_soft_pm_smooth_weight') || ...
                    ~isfinite(settings.stage2_soft_pm_smooth_weight) || settings.stage2_soft_pm_smooth_weight < 0
                settings.stage2_soft_pm_smooth_weight = settings.stage2_pm_smooth_weight;
            end
            settings.stage2_soft_pm_l2_weight = max(settings.stage2_soft_pm_l2_weight, 0);
            settings.stage2_soft_pm_stage1_prior_weight = max(settings.stage2_soft_pm_stage1_prior_weight, 0);
            settings.stage2_soft_be_weight = max(settings.stage2_soft_be_weight, 0);
            if ~isfield(settings, 'stage2_cost_smooth_points'), settings.stage2_cost_smooth_points = 5; end
            if ~isfield(settings, 'stage2_profile_cost_weight'), settings.stage2_profile_cost_weight = 1; end
            if ~isfield(settings, 'stage2_joint_rmse_cost_weight'), settings.stage2_joint_rmse_cost_weight = 1.0; end
            if ~isfield(settings, 'stage2_attention_cost_weight'), settings.stage2_attention_cost_weight = 0.35; end
            if ~isfield(settings, 'stage2_global_joint_weight_floor'), settings.stage2_global_joint_weight_floor = 1.0; end
            if ~isfield(settings, 'stage2_fit_alpha_min'), settings.stage2_fit_alpha_min = 0.02; end
            if ~isfield(settings, 'stage2_effective_pair_weight_floor'), settings.stage2_effective_pair_weight_floor = 0.04; end
        end

        function settings = apply_spatial_bin_alias(settings, bin_value)
            if isempty(bin_value)
                return;
            end
            vals = double(bin_value(:).');
            vals = vals(isfinite(vals));
            if isempty(vals)
                return;
            end
            if isscalar(vals)
                vals = [vals, vals];
            end
            settings.datacube_y_bin = max(1, round(vals(1)));
            settings.datacube_x_bin = max(1, round(vals(2)));
        end

        function result = solve_spectrum(ene, spectrum, settings)
            if nargin < 3, settings = struct(); end
            result = EELS_PearsonSparse_FullFit.solve_stage1(ene, spectrum, settings);
        end

        function result = promote_stage1_to_stage2(stage1_result, settings)
            if nargin < 2, settings = struct(); end
            if ~(isstruct(stage1_result) && isfield(stage1_result, 'stage1') && isfield(stage1_result, 'raw_data_norm'))
                error('promote_stage1_to_stage2 expects a Stage-1 result produced by EELS_PearsonSparse_FullFit.');
            end
            base_settings = struct();
            if isfield(stage1_result, 'settings') && isstruct(stage1_result.settings)
                base_settings = stage1_result.settings;
            elseif isfield(stage1_result, 'original_settings') && isstruct(stage1_result.original_settings)
                base_settings = stage1_result.original_settings;
            end
            names = fieldnames(settings);
            for ii = 1:numel(names)
                base_settings.(names{ii}) = settings.(names{ii});
            end
            base_settings = EELS_PearsonSparse_FullFit.migrate_stage2_legacy_defaults(base_settings, names);
            base_settings.run_stage2 = true;
            base_settings = EELS_PearsonSparse_FullFit.defaults(base_settings);
            base_settings.datacube_progress_callback = [];

            result = stage1_result;
            result.settings = base_settings;
            result.temperature_enabled = true;
            result.analysis_mode = 'pearson_sparse_stage2_reused_stage1';
            result.stage1_pm_ratio_reliability = EELS_PearsonSparse_FullFit.stage1_pm_ratio_reliability(result, base_settings);
            result.stage2 = EELS_PearsonSparse_FullFit.solve_stage2_temperature_window(result, base_settings);
            if isfield(result.stage2, 'T_fit') && result.stage2.enabled && strcmp(result.stage2.status, 'ok')
                result.T_fit = result.stage2.T_fit;
                result.T_uncertainty = result.stage2.T_uncertainty;
                result.reliability_score = result.stage2.reliability_score;
                if isfield(result.stage2, 'reliability_source')
                    result.reliability_source = result.stage2.reliability_source;
                else
                    result.reliability_source = EELS_PearsonSparse_FullFit.reliability_source_from_result(result);
                end
                if isfield(result.stage2, 'T_accuracy_K')
                    result.T_accuracy_K = result.stage2.T_accuracy_K;
                    result.T_accuracy_method = result.stage2.T_accuracy_method;
                end
                if isfield(result.stage2, 'dominant_error_source')
                    result.dominant_error_source = result.stage2.dominant_error_source;
                elseif ~isfield(result, 'dominant_error_source')
                    result.dominant_error_source = 'stage2_reused';
                end
                result = EELS_PearsonSparse_FullFit.apply_stage2_primary_model(result);
            else
                result.T_fit = NaN;
                result.T_uncertainty = NaN;
                result.T_accuracy_K = NaN;
                result.T_accuracy_method = 'stage2_cost_valley_10pct_width_unavailable';
                if isfield(result.stage2, 'reliability_score')
                    result.reliability_score = result.stage2.reliability_score;
                    if isfield(result.stage2, 'reliability_source')
                        result.reliability_source = result.stage2.reliability_source;
                    else
                        result.reliability_source = EELS_PearsonSparse_FullFit.reliability_source_from_result(result);
                    end
                else
                    result.reliability_score = result.stage1_pm_ratio_reliability.reliability_score;
                    result.reliability_source = 'stage2_v20_failed';
                end
                result.dominant_error_source = 'stage2_reuse_failed';
            end
            result.reused_stage1 = true;
            result = EELS_PearsonSparse_FullFit.enrich_result_diagnostics(result, result.stage1, base_settings);
        end

        function settings = migrate_stage2_legacy_defaults(settings, override_names)
            if nargin < 2 || isempty(override_names)
                override_names = {};
            end
            override_names = cellstr(override_names(:)); %#ok<NASGU>
            if isfield(settings, 'stage2_window_meV') && ...
                    isequal(double(settings.stage2_window_meV(:).'), [40, 150])
                settings = rmfield(settings, 'stage2_window_meV');
            end
            if isfield(settings, 'stage2_attention_window_meV') && ...
                    isequal(double(settings.stage2_attention_window_meV(:).'), [40, 180])
                settings = rmfield(settings, 'stage2_attention_window_meV');
            end
            if isfield(settings, 'stage2_far_tail_downweight_start_meV') && ...
                    isequal(double(settings.stage2_far_tail_downweight_start_meV), 175)
                settings = rmfield(settings, 'stage2_far_tail_downweight_start_meV');
            end
        end

        function cuberes = promote_datacube_stage1_to_stage2(stage1_cuberes, settings)
            if nargin < 2, settings = struct(); end
            if ~(isstruct(stage1_cuberes) && isfield(stage1_cuberes, 'results') && isfield(stage1_cuberes, 'sample_pixels'))
                error('promote_datacube_stage1_to_stage2 expects a datacube result with results and sample_pixels fields.');
            end
            settings = EELS_PearsonSparse_FullFit.defaults(settings);
            settings.run_stage2 = true;
            if ~isfolder(settings.test_output_dir), mkdir(settings.test_output_dir); end
            np = numel(stage1_cuberes.sample_pixels);
            cuberes = stage1_cuberes;
            if isfield(cuberes, 'stage1_sparse_pm_policy')
                cuberes.stage1_sparse_pm_policy = ...
                    EELS_PearsonSparse_FullFit.sanitize_stage1_sparse_pm_policy_for_stage2_reuse( ...
                    cuberes.stage1_sparse_pm_policy);
            end
            cuberes.results = stage1_cuberes.results;
            class_info = struct('enabled', false, 'mode', 'v20_disabled', 'n_classes', 1, ...
                'labels', zeros(1, np), 'class_score', NaN(1, np), 'class_counts', np);
            cuberes.stage2_attention_class_info = class_info;
            timer0 = tic;
            pixel_time = NaN(1, np);
            [use_parallel, parallel_workers, parallel_message] = EELS_PearsonSparse_FullFit.datacube_parallel_state(settings, np);
            if use_parallel
                fprintf('Stage2 reuse parallel enabled: %d workers, %d pixels\n', parallel_workers, np);
                EELS_PearsonSparse_FullFit.ensure_parallel_pool(parallel_workers);
                progress_settings = settings;
                progress_settings.datacube_parallel_effective_workers = max(1, round(parallel_workers));
                parallel_pool_cleanup = []; %#ok<NASGU>
                if isfield(settings, 'close_parallel_pool_after_run') && settings.close_parallel_pool_after_run
                    parallel_pool_cleanup = onCleanup(@()EELS_PearsonSparse_FullFit.close_parallel_pool());
                end
                worker_settings = settings;
                worker_settings.datacube_progress_callback = [];
                worker_settings.datacube_stop_callback = [];
                progress_queue = [];
                try
                    progress_queue = parallel.pool.DataQueue;
                    afterEach(progress_queue, @(t)EELS_PearsonSparse_FullFit.parallel_progress_update( ...
                        progress_settings, t, np, timer0, 0, []));
                catch
                    progress_queue = [];
                end
                sample_pixels = stage1_cuberes.sample_pixels;
                stage1_results = stage1_cuberes.results;
                stage2_results = cell(1, np);
                parallel_pixel_time = NaN(1, np);
                parfor ii = 1:np
                    t1 = tic;
                    old_res = stage1_results{ii};
                    if isstruct(old_res) && ~isfield(old_res, 'error')
                        try
                            pixel_settings = worker_settings;
                            stage2_results{ii} = EELS_PearsonSparse_FullFit.promote_stage1_to_stage2(old_res, pixel_settings);
                        catch ME
                            stage2_results{ii} = struct('error', ME.message, 'pixel_index', sample_pixels(ii));
                        end
                    else
                        stage2_results{ii} = old_res;
                    end
                    parallel_pixel_time(ii) = toc(t1);
                    if ~isempty(progress_queue)
                        send(progress_queue, parallel_pixel_time(ii));
                    end
                end
                cuberes.results = stage2_results;
                pixel_time = parallel_pixel_time;
            else
                if settings.datacube_parallel_enabled && ~isempty(parallel_message)
                    fprintf('Stage2 reuse parallel disabled: %s\n', parallel_message);
                end
                for ii = 1:np
                    t1 = tic;
                    old_res = stage1_cuberes.results{ii};
                    if isstruct(old_res) && ~isfield(old_res, 'error')
                        try
                            pixel_settings = settings;
                            cuberes.results{ii} = EELS_PearsonSparse_FullFit.promote_stage1_to_stage2(old_res, pixel_settings);
                        catch ME
                            cuberes.results{ii} = struct('error', ME.message, 'pixel_index', stage1_cuberes.sample_pixels(ii));
                        end
                    end
                    pixel_time(ii) = toc(t1);
                    EELS_PearsonSparse_FullFit.notify_progress(settings, ii, np, toc(timer0), pixel_time(1:ii), 'Stage2 reuse');
                end
            end
            cuberes = EELS_PearsonSparse_FullFit.rebuild_datacube_maps_from_results(cuberes, settings, pixel_time, timer0, true);
        end

        function cuberes = solve_datacube(filename, settings)
            if nargin < 2, settings = struct(); end
            cuberes = EELS_PearsonSparse_FullFit.run_test_datacube(filename, settings);
        end

        function settings = make_datacube_settings(varargin)
            settings = EELS_PearsonSparse_FullFit.defaults(struct());
            if nargin == 1 && isstruct(varargin{1})
                input = varargin{1};
                names = fieldnames(input);
                for ii = 1:numel(names)
                    settings.(names{ii}) = input.(names{ii});
                end
                settings = EELS_PearsonSparse_FullFit.defaults(settings);
                return;
            end
            if mod(nargin, 2) ~= 0
                error('make_datacube_settings expects name/value pairs or one settings struct.');
            end
            for ii = 1:2:nargin
                key = lower(char(varargin{ii}));
                val = varargin{ii + 1};
                switch key
                    case {'spatial_bin', 'binsize', 'bin_size'}
                        if isscalar(val)
                            settings.datacube_y_bin = val;
                            settings.datacube_x_bin = val;
                        else
                            settings.datacube_y_bin = val(1);
                            settings.datacube_x_bin = val(2);
                        end
                    case {'y_bin', 'datacube_y_bin'}
                        settings.datacube_y_bin = val;
                    case {'x_bin', 'datacube_x_bin'}
                        settings.datacube_x_bin = val;
                    case {'energy_bin', 'e_bin', 'datacube_energy_bin'}
                        settings.datacube_energy_bin = val;
                    case {'crop', 'crop_range', 'energy_range'}
                        settings.crop_range = double(val(:).');
                    case {'output_dir', 'test_output_dir'}
                        settings.test_output_dir = val;
                    case {'sample_strategy', 'datacube_sample_strategy'}
                        settings.datacube_sample_strategy = val;
                    case {'sample_pixels', 'datacube_sample_pixels'}
                        settings.datacube_sample_pixels = val;
                    case {'full', 'run_full'}
                        if logical(val)
                            settings.datacube_sample_strategy = 'all';
                            settings.datacube_sample_pixels = Inf;
                        end
                    case {'stage2', 'run_stage2'}
                        settings.run_stage2 = logical(val);
                    case {'spatial_refit', 'spatial_refit_enabled'}
                        settings.spatial_refit_enabled = logical(val);
                    otherwise
                        settings.(char(varargin{ii})) = val;
                end
            end
            settings = EELS_PearsonSparse_FullFit.defaults(settings);
        end

        function result = solve_stage1(ene_raw, y_raw, settings)
            if nargin < 3, settings = struct(); end
            settings = EELS_PearsonSparse_FullFit.defaults(settings);
            [ene, y, y_norm, norm_val] = EELS_PearsonSparse_FullFit.prepare_spectrum(ene_raw, y_raw, settings);
            [ene, zlp_align] = EELS_PearsonSparse_FullFit.apply_stage1_zlp_alignment(ene, y_norm, settings);
            opt_mask = EELS_PearsonSparse_FullFit.optimization_mask(ene, settings);
            ene_opt = ene(opt_mask);
            y_norm_opt = y_norm(opt_mask);
            rounds = settings.stage1_rounds;
            if settings.stage1_warm_start_skip_coarse && ~isempty(settings.initial_theta) && numel(rounds) > 1
                rounds = rounds(2:end);
            end
            history_cell = cell(numel(rounds), 1);
            prev_theta = settings.initial_theta;
            best = [];

            for rr = 1:numel(rounds)
                round_timer = tic;
                round_settings = settings;
                fields = fieldnames(rounds(rr));
                for kk = 1:numel(fields)
                    round_settings.(fields{kk}) = rounds(rr).(fields{kk});
                end
                round_settings = EELS_PearsonSparse_FullFit.defaults(round_settings);
                EM = EELS_PearsonSparse_FullFit.make_em_grid(ene, round_settings);
                theta0 = EELS_PearsonSparse_FullFit.initial_theta(ene_opt, y_norm_opt, round_settings, prev_theta);
                [lb, ub] = EELS_PearsonSparse_FullFit.theta_bounds(round_settings);
                obj_settings = round_settings;
                obj_settings.inner_solver = round_settings.nonlinear_inner_solver;
                obj_settings.tail_overshoot_constraint = round_settings.nonlinear_tail_overshoot_constraint;
                obj_settings.lambda_pm_l1 = 0;
                obj = @(theta) EELS_PearsonSparse_FullFit.stage1_objective(theta, ene_opt, y_norm_opt, EM, obj_settings);
                [theta_opt, resnorm, exitflag, output] = EELS_PearsonSparse_FullFit.run_lsq(obj, theta0, lb, ub, round_settings);
                fit_fast = EELS_PearsonSparse_FullFit.decode_stage1(theta_opt, ene, y_norm, EM, obj_settings);
                fit_candidate = EELS_PearsonSparse_FullFit.decode_stage1(theta_opt, ene, y_norm, EM, round_settings);
                fit = EELS_PearsonSparse_FullFit.select_stage1_candidate(fit_fast, fit_candidate, ene, y_norm, round_settings);
                fit.resnorm = resnorm;
                fit.exitflag = exitflag;
                fit.output = output;
                fit.round_name = round_settings.name;
                metrics = EELS_PearsonSparse_FullFit.metrics(y_norm, fit.model_norm, ene, round_settings);
                rec = metrics;
                rec.round_id = rr;
                rec.round_name = round_settings.name;
                rec.num_EM = numel(EM);
                rec.exitflag = exitflag;
                rec.resnorm = resnorm;
                rec.elapsed_seconds = toc(round_timer);
                if isstruct(output) && isfield(output, 'funcCount')
                    rec.funcCount = output.funcCount;
                else
                    rec.funcCount = NaN;
                end
                if isstruct(output) && isfield(output, 'iterations')
                    rec.iterations = output.iterations;
                else
                    rec.iterations = NaN;
                end
                rec.accepted = true;
                history_cell{rr} = rec;
                fit.metrics = metrics;
                fit.settings = round_settings;
                fit.EM = EM;
                fit.theta = theta_opt;
                if isempty(best) || fit.metrics.NRMSE < best.metrics.NRMSE
                    best = fit;
                end
                prev_theta = theta_opt;
                if EELS_PearsonSparse_FullFit.should_stop_stage1_early(metrics, rr, numel(rounds), round_settings)
                    history_cell = history_cell(1:rr);
                    break;
                end
            end

            if EELS_PearsonSparse_FullFit.should_run_stage1_strict_polish(settings, best)
                [best, history_cell] = EELS_PearsonSparse_FullFit.polish_stage1_solution( ...
                    best, ene, y_norm, history_cell, settings);
            end

            history = vertcat(history_cell{:});
            result = EELS_PearsonSparse_FullFit.package_result(best, ene, y, y_norm, norm_val, settings, history, zlp_align);
        end

        function tf = should_run_stage1_strict_polish(settings, best)
            tf = settings.stage1_strict_polish_enabled && ~isempty(best) && ...
                isfield(best, 'settings') && isfield(best.settings, 'nonlinear_inner_solver') && ...
                strcmpi(string(best.settings.nonlinear_inner_solver), "fast_nnls") && ...
                isfield(best, 'theta') && isfield(best, 'EM') && settings.stage1_strict_polish_iter > 0;
        end

        function [best, history_cell] = polish_stage1_solution(best, ene, y_norm, history_cell, settings)
            polish_settings = best.settings;
            polish_settings.name = 'strict_polish';
            polish_settings.max_nonlinear_iter = settings.stage1_strict_polish_iter;
            polish_settings.nonlinear_inner_solver = settings.stage1_strict_polish_inner_solver;
            obj_settings = polish_settings;
            obj_settings.inner_solver = settings.stage1_strict_polish_inner_solver;
            obj_settings.tail_overshoot_constraint = false;
            obj_settings.lambda_pm_l1 = 0;
            [lb, ub] = EELS_PearsonSparse_FullFit.theta_bounds(polish_settings);
            theta0 = min(max(best.theta, lb), ub);
            obj = @(theta) EELS_PearsonSparse_FullFit.stage1_objective(theta, ene, y_norm, best.EM, obj_settings);
            t_polish = tic;
            [theta_opt, resnorm, exitflag, output] = EELS_PearsonSparse_FullFit.run_lsq(obj, theta0, lb, ub, polish_settings);
            fit_fast = EELS_PearsonSparse_FullFit.decode_stage1(theta_opt, ene, y_norm, best.EM, obj_settings);
            fit_candidate = EELS_PearsonSparse_FullFit.decode_stage1(theta_opt, ene, y_norm, best.EM, polish_settings);
            fit = EELS_PearsonSparse_FullFit.select_stage1_candidate(fit_fast, fit_candidate, ene, y_norm, polish_settings);
            fit.resnorm = resnorm;
            fit.exitflag = exitflag;
            fit.output = output;
            fit.round_name = polish_settings.name;
            fit.settings = polish_settings;
            fit.EM = best.EM;
            fit.theta = theta_opt;
            metrics = EELS_PearsonSparse_FullFit.metrics(y_norm, fit.model_norm, ene, polish_settings);
            fit.metrics = metrics;
            rec = metrics;
            rec.round_id = numel(history_cell) + 1;
            rec.round_name = polish_settings.name;
            rec.num_EM = numel(best.EM);
            rec.exitflag = exitflag;
            rec.resnorm = resnorm;
            rec.elapsed_seconds = toc(t_polish);
            if isstruct(output) && isfield(output, 'funcCount')
                rec.funcCount = output.funcCount;
            else
                rec.funcCount = NaN;
            end
            if isstruct(output) && isfield(output, 'iterations')
                rec.iterations = output.iterations;
            else
                rec.iterations = NaN;
            end
            accept = metrics.weighted_NRMSE <= settings.stage1_strict_polish_accept_factor * max(best.metrics.weighted_NRMSE, eps) || ...
                (isfinite(metrics.highfreq_NRMSE) && isfinite(best.metrics.highfreq_NRMSE) && ...
                metrics.highfreq_NRMSE <= settings.stage1_strict_polish_accept_factor * max(best.metrics.highfreq_NRMSE, eps));
            rec.accepted = accept;
            history_cell{end + 1, 1} = rec;
            if accept
                best = fit;
            end
        end

        function tf = should_stop_stage1_early(metrics, round_id, num_rounds, settings)
            tf = false;
            if ~settings.stage1_early_stop_enabled || round_id < settings.stage1_early_stop_min_rounds || round_id >= num_rounds
                return;
            end
            high_ok = ~isfield(metrics, 'highfreq_NRMSE') || ~isfinite(metrics.highfreq_NRMSE) || ...
                metrics.highfreq_NRMSE <= settings.stage1_early_stop_highfreq_NRMSE;
            tf = isfield(metrics, 'R2') && isfield(metrics, 'weighted_NRMSE') && ...
                isfinite(metrics.R2) && isfinite(metrics.weighted_NRMSE) && ...
                metrics.R2 >= settings.stage1_early_stop_R2 && ...
                metrics.weighted_NRMSE <= settings.stage1_early_stop_weighted_NRMSE && high_ok;
        end

        function result = solve_stage2_temperature_window(stage1_result, settings)
            stage2_timer = tic;
            if nargin < 2, settings = struct(); end
            settings = EELS_PearsonSparse_FullFit.defaults(settings);
            result = struct();
            result.enabled = settings.run_stage2;
            result.status = 'disabled';
            result.message = 'Stage 2 disabled; v2.0 attention hard-BE + soft source-intensity fitting is available when settings.run_stage2 is true.';
            if isfield(stage1_result, 'T_soft')
                result.stage1_T_soft_seed = stage1_result.T_soft;
            else
                result.stage1_T_soft_seed = settings.T_soft_init;
            end
            if isfield(stage1_result, 'PM_stokes'), result.stage1_PM_stokes = stage1_result.PM_stokes; end
            if isfield(stage1_result, 'PM_antistokes'), result.stage1_PM_antistokes = stage1_result.PM_antistokes; end
            if isfield(stage1_result, 'stage1_zlp_alignment')
                result.stage1_zlp_alignment = stage1_result.stage1_zlp_alignment;
            end
            if isfield(stage1_result, 'stage1_zero_align_shift_meV')
                result.stage1_zero_align_shift_meV = stage1_result.stage1_zero_align_shift_meV;
            end
            result.settings = settings;
            if ~settings.run_stage2
                result.elapsed_seconds = toc(stage2_timer);
                return;
            end

            ene = stage1_result.ene_raw(:);
            ene_corr = stage1_result.ene_corrected(:);
            y_norm = stage1_result.raw_data_norm(:);
            params = stage1_result.stage1.params;
            zlp_norm = stage1_result.stage1.zlp_norm(:);
            win = sort(abs(settings.stage2_window_meV(:).'));
            if numel(win) < 2 || any(~isfinite(win(1:2)))
                win = [40, 190];
            end
            hard_window_mask = abs(ene_corr) >= win(1) & abs(ene_corr) <= win(2);
            attention = EELS_PearsonSparse_FullFit.build_stage2_attention(stage1_result, settings);
            attention = EELS_PearsonSparse_FullFit.update_stage2_attention_tail_masks(attention, ene_corr, settings);
            if ~isfield(attention, 'valid_mask') || numel(attention.valid_mask) ~= numel(ene_corr)
                attention.valid_mask = true(size(ene_corr));
            end
            if ~isfield(attention, 'feature_score') || numel(attention.feature_score) ~= numel(ene_corr)
                attention.feature_score = zeros(size(ene_corr));
            end
            if ~isfield(attention, 'active_mask') || numel(attention.active_mask) ~= numel(ene_corr)
                attention.active_mask = false(size(ene_corr));
            end
            mirrored_feature = interp1(ene_corr, min(max(attention.feature_score(:), 0), 1), -ene_corr, 'linear', 0);
            mirrored_active = interp1(ene_corr, double(attention.active_mask(:)), -ene_corr, 'nearest', 0) > 0;
            pair_feature = max(min(max(attention.feature_score(:), 0), 1), min(max(mirrored_feature(:), 0), 1));
            fit_threshold = min(max(settings.stage2_attention_fit_threshold, 0), 1);
            fit_support_seed = hard_window_mask(:) & attention.valid_mask(:) & ...
                (pair_feature(:) >= fit_threshold | attention.active_mask(:) | mirrored_active(:));
            if nnz(fit_support_seed) < 12
                fit_support_seed = hard_window_mask(:) & attention.valid_mask(:);
            end
            fit_support_alpha = EELS_PearsonSparse_FullFit.smooth_binary_mask_by_energy( ...
                hard_window_mask(:) & attention.valid_mask(:), ene_corr, settings.stage2_fit_sigmoid_taper_meV);
            attention_alpha = EELS_PearsonSparse_FullFit.smooth_binary_mask_by_energy( ...
                fit_support_seed, ene_corr, settings.stage2_fit_sigmoid_taper_meV);
            fit_alpha = min(max(fit_support_alpha(:) .* max(attention_alpha(:), pair_feature(:)), 0), 1);
            fit_mask = attention.valid_mask(:) & fit_alpha(:) >= min(max(settings.stage2_fit_alpha_min, 1e-5), 0.25);
            if nnz(fit_mask) < 12
                fit_alpha = min(max(fit_support_alpha(:), 0), 1);
                fit_mask = attention.valid_mask(:) & fit_alpha(:) > 1e-4;
            end
            if nnz(fit_mask) < 12
                result.status = 'insufficient_window_points';
                result.T_fit = NaN;
                result.T_uncertainty = NaN;
                result.T_accuracy_K = NaN;
                result.T_accuracy_method = 'stage2_v20_smoothed_cost_valley_unavailable';
                result.elapsed_seconds = toc(stage2_timer);
                return;
            end
            attention.fit_candidate_mask = fit_support_seed(:);
            attention.fit_alpha = fit_alpha(:);
            attention.pair_feature = pair_feature(:);
            soft_floor = min(max(settings.stage2_fit_sigmoid_min_weight, 0), 0.95);
            attention_weight = settings.stage2_attention_floor + (1 - settings.stage2_attention_floor) .* pair_feature(:);
            attention.weight = min(max(attention.weight(:), attention_weight(:)), 1);

            basis = EELS_PearsonSparse_FullFit.build_stage2_fixed_stage1_basis(stage1_result, attention, settings);
            if isempty(basis.EM_all)
                result.status = 'empty_stage1_em_grid';
                result.T_fit = NaN;
                result.T_uncertainty = NaN;
                result.T_accuracy_K = NaN;
                result.T_accuracy_method = 'stage2_v20_empty_stage1_em_grid';
                result.elapsed_seconds = toc(stage2_timer);
                return;
            end
            [A, layout, ~, kernels] = EELS_PearsonSparse_FullFit.build_design_matrix(ene, basis.EM_all, params, settings);
            if EELS_PearsonSparse_FullFit.is_fitted_zlp_kernel_mode(settings.kernel_mode)
                [A, kernels] = EELS_PearsonSparse_FullFit.replace_signal_kernel_columns( ...
                    A, layout, ene_corr(:), basis.EM_all(:).', zlp_norm(:));
            end
            base_sqrtW = EELS_PearsonSparse_FullFit.build_weights(ene_corr, y_norm, settings);
            joint_floor = min(max(settings.stage2_global_joint_weight_floor, 0.01), 1);
            local_profile = (soft_floor + (1 - soft_floor) .* fit_alpha(:)) .* attention.weight(:);
            joint_profile = max(joint_floor, min(max(local_profile(:), 0), 1));
            sqrtW = base_sqrtW(:) .* sqrt(joint_profile(:));
            anti_fit = ene_corr < -win(1) & ene_corr >= -win(2);
            if settings.stage2_antistokes_fit_weight_boost > 1
                sqrtW(anti_fit) = sqrtW(anti_fit) .* sqrt(settings.stage2_antistokes_fit_weight_boost);
            end
            sqrtW = EELS_PearsonSparse_FullFit.apply_stage2_far_tail_downweight(ene_corr, sqrtW, settings, attention);
            attention.joint_weight_profile = joint_profile(:);

            scan_settings = settings;
            scan_settings.inner_solver = settings.stage2_scan_inner_solver;
            scan_settings.tail_overshoot_constraint = settings.stage2_scan_tail_overshoot_constraint;
            T_grid = settings.stage2_T_grid(:).';
            scan = EELS_PearsonSparse_FullFit.stage2_scan_temperatures_adaptive(T_grid, basis, kernels, A, layout, ...
                fit_mask, y_norm, zlp_norm, sqrtW, stage1_result, scan_settings);
            scan = EELS_PearsonSparse_FullFit.stage2_refine_temperature_local(scan, T_grid, basis, kernels, A, layout, ...
                fit_mask, y_norm, zlp_norm, sqrtW, stage1_result, scan_settings);
            best = scan.best;
            if ~isfield(best, 'T_fit') || ~isfinite(best.T_fit)
                result.status = 'temperature_scan_failed';
                result.T_fit = NaN;
                result.T_grid = scan.T_grid;
                result.cost_grid = scan.cost_grid;
                if isfield(scan, 'cost_grid_smooth'), result.cost_grid_smooth = scan.cost_grid_smooth; end
                result.fit_mask = fit_mask;
                result.attention = attention;
                result.EM_final = basis.EM_all;
                result.elapsed_seconds = toc(stage2_timer);
                return;
            end

            stage2_valley = EELS_PearsonSparse_FullFit.stage2_cost_valley(scan.T_grid, scan.cost_grid_smooth);
            result.status = 'ok';
            result.message = 'Stage 2 v2.0 attention hard-BE + soft source intensity: only attention-paired mode energies constrain temperature; remaining centers fit independent source intensities.';
            result.temperature_algorithm = 'grid_scan';
            result.temperature_policy = 'v20_attention_hardBE_softPM';
            result.T_fit = best.T_fit;
            result.cost = best.cost;
            result.inner_solver = best.inner_solver;
            result.model_norm = best.model_norm;
            result.zlp_norm = best.zlp_norm;
            result.zlp_delta_norm = best.zlp_delta_norm;
            result.zlp_refine = best.zlp_refine;
            result.vib_norm = best.vib_norm;
            result.vib_hard_norm = best.vib_hard_norm;
            result.vib_soft_norm = best.vib_soft_norm;
            result.bg_norm = zeros(size(best.vib_norm));
            result.metrics = EELS_PearsonSparse_FullFit.metrics(y_norm, best.model_norm, ene_corr, settings);
            result.joint_RMSE_norm = sqrt(mean((best.model_norm(:) - y_norm(:)) .^ 2, 'omitnan'));
            result.joint_NRMSE = result.metrics.NRMSE;
            result.joint_R2 = result.metrics.R2;
            result.PM_stokes = best.PM_stokes;
            result.PM_antistokes = best.PM_antistokes;
            result.PM_stokes_hard = best.PM_stokes_hard;
            result.PM_antistokes_hard = best.PM_antistokes_hard;
            result.PM_stokes_soft = best.PM_stokes_soft;
            result.PM_antistokes_soft = best.PM_antistokes_soft;
            result.EM_final = basis.EM_all;
            result.EM_hard = basis.EM_hard;
            result.EM_soft = basis.EM_soft;
            result.hard_idx = basis.hard_idx;
            result.soft_idx = basis.soft_idx;
            result.hard_fraction = basis.hard_fraction;
            result.hard_pair_score = basis.hard_pair_score;
            result.hard_pair_score_hard = basis.hard_pair_score_hard;
            result.soft_pair_score = basis.soft_pair_score;
            result.hard_be_attention_status = basis.hard_selection_status;
            result.hard_pair_threshold = basis.hard_pair_threshold;
            result.hard_pair_floor = basis.hard_pair_floor;
            result.min_hard_pairs = basis.min_hard_pairs;
            result.num_hard_pairs = basis.num_hard_pairs;
            result.num_soft_pairs = basis.num_soft_pairs;
            result.num_stage1_em = basis.num_stage1_em;
            result.num_stage2_em = basis.num_stage2_em;
            result.num_hard_dense_added = basis.num_hard_dense_added;
            result.num_soft_dense_added = basis.num_soft_dense_added;
            result.EM_stage1 = basis.EM_stage1;
            result.EM_hard_dense_added = basis.EM_hard_dense_added;
            result.EM_soft_dense_added = basis.EM_soft_dense_added;
            result.temperature_reliability_flag = '';
            if strcmp(basis.hard_selection_status, 'insufficient_hard_be_attention_pairs')
                result.temperature_reliability_flag = 'insufficient_hard_be_attention_pairs';
            end
            result.stage2_EM_matches_stage1 = EELS_PearsonSparse_FullFit.stage2_fixed_grid_matches_stage1(result.EM_final, stage1_result.EM_final);
            result.window_meV = win;
            result.basis = basis;
            result.fit_mask = fit_mask;
            result.fit_alpha = fit_alpha;
            result.hard_window_mask = hard_window_mask;
            result.transition_mask = false(size(fit_mask));
            result.attention = attention;
            result.T_grid = scan.T_grid;
            result.cost_grid = scan.cost_grid;
            result.cost_grid_smooth = scan.cost_grid_smooth;
            result.temperature_profile = scan.temperature_profile;
            if isfield(scan, 'scan_strategy')
                result.scan_strategy = scan.scan_strategy;
            end
            if isfield(scan, 'temperature_refine')
                result.temperature_refine = scan.temperature_refine;
            end
            result.T_cost_min = stage2_valley.best_T;
            result.T_profile_offset_K = result.T_fit - stage2_valley.best_T;
            result.T_cost_valley_lower_K = stage2_valley.lower_K;
            result.T_cost_valley_upper_K = stage2_valley.upper_K;
            result.T_uncertainty = stage2_valley.uncertainty_K;
            result.T_accuracy_K = result.T_uncertainty;
            result.T_accuracy_method = 'stage2_v20_smoothed_cost_valley_10pct_width';
            result.T_uncertainty_method = result.T_accuracy_method;
            result.reliability_score = EELS_PearsonSparse_FullFit.stage2_reliability(result, stage1_result);
            result.stage2_reliability_score = result.reliability_score;
            result.reliability_source = 'stage2_v20_smoothed_cost_attention_linearity';
            result.elapsed_seconds = toc(stage2_timer);
        end

        function testres = run_test_simu(settings)
            if nargin < 1, settings = struct(); end
            settings = EELS_PearsonSparse_FullFit.defaults(settings);
            if ~isfolder(settings.test_output_dir), mkdir(settings.test_output_dir); end
            if settings.save_plots
                plot_dir = fullfile(settings.test_output_dir, 'plots');
                if ~isfolder(plot_dir), mkdir(plot_dir); end
            end

            S = load('simu2.mat', 'simu2', 'temp2');
            E = load('eneref.mat', 'eneref');
            ene = double(E.eneref(:));
            T_true = double(S.temp2(:));
            n_spec = size(S.simu2, 2);
            if isempty(settings.test_indices)
                idx = EELS_PearsonSparse_FullFit.pick_test_indices(T_true, n_spec, settings.test_max_spectra);
            else
                idx = settings.test_indices(:).';
            end
            n = numel(idx);
            results = cell(1, n);
            stage2 = cell(1, n);
            R2_all = NaN(1, n);
            NRMSE_all = NaN(1, n);
            weighted_NRMSE_all = NaN(1, n);
            T_fit_all = NaN(1, n);
            T_error_all = NaN(1, n);
            fit_status = zeros(1, n);
            for ii = 1:n
                try
                    spectrum = double(S.simu2(:, idx(ii)));
                    res = EELS_PearsonSparse_FullFit.solve_stage1(ene, spectrum, settings);
                    st2 = res.stage2;
                    results{ii} = res;
                    stage2{ii} = st2;
                    R2_all(ii) = res.metrics.R2;
                    NRMSE_all(ii) = res.metrics.NRMSE;
                    weighted_NRMSE_all(ii) = res.metrics.weighted_NRMSE;
                    if isfield(res.stage2, 'T_fit') && res.stage2.enabled && strcmp(res.stage2.status, 'ok')
                        T_fit_all(ii) = res.stage2.T_fit;
                        T_error_all(ii) = T_fit_all(ii) - T_true(idx(ii));
                    end
                    fit_status(ii) = 1;
                    if settings.save_plots
                        fig = EELS_PearsonSparse_FullFit.plot_result(res, sprintf('simu2 index %d, T true %.2f K', idx(ii), T_true(idx(ii))));
                        saveas(fig, fullfile(settings.test_output_dir, 'plots', sprintf('simu2_%06d_stage1.png', idx(ii))));
                        close(fig);
                    end
                catch ME
                    fit_status(ii) = -1;
                    results{ii} = struct('error', ME.message, 'index', idx(ii));
                end
            end

            testres = struct();
            testres.indices = idx;
            testres.T_true = T_true(idx);
            testres.results = results;
            testres.stage2 = stage2;
            testres.R2_all = R2_all;
            testres.NRMSE_all = NRMSE_all;
            testres.weighted_NRMSE_all = weighted_NRMSE_all;
            testres.T_fit_all = T_fit_all;
            testres.T_error_all = T_error_all;
            testres.fit_status = fit_status;
            testres.settings = settings;
            save(fullfile(settings.test_output_dir, 'simu_stage1_result.mat'), 'testres', '-v7.3');
            EELS_PearsonSparse_FullFit.write_summary(fullfile(settings.test_output_dir, 'simu_stage1_summary.txt'), testres);
        end

        function cuberes = run_test_datacube(filename, settings)
            if nargin < 2, settings = struct(); end
            settings = EELS_PearsonSparse_FullFit.defaults(settings);
            figure_visibility_cleanup = EELS_PearsonSparse_FullFit.hidden_figure_cleanup(settings); %#ok<NASGU>
            [~, cube_tag] = fileparts(char(filename));
            if isempty(cube_tag), cube_tag = 'datacube'; end
            cube_tag = regexprep(cube_tag, '[<>:"/\\|?*\s]+', '_');
            if ~isfolder(settings.test_output_dir), mkdir(settings.test_output_dir); end
            if settings.save_plots
                plot_dir = fullfile(settings.test_output_dir, 'plots');
                if ~isfolder(plot_dir), mkdir(plot_dir); end
            end
            [cube, ene, cube_info] = EELS_PearsonSparse_FullFit.load_datacube_mat(filename);
            [cube, bin_info] = EELS_PearsonSparse_FullFit.bin_datacube_spatial(cube, settings.datacube_y_bin, settings.datacube_x_bin);
            [cube, ene, energy_bin_info] = EELS_PearsonSparse_FullFit.bin_datacube_energy(cube, ene, settings.datacube_energy_bin);
            cube_info.original_size = cube_info.size;
            cube_info.size = size(cube);
            cube_info.spatial_bin = bin_info;
            cube_info.energy_bin = energy_bin_info;
            [ny, nx, ~] = size(cube);
            pixels = EELS_PearsonSparse_FullFit.pick_datacube_pixels(ny, nx, settings);
            np = numel(pixels);
            settings = EELS_PearsonSparse_FullFit.calibrate_datacube_stage1_sparse_pm(cube, ene, pixels, settings);
            results = cell(1, np);
            R2_map = NaN(ny, nx);
            NRMSE_map = NaN(ny, nx);
            T_map_raw = NaN(ny, nx);
            T_uncertainty_map = NaN(ny, nx);
            reliability_score_map = NaN(ny, nx);
            reliability_source_map = repmat({''}, ny, nx);
            stage1_PM_ratio_T_map = NaN(ny, nx);
            stage1_PM_ratio_reliability_map = NaN(ny, nx);
            stage1_PM_ratio_pairs_map = NaN(ny, nx);
            stage1_PM_ratio_bestline_R2_map = NaN(ny, nx);
            stage1_PM_ratio_selected_pairs_map = NaN(ny, nx);
            stage1_PM_ratio_selected_EM_lo_map = NaN(ny, nx);
            stage1_PM_ratio_selected_EM_hi_map = NaN(ny, nx);
            stage2_zlp_refine_max_delta_map = NaN(ny, nx);
            stage2_zlp_refine_fit_rms_map = NaN(ny, nx);
            attention_score_map = NaN(ny, nx);
            zlp_guard_start_map = NaN(ny, nx);
            zlp_guard_fraction_map = NaN(ny, nx);
            zlp_core_width_map = NaN(ny, nx);
            nbins = settings.spatial_attention_profile_bins;
            attention_profile_cube = NaN(ny, nx, nbins);
            signal_profile_cube = NaN(ny, nx, nbins);
            fit_status_map = zeros(ny, nx);
            run_timer = tic;
            [use_parallel, parallel_workers, parallel_message] = EELS_PearsonSparse_FullFit.datacube_parallel_state(settings, np);
            if use_parallel
                if ~settings.silent
                    fprintf('datacube parallel enabled: %d workers, %d pixels\n', parallel_workers, np);
                end
                [results, pixel_time] = EELS_PearsonSparse_FullFit.fit_datacube_pixels_parallel(cube, ene, pixels, settings, parallel_workers);
                progress_settings = settings;
                progress_settings.datacube_parallel_effective_workers = parallel_workers;
                EELS_PearsonSparse_FullFit.notify_progress(progress_settings, np, np, toc(run_timer), pixel_time, 'datacube parallel');
            else
                if settings.datacube_parallel_enabled && ~isempty(parallel_message)
                    fprintf('datacube parallel disabled: %s\n', parallel_message);
                end
                [results, pixel_time] = EELS_PearsonSparse_FullFit.fit_datacube_pixels_serial(cube, ene, pixels, settings, run_timer);
            end
            EELS_PearsonSparse_FullFit.notify_progress(settings, np, np, toc(run_timer), pixel_time, 'building maps and diagnostics');
            for ii = 1:np
                res = results{ii};
                if ~(isstruct(res) && isfield(res, 'error'))
                    [yy, xx] = ind2sub([ny, nx], pixels(ii));
                    R2_map(yy, xx) = res.metrics.R2;
                    NRMSE_map(yy, xx) = res.metrics.NRMSE;
                    T_map_raw(yy, xx) = res.T_fit;
                    T_uncertainty_map(yy, xx) = res.T_uncertainty;
                    reliability_score_map(yy, xx) = res.reliability_score;
                    reliability_source_map{yy, xx} = EELS_PearsonSparse_FullFit.reliability_source_from_result(res);
                    [stage1_PM_ratio_T_map(yy, xx), stage1_PM_ratio_reliability_map(yy, xx), ...
                        stage1_PM_ratio_pairs_map(yy, xx), stage1_PM_ratio_bestline_R2_map(yy, xx), ...
                        stage1_PM_ratio_selected_pairs_map(yy, xx), stage1_PM_ratio_selected_EM_lo_map(yy, xx), ...
                        stage1_PM_ratio_selected_EM_hi_map(yy, xx)] = EELS_PearsonSparse_FullFit.stage1_ratio_summary_from_result(res);
                    [stage2_zlp_refine_max_delta_map(yy, xx), stage2_zlp_refine_fit_rms_map(yy, xx)] = ...
                        EELS_PearsonSparse_FullFit.stage2_zlp_refine_summary_from_result(res);
                    if isfield(res, 'zlp_signal_guard')
                        zlp_guard_start_map(yy, xx) = res.zlp_signal_guard.energy_range_meV(1);
                        zlp_guard_fraction_map(yy, xx) = res.zlp_signal_guard.zlp_area_fraction;
                    end
                    if isfield(res, 'zlp_region_policy') && isfield(res.zlp_region_policy, 'core_downweight_width_meV')
                        zlp_core_width_map(yy, xx) = res.zlp_region_policy.core_downweight_width_meV;
                    end
                    if isfield(res.stage2, 'attention') && isfield(res.stage2.attention, 'selected_fraction')
                        attention_score_map(yy, xx) = res.stage2.attention.selected_fraction;
                        attention_profile_cube(yy, xx, :) = EELS_PearsonSparse_FullFit.bin_energy_profile( ...
                            res.ene_corrected, res.stage2.attention.feature_score, nbins);
                    elseif isfield(res, 'attention_summary')
                        attention_profile_cube(yy, xx, :) = EELS_PearsonSparse_FullFit.bin_energy_profile( ...
                            res.ene_corrected, res.attention_summary.attention_profile, nbins);
                    end
                    signal_profile_cube(yy, xx, :) = EELS_PearsonSparse_FullFit.bin_energy_profile( ...
                        res.ene_corrected, max(res.raw_data_norm(:) - res.stage1.zlp_norm(:), 0), nbins);
                    fit_status_map(yy, xx) = 1;
                    if settings.save_plots && ii <= min(np, settings.max_saved_spectrum_plots)
                        EELS_PearsonSparse_FullFit.notify_progress(settings, np, np, toc(run_timer), pixel_time, ...
                            sprintf('saving spectrum diagnostic %d/%d', ii, min(np, settings.max_saved_spectrum_plots)));
                        fig = EELS_PearsonSparse_FullFit.plot_result(res, sprintf('%s pixel %d (%d,%d)', cube_tag, pixels(ii), yy, xx));
                        saveas(fig, fullfile(settings.test_output_dir, 'plots', sprintf('%s_pixel_%06d_fit.png', ...
                            cube_tag, pixels(ii))));
                        close(fig);
                    end
                else
                    [yy, xx] = ind2sub([ny, nx], pixels(ii));
                    fit_status_map(yy, xx) = -1;
                end
            end
            maps = struct();
            maps.R2_map = R2_map;
            maps.NRMSE_map = NRMSE_map;
            maps.T_map_raw = T_map_raw;
            maps.T_uncertainty_map = T_uncertainty_map;
            maps.reliability_score_map = reliability_score_map;
            maps.reliability_source_map = reliability_source_map;
            maps.stage1_PM_ratio_T_map = stage1_PM_ratio_T_map;
            maps.stage1_PM_ratio_reliability_map = stage1_PM_ratio_reliability_map;
            maps.stage1_PM_ratio_pairs_map = stage1_PM_ratio_pairs_map;
            maps.stage1_PM_ratio_bestline_R2_map = stage1_PM_ratio_bestline_R2_map;
            maps.stage1_PM_ratio_selected_pairs_map = stage1_PM_ratio_selected_pairs_map;
            maps.stage1_PM_ratio_selected_EM_lo_map = stage1_PM_ratio_selected_EM_lo_map;
            maps.stage1_PM_ratio_selected_EM_hi_map = stage1_PM_ratio_selected_EM_hi_map;
            maps.stage2_zlp_refine_max_delta_map = stage2_zlp_refine_max_delta_map;
            maps.stage2_zlp_refine_fit_rms_map = stage2_zlp_refine_fit_rms_map;
            maps.attention_score_map = attention_score_map;
            maps.zlp_guard_start_map = zlp_guard_start_map;
            maps.zlp_guard_fraction_map = zlp_guard_fraction_map;
            maps.zlp_core_width_map = zlp_core_width_map;
            maps.attention_profile_cube = attention_profile_cube;
            maps.signal_profile_cube = signal_profile_cube;
            maps.fit_status_map = fit_status_map;
            spatial = EELS_PearsonSparse_FullFit.evaluate_spatial_consistency(maps.T_map_raw, maps.R2_map, ...
                maps.NRMSE_map, maps.reliability_score_map, maps.attention_score_map, maps.fit_status_map, settings, ...
                maps.attention_profile_cube, maps.signal_profile_cube);
            spatial_refit_history = struct('iteration', {}, 'candidate_count', {}, 'attempted', {}, 'accepted', {});
            if settings.spatial_refit_enabled
                [results, maps, spatial, spatial_refit_history] = EELS_PearsonSparse_FullFit.spatial_refit_datacube( ...
                    ene, cube, pixels, results, maps, settings, spatial);
            end
            R2_map = maps.R2_map;
            NRMSE_map = maps.NRMSE_map;
            T_map_raw = maps.T_map_raw;
            T_uncertainty_map = maps.T_uncertainty_map;
            reliability_score_map = maps.reliability_score_map;
            reliability_source_map = maps.reliability_source_map;
            stage1_PM_ratio_T_map = maps.stage1_PM_ratio_T_map;
            stage1_PM_ratio_reliability_map = maps.stage1_PM_ratio_reliability_map;
            stage1_PM_ratio_pairs_map = maps.stage1_PM_ratio_pairs_map;
            stage1_PM_ratio_bestline_R2_map = maps.stage1_PM_ratio_bestline_R2_map;
            stage1_PM_ratio_selected_pairs_map = maps.stage1_PM_ratio_selected_pairs_map;
            stage1_PM_ratio_selected_EM_lo_map = maps.stage1_PM_ratio_selected_EM_lo_map;
            stage1_PM_ratio_selected_EM_hi_map = maps.stage1_PM_ratio_selected_EM_hi_map;
            stage2_zlp_refine_max_delta_map = maps.stage2_zlp_refine_max_delta_map;
            stage2_zlp_refine_fit_rms_map = maps.stage2_zlp_refine_fit_rms_map;
            attention_score_map = maps.attention_score_map;
            zlp_guard_start_map = maps.zlp_guard_start_map;
            zlp_guard_fraction_map = maps.zlp_guard_fraction_map;
            zlp_core_width_map = maps.zlp_core_width_map;
            attention_profile_cube = maps.attention_profile_cube;
            signal_profile_cube = maps.signal_profile_cube;
            fit_status_map = maps.fit_status_map;
            cuberes = struct();
            cuberes.filename = filename;
            cuberes.cube_info = cube_info;
            cuberes.sample_pixels = pixels;
            cuberes.results = results;
            cuberes.R2_map = R2_map;
            cuberes.NRMSE_map = NRMSE_map;
            cuberes.T_map_raw = T_map_raw;
            cuberes.T_map_repaired = spatial.T_map_repaired;
            cuberes.T_map_spatial_refined = spatial.T_map_spatial_refined;
            cuberes.spatial_refine_delta_map = spatial.spatial_refine_delta_map;
            cuberes.spatial_refine_mask = spatial.spatial_refine_mask;
            cuberes.T_uncertainty_map = T_uncertainty_map;
            cuberes.reliability_score_map = reliability_score_map;
            cuberes.reliability_source_map = reliability_source_map;
            cuberes.stage1_PM_ratio_T_map = stage1_PM_ratio_T_map;
            cuberes.stage1_PM_ratio_reliability_map = stage1_PM_ratio_reliability_map;
            cuberes.stage1_PM_ratio_pairs_map = stage1_PM_ratio_pairs_map;
            cuberes.stage1_PM_ratio_bestline_R2_map = stage1_PM_ratio_bestline_R2_map;
            cuberes.stage1_PM_ratio_selected_pairs_map = stage1_PM_ratio_selected_pairs_map;
            cuberes.stage1_PM_ratio_selected_EM_lo_map = stage1_PM_ratio_selected_EM_lo_map;
            cuberes.stage1_PM_ratio_selected_EM_hi_map = stage1_PM_ratio_selected_EM_hi_map;
            cuberes.stage2_zlp_refine_max_delta_map = stage2_zlp_refine_max_delta_map;
            cuberes.stage2_zlp_refine_fit_rms_map = stage2_zlp_refine_fit_rms_map;
            cuberes.attention_score_map = attention_score_map;
            cuberes.zlp_guard_start_map = zlp_guard_start_map;
            cuberes.zlp_guard_fraction_map = zlp_guard_fraction_map;
            cuberes.zlp_core_width_map = zlp_core_width_map;
            cuberes.attention_profile_cube = attention_profile_cube;
            cuberes.signal_profile_cube = signal_profile_cube;
            cuberes.warning_map = spatial.warning_map;
            cuberes.bad_pixel_mask = spatial.bad_pixel_mask;
            cuberes.interface_mask = spatial.interface_mask;
            cuberes.neighbor_loss_map = spatial.neighbor_loss_map;
            cuberes.spectrum_neighbor_loss_map = spatial.spectrum_neighbor_loss_map;
            cuberes.temperature_neighbor_z_map = spatial.temperature_neighbor_z_map;
            cuberes.attention_peak_delta_map = spatial.attention_peak_delta_map;
            cuberes.signal_peak_delta_map = spatial.signal_peak_delta_map;
            cuberes.spatial_refit_candidate_mask = spatial.spatial_refit_candidate_mask;
            cuberes.low_R2_candidate_mask = spatial.low_R2_candidate_mask;
            cuberes.region_label_map = spatial.region_label_map;
            cuberes.region_median_map = spatial.region_median_map;
            cuberes.region_delta_map = spatial.region_delta_map;
            cuberes.region_z_map = spatial.region_z_map;
            cuberes.region_size_map = spatial.region_size_map;
            cuberes.region_post_smooth_mask = spatial.region_post_smooth_mask;
            cuberes.spatial_refit_history = spatial_refit_history;
            cuberes.fit_status_map = fit_status_map;
            if isfield(settings, 'stage1_sparse_pm_policy')
                cuberes.stage1_sparse_pm_policy = settings.stage1_sparse_pm_policy;
            end
            cuberes.spatial = spatial;
            cuberes.benchmark = struct('total_seconds', toc(run_timer), 'median_seconds_per_pixel', median(pixel_time, 'omitnan'), ...
                'mean_seconds_per_pixel', mean(pixel_time, 'omitnan'), 'pixel_seconds', pixel_time, ...
                'pixels_per_second', np ./ max(toc(run_timer), eps), 'parallel_enabled', use_parallel, ...
                'parallel_workers', parallel_workers, 'parallel_message', parallel_message);
            cuberes.settings = EELS_PearsonSparse_FullFit.strip_runtime_settings(settings);
            EELS_PearsonSparse_FullFit.notify_progress(settings, np, np, toc(run_timer), pixel_time, 'saving result MAT');
            cuberes = EELS_PearsonSparse_FullFit.sanitize_for_save(cuberes);
            save(fullfile(settings.test_output_dir, sprintf('%s_datacube_result.mat', cube_tag)), 'cuberes', '-v7.3');
            EELS_PearsonSparse_FullFit.notify_progress(settings, np, np, toc(run_timer), pixel_time, 'writing summary');
            EELS_PearsonSparse_FullFit.write_datacube_summary(fullfile(settings.test_output_dir, ...
                sprintf('%s_datacube_summary.txt', cube_tag)), cuberes);
            if settings.save_plots
                EELS_PearsonSparse_FullFit.notify_progress(settings, np, np, toc(run_timer), pixel_time, 'saving datacube map plot');
                fig = EELS_PearsonSparse_FullFit.plot_datacube_maps(cuberes);
                saveas(fig, fullfile(settings.test_output_dir, 'plots', sprintf('%s_maps.png', cube_tag)));
                close(fig);
            end
            EELS_PearsonSparse_FullFit.notify_progress(settings, np, np, toc(run_timer), pixel_time, 'complete');
        end

        function cuberes = rebuild_datacube_maps_from_results(cuberes, settings, pixel_time, run_timer, reused_stage1)
            if nargin < 5, reused_stage1 = false; end
            pixels = cuberes.sample_pixels(:).';
            np = numel(pixels);
            if isfield(cuberes, 'R2_map')
                [ny, nx] = size(cuberes.R2_map);
            elseif isfield(cuberes, 'cube_info') && isfield(cuberes.cube_info, 'size')
                ny = cuberes.cube_info.size(1);
                nx = cuberes.cube_info.size(2);
            else
                error('Cannot infer datacube map size from stage1 cuberes.');
            end
            nbins = settings.spatial_attention_profile_bins;
            maps = struct();
            maps.R2_map = NaN(ny, nx);
            maps.NRMSE_map = NaN(ny, nx);
            maps.T_map_raw = NaN(ny, nx);
            maps.T_uncertainty_map = NaN(ny, nx);
            maps.reliability_score_map = NaN(ny, nx);
            maps.reliability_source_map = repmat({''}, ny, nx);
            maps.stage1_PM_ratio_T_map = NaN(ny, nx);
            maps.stage1_PM_ratio_reliability_map = NaN(ny, nx);
            maps.stage1_PM_ratio_pairs_map = NaN(ny, nx);
            maps.stage1_PM_ratio_bestline_R2_map = NaN(ny, nx);
            maps.stage1_PM_ratio_selected_pairs_map = NaN(ny, nx);
            maps.stage1_PM_ratio_selected_EM_lo_map = NaN(ny, nx);
            maps.stage1_PM_ratio_selected_EM_hi_map = NaN(ny, nx);
            maps.stage2_zlp_refine_max_delta_map = NaN(ny, nx);
            maps.stage2_zlp_refine_fit_rms_map = NaN(ny, nx);
            maps.attention_score_map = NaN(ny, nx);
            maps.zlp_guard_start_map = NaN(ny, nx);
            maps.zlp_guard_fraction_map = NaN(ny, nx);
            maps.zlp_core_width_map = NaN(ny, nx);
            maps.attention_profile_cube = NaN(ny, nx, nbins);
            maps.signal_profile_cube = NaN(ny, nx, nbins);
            maps.fit_status_map = zeros(ny, nx);
            for ii = 1:np
                res = cuberes.results{ii};
                [yy, xx] = ind2sub([ny, nx], pixels(ii));
                if isstruct(res) && ~isfield(res, 'error')
                    maps = EELS_PearsonSparse_FullFit.update_datacube_maps_from_result( ...
                        maps, res, yy, xx, settings, nbins, 1);
                else
                    maps.fit_status_map(yy, xx) = -1;
                end
            end
            spatial = EELS_PearsonSparse_FullFit.evaluate_spatial_consistency(maps.T_map_raw, maps.R2_map, ...
                maps.NRMSE_map, maps.reliability_score_map, maps.attention_score_map, maps.fit_status_map, settings, ...
                maps.attention_profile_cube, maps.signal_profile_cube);
            cuberes.R2_map = maps.R2_map;
            cuberes.NRMSE_map = maps.NRMSE_map;
            cuberes.T_map_raw = maps.T_map_raw;
            cuberes.T_map_repaired = spatial.T_map_repaired;
            cuberes.T_map_spatial_refined = spatial.T_map_spatial_refined;
            cuberes.spatial_refine_delta_map = spatial.spatial_refine_delta_map;
            cuberes.spatial_refine_mask = spatial.spatial_refine_mask;
            cuberes.T_uncertainty_map = maps.T_uncertainty_map;
            cuberes.reliability_score_map = maps.reliability_score_map;
            cuberes.reliability_source_map = maps.reliability_source_map;
            cuberes.stage1_PM_ratio_T_map = maps.stage1_PM_ratio_T_map;
            cuberes.stage1_PM_ratio_reliability_map = maps.stage1_PM_ratio_reliability_map;
            cuberes.stage1_PM_ratio_pairs_map = maps.stage1_PM_ratio_pairs_map;
            cuberes.stage1_PM_ratio_bestline_R2_map = maps.stage1_PM_ratio_bestline_R2_map;
            cuberes.stage1_PM_ratio_selected_pairs_map = maps.stage1_PM_ratio_selected_pairs_map;
            cuberes.stage1_PM_ratio_selected_EM_lo_map = maps.stage1_PM_ratio_selected_EM_lo_map;
            cuberes.stage1_PM_ratio_selected_EM_hi_map = maps.stage1_PM_ratio_selected_EM_hi_map;
            cuberes.stage2_zlp_refine_max_delta_map = maps.stage2_zlp_refine_max_delta_map;
            cuberes.stage2_zlp_refine_fit_rms_map = maps.stage2_zlp_refine_fit_rms_map;
            cuberes.attention_score_map = maps.attention_score_map;
            cuberes.zlp_guard_start_map = maps.zlp_guard_start_map;
            cuberes.zlp_guard_fraction_map = maps.zlp_guard_fraction_map;
            cuberes.zlp_core_width_map = maps.zlp_core_width_map;
            cuberes.attention_profile_cube = maps.attention_profile_cube;
            cuberes.signal_profile_cube = maps.signal_profile_cube;
            cuberes.warning_map = spatial.warning_map;
            cuberes.bad_pixel_mask = spatial.bad_pixel_mask;
            cuberes.interface_mask = spatial.interface_mask;
            cuberes.neighbor_loss_map = spatial.neighbor_loss_map;
            cuberes.spectrum_neighbor_loss_map = spatial.spectrum_neighbor_loss_map;
            cuberes.temperature_neighbor_z_map = spatial.temperature_neighbor_z_map;
            cuberes.attention_peak_delta_map = spatial.attention_peak_delta_map;
            cuberes.signal_peak_delta_map = spatial.signal_peak_delta_map;
            cuberes.spatial_refit_candidate_mask = spatial.spatial_refit_candidate_mask;
            cuberes.low_R2_candidate_mask = spatial.low_R2_candidate_mask;
            cuberes.region_label_map = spatial.region_label_map;
            cuberes.region_median_map = spatial.region_median_map;
            cuberes.region_delta_map = spatial.region_delta_map;
            cuberes.region_z_map = spatial.region_z_map;
            cuberes.region_size_map = spatial.region_size_map;
            cuberes.region_post_smooth_mask = spatial.region_post_smooth_mask;
            cuberes.fit_status_map = maps.fit_status_map;
            cuberes.spatial = spatial;
            cuberes.reused_stage1_for_stage2 = reused_stage1;
            elapsed = toc(run_timer);
            cuberes.benchmark = struct('total_seconds', elapsed, 'median_seconds_per_pixel', median(pixel_time, 'omitnan'), ...
                'mean_seconds_per_pixel', mean(pixel_time, 'omitnan'), 'pixel_seconds', pixel_time, ...
                'pixels_per_second', np ./ max(elapsed, eps), 'parallel_enabled', false, ...
                'parallel_workers', 0, 'parallel_message', 'stage2 reused existing stage1 results');
            cuberes.settings = EELS_PearsonSparse_FullFit.strip_runtime_settings(settings);
            EELS_PearsonSparse_FullFit.notify_progress(settings, numel(pixels), numel(pixels), elapsed, pixel_time, 'saving Stage2 result MAT');
            cuberes = EELS_PearsonSparse_FullFit.sanitize_for_save(cuberes);
            save(fullfile(settings.test_output_dir, 'stage2_reused_datacube_result.mat'), 'cuberes', '-v7.3');
            EELS_PearsonSparse_FullFit.notify_progress(settings, numel(pixels), numel(pixels), elapsed, pixel_time, 'writing Stage2 summary');
            EELS_PearsonSparse_FullFit.write_datacube_summary(fullfile(settings.test_output_dir, ...
                'stage2_reused_datacube_summary.txt'), cuberes);
            if settings.save_plots
                plot_dir = fullfile(settings.test_output_dir, 'plots');
                if ~isfolder(plot_dir), mkdir(plot_dir); end
                EELS_PearsonSparse_FullFit.notify_progress(settings, numel(pixels), numel(pixels), elapsed, pixel_time, 'saving Stage2 map plot');
                fig = EELS_PearsonSparse_FullFit.plot_datacube_maps(cuberes);
                saveas(fig, fullfile(plot_dir, 'stage2_reused_maps.png'));
                close(fig);
            end
            EELS_PearsonSparse_FullFit.notify_progress(settings, numel(pixels), numel(pixels), elapsed, pixel_time, 'complete');
        end

        function settings_out = strip_runtime_settings(settings)
            settings_out = EELS_PearsonSparse_FullFit.sanitize_for_save(settings);
            if isfield(settings_out, 'datacube_progress_callback')
                settings_out.datacube_progress_callback = [];
            end
            if isfield(settings_out, 'datacube_stop_callback')
                settings_out.datacube_stop_callback = [];
            end
        end

        function out = sanitize_for_save(in, depth)
            if nargin < 2, depth = 0; end
            if depth > 12
                out = [];
                return;
            end
            if isa(in, 'function_handle')
                out = [];
                return;
            end
            if isobject(in)
                if isa(in, 'datetime') || isa(in, 'duration') || isa(in, 'calendarDuration')
                    out = in;
                else
                    out = [];
                end
                return;
            end
            if isstruct(in)
                out = in;
                f = fieldnames(in);
                for kk = 1:numel(in)
                    for ff = 1:numel(f)
                        out(kk).(f{ff}) = EELS_PearsonSparse_FullFit.sanitize_for_save(in(kk).(f{ff}), depth + 1);
                    end
                end
                return;
            end
            if iscell(in)
                out = cell(size(in));
                for kk = 1:numel(in)
                    out{kk} = EELS_PearsonSparse_FullFit.sanitize_for_save(in{kk}, depth + 1);
                end
                return;
            end
            out = in;
        end

        function tf = datacube_stop_requested(settings)
            tf = false;
            if isfield(settings, 'datacube_stop_callback') && isa(settings.datacube_stop_callback, 'function_handle')
                try
                    tf = logical(settings.datacube_stop_callback());
                catch
                    tf = false;
                end
            end
        end

        function save_datacube_checkpoint(settings, np, results, pixel_time, pixels, completed)
            completed = max(0, min(round(completed), np));
            checkpoint_settings = EELS_PearsonSparse_FullFit.strip_runtime_settings(settings);
            checkpoint = struct('completed', completed, 'total', np, ...
                'results', {results(1:completed)}, ...
                'pixel_time', pixel_time(1:completed), ...
                'pixels', pixels(1:completed), 'settings', checkpoint_settings);
            save(fullfile(settings.test_output_dir, 'datacube_fit_checkpoint.mat'), 'checkpoint', '-v7.3');
        end

        function notify_progress(settings, done, total, elapsed, pixel_times, label)
            if nargin < 6, label = 'datacube'; end
            med_sec = median(pixel_times(:), 'omitnan');
            worker_factor = 1;
            if isfield(settings, 'datacube_parallel_effective_workers') && ...
                    ~isempty(settings.datacube_parallel_effective_workers) && ...
                    all(isfinite(settings.datacube_parallel_effective_workers))
                worker_factor = max(1, round(settings.datacube_parallel_effective_workers));
            elseif contains(lower(char(string(label))), 'parallel') && ...
                    isfield(settings, 'datacube_parallel_workers') && ~isempty(settings.datacube_parallel_workers) && ...
                    all(isfinite(settings.datacube_parallel_workers)) && settings.datacube_parallel_workers > 0
                worker_factor = max(1, round(settings.datacube_parallel_workers));
            end
            if isfinite(med_sec)
                eta = med_sec .* max(total - done, 0) ./ worker_factor;
            else
                eta = NaN;
            end
            if isfield(settings, 'datacube_progress_callback') && isa(settings.datacube_progress_callback, 'function_handle')
                try
                    settings.datacube_progress_callback(done, total, elapsed, eta, med_sec, label);
                catch
                end
            end
            if settings.datacube_progress_every > 0 && (mod(done, settings.datacube_progress_every) == 0 || done == total)
                if worker_factor > 1
                    fprintf('%s progress %d/%d pixels | elapsed %.1f s | median %.2f s/pixel | workers %d | ETA %.1f s\n', ...
                        label, done, total, elapsed, med_sec, worker_factor, eta);
                else
                    fprintf('%s progress %d/%d pixels | elapsed %.1f s | median %.2f s/pixel | ETA %.1f s\n', ...
                        label, done, total, elapsed, med_sec, eta);
                end
            end
        end

        function [use_parallel, workers, message] = datacube_parallel_state(settings, np)
            use_parallel = false;
            workers = 0;
            message = '';
            if ~settings.datacube_parallel_enabled
                message = 'parallel flag is off';
                return;
            end
            if np < settings.datacube_parallel_min_pixels
                message = sprintf('pixel count %d < datacube_parallel_min_pixels %d', np, settings.datacube_parallel_min_pixels);
                return;
            end
            if exist('parpool', 'file') ~= 2 || exist('gcp', 'file') ~= 2
                message = 'Parallel Computing Toolbox functions are unavailable';
                return;
            end
            try
                if ~license('test', 'Distrib_Computing_Toolbox')
                    message = 'Parallel Computing Toolbox license is unavailable';
                    return;
                end
            catch
            end
            if settings.datacube_parallel_workers > 0
                workers = round(settings.datacube_parallel_workers);
            else
                try
                    workers = max(1, feature('numcores') - 1);
                catch
                    workers = 4;
                end
                workers = min(workers, 8);
            end
            workers = max(1, min([workers, np]));
            use_parallel = workers > 1;
            if ~use_parallel
                message = 'parallel worker count resolved to 1';
            else
                message = sprintf('using %d workers', workers);
            end
        end

        function settings = calibrate_datacube_stage1_sparse_pm(cube, ene, pixels, settings)
            if ~settings.datacube_stage1_auto_sparse_pm_enabled || ~isempty(settings.stage1_fixed_EM) || isempty(pixels)
                return;
            end
            policy = struct('enabled', true, 'status', 'not_run', 'probe_pixel', NaN, ...
                'selected_count', NaN, 'reference_count', NaN, 'selected_R2', NaN, ...
                'reference_R2', NaN, 'message', '');
            try
                [ny, nx, ~] = size(cube);
                probe_pixel = EELS_PearsonSparse_FullFit.pick_stage1_sparse_probe_pixel(cube, pixels, settings);
                [yy, xx] = ind2sub([ny, nx], probe_pixel);
                probe_settings = settings;
                probe_settings.datacube_stage1_auto_sparse_pm_enabled = false;
                probe_settings.stage1_fixed_EM = [];
                probe_settings.run_stage2 = false;
                probe_settings.save_plots = false;
                probe_settings.save_diagnostic_plots = false;
                probe_settings.datacube_progress_callback = [];
                probe_settings.datacube_stop_callback = [];
                spectrum = squeeze(cube(yy, xx, :));
                probe_res = EELS_PearsonSparse_FullFit.solve_stage1(ene, spectrum, probe_settings);
                [selected_EM, sparse_policy] = EELS_PearsonSparse_FullFit.estimate_stage1_sparse_em_from_probe(probe_res, settings);
                if ~isempty(selected_EM)
                    settings.stage1_fixed_EM = selected_EM(:).';
                    policy = sparse_policy;
                    policy.enabled = true;
                    policy.status = 'ok';
                    policy.probe_pixel = probe_pixel;
                    policy.probe_y = yy;
                    policy.probe_x = xx;
                    policy.message = sprintf('Stage1 sparse source-intensity grid calibrated from pixel %d (%d,%d): %d -> %d mode-energy centers.', ...
                        probe_pixel, yy, xx, policy.reference_count, policy.selected_count);
                    if ~settings.silent
                        fprintf('%s\n', policy.message);
                    end
                else
                    policy.status = 'skipped';
                    policy.probe_pixel = probe_pixel;
                    policy.message = 'No sparse EM candidate preserved the probe Stage1 vib shape.';
                end
            catch ME
                policy.status = 'failed';
                policy.message = ME.message;
            end
            settings.stage1_sparse_pm_policy = policy;
        end

        function probe_pixel = pick_stage1_sparse_probe_pixel(cube, pixels, settings)
            if ~isempty(settings.stage1_auto_sparse_pm_probe_pixel) && ...
                    any(pixels(:) == settings.stage1_auto_sparse_pm_probe_pixel)
                probe_pixel = settings.stage1_auto_sparse_pm_probe_pixel;
                return;
            end
            [ny, nx, ~] = size(cube);
            [yy, xx] = ind2sub([ny, nx], pixels(:));
            cy = (ny + 1) ./ 2;
            cx = (nx + 1) ./ 2;
            center_dist = hypot(double(yy) - cy, double(xx) - cx);
            score = -center_dist(:).';
            try
                eabs = abs(settings.crop_range);
                use_energy = true;
                if any(~isfinite(eabs)), use_energy = false; end
                for ii = 1:numel(pixels)
                    spec = squeeze(cube(yy(ii), xx(ii), :));
                    if use_energy
                        amp = max(abs(spec), [], 'omitnan');
                    else
                        amp = max(abs(spec), [], 'omitnan');
                    end
                    if isfinite(amp)
                        score(ii) = score(ii) + 0.01 .* amp ./ max(abs(amp), eps);
                    end
                end
            catch
            end
            [~, pick] = max(score);
            if isempty(pick) || ~isfinite(pick)
                pick = ceil(numel(pixels) ./ 2);
            end
            probe_pixel = pixels(pick);
        end

        function [selected_EM, policy] = estimate_stage1_sparse_em_from_probe(probe_res, settings)
            selected_EM = [];
            policy = struct('status', 'not_run', 'reference_count', NaN, 'selected_count', NaN, ...
                'reference_R2', NaN, 'selected_R2', NaN, 'candidate_counts', [], ...
                'candidate_R2', [], 'R2_threshold', NaN);
            if ~(isstruct(probe_res) && isfield(probe_res, 'stage1') && isfield(probe_res.stage1, 'vib_norm'))
                policy.status = 'missing_probe_result';
                return;
            end
            EM = probe_res.EM_final(:).';
            EM = unique(EM(isfinite(EM) & EM > 0));
            if numel(EM) <= settings.stage1_auto_sparse_pm_min_count
                policy.status = 'already_sparse';
                selected_EM = EM;
                policy.reference_count = numel(EM);
                policy.selected_count = numel(EM);
                return;
            end
            ene = probe_res.ene_raw(:);
            ene_corr = probe_res.ene_corrected(:);
            target = probe_res.stage1.vib_norm(:);
            params = probe_res.stage1.params;
            [~, ~, ~, kernels] = EELS_PearsonSparse_FullFit.build_design_matrix(ene, EM, params, settings);
            pm_s = EELS_PearsonSparse_FullFit.interpolate_pm(probe_res.EM_final(:), probe_res.PM_stokes(:), EM(:));
            pm_a = EELS_PearsonSparse_FullFit.interpolate_pm(probe_res.EM_final(:), probe_res.PM_antistokes(:), EM(:));
            pm_mass = max(abs(pm_s(:).'), abs(pm_a(:).'));
            signal = abs(probe_res.raw_data_norm(:) - probe_res.stage1.zlp_norm(:));
            signal = movmedian(signal, 7, 'omitnan');
            feature_pos = interp1(ene_corr, signal, EM, 'linear', 0);
            feature_neg = interp1(ene_corr, signal, -EM, 'linear', 0);
            feature = max(feature_pos, feature_neg);
            feature = EELS_PearsonSparse_FullFit.unit_scale(feature(:), true(numel(feature), 1)).';
            score = EELS_PearsonSparse_FullFit.unit_scale(pm_mass(:), true(numel(pm_mass), 1)).' + 0.45 .* feature;
            protect = score >= settings.stage1_auto_sparse_pm_feature_threshold | ...
                pm_mass >= 0.02 .* max(pm_mass, [], 'omitnan');
            protected_EM = zeros(1, 0);
            if any(protect)
                protect_score = score(protect) + 0.25 .* feature(protect);
                protected_EM = EELS_PearsonSparse_FullFit.select_spaced_em(EM(protect), protect_score, ...
                    min(settings.stage1_auto_sparse_pm_protect_max_centers, nnz(protect)), ...
                    max(0.75, 0.75 .* settings.stage1_auto_sparse_pm_min_spacing_meV));
            end
            sqrtW = EELS_PearsonSparse_FullFit.build_weights(ene_corr, probe_res.raw_data_norm(:), settings);
            fit_mask = abs(ene_corr) >= max(4, settings.em_min_meV) & isfinite(target);
            sqrtW(~fit_mask) = 0;
            reference_R2 = EELS_PearsonSparse_FullFit.stage2_pruned_basis_R2(kernels, true(size(EM)), target, sqrtW);
            max_count = min(numel(EM), settings.stage1_auto_sparse_pm_max_count);
            min_count = min(max(4, settings.stage1_auto_sparse_pm_min_count), max_count);
            counts = unique([min_count:settings.stage1_auto_sparse_pm_step:max_count, max_count]);
            trials = struct('count', {}, 'selected_count', {}, 'R2', {}, 'EM', {});
            for ii = 1:numel(counts)
                cand = EELS_PearsonSparse_FullFit.select_spaced_em(EM, score, counts(ii), ...
                    settings.stage1_auto_sparse_pm_min_spacing_meV);
                cand = sort(unique([cand(:); protected_EM(:)]).');
                keep = ismembertol(EM, cand, 1e-8);
                R2 = EELS_PearsonSparse_FullFit.stage2_pruned_basis_R2(kernels, keep, target, sqrtW);
                trials(end + 1).count = counts(ii); %#ok<AGROW>
                trials(end).selected_count = nnz(keep);
                trials(end).R2 = R2;
                trials(end).EM = EM(keep);
            end
            if isempty(trials)
                policy.status = 'no_trials';
                return;
            end
            R2s = [trials.R2];
            maxR2 = max(R2s, [], 'omitnan');
            if isfinite(reference_R2)
                threshold = max(settings.stage1_auto_sparse_pm_min_R2, reference_R2 - settings.stage1_auto_sparse_pm_R2_tolerance);
            else
                threshold = max(settings.stage1_auto_sparse_pm_min_R2, maxR2 - settings.stage1_auto_sparse_pm_R2_tolerance);
            end
            eligible = find(R2s >= threshold);
            if isempty(eligible)
                [~, pick] = max(R2s);
                policy.status = 'best_available';
            else
                [~, ord] = sort([trials(eligible).selected_count], 'ascend');
                pick = eligible(ord(1));
                policy.status = 'sparsest_within_R2_tolerance';
            end
            selected_EM = trials(pick).EM;
            policy.reference_count = numel(EM);
            policy.selected_count = numel(selected_EM);
            policy.reference_R2 = reference_R2;
            policy.selected_R2 = trials(pick).R2;
            policy.candidate_counts = [trials.selected_count];
            policy.candidate_R2 = R2s;
            policy.R2_threshold = threshold;
        end

        function R2 = stage2_pruned_basis_R2(kernels, keep, target, sqrtW)
            R2 = NaN;
            if ~(isstruct(kernels) && isfield(kernels, 'stokes') && isfield(kernels, 'anti'))
                return;
            end
            target = target(:);
            sqrtW = sqrtW(:);
            n = min([numel(keep), size(kernels.stokes, 2), size(kernels.anti, 2)]);
            if n <= 0 || numel(target) ~= size(kernels.stokes, 1)
                return;
            end
            keep = logical(keep(:));
            keep = keep(1:n);
            if ~any(keep)
                return;
            end
            stokes = kernels.stokes(:, 1:n);
            anti = kernels.anti(:, 1:n);
            X = [stokes(:, keep), anti(:, keep)];
            X(~isfinite(X)) = 0;
            if numel(sqrtW) ~= numel(target)
                sqrtW = ones(size(target));
            end
            good = isfinite(target) & isfinite(sqrtW) & sqrtW > 0;
            if nnz(good) < 3 || isempty(X)
                return;
            end
            C = X(good, :) .* sqrtW(good);
            d = target(good) .* sqrtW(good);
            if isempty(C) || size(C, 2) == 0
                return;
            end
            coef = lsqnonneg(C, d);
            model = X(good, :) * coef;
            y = target(good);
            w = sqrtW(good) .^ 2;
            wsum = sum(w, 'omitnan');
            if ~(isfinite(wsum) && wsum > 0)
                return;
            end
            ybar = sum(w .* y, 'omitnan') ./ wsum;
            ss_tot = sum(w .* (y - ybar) .^ 2, 'omitnan');
            ss_res = sum(w .* (model - y) .^ 2, 'omitnan');
            if ss_tot <= eps
                R2 = double(ss_res <= 1e-12);
            else
                R2 = 1 - ss_res ./ ss_tot;
            end
        end

        function policy = sanitize_stage1_sparse_pm_policy_for_stage2_reuse(policy)
            if ~isstruct(policy) || ~isfield(policy, 'status')
                return;
            end
            msg = '';
            if isfield(policy, 'message')
                msg = char(string(policy.message));
            end
            failed = strcmpi(char(string(policy.status)), 'failed');
            helper_missing = contains(msg, 'stage2_pruned_basis_R2') || ...
                contains(msg, '没有名为 ''stage2_pruned_basis_R2''');
            if failed && helper_missing
                policy.status = 'not_rerun_stage2_reuse_helper_restored';
                policy.message = ['Previous Stage1 sparse source-intensity calibration failed because ', ...
                    'stage2_pruned_basis_R2 was missing; the helper is now available, ', ...
                    'but this Stage2 reuse run did not rerun Stage1 sparse calibration.'];
            end
        end

        function [results, pixel_time] = fit_datacube_pixels_serial(cube, ene, pixels, settings, run_timer)
            np = numel(pixels);
            results = cell(1, np);
            pixel_time = NaN(1, np);
            prev_theta = [];
            [ny, nx, ~] = size(cube);
            start_index = 1;
            checkpoint_file = fullfile(settings.test_output_dir, 'datacube_fit_checkpoint.mat');
            if settings.datacube_resume_checkpoint && isfile(checkpoint_file)
                try
                    C = load(checkpoint_file, 'checkpoint');
                    cp = C.checkpoint;
                    same_total = isfield(cp, 'total') && isequal(double(cp.total), double(np));
                    same_pixels = false;
                    if isfield(cp, 'pixels') && numel(cp.pixels) <= numel(pixels)
                        pixel_prefix = pixels(1:numel(cp.pixels));
                        same_pixels = isequal(cp.pixels(:), pixel_prefix(:));
                    end
                    if same_total && same_pixels && isfield(cp, 'results')
                        done = min([numel(cp.results), np]);
                        results(1:done) = cp.results(1:done);
                        if isfield(cp, 'pixel_time')
                            pixel_time(1:min(numel(cp.pixel_time), done)) = cp.pixel_time(1:min(numel(cp.pixel_time), done));
                        end
                        start_index = done + 1;
                        if done > 0 && isstruct(results{done}) && ~isfield(results{done}, 'error') && isfield(results{done}, 'theta')
                            prev_theta = results{done}.theta;
                        end
                        fprintf('datacube resume checkpoint: continuing from pixel %d/%d\n', start_index, np);
                    end
                catch ME
                    fprintf('datacube resume checkpoint ignored: %s\n', ME.message);
                end
            end
            for ii = start_index:np
                if EELS_PearsonSparse_FullFit.datacube_stop_requested(settings)
                    EELS_PearsonSparse_FullFit.save_datacube_checkpoint(settings, np, results, pixel_time, pixels, ii - 1);
                    error('EELS:DatacubePaused', 'Datacube fitting paused after pixel %d/%d. Checkpoint saved in %s.', ...
                        ii - 1, np, settings.test_output_dir);
                end
                [yy, xx] = ind2sub([ny, nx], pixels(ii));
                spectrum = squeeze(cube(yy, xx, :));
                [results{ii}, pixel_time(ii)] = EELS_PearsonSparse_FullFit.fit_one_datacube_pixel( ...
                    ene, spectrum, pixels(ii), yy, xx, settings, prev_theta);
                if isstruct(results{ii}) && ~isfield(results{ii}, 'error') && isfield(results{ii}, 'theta')
                    prev_theta = results{ii}.theta;
                end
                if settings.datacube_checkpoint_enabled && settings.datacube_checkpoint_every > 0 && ...
                        (mod(ii, settings.datacube_checkpoint_every) == 0 || ii == np)
                    EELS_PearsonSparse_FullFit.save_datacube_checkpoint(settings, np, results, pixel_time, pixels, ii);
                end
                EELS_PearsonSparse_FullFit.notify_progress(settings, ii, np, toc(run_timer), pixel_time(1:ii), 'datacube');
                if EELS_PearsonSparse_FullFit.datacube_stop_requested(settings)
                    EELS_PearsonSparse_FullFit.save_datacube_checkpoint(settings, np, results, pixel_time, pixels, ii);
                    error('EELS:DatacubePaused', 'Datacube fitting paused after pixel %d/%d. Checkpoint saved in %s.', ...
                        ii, np, settings.test_output_dir);
                end
            end
        end

        function [results, pixel_time] = fit_datacube_pixels_parallel(cube, ene, pixels, settings, workers)
            np = numel(pixels);
            EELS_PearsonSparse_FullFit.ensure_parallel_pool(workers);
            progress_settings = settings;
            progress_settings.datacube_parallel_effective_workers = max(1, round(workers));
            parallel_pool_cleanup = []; %#ok<NASGU>
            if isfield(settings, 'close_parallel_pool_after_run') && settings.close_parallel_pool_after_run
                parallel_pool_cleanup = onCleanup(@()EELS_PearsonSparse_FullFit.close_parallel_pool());
            end
            chunks = EELS_PearsonSparse_FullFit.datacube_pixel_chunks(np, workers, settings.datacube_parallel_chunks);
            n_chunks = numel(chunks);
            chunk_pixels = cell(1, n_chunks);
            for cc = 1:n_chunks
                chunk_pixels{cc} = pixels(chunks{cc});
            end
            chunk_results = cell(1, n_chunks);
            chunk_times = cell(1, n_chunks);
            [ny, nx, ~] = size(cube);
            worker_settings = settings;
            worker_settings.datacube_progress_callback = [];
            worker_settings.datacube_stop_callback = [];
            progress_queue = [];
            progress_timer = tic;
            try
                progress_queue = parallel.pool.DataQueue;
                completed_pixels = 0;
                observed_times = [];
                afterEach(progress_queue, @(t)EELS_PearsonSparse_FullFit.parallel_progress_update( ...
                    progress_settings, t, np, progress_timer, completed_pixels, observed_times));
            catch
                progress_queue = [];
            end
            parfor cc = 1:n_chunks
                [chunk_results{cc}, chunk_times{cc}] = EELS_PearsonSparse_FullFit.fit_datacube_pixel_chunk( ...
                    cube, ene, chunk_pixels{cc}, ny, nx, worker_settings, progress_queue);
            end
            results = cell(1, np);
            pixel_time = NaN(1, np);
            for cc = 1:n_chunks
                idx = chunks{cc};
                results(idx) = chunk_results{cc};
                pixel_time(idx) = chunk_times{cc};
            end
        end

        function parallel_progress_update(settings, new_times, total, timer0, completed_pixels, observed_times)
            persistent par_done par_times par_total par_timer_marker
            if isempty(par_timer_marker) || ~isequal(par_total, total) || ~isequal(par_timer_marker, timer0)
                par_done = 0;
                par_times = [];
                par_total = total;
                par_timer_marker = timer0;
            end
            if ~isempty(completed_pixels) && completed_pixels > par_done
                par_done = completed_pixels;
            end
            if ~isempty(observed_times)
                par_times = observed_times;
            end
            new_times = double(new_times(:).');
            new_times = new_times(isfinite(new_times));
            par_done = min(total, par_done + numel(new_times));
            par_times = [par_times, new_times]; %#ok<AGROW>
            EELS_PearsonSparse_FullFit.notify_progress(settings, par_done, total, toc(timer0), par_times, 'datacube parallel');
            if par_done >= total
                par_done = [];
                par_times = [];
                par_total = [];
                par_timer_marker = [];
            end
        end

        function ensure_parallel_pool(workers)
            pool = gcp('nocreate');
            if ~isempty(pool)
                try
                    if pool.NumWorkers ~= workers
                        delete(pool);
                        pool = [];
                    end
                catch
                    pool = [];
                end
            end
            if isempty(pool)
                parpool('local', workers);
            end
        end

        function close_parallel_pool()
            pool = gcp('nocreate');
            if ~isempty(pool)
                delete(pool);
            end
        end

        function cleanup = hidden_figure_cleanup(settings)
            cleanup = [];
            if isfield(settings, 'force_hidden_figures') && settings.force_hidden_figures
                old_visible = get(groot, 'DefaultFigureVisible');
                set(groot, 'DefaultFigureVisible', 'off');
                cleanup = onCleanup(@()set(groot, 'DefaultFigureVisible', old_visible));
            end
        end

        function chunks = datacube_pixel_chunks(np, workers, requested_chunks)
            if requested_chunks > 0
                n_chunks = max(1, min(np, round(requested_chunks)));
            else
                n_chunks = max(workers, min(np, workers .* 3));
            end
            edges = unique(round(linspace(1, np + 1, n_chunks + 1)));
            chunks = {};
            for ii = 1:(numel(edges) - 1)
                idx = edges(ii):(edges(ii + 1) - 1);
                if ~isempty(idx)
                    chunks{end + 1} = idx; %#ok<AGROW>
                end
            end
        end

        function [results, pixel_time] = fit_datacube_pixel_chunk(cube, ene, chunk_pixels, ny, nx, settings, progress_queue)
            if nargin < 7
                progress_queue = [];
            end
            n = numel(chunk_pixels);
            results = cell(1, n);
            pixel_time = NaN(1, n);
            prev_theta = [];
            for ii = 1:n
                [yy, xx] = ind2sub([ny, nx], chunk_pixels(ii));
                spectrum = squeeze(cube(yy, xx, :));
                warm_theta = [];
                if settings.datacube_parallel_warm_start_within_chunk
                    warm_theta = prev_theta;
                end
                [results{ii}, pixel_time(ii)] = EELS_PearsonSparse_FullFit.fit_one_datacube_pixel( ...
                    ene, spectrum, chunk_pixels(ii), yy, xx, settings, warm_theta);
                if ~isempty(progress_queue)
                    send(progress_queue, pixel_time(ii));
                end
                if isstruct(results{ii}) && ~isfield(results{ii}, 'error') && isfield(results{ii}, 'theta')
                    prev_theta = results{ii}.theta;
                end
            end
        end

        function [res, elapsed] = fit_one_datacube_pixel(ene, spectrum, pixel_index, yy, xx, settings, warm_theta)
            if nargin < 7, warm_theta = []; end
            n_attempts = max(1, 1 + round(settings.datacube_fit_retries));
            elapsed_timer = tic;
            last_error = [];
            for attempt = 1:n_attempts
                pixel_settings = settings;
                pixel_settings.datacube_progress_callback = [];
                pixel_settings.datacube_stop_callback = [];
                if settings.datacube_use_warm_start && attempt == 1 && ~isempty(warm_theta)
                    pixel_settings.initial_theta = warm_theta;
                else
                    pixel_settings.initial_theta = [];
                end
                if attempt > 1 && settings.datacube_retry_with_strict_nnls
                    pixel_settings.nonlinear_inner_solver = 'lsqnonneg';
                    pixel_settings.stage2_scan_inner_solver = 'lsqlin';
                end
                try
                    res = EELS_PearsonSparse_FullFit.solve_stage1(ene, spectrum, pixel_settings);
                    res = EELS_PearsonSparse_FullFit.sanitize_for_save(res);
                    res.pixel_y = yy;
                    res.pixel_x = xx;
                    res.pixel_index = pixel_index;
                    res.fit_attempts = attempt;
                    elapsed = toc(elapsed_timer);
                    return;
                catch ME
                    last_error = ME;
                end
            end
            elapsed = toc(elapsed_timer);
            if isempty(last_error)
                msg = 'unknown fitting failure';
            else
                msg = last_error.message;
            end
            res = struct('error', msg, 'pixel_y', yy, 'pixel_x', xx, ...
                'pixel_index', pixel_index, 'fit_attempts', n_attempts);
        end

        function [ene, y, y_norm, norm_val] = prepare_spectrum(ene_raw, y_raw, settings)
            ene_raw = double(ene_raw(:));
            y_raw = double(y_raw(:));
            mask = isfinite(ene_raw) & isfinite(y_raw) & ene_raw >= settings.crop_range(1) & ene_raw <= settings.crop_range(2);
            ene = ene_raw(mask);
            y = y_raw(mask);
            y(y < 0 & abs(y) < 100 * eps(max(abs(y)))) = 0;
            norm_val = max(abs(y));
            if ~isfinite(norm_val) || norm_val <= 0, norm_val = 1; end
            y_norm = y ./ norm_val;
        end

        function [ene_aligned, info] = apply_stage1_zlp_alignment(ene, y_norm, settings)
            ene = double(ene(:));
            y_norm = double(y_norm(:));
            ene_aligned = ene;
            info = struct('enabled', logical(settings.stage1_zlp_align_enabled), ...
                'status', 'disabled', 'estimated_peak_meV', NaN, ...
                'applied_shift_meV', 0, ...
                'max_shift_meV', settings.stage1_zlp_align_max_shift_meV, ...
                'search_meV', settings.stage1_zlp_align_search_meV, ...
                'method', 'quadratic_peak_near_zero');
            if ~settings.stage1_zlp_align_enabled
                return;
            end
            max_shift = max(0, abs(settings.stage1_zlp_align_max_shift_meV));
            if max_shift <= 0 || numel(ene) < 3 || numel(ene) ~= numel(y_norm)
                info.status = 'insufficient_data';
                return;
            end
            search_meV = max(max_shift, abs(settings.stage1_zlp_align_search_meV));
            mask = isfinite(ene) & isfinite(y_norm) & abs(ene) <= search_meV;
            if nnz(mask) < 3
                info.status = 'insufficient_search_points';
                return;
            end
            ys = y_norm(:);
            smooth_pts = max(1, round(settings.stage1_zlp_align_smooth_points));
            if smooth_pts > 1 && numel(ys) >= smooth_pts
                ys = movmean(ys, smooth_pts, 'omitnan');
            end
            idx_all = find(mask);
            [peak_height, rel_idx] = max(ys(idx_all), [], 'omitnan');
            if isempty(rel_idx) || ~isfinite(peak_height)
                info.status = 'no_finite_peak';
                return;
            end
            global_height = max(ys(isfinite(ys)), [], 'omitnan');
            if ~isfinite(global_height) || global_height <= 0 || ...
                    peak_height < settings.stage1_zlp_align_min_peak_fraction .* global_height
                info.status = 'weak_peak_near_zero';
                return;
            end
            idx = idx_all(rel_idx);
            peak_center = ene(idx);
            if idx > 1 && idx < numel(ene)
                x = ene((idx - 1):(idx + 1));
                yy = y_norm((idx - 1):(idx + 1));
                if all(isfinite(x)) && all(isfinite(yy)) && numel(unique(x)) == 3
                    p = polyfit(x(:), yy(:), 2);
                    denom = 2 .* p(1);
                    candidate = NaN;
                    if abs(denom) > eps
                        candidate = -p(2) ./ denom;
                    end
                    if p(1) < 0 && isfinite(candidate) && candidate >= min(x) && candidate <= max(x)
                        peak_center = candidate;
                    end
                end
            end
            shift = min(max(peak_center, -max_shift), max_shift);
            if ~isfinite(shift)
                info.status = 'invalid_peak_center';
                return;
            end
            ene_aligned = ene - shift;
            info.status = 'ok';
            info.estimated_peak_meV = peak_center;
            info.applied_shift_meV = shift;
        end

        function EM = make_em_grid(ene, settings)
            if isfield(settings, 'stage1_fixed_EM') && ~isempty(settings.stage1_fixed_EM)
                max_abs_e = max(abs(ene(:)), [], 'omitnan');
                EM = double(settings.stage1_fixed_EM(:).');
                EM = unique(EM(isfinite(EM) & EM >= settings.em_min_meV & EM <= max_abs_e));
                if ~isempty(settings.em_max_meV)
                    EM = EM(EM <= settings.em_max_meV);
                end
                if ~isempty(EM)
                    return;
                end
            end
            emax = max(abs(ene));
            if isempty(settings.em_max_meV)
                max_em = max(settings.em_min_meV, emax - settings.em_exclude_core_meV);
            else
                max_em = min(settings.em_max_meV, emax);
            end
            EM = settings.em_min_meV:settings.em_step_meV:max_em;
            if settings.em_lowfreq_dense_enabled
                dense_range = sort(abs(settings.em_lowfreq_dense_range_meV(:).'));
                if numel(dense_range) >= 2
                    dense_start = max(settings.em_min_meV, dense_range(1));
                    dense_end = min(max_em, dense_range(2));
                    dense_step = min(settings.em_lowfreq_dense_max_step_meV, ...
                        max(0.75, settings.em_step_meV .* settings.em_lowfreq_dense_step_factor));
                    if dense_end > dense_start + dense_step
                        EM = [EM, dense_start:dense_step:dense_end];
                    end
                end
            end
            EM = unique(EM(EM >= settings.em_min_meV & EM <= max_em));
        end

        function mask = optimization_mask(ene, settings)
            n = numel(ene);
            mask = true(n, 1);
            max_points = settings.optimization_max_points;
            if ~isfinite(max_points) || max_points <= 0 || n <= max_points
                return;
            end
            max_points = max(20, round(max_points));
            core = abs(ene(:)) <= settings.optimization_keep_core_meV;
            keep = false(n, 1);
            keep(core) = true;
            remaining = find(~core);
            budget = max(max_points - nnz(keep), 0);
            if budget > 0 && ~isempty(remaining)
                pos = unique(round(linspace(1, numel(remaining), min(budget, numel(remaining)))));
                keep(remaining(pos)) = true;
            end
            if nnz(keep) < min(max_points, n)
                pos = unique(round(linspace(1, n, min(max_points, n))));
                keep(pos) = true;
            end
            mask = keep;
        end

        function theta0 = initial_theta(ene, y_norm, settings, prev_theta)
            nz = settings.num_zlp_components;
            if ~isempty(prev_theta) && numel(prev_theta) == 2 + 5 * nz
                theta0 = prev_theta;
                [lb, ub] = EELS_PearsonSparse_FullFit.theta_bounds(settings);
                theta0 = min(max(theta0, lb), ub);
                return;
            end
            [~, imax] = max(y_norm);
            shift0 = ene(imax);
            centers = zeros(1, nz);
            width_max = EELS_PearsonSparse_FullFit.expand_component_setting(settings.zlp_component_width_max_meV, nz, settings.zlp_width_bounds(2));
            widths = min(linspace(0.8, 26, nz), 0.65 .* width_max);
            if settings.stage1_zlp_gain_heavy_init_enabled
                widthsL = widths .* linspace(1.02, max(settings.stage1_zlp_gain_width_boost, 1.02), nz);
                widthsR = widths .* linspace(1.00, min(max(settings.stage1_zlp_loss_width_scale, 0.55), 1.25), nz);
            else
                widthsL = widths;
                widthsR = widths .* linspace(1.05, 1.8, nz);
            end
            if EELS_PearsonSparse_FullFit.is_pvoigt_model(settings.zlp_model_type)
                eta0 = linspace(0.12, 0.36, nz);
                if settings.stage1_zlp_gain_heavy_init_enabled
                    mL = min(eta0 .* linspace(1.05, 1.35, nz), settings.pvoigt_eta_bounds(2));
                    mR = max(eta0 .* linspace(1.00, 0.82, nz), settings.pvoigt_eta_bounds(1));
                else
                    mL = eta0;
                    mR = min(eta0 .* linspace(1.0, 1.35, nz), settings.pvoigt_eta_bounds(2));
                end
            else
                if settings.stage1_zlp_gain_heavy_init_enabled
                    mbase = linspace(4.4, 1.8, nz);
                    gain_shape = linspace(1.00, min(max(settings.stage1_zlp_gain_shape_heavy_factor, 0.55), 1.05), nz);
                    loss_shape = linspace(1.04, max(settings.stage1_zlp_loss_shape_light_factor, 1.04), nz);
                    mL = mbase .* gain_shape;
                    mR = mbase .* loss_shape;
                else
                    mL = linspace(4.5, 1.8, nz);
                    mR = linspace(3.8, 1.6, nz);
                end
            end
            Tsoft = settings.T_soft_init;
            theta0 = [shift0, centers, widthsL, widthsR, mL, mR, log(Tsoft)];
            [lb, ub] = EELS_PearsonSparse_FullFit.theta_bounds(settings);
            theta0 = min(max(theta0, lb), ub);
        end

        function [lb, ub] = theta_bounds(settings)
            nz = settings.num_zlp_components;
            width_min = EELS_PearsonSparse_FullFit.expand_component_setting(settings.zlp_component_width_min_meV, nz, settings.zlp_width_bounds(1));
            width_max = EELS_PearsonSparse_FullFit.expand_component_setting(settings.zlp_component_width_max_meV, nz, settings.zlp_width_bounds(2));
            width_min = max(width_min, settings.zlp_width_bounds(1));
            width_max = min(max(width_max, width_min + eps), settings.zlp_width_bounds(2));
            if EELS_PearsonSparse_FullFit.is_pvoigt_model(settings.zlp_model_type)
                shape_bounds = settings.pvoigt_eta_bounds;
            else
                shape_bounds = settings.pearson_m_bounds;
            end
            lb = [settings.eshift_range(1), ...
                -settings.zlp_center_limit .* ones(1, nz), ...
                width_min, ...
                width_min, ...
                shape_bounds(1) .* ones(1, nz), ...
                shape_bounds(1) .* ones(1, nz), ...
                log(settings.T_soft_bounds(1))];
            ub = [settings.eshift_range(2), ...
                settings.zlp_center_limit .* ones(1, nz), ...
                width_max, ...
                width_max, ...
                shape_bounds(2) .* ones(1, nz), ...
                shape_bounds(2) .* ones(1, nz), ...
                log(settings.T_soft_bounds(2))];
        end

        function values = expand_component_setting(raw_values, n, fallback)
            if isempty(raw_values)
                values = fallback .* ones(1, n);
                return;
            end
            raw_values = double(raw_values(:).');
            if isscalar(raw_values)
                values = raw_values .* ones(1, n);
            elseif numel(raw_values) >= n
                values = raw_values(1:n);
            else
                values = [raw_values, raw_values(end) .* ones(1, n - numel(raw_values))];
            end
            values(~isfinite(values)) = fallback;
        end

        function tf = is_pvoigt_model(model_type)
            tf = ismember(lower(string(model_type)), ...
                ["asym_pvoigt", "pvoigt", "asym_pseudo_voigt", "pseudo_voigt"]);
        end

        function residual = stage1_objective(theta, ene, y_norm, EM, settings)
            fit = EELS_PearsonSparse_FullFit.decode_stage1(theta, ene, y_norm, EM, settings);
            residual = [fit.optimizer_residual; ...
                EELS_PearsonSparse_FullFit.stage1_zlp_gain_asym_shape_residual(fit.params, settings)];
        end

        function fit = decode_stage1(theta, ene, y_norm, EM, settings)
            params = EELS_PearsonSparse_FullFit.decode_theta(theta, settings);
            [fit_settings, zlp_region_policy] = EELS_PearsonSparse_FullFit.adaptive_zlp_region_settings(ene, y_norm, settings);
            [A, layout, zlp_shapes, kernel_columns] = EELS_PearsonSparse_FullFit.build_design_matrix(ene, EM, params, fit_settings);
            sqrtW = EELS_PearsonSparse_FullFit.build_weights(ene, y_norm, fit_settings);
            [R, r0] = EELS_PearsonSparse_FullFit.regularization_rows(EM, layout, params.T_soft, size(A, 2), fit_settings, A, ene, y_norm);
            l1 = EELS_PearsonSparse_FullFit.stage1_l1_weights(layout, size(A, 2), fit_settings);
            [coef, inner_info] = EELS_PearsonSparse_FullFit.solve_inner_linear(A, y_norm(:), sqrtW(:), R, r0, ene(:), fit_settings, l1);
            if EELS_PearsonSparse_FullFit.is_fitted_zlp_kernel_mode(fit_settings.kernel_mode) && ...
                    ~isempty(coef) && all(isfinite(coef))
                zlp_seed = zlp_shapes * coef(layout.zlp);
                [A, kernel_columns] = EELS_PearsonSparse_FullFit.replace_signal_kernel_columns( ...
                    A, layout, ene(:), EM(:).', zlp_seed);
                [coef_refit, inner_info_refit] = EELS_PearsonSparse_FullFit.solve_inner_linear( ...
                    A, y_norm(:), sqrtW(:), R, r0, ene(:), fit_settings, l1);
                if ~isempty(coef_refit) && all(isfinite(coef_refit))
                    coef = coef_refit;
                    inner_info.fitted_zlp_kernel_refit = inner_info_refit;
                end
            end
            model = A * coef;
            zlp = zlp_shapes * coef(layout.zlp);
            pm_s = coef(layout.pm_s);
            pm_a = coef(layout.pm_a);
            vib = kernel_columns.stokes * pm_s + kernel_columns.anti * pm_a;
            if isempty(layout.bg)
                bg = zeros(size(y_norm(:)));
            else
                bg = A(:, layout.bg) * coef(layout.bg);
            end
            zlp_component_norm = zlp_shapes .* reshape(coef(layout.zlp), 1, []);
            fit = struct();
            fit.params = params;
            fit.coef = coef;
            fit.model_norm = model;
            fit.zlp_norm = zlp;
            fit.zlp_component_norm = zlp_component_norm;
            fit.vib_norm = vib;
            fit.bg_norm = bg;
            fit.PM_stokes = pm_s(:).';
            fit.PM_antistokes = pm_a(:).';
            fit.optimizer_residual = [(model - y_norm(:)) .* sqrtW(:); R * coef - r0];
            fit.weights = sqrtW(:) .^ 2;
            fit.layout = layout;
            fit.EM = EM;
            fit.T_soft = params.T_soft;
            fit.inner_solver = inner_info;
            fit.zlp_region_policy = zlp_region_policy;
            fit.effective_settings = fit_settings;
            fit.zlp_core_fartail = EELS_PearsonSparse_FullFit.stage1_zlp_core_fartail_summary(ene(:), y_norm(:), zlp(:), fit_settings);
        end

        function fit = select_stage1_candidate(fit_fast, fit_candidate, ene, y_norm, settings)
            fast_metrics = EELS_PearsonSparse_FullFit.metrics(y_norm, fit_fast.model_norm, ene, settings);
            cand_metrics = EELS_PearsonSparse_FullFit.metrics(y_norm, fit_candidate.model_norm, ene, settings);
            tail_fast = EELS_PearsonSparse_FullFit.tail_overshoot_score(fit_fast.model_norm, y_norm, ene, settings);
            tail_cand = EELS_PearsonSparse_FullFit.tail_overshoot_score(fit_candidate.model_norm, y_norm, ene, settings);
            solver_name = "";
            if isfield(fit_candidate, 'inner_solver') && isfield(fit_candidate.inner_solver, 'solver')
                solver_name = string(fit_candidate.inner_solver.solver);
            end
            is_plain = contains(lower(solver_name), "lsqnonneg");
            shape_ok = cand_metrics.weighted_NRMSE <= 1.03 * max(fast_metrics.weighted_NRMSE, eps);
            tail_better = isfinite(tail_fast) && tail_fast > 0 && isfinite(tail_cand) && tail_cand < 0.55 * tail_fast;
            high_better = isfinite(cand_metrics.highfreq_NRMSE) && isfinite(fast_metrics.highfreq_NRMSE) && ...
                cand_metrics.highfreq_NRMSE < 0.95 * fast_metrics.highfreq_NRMSE;
            accept_candidate = is_plain || (shape_ok && (tail_better || high_better));
            if accept_candidate
                fit = fit_candidate;
                fit.inner_solver.accepted = true;
            else
                fit = fit_fast;
                fit.inner_solver.accepted = false;
            end
            fit.inner_solver.selection = struct('candidate_solver', char(solver_name), ...
                'accepted_candidate', accept_candidate, 'fast_weighted_NRMSE', fast_metrics.weighted_NRMSE, ...
                'candidate_weighted_NRMSE', cand_metrics.weighted_NRMSE, ...
                'fast_highfreq_NRMSE', fast_metrics.highfreq_NRMSE, ...
                'candidate_highfreq_NRMSE', cand_metrics.highfreq_NRMSE, ...
                'fast_tail_overshoot', tail_fast, 'candidate_tail_overshoot', tail_cand, ...
                'reason', EELS_PearsonSparse_FullFit.selection_reason(accept_candidate, tail_better, high_better, shape_ok));
        end

        function reason = selection_reason(accepted, tail_better, high_better, shape_ok)
            if accepted && tail_better
                reason = 'accepted_constrained_tail_improved';
            elseif accepted && high_better
                reason = 'accepted_constrained_highfreq_improved';
            elseif accepted && shape_ok
                reason = 'accepted_plain_or_shape_safe';
            elseif ~shape_ok
                reason = 'rejected_shape_worse';
            else
                reason = 'rejected_no_tail_or_highfreq_gain';
            end
        end

        function score = tail_overshoot_score(model, target, ene, settings)
            target = target(:);
            model = model(:);
            ene = ene(:);
            if numel(target) ~= numel(model) || numel(ene) ~= numel(model)
                score = NaN;
                return;
            end
            mask = abs(ene) >= settings.tail_overshoot_start_meV & isfinite(target) & isfinite(model);
            if nnz(mask) < 5
                score = NaN;
                return;
            end
            noise = EELS_PearsonSparse_FullFit.robust_noise(target, mask);
            amp = max(abs(target(mask)), [], 'omitnan');
            margin = max([settings.tail_overshoot_abs_margin, ...
                settings.tail_overshoot_noise_sigma .* max(noise, 0), ...
                settings.tail_overshoot_rel_margin .* max(amp, eps)]);
            over = model(mask) - target(mask) - margin;
            score = max(over, [], 'omitnan');
        end

        function params = decode_theta(theta, settings)
            nz = settings.num_zlp_components;
            pos = 1;
            params.E_shift = theta(pos); pos = pos + 1;
            params.center = theta(pos:(pos + nz - 1)); pos = pos + nz;
            params.sL = theta(pos:(pos + nz - 1)); pos = pos + nz;
            params.sR = theta(pos:(pos + nz - 1)); pos = pos + nz;
            params.mL = theta(pos:(pos + nz - 1)); pos = pos + nz;
            params.mR = theta(pos:(pos + nz - 1)); pos = pos + nz;
            params.T_soft = exp(theta(pos));
        end

        function [A, layout, zlp_shapes, kernels] = build_design_matrix(ene, EM, params, settings)
            e = ene(:) - params.E_shift;
            nz = settings.num_zlp_components;
            n = numel(e);
            zlp_shapes = zeros(n, nz);
            for kk = 1:nz
                zlp_shapes(:, kk) = EELS_PearsonSparse_FullFit.zlp_shape(e - params.center(kk), ...
                    params.sL(kk), params.sR(kk), params.mL(kk), params.mR(kk), settings);
            end
            kernel_base = EELS_PearsonSparse_FullFit.signal_kernel(e, params, settings);
            ns = numel(EM);
            S = EELS_PearsonSparse_FullFit.shift_kernel_matrix(ene(:), kernel_base, EM(:).');
            Aanti = EELS_PearsonSparse_FullFit.shift_kernel_matrix(ene(:), kernel_base, -EM(:).');
            bg = zeros(n, 0);
            A = [zlp_shapes, S, Aanti, bg];
            layout = struct();
            layout.zlp = 1:nz;
            layout.pm_s = (nz + 1):(nz + ns);
            layout.pm_a = (nz + ns + 1):(nz + 2 * ns);
            layout.bg = (nz + 2 * ns + 1):size(A, 2);
            kernels = struct('stokes', S, 'anti', Aanti);
        end

        function k = signal_kernel(e, params, settings)
            switch lower(settings.kernel_mode)
                case {'full_zlp', 'fitted_zlp'}
                    k = zeros(size(e(:)));
                    for kk = 1:numel(params.center)
                        k = k + EELS_PearsonSparse_FullFit.zlp_shape(e(:) - params.center(kk), ...
                            params.sL(kk), params.sR(kk), params.mL(kk), params.mR(kk), settings);
                    end
                otherwise
                    [~, kk] = min(0.5 .* (params.sL + params.sR));
                    k = EELS_PearsonSparse_FullFit.zlp_shape(e(:) - params.center(kk), ...
                        params.sL(kk), params.sR(kk), params.mL(kk), params.mR(kk), settings);
            end
            k(~isfinite(k) | k < 0) = 0;
            s = sum(k);
            if s <= 0, s = 1; end
            k = k ./ s;
        end

        function tf = is_fitted_zlp_kernel_mode(kernel_mode)
            tf = strcmpi(char(string(kernel_mode)), 'fitted_zlp');
        end

        function [A, kernels] = replace_signal_kernel_columns(A, layout, ene, EM, kernel_base)
            kernel_base = EELS_PearsonSparse_FullFit.normalize_signal_kernel(kernel_base);
            S = EELS_PearsonSparse_FullFit.shift_kernel_matrix(ene(:), kernel_base(:), EM(:).');
            Aanti = EELS_PearsonSparse_FullFit.shift_kernel_matrix(ene(:), kernel_base(:), -EM(:).');
            if isfield(layout, 'pm_s') && ~isempty(layout.pm_s)
                A(:, layout.pm_s) = S;
            end
            if isfield(layout, 'pm_a') && ~isempty(layout.pm_a)
                A(:, layout.pm_a) = Aanti;
            end
            kernels = struct('stokes', S, 'anti', Aanti, 'base', kernel_base(:));
        end

        function k = normalize_signal_kernel(k)
            k = k(:);
            k(~isfinite(k) | k < 0) = 0;
            s = sum(k);
            if s <= 0
                k = zeros(size(k));
                if ~isempty(k)
                    k(ceil(numel(k) ./ 2)) = 1;
                end
                return;
            end
            k = k ./ s;
        end

        function y = zlp_shape(x, sL, sR, pL, pR, settings)
            switch lower(settings.zlp_model_type)
                case {'asym_pvoigt', 'pvoigt', 'asym_pseudo_voigt', 'pseudo_voigt'}
                    y = EELS_PearsonSparse_FullFit.asym_pvoigt(x, sL, sR, pL, pR);
                otherwise
                    y = EELS_PearsonSparse_FullFit.asym_pearson(x, sL, sR, pL, pR);
            end
        end

        function y = asym_pearson(x, sL, sR, mL, mR)
            x = x(:);
            s = sR .* ones(size(x));
            m = mR .* ones(size(x));
            left = x < 0;
            s(left) = sL;
            m(left) = mL;
            s = max(s, eps);
            m = max(m, 0.51);
            y = (1 + (x ./ s) .^ 2 ./ m) .^ (-m);
            y(~isfinite(y)) = 0;
        end

        function y = asym_pvoigt(x, sL, sR, etaL, etaR)
            x = x(:);
            s = sR .* ones(size(x));
            eta = etaR .* ones(size(x));
            left = x < 0;
            s(left) = sL;
            eta(left) = etaL;
            s = max(s, eps);
            eta = min(max(eta, 0), 1);
            gaussian = exp(-0.5 .* (x ./ s) .^ 2);
            lorentz = 1 ./ (1 + (x ./ s) .^ 2);
            y = eta .* lorentz + (1 - eta) .* gaussian;
            y(~isfinite(y)) = 0;
        end

        function residual = stage1_zlp_gain_asym_shape_residual(params, settings)
            residual = zeros(0, 1);
            if ~settings.stage1_zlp_gain_asym_shape_prior_enabled || ...
                    settings.stage1_zlp_gain_asym_shape_prior_weight <= 0 || ...
                    ~isfield(params, 'sL') || ~isfield(params, 'sR')
                return;
            end
            nz = numel(params.sL);
            if nz == 0
                return;
            end
            start_idx = max(1, min(nz, round(settings.stage1_zlp_gain_asym_prior_start_index)));
            idx = start_idx:nz;
            width_ratio = max(settings.stage1_zlp_gain_asym_width_ratio, 1);
            width_penalty = max(0, log(width_ratio .* max(params.sR(idx), eps) ./ max(params.sL(idx), eps)));
            shape_ratio = max(settings.stage1_zlp_gain_asym_shape_ratio, 1);
            if EELS_PearsonSparse_FullFit.is_pvoigt_model(settings.zlp_model_type)
                shape_penalty = max(0, log(shape_ratio .* max(params.mR(idx), eps) ./ max(params.mL(idx), eps)));
            else
                shape_penalty = max(0, log(shape_ratio .* max(params.mL(idx), eps) ./ max(params.mR(idx), eps)));
            end
            residual = sqrt(settings.stage1_zlp_gain_asym_shape_prior_weight) .* ...
                [width_penalty(:); shape_penalty(:)];
            residual(~isfinite(residual)) = 0;
        end

        function shifted = shift_kernel(ene, kernel, delta)
            shifted = interp1(ene(:), kernel(:), ene(:) - delta, 'linear', 0);
            s = sum(shifted);
            if s > 0
                shifted = shifted ./ s;
            end
        end

        function M = shift_kernel_matrix(ene, kernel, deltas)
            ene = ene(:);
            kernel = kernel(:);
            deltas = deltas(:).';
            n = numel(ene);
            nd = numel(deltas);
            if nd == 0
                M = zeros(n, 0);
                return;
            end
            if n < 2 || numel(kernel) ~= n
                M = zeros(n, nd);
                for ii = 1:nd
                    M(:, ii) = EELS_PearsonSparse_FullFit.shift_kernel(ene, kernel, deltas(ii));
                end
                return;
            end
            de = diff(ene);
            step = median(de, 'omitnan');
            uniform_tol = max(1e-9, 1e-7 .* max(abs(step), eps));
            if ~isfinite(step) || step == 0 || max(abs(de - step), [], 'omitnan') > uniform_tol
                M = zeros(n, nd);
                for ii = 1:nd
                    M(:, ii) = EELS_PearsonSparse_FullFit.shift_kernel(ene, kernel, deltas(ii));
                end
                return;
            end
            frac_idx = (ene - deltas - ene(1)) ./ step + 1;
            lo = floor(frac_idx);
            alpha = frac_idx - lo;
            valid = lo >= 1 & lo < n;
            lo(~valid) = 1;
            hi = lo + 1;
            M = (1 - alpha) .* reshape(kernel(lo), n, nd) + alpha .* reshape(kernel(hi), n, nd);
            M(~valid) = 0;
            M(~isfinite(M) | M < 0) = 0;
            colsum = sum(M, 1);
            good = colsum > 0;
            if any(good)
                M(:, good) = M(:, good) ./ colsum(good);
            end
        end

        function [R, r0] = regularization_rows(EM, layout, Tsoft, ncoef, settings, A, ene, y_norm)
            if nargin < 6, A = []; end
            if nargin < 7, ene = []; end
            if nargin < 8, y_norm = []; end
            rows = cell(18, 1);
            rhs = cell(18, 1);
            n_rows = 0;
            ns = numel(EM);
            if settings.lambda_be_soft > 0 && ns > 0
                kb = 0.08617333262;
                q = exp(-EM(:) ./ max(kb * Tsoft, eps));
                Rbe = zeros(ns, ncoef);
                for ii = 1:ns
                    Rbe(ii, layout.pm_a(ii)) = 1;
                    Rbe(ii, layout.pm_s(ii)) = -q(ii);
                end
                n_rows = n_rows + 1;
                rows{n_rows} = sqrt(settings.lambda_be_soft) .* Rbe;
                rhs{n_rows} = zeros(ns, 1);
            end
            if settings.lambda_pm_smooth > 0 && ns >= 3
                D2 = diff(eye(ns), 2, 1);
                Rs = zeros(size(D2, 1), ncoef);
                Ra = Rs;
                Rs(:, layout.pm_s) = D2;
                Ra(:, layout.pm_a) = D2;
                n_rows = n_rows + 1;
                rows{n_rows} = sqrt(settings.lambda_pm_smooth) .* Rs;
                rhs{n_rows} = zeros(size(D2, 1), 1);
                n_rows = n_rows + 1;
                rows{n_rows} = sqrt(settings.lambda_pm_smooth) .* Ra;
                rhs{n_rows} = zeros(size(D2, 1), 1);
            end
            if settings.lambda_pm_l2 > 0 && ns > 0
                Rpm = zeros(2 * ns, ncoef);
                Rpm(1:ns, layout.pm_s) = eye(ns);
                Rpm(ns + (1:ns), layout.pm_a) = eye(ns);
                n_rows = n_rows + 1;
                rows{n_rows} = sqrt(settings.lambda_pm_l2) .* Rpm;
                rhs{n_rows} = zeros(2 * ns, 1);
            end
            if settings.lambda_bg_l2 > 0 && isfield(layout, 'bg') && ~isempty(layout.bg)
                Rbg = zeros(numel(layout.bg), ncoef);
                Rbg(:, layout.bg) = eye(numel(layout.bg));
                n_rows = n_rows + 1;
                rows{n_rows} = sqrt(settings.lambda_bg_l2) .* Rbg;
                rhs{n_rows} = zeros(numel(layout.bg), 1);
            end
            if settings.zlp_wide_component_l2 > 0 && isfield(layout, 'zlp') && ~isempty(layout.zlp)
                zidx = layout.zlp(:).';
                start_idx = max(1, min(numel(zidx), round(settings.zlp_wide_component_start_index)));
                penalized = zidx(start_idx:end);
                if ~isempty(penalized)
                    Rwide = zeros(numel(penalized), ncoef);
                    Rwide(:, penalized) = eye(numel(penalized));
                    n_rows = n_rows + 1;
                    rows{n_rows} = sqrt(settings.zlp_wide_component_l2) .* Rwide;
                    rhs{n_rows} = zeros(numel(penalized), 1);
                end
            end
            if settings.zlp_signal_guard_enabled && settings.zlp_signal_guard_weight > 0 && ...
                    ~isempty(A) && ~isempty(ene) && isfield(layout, 'zlp') && ~isempty(layout.zlp)
                eabs = abs(ene(:));
                guard_start = settings.zlp_signal_guard_start_meV;
                guard_end = settings.zlp_signal_guard_end_meV;
                guard = eabs >= guard_start & eabs <= guard_end;
                if guard_end > guard_start && nnz(guard) >= 3
                    taper_width = max(settings.zlp_signal_guard_taper_meV, eps);
                    low_taper = 1 ./ (1 + exp(-(eabs(guard) - guard_start) ./ taper_width));
                    high_taper = 1 ./ (1 + exp((eabs(guard) - guard_end) ./ taper_width));
                    release_width = max(settings.zlp_signal_guard_tail_release_width_meV, eps);
                    release_start = settings.zlp_signal_guard_tail_release_start_meV;
                    release_floor = min(max(settings.zlp_signal_guard_tail_floor, 0), 1);
                    tail_release = release_floor + (1 - release_floor) ./ ...
                        (1 + exp((eabs(guard) - release_start) ./ release_width));
                    guard_profile = EELS_PearsonSparse_FullFit.zlp_signal_guard_profile(ene, y_norm, settings);
                    taper = low_taper .* high_taper .* tail_release .* guard_profile(guard);
                    Rz = zeros(nnz(guard), ncoef);
                    Rz(:, layout.zlp) = A(guard, layout.zlp) .* taper(:);
                    n_rows = n_rows + 1;
                    rows{n_rows} = sqrt(settings.zlp_signal_guard_weight) .* Rz;
                    rhs{n_rows} = zeros(nnz(guard), 1);
                end
            end
            if settings.zlp_lowfreq_guard_enabled && settings.zlp_lowfreq_guard_weight > 0 && ...
                    ~isempty(A) && ~isempty(ene) && isfield(layout, 'zlp') && ~isempty(layout.zlp)
                eabs = abs(ene(:));
                low_start = settings.zlp_lowfreq_guard_start_meV;
                low_end = settings.zlp_lowfreq_guard_end_meV;
                low_guard = eabs >= low_start & eabs <= low_end;
                if low_end > low_start && nnz(low_guard) >= 3
                    taper_width = max(settings.zlp_lowfreq_guard_taper_meV, eps);
                    low_taper = 1 ./ (1 + exp(-(eabs(low_guard) - low_start) ./ taper_width));
                    high_taper = 1 ./ (1 + exp((eabs(low_guard) - low_end) ./ taper_width));
                    guard_profile = EELS_PearsonSparse_FullFit.zlp_signal_guard_profile(ene, y_norm, settings);
                    taper = low_taper .* high_taper .* (0.35 + 0.65 .* guard_profile(low_guard));
                    Rlow = zeros(nnz(low_guard), ncoef);
                    Rlow(:, layout.zlp) = A(low_guard, layout.zlp) .* taper(:);
                    n_rows = n_rows + 1;
                    rows{n_rows} = sqrt(settings.zlp_lowfreq_guard_weight) .* Rlow;
                    rhs{n_rows} = zeros(nnz(low_guard), 1);
                end
            end
            if settings.stage1_zlp_core_fartail_enabled && ...
                    ~isempty(A) && ~isempty(ene) && ~isempty(y_norm) && ...
                    isfield(layout, 'zlp') && ~isempty(layout.zlp)
                ene_col = ene(:);
                y_col = y_norm(:);
                eabs = abs(ene_col);
                good = isfinite(ene_col) & isfinite(y_col);
                if settings.stage1_zlp_core_anchor_weight > 0
                    core_width = max(settings.stage1_zlp_core_anchor_width_meV, eps);
                    taper_width = max(0.25 .* core_width, 1);
                    core_mask = good & eabs <= core_width + 3 .* taper_width;
                    if nnz(core_mask) >= 3
                        core_taper = 1 ./ (1 + exp((eabs(core_mask) - core_width) ./ taper_width));
                        Rcore = zeros(nnz(core_mask), ncoef);
                        Rcore(:, layout.zlp) = A(core_mask, layout.zlp) .* core_taper(:);
                        n_rows = n_rows + 1;
                        rows{n_rows} = sqrt(settings.stage1_zlp_core_anchor_weight) .* Rcore;
                        rhs{n_rows} = sqrt(settings.stage1_zlp_core_anchor_weight) .* y_col(core_mask) .* core_taper(:);
                    end
                end
                if settings.stage1_zlp_fartail_anchor_weight > 0
                    [tail_profile, tail_mask] = EELS_PearsonSparse_FullFit.stage1_zlp_fartail_profile(ene_col, settings);
                    tail_mask = tail_mask & good;
                    anchor_side = lower(string(settings.stage1_zlp_fartail_anchor_side));
                    if anchor_side == "gain" || anchor_side == "left" || anchor_side == "negative"
                        tail_mask = tail_mask & ene_col < 0;
                    elseif anchor_side == "loss" || anchor_side == "right" || anchor_side == "positive"
                        tail_mask = tail_mask & ene_col > 0;
                    end
                    if nnz(tail_mask) >= 4
                        profile = tail_profile(tail_mask);
                        Rtail = zeros(nnz(tail_mask), ncoef);
                        Rtail(:, layout.zlp) = A(tail_mask, layout.zlp) .* profile(:);
                        n_rows = n_rows + 1;
                        rows{n_rows} = sqrt(settings.stage1_zlp_fartail_anchor_weight) .* Rtail;
                        rhs{n_rows} = sqrt(settings.stage1_zlp_fartail_anchor_weight) .* y_col(tail_mask) .* profile(:);
                    end
                end
                if settings.stage1_zlp_gain_tail_prior_enabled && settings.stage1_zlp_gain_tail_prior_weight > 0
                    Rpair = EELS_PearsonSparse_FullFit.stage1_zlp_gain_tail_prior_rows(A, layout, ene_col, ncoef, settings);
                    if ~isempty(Rpair)
                        n_rows = n_rows + 1;
                        rows{n_rows} = sqrt(settings.stage1_zlp_gain_tail_prior_weight) .* Rpair;
                        rhs{n_rows} = zeros(size(Rpair, 1), 1);
                    end
                end
            end
            rows = rows(1:n_rows);
            rhs = rhs(1:n_rows);
            if isempty(rows)
                R = zeros(0, ncoef);
                r0 = zeros(0, 1);
            else
                R = vertcat(rows{:});
                r0 = vertcat(rhs{:});
            end
        end

        function l1 = stage1_l1_weights(layout, ncoef, settings)
            l1 = zeros(ncoef, 1);
            if settings.lambda_pm_l1 > 0
                l1(layout.pm_s) = settings.lambda_pm_l1;
                l1(layout.pm_a) = settings.lambda_pm_l1;
            end
        end

        function info = stage2_zlp_refine_basis(stage1_result, zlp_norm, fit_mask, transition_mask, settings)
            n = numel(zlp_norm);
            info = struct('enabled', settings.stage2_zlp_refine_enabled, ...
                'B', zeros(n, 0), 'mode_names', {cell(1, 0)}, 'lb', zeros(0, 1), ...
                'ub', zeros(0, 1), 'l2', zeros(0, 1), 'energy_axis', [], ...
                'fit_points', nnz(fit_mask), 'transition_points', nnz(transition_mask));
            if ~settings.stage2_zlp_refine_enabled || ...
                    ~isfield(stage1_result, 'ene_corrected') || numel(stage1_result.ene_corrected) ~= n
                info.enabled = false;
                return;
            end
            ene = stage1_result.ene_corrected(:);
            z = zlp_norm(:);
            valid_mask = (fit_mask(:) | transition_mask(:)) & isfinite(z) & isfinite(ene);
            if nnz(valid_mask) < 8
                info.enabled = false;
                return;
            end
            B = {};
            names = {};
            lb = [];
            ub = [];
            l2 = [];
            if settings.stage2_zlp_refine_scale_enabled
                B{end + 1} = z; %#ok<AGROW>
                names{end + 1} = 'zlp_scale'; %#ok<AGROW>
                b = abs(settings.stage2_zlp_refine_scale_bound);
                lb(end + 1, 1) = -b; %#ok<AGROW>
                ub(end + 1, 1) = b; %#ok<AGROW>
                l2(end + 1, 1) = max(settings.stage2_zlp_refine_scale_l2, 0); %#ok<AGROW>
            end
            if settings.stage2_zlp_refine_tilt_enabled
                tilt = z .* tanh(ene ./ max(settings.stage2_zlp_refine_tilt_scale_meV, eps));
                B{end + 1} = tilt; %#ok<AGROW>
                names{end + 1} = 'zlp_left_right_tilt'; %#ok<AGROW>
                b = abs(settings.stage2_zlp_refine_tilt_bound);
                lb(end + 1, 1) = -b; %#ok<AGROW>
                ub(end + 1, 1) = b; %#ok<AGROW>
                l2(end + 1, 1) = max(settings.stage2_zlp_refine_tilt_l2, 0); %#ok<AGROW>
            end
            if settings.stage2_zlp_refine_shift_enabled
                dE = gradient(ene);
                dzdE = gradient(z) ./ max(abs(dE), eps);
                shift_col = -dzdE;
                B{end + 1} = shift_col; %#ok<AGROW>
                names{end + 1} = 'zlp_shift_meV'; %#ok<AGROW>
                b = abs(settings.stage2_zlp_refine_shift_bound_meV);
                lb(end + 1, 1) = -b; %#ok<AGROW>
                ub(end + 1, 1) = b; %#ok<AGROW>
                l2(end + 1, 1) = max(settings.stage2_zlp_refine_shift_l2, 0); %#ok<AGROW>
            end
            dE = gradient(ene);
            dzdE = gradient(z) ./ max(abs(dE), eps);
            if settings.stage2_zlp_refine_width_enabled
                width_col = -ene .* dzdE;
                B{end + 1} = width_col; %#ok<AGROW>
                names{end + 1} = 'zlp_width_stretch'; %#ok<AGROW>
                b = abs(settings.stage2_zlp_refine_width_bound);
                lb(end + 1, 1) = -b; %#ok<AGROW>
                ub(end + 1, 1) = b; %#ok<AGROW>
                l2(end + 1, 1) = max(settings.stage2_zlp_refine_width_l2, 0); %#ok<AGROW>
            end
            tail_scale = max(settings.stage2_zlp_refine_tail_scale_meV, eps);
            tail_gain = log1p((abs(ene) ./ tail_scale) .^ 2);
            valid_tail_gain = tail_gain(valid_mask & isfinite(tail_gain));
            if ~isempty(valid_tail_gain)
                tail_gain = tail_gain ./ max(max(valid_tail_gain), eps);
            end
            if settings.stage2_zlp_refine_tail_enabled
                tail_col = z .* tail_gain;
                B{end + 1} = tail_col; %#ok<AGROW>
                names{end + 1} = 'zlp_tail_symmetric'; %#ok<AGROW>
                b = abs(settings.stage2_zlp_refine_tail_bound);
                lb(end + 1, 1) = -b; %#ok<AGROW>
                ub(end + 1, 1) = b; %#ok<AGROW>
                l2(end + 1, 1) = max(settings.stage2_zlp_refine_tail_l2, 0); %#ok<AGROW>
            end
            if settings.stage2_zlp_refine_tail_asym_enabled
                tail_asym_col = z .* tail_gain .* tanh(ene ./ max(settings.stage2_zlp_refine_tilt_scale_meV, eps));
                B{end + 1} = tail_asym_col; %#ok<AGROW>
                names{end + 1} = 'zlp_tail_asymmetry'; %#ok<AGROW>
                b = abs(settings.stage2_zlp_refine_tail_asym_bound);
                lb(end + 1, 1) = -b; %#ok<AGROW>
                ub(end + 1, 1) = b; %#ok<AGROW>
                l2(end + 1, 1) = max(settings.stage2_zlp_refine_tail_asym_l2, 0); %#ok<AGROW>
            end
            if isempty(B)
                info.enabled = false;
                return;
            end
            M = [B{:}];
            M(~isfinite(M)) = 0;
            keep = any(abs(M(valid_mask, :)) > 0, 1);
            M = M(:, keep);
            names = names(keep);
            lb = lb(keep);
            ub = ub(keep);
            l2 = l2(keep);
            info.enabled = ~isempty(names);
            info.B = M;
            info.mode_names = names;
            info.lb = lb(:);
            info.ub = ub(:);
            info.l2 = l2(:);
            info.energy_axis = ene;
        end

        function [Reg, rhs] = stage2_zlp_refine_regularization(zlp_refine, base_ncoef)
            n = size(zlp_refine.B, 2);
            if n <= 0
                Reg = zeros(0, base_ncoef);
                rhs = zeros(0, 1);
                return;
            end
            Reg = zeros(n, base_ncoef + n);
            Reg(:, base_ncoef + (1:n)) = diag(sqrt(max(zlp_refine.l2(:), 0)));
            rhs = zeros(n, 1);
        end

        function bounds = stage2_append_zlp_refine_bounds(bounds, zlp_refine)
            n = size(zlp_refine.B, 2);
            if n <= 0
                return;
            end
            bounds.lb = [bounds.lb(:); zlp_refine.lb(:)];
            bounds.ub = [bounds.ub(:); zlp_refine.ub(:)];
            bounds.enabled = isfield(bounds, 'enabled') && bounds.enabled;
            bounds.zlp_refine_enabled = true;
            bounds.zlp_refine_modes = zlp_refine.mode_names;
        end

        function out = stage2_zlp_refine_result(zlp_refine, coef, delta_norm, fit_mask)
            coef = coef(:);
            out = struct('enabled', zlp_refine.enabled && ~isempty(coef), ...
                'mode_names', {zlp_refine.mode_names}, 'coef', coef(:).', ...
                'max_abs_delta_norm', NaN, 'rms_delta_norm', NaN, ...
                'fit_rms_delta_norm', NaN, 'fit_points', nnz(fit_mask));
            if isempty(coef) || isempty(delta_norm)
                out.enabled = false;
                return;
            end
            d = delta_norm(:);
            out.max_abs_delta_norm = max(abs(d), [], 'omitnan');
            out.rms_delta_norm = sqrt(mean(d .^ 2, 'omitnan'));
            m = fit_mask(:) & isfinite(d);
            if any(m)
                out.fit_rms_delta_norm = sqrt(mean(d(m) .^ 2, 'omitnan'));
            end
        end

        function [coef, info] = solve_inner_linear(A, y, sqrtW, R, r0, ene, settings, l1, bounds)
            ncoef = size(A, 2);
            if nargin < 8 || isempty(l1), l1 = zeros(ncoef, 1); end
            if nargin < 9 || isempty(bounds), bounds = struct(); end
            l1 = l1(:);
            if numel(l1) ~= ncoef, l1 = zeros(ncoef, 1); end
            y = y(:);
            sqrtW = sqrtW(:);
            C = [A .* sqrtW; R];
            d = [y .* sqrtW; r0(:)];
            lb = zeros(ncoef, 1);
            ub = Inf(ncoef, 1);
            if isfield(bounds, 'lb') && numel(bounds.lb) == ncoef
                lb = bounds.lb(:);
                lb(~isfinite(lb)) = -Inf;
            end
            if isfield(bounds, 'ub') && numel(bounds.ub) == ncoef
                ub = bounds.ub(:);
                ub(~isfinite(ub)) = Inf;
                ub = max(ub, lb + eps);
            end
            [Aineq, bineq, tail_mask] = EELS_PearsonSparse_FullFit.tail_overshoot_inequality(A, y, ene, settings);
            methods = EELS_PearsonSparse_FullFit.inner_solver_order(settings.inner_solver);
            info = struct('solver', '', 'exitflag', NaN, 'message', '', ...
                'num_tail_constraints', nnz(tail_mask), 'max_tail_overshoot', NaN, ...
                'finite_bounds_enabled', isfield(bounds, 'enabled') && bounds.enabled, ...
                'finite_bounds_lower_count', nnz(lb > 0), ...
                'finite_bounds_upper_count', nnz(isfinite(ub)));
            coef = [];
            has_finite_bounds = any(lb > 0) || any(isfinite(ub));
            for mm = 1:numel(methods)
                method = methods{mm};
                try
                    switch method
                        case 'fast_nnls'
                            if ~isempty(Aineq) || has_finite_bounds
                                continue;
                            end
                            tmp = EELS_PearsonSparse_FullFit.fast_projected_nnls(C, d, settings.fast_nnls_refine_passes);
                            exitflag = 1;
                            output = struct('message', 'fast projected NNLS approximation');
                        case 'quadprog'
                            if exist('quadprog', 'file') ~= 2, continue; end
                            H = C' * C;
                            H = 0.5 .* (H + H') + 1e-12 .* speye(ncoef);
                            f = -C' * d + l1;
                            opts = optimoptions('quadprog', 'Display', 'off', ...
                                'MaxIterations', settings.inner_solver_max_iter, ...
                                'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-8);
                            [tmp, ~, exitflag, output] = quadprog(H, f, Aineq, bineq, [], [], lb, ub, [], opts);
                        case 'lsqlin'
                            if exist('lsqlin', 'file') ~= 2, continue; end
                            opts = optimoptions('lsqlin', 'Display', 'off', ...
                                'MaxIterations', settings.inner_solver_max_iter, ...
                                'OptimalityTolerance', 1e-8, 'ConstraintTolerance', 1e-8);
                            [tmp, ~, ~, exitflag, output] = lsqlin(C, d, Aineq, bineq, [], [], lb, ub, [], opts);
                        otherwise
                            tmp = lsqnonneg(C, d);
                            tmp = min(max(tmp, lb), ub);
                            exitflag = 1;
                            output = struct('message', 'lsqnonneg fallback');
                    end
                    if ~isempty(tmp) && all(isfinite(tmp)) && exitflag > 0
                        coef = tmp(:);
                        info.solver = method;
                        info.exitflag = exitflag;
                        if isstruct(output) && isfield(output, 'message')
                            info.message = output.message;
                        end
                        break;
                    end
                catch ME
                    info.message = ME.message;
                end
            end
            if isempty(coef)
                coef = lsqnonneg(C, d);
                coef = min(max(coef, lb), ub);
                info.solver = 'lsqnonneg_final_fallback';
                info.exitflag = 1;
            end
            if ~isempty(Aineq)
                info.max_tail_overshoot = max(Aineq * coef - bineq, [], 'omitnan');
            else
                info.max_tail_overshoot = NaN;
            end
        end

        function coef = fast_projected_nnls(C, d, refine_passes)
            if nargin < 3 || isempty(refine_passes), refine_passes = 2; end
            refine_passes = max(0, round(refine_passes));
            ncoef = size(C, 2);
            if isempty(C) || ncoef == 0
                coef = zeros(ncoef, 1);
                return;
            end
            coef = EELS_PearsonSparse_FullFit.fast_lstsq(C, d);
            coef(~isfinite(coef)) = 0;
            coef(coef < 0) = 0;
            active = coef > max(1e-12, 1e-9 .* max(coef));
            if nnz(active) == 0
                [~, imax] = max(abs(C' * d));
                active(imax) = true;
            end
            for pp = 1:refine_passes
                tmp = zeros(ncoef, 1);
                tmp(active) = EELS_PearsonSparse_FullFit.fast_lstsq(C(:, active), d);
                tmp(~isfinite(tmp)) = 0;
                bad = active & tmp < 0;
                tmp(tmp < 0) = 0;
                coef = tmp;
                if ~any(bad)
                    corr = C' * (d - C * coef);
                    add = corr > max(1e-10, 1e-6 .* max(abs(corr)));
                    add = add & ~active;
                    if ~any(add)
                        break;
                    end
                    [~, order] = sort(corr(add), 'descend');
                    add_idx = find(add);
                    add_idx = add_idx(order(1:min(4, numel(order))));
                    active(add_idx) = true;
                else
                    active(bad) = false;
                    if nnz(active) == 0
                        [~, imax] = max(abs(C' * d));
                        active(imax) = true;
                    end
                end
            end
            coef(~isfinite(coef)) = 0;
            coef(coef < 0) = 0;
        end

        function x = fast_lstsq(C, d)
            if isempty(C) || size(C, 2) == 0
                x = zeros(size(C, 2), 1);
                return;
            end
            if exist('lsqminnorm', 'file') == 2
                x = lsqminnorm(C, d);
            else
                warn_state = warning('off', 'MATLAB:rankDeficientMatrix');
                cleanup = onCleanup(@() warning(warn_state));
                x = C \ d;
            end
            x(~isfinite(x)) = 0;
        end

        function methods = inner_solver_order(requested)
            requested = lower(string(requested));
            switch requested
                case "fast_nnls"
                    methods = {'fast_nnls'};
                case "quadprog"
                    methods = {'quadprog', 'lsqlin', 'lsqnonneg'};
                case "lsqlin"
                    methods = {'lsqlin', 'quadprog', 'lsqnonneg'};
                case "lsqnonneg"
                    methods = {'lsqnonneg'};
                otherwise
                    methods = {'lsqlin', 'lsqnonneg', 'quadprog'};
            end
        end

        function [Aineq, bineq, mask] = tail_overshoot_inequality(A, target, ene, settings)
            n = size(A, 1);
            mask = false(n, 1);
            Aineq = [];
            bineq = [];
            if ~settings.tail_overshoot_constraint || numel(ene) ~= n || numel(target) ~= n
                return;
            end
            ene = ene(:);
            target = target(:);
            eabs = abs(ene);
            mask = eabs >= settings.tail_overshoot_start_meV & isfinite(target);
            if nnz(mask) < 5
                return;
            end
            noise = EELS_PearsonSparse_FullFit.robust_noise(target, mask);
            amp = max(abs(target(mask)), [], 'omitnan');
            if ~isfinite(amp) || amp <= 0, amp = max(abs(target), [], 'omitnan'); end
            if ~isfinite(noise), noise = 0; end
            margin = max([settings.tail_overshoot_abs_margin, ...
                settings.tail_overshoot_noise_sigma .* noise, ...
                settings.tail_overshoot_rel_margin .* max(amp, eps)]);
            Aineq = A(mask, :);
            bineq = max(target(mask) + margin, margin);
        end

        function [profile, mask] = stage1_zlp_fartail_profile(ene, settings)
            ene = ene(:);
            eabs = abs(ene);
            start_e = max(settings.stage1_zlp_fartail_start_meV, 0);
            end_e = settings.stage1_zlp_fartail_end_meV;
            if isempty(end_e) || ~isfinite(end_e)
                end_e = Inf;
            end
            taper_width = max(settings.stage1_zlp_fartail_taper_meV, eps);
            low_taper = 1 ./ (1 + exp(-(eabs - start_e) ./ taper_width));
            if isfinite(end_e)
                high_taper = 1 ./ (1 + exp((eabs - end_e) ./ taper_width));
            else
                high_taper = ones(size(eabs));
            end
            side_boost = ones(size(eabs));
            gain_boost = max(settings.stage1_zlp_fartail_gain_boost, 1);
            side_boost(ene < 0) = gain_boost;
            profile = low_taper .* high_taper .* side_boost;
            mask = isfinite(ene) & eabs >= start_e & eabs <= end_e;
            profile(~mask | ~isfinite(profile)) = 0;
        end

        function Rpair = stage1_zlp_gain_tail_prior_rows(A, layout, ene, ncoef, settings)
            Rpair = [];
            if isempty(A) || isempty(ene) || ~isfield(layout, 'zlp') || isempty(layout.zlp)
                return;
            end
            ene = ene(:);
            start_e = max(settings.stage1_zlp_fartail_start_meV, 0);
            end_e = settings.stage1_zlp_fartail_end_meV;
            if isempty(end_e) || ~isfinite(end_e)
                end_e = max(abs(ene(isfinite(ene))), [], 'omitnan');
            end
            pos = ene > 0 & ene >= start_e & ene <= end_e & isfinite(ene);
            neg = ene < 0 & -ene >= start_e & -ene <= end_e & isfinite(ene);
            if nnz(pos) < 3 || nnz(neg) < 3
                return;
            end
            max_pair_e = min([max(ene(pos), [], 'omitnan'), max(-ene(neg), [], 'omitnan'), end_e]);
            if ~isfinite(max_pair_e) || max_pair_e <= start_e
                return;
            end
            step = max(settings.stage1_zlp_gain_tail_pair_step_meV, eps);
            pair_e = start_e:step:max_pair_e;
            if numel(pair_e) < 2
                pair_e = linspace(start_e, max_pair_e, 3);
            end
            ratio = max(settings.stage1_zlp_gain_tail_prior_ratio, 1);
            grid_step = median(abs(diff(sort(ene))), 'omitnan');
            if ~isfinite(grid_step) || grid_step <= 0
                grid_step = step;
            end
            tol = max(step, grid_step .* 2);
            Rtmp = zeros(numel(pair_e), ncoef);
            n = 0;
            for ii = 1:numel(pair_e)
                e0 = pair_e(ii);
                [dpos, ip] = min(abs(ene - e0));
                [dneg, ineg] = min(abs(ene + e0));
                if isfinite(dpos) && isfinite(dneg) && dpos <= tol && dneg <= tol
                    n = n + 1;
                    Rtmp(n, layout.zlp) = A(ineg, layout.zlp) - ratio .* A(ip, layout.zlp);
                end
            end
            if n > 0
                Rpair = Rtmp(1:n, :);
            end
        end

        function sqrtW = build_weights(ene, y_norm, settings)
            eabs = abs(ene(:));
            tw = 4;
            core = 1 ./ (1 + exp((eabs - settings.core_downweight_width) ./ tw));
            shoulder = 1 ./ (1 + exp(-(eabs - settings.core_downweight_width) ./ tw)) .* ...
                1 ./ (1 + exp((eabs - settings.highfreq_start_meV) ./ 12));
            high = 1 ./ (1 + exp(-(eabs - settings.highfreq_start_meV) ./ 12));
            tail_start = max(settings.highfreq_start_meV, settings.tail_start_fraction .* max(eabs));
            tail = 1 ./ (1 + exp(-(eabs - tail_start) ./ 12));
            w = settings.core_weight .* core + settings.shoulder_weight .* shoulder + ...
                settings.signal_weight .* (1 - core) + settings.highfreq_weight .* high + settings.tail_weight .* tail;
            if settings.stage1_zlp_core_fartail_enabled
                core_width = max(settings.stage1_zlp_core_anchor_width_meV, eps);
                core_taper = 1 ./ (1 + exp((eabs - core_width) ./ max(0.25 .* core_width, 1)));
                [tail_profile, tail_mask] = EELS_PearsonSparse_FullFit.stage1_zlp_fartail_profile(ene(:), settings);
                w = w + settings.stage1_zlp_core_data_weight .* core_taper + ...
                    settings.stage1_zlp_fartail_data_weight .* tail_profile .* double(tail_mask);
            end
            noise = EELS_PearsonSparse_FullFit.robust_noise(y_norm, eabs > tail_start);
            local = movmedian(max(y_norm(:), 0), 15, 'omitnan');
            ref = median(local(eabs > settings.core_downweight_width), 'omitnan');
            if ~isfinite(ref) || ref <= 0, ref = max(local); end
            floor_signal = max(settings.dynamic_weight_noise_sigma .* noise, 0.01 .* ref);
            balance = (ref ./ max(local, floor_signal)) .^ 0.45;
            balance = min(max(balance, 1), settings.dynamic_weight_max_gain);
            w = w .* (1 + settings.dynamic_weight_strength .* (1 - core) .* (balance - 1));
            w(~isfinite(w) | w <= 0) = 1;
            sqrtW = sqrt(w);
        end

        function [settings, policy] = adaptive_zlp_region_settings(ene, y_norm, settings)
            policy = struct();
            policy.enabled = settings.adaptive_zlp_region_enabled;
            policy.method = 'static';
            policy.peak_center_meV = NaN;
            policy.core_half_width_meV = NaN;
            policy.guard_start_meV = settings.zlp_signal_guard_start_meV;
            policy.guard_end_meV = settings.zlp_signal_guard_end_meV;
            policy.core_downweight_width_meV = settings.core_downweight_width;
            if ~settings.adaptive_zlp_region_enabled
                return;
            end

            ene = ene(:);
            y = y_norm(:);
            good = isfinite(ene) & isfinite(y);
            search = good & abs(ene) <= settings.adaptive_zlp_region_search_meV;
            if nnz(search) < 7
                return;
            end

            y0 = max(y, 0);
            y_smooth = movmean(movmedian(y0, 5, 'omitnan'), 5, 'omitnan');
            search_idx = find(search);
            [peak, rel_idx] = max(y_smooth(search_idx));
            if ~isfinite(peak) || peak <= 0
                return;
            end
            peak_idx = search_idx(rel_idx);
            peak_center = ene(peak_idx);
            centered_e = ene - peak_center;
            core_threshold = settings.adaptive_zlp_core_fraction .* peak;
            core_mask = search & y_smooth >= core_threshold;
            if nnz(core_mask) < 3
                return;
            end

            core_half_width = max(abs(centered_e(core_mask)), [], 'omitnan');
            if ~isfinite(core_half_width) || core_half_width <= 0
                return;
            end
            guard_start = core_half_width + settings.adaptive_zlp_guard_margin_meV;
            guard_start = min(max(guard_start, settings.adaptive_zlp_guard_min_start_meV), ...
                settings.adaptive_zlp_guard_max_start_meV);
            max_abs_e = max(abs(ene(good)), [], 'omitnan');
            guard_end = min(settings.zlp_signal_guard_end_meV, ...
                settings.adaptive_zlp_guard_end_fraction .* max_abs_e);
            if ~isfinite(guard_end) || guard_end <= guard_start + 5
                guard_end = settings.zlp_signal_guard_end_meV;
            end
            if guard_end <= guard_start + 5
                return;
            end

            core_width = core_half_width + 0.5 .* settings.adaptive_zlp_guard_margin_meV;
            core_width = min(max(core_width, settings.adaptive_zlp_core_min_meV), ...
                settings.adaptive_zlp_core_max_meV);

            settings.zlp_signal_guard_start_meV = guard_start;
            settings.zlp_signal_guard_end_meV = guard_end;
            if settings.adaptive_zlp_region_adjust_weights
                settings.core_downweight_width = core_width;
            end
            policy.method = 'adaptive_peak_fraction';
            policy.peak_center_meV = peak_center;
            policy.core_half_width_meV = core_half_width;
            policy.guard_start_meV = guard_start;
            policy.guard_end_meV = guard_end;
            policy.core_downweight_width_meV = settings.core_downweight_width;
            policy.peak_fraction = settings.adaptive_zlp_core_fraction;
        end

        function [theta_opt, resnorm, exitflag, output] = run_lsq(obj, theta0, lb, ub, settings)
            if exist('lsqnonlin', 'file') == 2
                display_mode = 'off';
                if ~settings.silent, display_mode = 'iter'; end
                opts = optimoptions('lsqnonlin', 'Display', display_mode, ...
                    'MaxIterations', settings.max_nonlinear_iter, ...
                    'FunctionTolerance', 1e-6, 'StepTolerance', 1e-6, ...
                    'FiniteDifferenceStepSize', 1e-4);
                [theta_opt, resnorm, ~, exitflag, output] = lsqnonlin(obj, theta0, lb, ub, opts);
            else
                penalty = @(x) sum(obj(min(max(x, lb), ub)).^2);
                opts = optimset('Display', 'off', 'MaxIter', settings.max_nonlinear_iter * 30, 'TolX', 1e-5, 'TolFun', 1e-5);
                [tmp, resnorm, exitflag, output] = fminsearch(penalty, theta0, opts);
                theta_opt = min(max(tmp, lb), ub);
            end
        end

        function result = package_result(fit, ene, y, y_norm, norm_val, settings, history, zlp_align)
            if nargin < 8 || isempty(zlp_align)
                zlp_align = struct('enabled', false, 'status', 'not_recorded', ...
                    'estimated_peak_meV', NaN, 'applied_shift_meV', 0, ...
                    'max_shift_meV', NaN, 'search_meV', NaN, 'method', 'none');
            end
            effective_settings = settings;
            if isfield(fit, 'effective_settings') && isstruct(fit.effective_settings)
                effective_settings = fit.effective_settings;
            end
            result = struct();
            result.analysis_mode = 'pearson_sparse_stage1';
            result.temperature_enabled = false;
            result.ene_corrected = ene(:) - fit.params.E_shift;
            result.ene_raw = ene(:);
            result.ene_stage1_aligned = ene(:);
            result.raw_data = y(:);
            result.raw_data_norm = y_norm(:);
            result.norm_val = norm_val;
            result.fitted_spectrum_norm = fit.model_norm(:);
            result.fitted_spectrum = fit.model_norm(:) .* norm_val;
            result.zlp_only = fit.zlp_norm(:) .* norm_val;
            result.zlp_components = fit.zlp_component_norm .* norm_val;
            result.vib_only = fit.vib_norm(:) .* norm_val;
            result.bg_only = fit.bg_norm(:) .* norm_val;
            result.residual = result.raw_data - result.fitted_spectrum;
            result.residual_norm = y_norm(:) - fit.model_norm(:);
            result.metrics = fit.metrics;
            result.EM_final = fit.EM;
            result.PM_stokes = fit.PM_stokes;
            result.PM_antistokes = fit.PM_antistokes;
            result.PM_final = [flip(fit.PM_antistokes), fit.PM_stokes];
            result.vib_EM_full = [-flip(fit.EM), fit.EM];
            result.T_soft = fit.T_soft;
            result.theta = fit.theta;
            result.zlp_params = fit.params;
            result.stage1_zlp_alignment = zlp_align;
            result.stage1_zero_align_shift_meV = zlp_align.applied_shift_meV;
            result.total_energy_shift_meV = zlp_align.applied_shift_meV + fit.params.E_shift;
            result.weights_used = fit.weights;
            result.iter_history = history;
            if ~isempty(history) && isfield(history, 'elapsed_seconds')
                result.stage1_elapsed_seconds = sum([history.elapsed_seconds], 'omitnan');
            else
                result.stage1_elapsed_seconds = NaN;
            end
            result.stage1 = fit;
            if isfield(fit, 'zlp_core_fartail')
                result.zlp_core_fartail = fit.zlp_core_fartail;
            else
                result.zlp_core_fartail = EELS_PearsonSparse_FullFit.stage1_zlp_core_fartail_summary( ...
                    result.ene_corrected(:), result.raw_data_norm(:), fit.zlp_norm(:), effective_settings);
            end
            if isfield(fit, 'zlp_region_policy')
                result.zlp_region_policy = fit.zlp_region_policy;
            else
                result.zlp_region_policy = struct();
            end
            result.original_settings = settings;
            result.settings = effective_settings;
            result.stage1_pm_ratio_reliability = EELS_PearsonSparse_FullFit.stage1_pm_ratio_reliability(result, effective_settings);
            result.stage2 = EELS_PearsonSparse_FullFit.solve_stage2_temperature_window(result, effective_settings);
            if isfield(result.stage2, 'T_fit') && result.stage2.enabled && strcmp(result.stage2.status, 'ok')
                result.T_fit = result.stage2.T_fit;
                result.T_uncertainty = result.stage2.T_uncertainty;
                result.reliability_score = result.stage2.reliability_score;
                if isfield(result.stage2, 'reliability_source')
                    result.reliability_source = result.stage2.reliability_source;
                else
                    result.reliability_source = EELS_PearsonSparse_FullFit.reliability_source_from_result(result);
                end
                if isfield(result.stage2, 'T_accuracy_K')
                    result.T_accuracy_K = result.stage2.T_accuracy_K;
                    result.T_accuracy_method = result.stage2.T_accuracy_method;
                else
                    result.T_accuracy_K = result.T_uncertainty;
                    result.T_accuracy_method = 'stage2_cost_valley_10pct_width';
                end
                result = EELS_PearsonSparse_FullFit.apply_stage2_primary_model(result);
            else
                result.T_fit = NaN;
                result.T_uncertainty = NaN;
                result.T_accuracy_K = NaN;
                result.T_accuracy_method = 'stage2_cost_valley_10pct_width_unavailable';
                if isfield(result.stage2, 'reliability_score')
                    result.reliability_score = result.stage2.reliability_score;
                    if isfield(result.stage2, 'reliability_source')
                        result.reliability_source = result.stage2.reliability_source;
                    else
                        result.reliability_source = EELS_PearsonSparse_FullFit.reliability_source_from_result(result);
                    end
                else
                    result.reliability_score = result.stage1_pm_ratio_reliability.reliability_score;
                    result.reliability_source = 'stage2_v20_unavailable';
                end
            end
            if isfield(result, 'stage2_primary_model') && result.stage2_primary_model
                result.diagnostic_flags = {'stage2_v20_strict_hardBE_joint_fit'};
                result.dominant_error_source = 'stage2_joint_fit';
            else
                result.diagnostic_flags = {'stage1_shape_only', 'BE_soft_penalty_only'};
                result.dominant_error_source = 'temperature_not_requested';
            end
            result = EELS_PearsonSparse_FullFit.enrich_result_diagnostics(result, fit, effective_settings);
        end

        function m = metrics(y, yfit, ene, settings)
            y = y(:);
            yfit = yfit(:);
            r = y - yfit;
            ss_res = sum(r.^2, 'omitnan');
            ss_tot = sum((y - mean(y, 'omitnan')).^2, 'omitnan');
            if ss_tot <= 0, ss_tot = eps; end
            sqrtW = EELS_PearsonSparse_FullFit.build_weights(ene, y, settings);
            rw = sqrtW(:) .* r;
            yw = sqrtW(:) .* y;
            m = struct();
            m.R2 = 1 - ss_res ./ ss_tot;
            m.mean_R2 = m.R2;
            m.median_R2 = m.R2;
            m.NRMSE = sqrt(mean(r.^2, 'omitnan')) ./ max(max(y) - min(y), eps);
            m.weighted_NRMSE = sqrt(mean(rw.^2, 'omitnan')) ./ max(max(yw) - min(yw), eps);
            m.residual_RMS = sqrt(mean(r.^2, 'omitnan'));
            m.max_abs_residual = max(abs(r), [], 'omitnan');
            eabs = abs(ene(:));
            mask_high = eabs >= settings.highfreq_start_meV;
            if nnz(mask_high) >= 5
                m.highfreq_NRMSE = sqrt(mean(r(mask_high).^2, 'omitnan')) ./ ...
                    max(max(y(mask_high)) - min(y(mask_high)), eps);
            else
                m.highfreq_NRMSE = NaN;
            end
        end

        function result = apply_stage2_primary_model(result)
            if ~(isfield(result, 'stage2') && isstruct(result.stage2) && ...
                    isfield(result.stage2, 'status') && strcmp(result.stage2.status, 'ok') && ...
                    isfield(result.stage2, 'model_norm'))
                return;
            end
            st2 = result.stage2;
            if isfield(result, 'metrics') && ~isfield(result, 'stage1_metrics')
                result.stage1_metrics = result.metrics;
            end
            if isfield(result, 'fitted_spectrum_norm') && ~isfield(result, 'stage1_fitted_spectrum_norm')
                result.stage1_fitted_spectrum_norm = result.fitted_spectrum_norm;
            end
            if isfield(result, 'fitted_spectrum') && ~isfield(result, 'stage1_fitted_spectrum')
                result.stage1_fitted_spectrum = result.fitted_spectrum;
            end
            if isfield(result, 'residual_norm') && ~isfield(result, 'stage1_residual_norm')
                result.stage1_residual_norm = result.residual_norm;
            end
            if isfield(result, 'residual') && ~isfield(result, 'stage1_residual')
                result.stage1_residual = result.residual;
            end
            result.fitted_spectrum_norm = st2.model_norm(:);
            result.residual_norm = result.raw_data_norm(:) - result.fitted_spectrum_norm(:);
            if isfield(st2, 'metrics') && isstruct(st2.metrics)
                result.metrics = st2.metrics;
            else
                metric_settings = struct();
                if isfield(result, 'settings') && isstruct(result.settings)
                    metric_settings = result.settings;
                end
                metric_settings = EELS_PearsonSparse_FullFit.defaults(metric_settings);
                result.metrics = EELS_PearsonSparse_FullFit.metrics( ...
                    result.raw_data_norm(:), result.fitted_spectrum_norm(:), result.ene_corrected(:), metric_settings);
            end
            if isfield(result, 'norm_val') && isfinite(result.norm_val) && result.norm_val > 0
                result.fitted_spectrum = result.fitted_spectrum_norm(:) .* result.norm_val;
                result.residual = result.raw_data(:) - result.fitted_spectrum(:);
                if isfield(st2, 'zlp_norm')
                    result.zlp_only = st2.zlp_norm(:) .* result.norm_val;
                end
                if isfield(st2, 'vib_norm')
                    result.vib_only = st2.vib_norm(:) .* result.norm_val;
                end
                if isfield(st2, 'bg_norm')
                    result.bg_only = st2.bg_norm(:) .* result.norm_val;
                end
            end
            if isfield(st2, 'EM_final')
                result.EM_final = st2.EM_final(:).';
                result.vib_EM_full = [-flip(result.EM_final), result.EM_final];
            end
            if isfield(st2, 'PM_stokes')
                result.PM_stokes = st2.PM_stokes(:).';
            end
            if isfield(st2, 'PM_antistokes')
                result.PM_antistokes = st2.PM_antistokes(:).';
            end
            if isfield(result, 'PM_stokes') && isfield(result, 'PM_antistokes')
                result.PM_final = [flip(result.PM_antistokes), result.PM_stokes];
            end
            result.stage2_primary_model = true;
        end

        function result = enrich_result_diagnostics(result, fit, settings)
            if ~isfield(result, 'zlp_core_fartail') || ~isstruct(result.zlp_core_fartail)
                result.zlp_core_fartail = EELS_PearsonSparse_FullFit.stage1_zlp_core_fartail_summary( ...
                    result.ene_corrected(:), result.raw_data_norm(:), fit.zlp_norm(:), settings);
            end
            result.zlp_signal_guard = EELS_PearsonSparse_FullFit.zlp_signal_guard_summary(result, settings);
            result.attention_summary = EELS_PearsonSparse_FullFit.summarize_attention_regions(result, settings);
            result.loss_breakdown = EELS_PearsonSparse_FullFit.compute_loss_breakdown(result, fit, settings);
            result.confidence_interval = EELS_PearsonSparse_FullFit.energy_confidence_interval(result, settings);
            result.error_budget = EELS_PearsonSparse_FullFit.compute_error_budget(result, settings);
        end

        function summary = stage1_zlp_core_fartail_summary(ene, y_norm, zlp_norm, settings)
            ene = ene(:);
            y = y_norm(:);
            z = zlp_norm(:);
            n = min([numel(ene), numel(y), numel(z)]);
            ene = ene(1:n);
            y = y(1:n);
            z = z(1:n);
            eabs = abs(ene);
            core_width = max(settings.stage1_zlp_core_anchor_width_meV, eps);
            core = isfinite(ene) & isfinite(y) & isfinite(z) & eabs <= core_width;
            [tail_profile, tail] = EELS_PearsonSparse_FullFit.stage1_zlp_fartail_profile(ene, settings);
            tail = tail & isfinite(y) & isfinite(z);
            gain_tail = tail & ene < 0;
            loss_tail = tail & ene > 0;
            summary = struct();
            summary.enabled = settings.stage1_zlp_core_fartail_enabled;
            summary.core_width_meV = core_width;
            summary.fartail_range_meV = [settings.stage1_zlp_fartail_start_meV, settings.stage1_zlp_fartail_end_meV];
            summary.fartail_gain_boost = settings.stage1_zlp_fartail_gain_boost;
            summary.gain_tail_prior_enabled = settings.stage1_zlp_gain_tail_prior_enabled;
            summary.gain_tail_prior_ratio = settings.stage1_zlp_gain_tail_prior_ratio;
            summary.core_points = nnz(core);
            summary.fartail_points = nnz(tail);
            summary.gain_tail_points = nnz(gain_tail);
            summary.loss_tail_points = nnz(loss_tail);
            summary.core_zlp_to_data_area = sum(max(z(core), 0), 'omitnan') ./ max(sum(max(y(core), 0), 'omitnan'), eps);
            summary.fartail_zlp_to_data_area = sum(max(z(tail), 0), 'omitnan') ./ max(sum(max(y(tail), 0), 'omitnan'), eps);
            summary.gain_tail_zlp_mean = EELS_PearsonSparse_FullFit.safe_mean(max(z(gain_tail), 0));
            summary.loss_tail_zlp_mean = EELS_PearsonSparse_FullFit.safe_mean(max(z(loss_tail), 0));
            summary.gain_tail_data_mean = EELS_PearsonSparse_FullFit.safe_mean(max(y(gain_tail), 0));
            summary.loss_tail_data_mean = EELS_PearsonSparse_FullFit.safe_mean(max(y(loss_tail), 0));
            summary.gain_to_loss_zlp_ratio = summary.gain_tail_zlp_mean ./ max(summary.loss_tail_zlp_mean, eps);
            summary.gain_to_loss_data_ratio = summary.gain_tail_data_mean ./ max(summary.loss_tail_data_mean, eps);
            summary.mean_fartail_profile = EELS_PearsonSparse_FullFit.safe_mean(tail_profile(tail));
        end

        function summary = zlp_signal_guard_summary(result, settings)
            ene = result.ene_corrected(:);
            eabs = abs(ene);
            guard = settings.zlp_signal_guard_enabled & ...
                eabs >= settings.zlp_signal_guard_start_meV & ...
                eabs <= settings.zlp_signal_guard_end_meV;
            zlp = max(result.zlp_only(:), 0);
            total = max(result.fitted_spectrum(:), 0);
            zlp_area = sum(zlp(guard), 'omitnan');
            total_area = sum(total(guard), 'omitnan');
            summary = struct();
            summary.enabled = settings.zlp_signal_guard_enabled;
            summary.energy_range_meV = [settings.zlp_signal_guard_start_meV, settings.zlp_signal_guard_end_meV];
            summary.weight = settings.zlp_signal_guard_weight;
            summary.tail_release_start_meV = settings.zlp_signal_guard_tail_release_start_meV;
            summary.tail_release_floor = settings.zlp_signal_guard_tail_floor;
            summary.lowfreq_guard_enabled = settings.zlp_lowfreq_guard_enabled;
            summary.lowfreq_guard_range_meV = [settings.zlp_lowfreq_guard_start_meV, settings.zlp_lowfreq_guard_end_meV];
            summary.lowfreq_guard_weight = settings.zlp_lowfreq_guard_weight;
            summary.zlp_area_fraction = zlp_area ./ max(total_area, eps);
            summary.guard_points = nnz(guard);
        end

        function summary = summarize_attention_regions(result, settings)
            ene = result.ene_corrected(:);
            if isfield(result, 'stage2') && isfield(result.stage2, 'attention') && ...
                    isfield(result.stage2.attention, 'feature_score') && strcmp(result.stage2.status, 'ok')
                score = result.stage2.attention.feature_score(:);
                weight = result.stage2.attention.weight(:);
            else
                weight = result.weights_used(:);
                score = EELS_PearsonSparse_FullFit.unit_scale(weight, isfinite(weight));
            end
            regions = EELS_PearsonSparse_FullFit.energy_region_masks(ene, settings);
            names = fieldnames(regions);
            mean_score = NaN(numel(names), 1);
            max_score = NaN(numel(names), 1);
            mean_weight = NaN(numel(names), 1);
            active_fraction = NaN(numel(names), 1);
            for ii = 1:numel(names)
                mask = regions.(names{ii});
                if nnz(mask) > 0
                    vals = score(mask);
                    wvals = weight(mask);
                    mean_score(ii) = mean(vals, 'omitnan');
                    max_score(ii) = max(vals, [], 'omitnan');
                    mean_weight(ii) = mean(wvals, 'omitnan');
                    active_fraction(ii) = nnz(vals >= settings.stage2_attention_threshold) ./ max(nnz(mask), 1);
                end
            end
            summary = struct();
            summary.region_names = names(:).';
            summary.mean_attention = mean_score(:).';
            summary.max_attention = max_score(:).';
            summary.mean_weight = mean_weight(:).';
            summary.active_fraction = active_fraction(:).';
            summary.energy = ene;
            summary.attention_profile = score;
            summary.weight_profile = weight;
        end

        function breakdown = compute_loss_breakdown(result, fit, settings)
            ene = result.ene_corrected(:);
            residual = result.raw_data_norm(:) - result.fitted_spectrum_norm(:);
            sqrtW = sqrt(max(result.weights_used(:), eps));
            if isfield(settings, 'stage2_far_tail_loss_budget_downweight_enabled') && ...
                    settings.stage2_far_tail_loss_budget_downweight_enabled
                attention = [];
                if isfield(result, 'stage2') && isstruct(result.stage2) && ...
                        isfield(result.stage2, 'attention')
                    attention = result.stage2.attention;
                end
                sqrtW = EELS_PearsonSparse_FullFit.apply_stage2_far_tail_downweight(ene, sqrtW, settings, attention);
            end
            rw2 = (sqrtW .* residual) .^ 2;
            regions = EELS_PearsonSparse_FullFit.energy_region_masks(ene, settings);
            names = fieldnames(regions);
            values = zeros(1, numel(names));
            for ii = 1:numel(names)
                values(ii) = sum(rw2(regions.(names{ii})), 'omitnan');
            end
            data_sse = sum(values, 'omitnan');
            data_n = numel(result.ene_raw);
            reg_sse = 0;
            if isfield(fit, 'optimizer_residual') && numel(fit.optimizer_residual) > data_n
                reg_sse = sum(fit.optimizer_residual((data_n + 1):end) .^ 2, 'omitnan');
            end
            zlp_guard_sse = EELS_PearsonSparse_FullFit.zlp_signal_guard_sse(result, settings);
            zlp_lowfreq_sse = EELS_PearsonSparse_FullFit.zlp_lowfreq_guard_sse(result, settings);
            reg_sse = max(reg_sse - zlp_guard_sse - zlp_lowfreq_sse, 0);
            be_name = 'soft_BE';
            if isfield(result, 'stage2_primary_model') && result.stage2_primary_model
                be_sse = 0;
                be_name = 'hard_BE_consistency';
            else
                be_sse = EELS_PearsonSparse_FullFit.soft_be_penalty_value(result, settings);
            end
            tail_over = EELS_PearsonSparse_FullFit.tail_overshoot_score(result.fitted_spectrum_norm, ...
                result.raw_data_norm, ene, settings);
            tail_budget_scale = 1;
            if isfield(settings, 'stage2_far_tail_loss_budget_scale')
                tail_budget_scale = min(max(settings.stage2_far_tail_loss_budget_scale, 0), 1);
            end
            extra_names = {'regularization', 'zlp_signal_guard', 'zlp_lowfreq_guard', be_name, 'tail_overshoot'};
            extra_values = [reg_sse, zlp_guard_sse, zlp_lowfreq_sse, be_sse, ...
                tail_budget_scale .* max(tail_over, 0) .^ 2];
            all_names = [names(:).', extra_names];
            all_values = [values, extra_values];
            total = sum(all_values, 'omitnan');
            if total <= 0, total = eps; end
            breakdown = struct();
            breakdown.names = all_names;
            breakdown.values = all_values;
            breakdown.fraction = all_values ./ total;
            breakdown.data_sse = data_sse;
            breakdown.regularization_sse = reg_sse;
            breakdown.zlp_signal_guard_sse = zlp_guard_sse;
            breakdown.zlp_lowfreq_guard_sse = zlp_lowfreq_sse;
            breakdown.soft_BE_sse = be_sse;
            breakdown.tail_overshoot_score = tail_over;
            [~, imax] = max(breakdown.fraction);
            breakdown.dominant = all_names{imax};
        end

        function sse = zlp_signal_guard_sse(result, settings)
            sse = 0;
            if ~settings.zlp_signal_guard_enabled || settings.zlp_signal_guard_weight <= 0
                return;
            end
            ene = result.ene_corrected(:);
            eabs = abs(ene);
            guard_start = settings.zlp_signal_guard_start_meV;
            guard_end = settings.zlp_signal_guard_end_meV;
            guard = eabs >= guard_start & eabs <= guard_end;
            if guard_end <= guard_start || nnz(guard) < 3
                return;
            end
            taper_width = max(settings.zlp_signal_guard_taper_meV, eps);
            low_taper = 1 ./ (1 + exp(-(eabs(guard) - guard_start) ./ taper_width));
            high_taper = 1 ./ (1 + exp((eabs(guard) - guard_end) ./ taper_width));
            release_width = max(settings.zlp_signal_guard_tail_release_width_meV, eps);
            release_start = settings.zlp_signal_guard_tail_release_start_meV;
            release_floor = min(max(settings.zlp_signal_guard_tail_floor, 0), 1);
            tail_release = release_floor + (1 - release_floor) ./ ...
                (1 + exp((eabs(guard) - release_start) ./ release_width));
            guard_profile = EELS_PearsonSparse_FullFit.zlp_signal_guard_profile(ene, result.raw_data_norm(:), settings);
            zlp_norm = result.zlp_only(:) ./ max(result.norm_val, eps);
            vals = zlp_norm(guard) .* low_taper .* high_taper .* tail_release .* guard_profile(guard);
            sse = settings.zlp_signal_guard_weight .* sum(vals .^ 2, 'omitnan');
        end

        function sse = zlp_lowfreq_guard_sse(result, settings)
            sse = 0;
            if ~settings.zlp_lowfreq_guard_enabled || settings.zlp_lowfreq_guard_weight <= 0
                return;
            end
            ene = result.ene_corrected(:);
            eabs = abs(ene);
            low_start = settings.zlp_lowfreq_guard_start_meV;
            low_end = settings.zlp_lowfreq_guard_end_meV;
            guard = eabs >= low_start & eabs <= low_end;
            if low_end <= low_start || nnz(guard) < 3
                return;
            end
            taper_width = max(settings.zlp_lowfreq_guard_taper_meV, eps);
            low_taper = 1 ./ (1 + exp(-(eabs(guard) - low_start) ./ taper_width));
            high_taper = 1 ./ (1 + exp((eabs(guard) - low_end) ./ taper_width));
            guard_profile = EELS_PearsonSparse_FullFit.zlp_signal_guard_profile(ene, result.raw_data_norm(:), settings);
            zlp_norm = result.zlp_only(:) ./ max(result.norm_val, eps);
            vals = zlp_norm(guard) .* low_taper .* high_taper .* (0.35 + 0.65 .* guard_profile(guard));
            sse = settings.zlp_lowfreq_guard_weight .* sum(vals .^ 2, 'omitnan');
        end

        function profile = zlp_signal_guard_profile(ene, y_norm, settings)
            ene = ene(:);
            if ~settings.zlp_signal_guard_feature_weighted || isempty(y_norm)
                profile = ones(size(ene));
                return;
            end
            y = double(y_norm(:));
            if numel(y) ~= numel(ene)
                profile = ones(size(ene));
                return;
            end
            y(~isfinite(y)) = 0;
            y = max(y, 0);
            eabs = abs(ene);
            guard_mask = eabs >= settings.zlp_signal_guard_start_meV & eabs <= settings.zlp_signal_guard_end_meV;
            if nnz(guard_mask) < 8
                profile = ones(size(ene));
                return;
            end
            smooth_y = movmedian(y, 7, 'omitnan');
            smooth_y = movmean(smooth_y, 7, 'omitnan');
            noise = EELS_PearsonSparse_FullFit.robust_noise(smooth_y, guard_mask);
            local_base = movmedian(smooth_y, 35, 'omitnan');
            contrast = max(smooth_y - local_base, 0) ./ max(local_base + 4 .* noise, eps);
            log_y = log(max(smooth_y, max(noise, eps)));
            slope = [0; abs(diff(movmean(log_y, 7, 'omitnan')))];
            curvature = [0; abs(diff(movmean(slope, 7, 'omitnan')))];
            feature = 0.58 .* EELS_PearsonSparse_FullFit.unit_scale(contrast, guard_mask) + ...
                0.27 .* EELS_PearsonSparse_FullFit.unit_scale(slope, guard_mask) + ...
                0.15 .* EELS_PearsonSparse_FullFit.unit_scale(curvature, guard_mask);
            n_smooth = max(3, round(settings.zlp_signal_guard_feature_smooth_points));
            feature = movmean(feature, n_smooth, 'omitnan');
            feature(~guard_mask | ~isfinite(feature)) = 0;
            floor_val = min(max(settings.zlp_signal_guard_min_profile, 0), 1);
            profile = floor_val + (1 - floor_val) .* min(max(feature, 0), 1);
        end

        function ci = energy_confidence_interval(result, settings)
            residual = result.residual(:);
            n = numel(residual);
            win = min(max(5, settings.ci_window_points), max(5, n));
            if mod(win, 2) == 0, win = win + 1; end
            local_center = movmedian(residual, win, 'omitnan');
            local_mad = movmedian(abs(residual - local_center), win, 'omitnan');
            sigma = 1.4826 .* local_mad;
            floor_sigma = max(settings.ci_min_noise_fraction .* max(abs(result.raw_data), [], 'omitnan'), eps);
            sigma(~isfinite(sigma) | sigma < floor_sigma) = floor_sigma;
            z = EELS_PearsonSparse_FullFit.normal_z_from_level(settings.ci_level);
            half_width = z .* sigma;
            ci = struct();
            ci.level = settings.ci_level;
            ci.energy = result.ene_corrected(:);
            ci.sigma = sigma(:);
            ci.half_width = half_width(:);
            ci.fit_lower = result.fitted_spectrum(:) - half_width(:);
            ci.fit_upper = result.fitted_spectrum(:) + half_width(:);
            if isfield(result, 'stage2') && isfield(result.stage2, 'model_norm') && strcmp(result.stage2.status, 'ok')
                stage2_fit = result.stage2.model_norm(:) .* result.norm_val;
                stage2_res = result.raw_data(:) - stage2_fit;
                stage2_center = movmedian(stage2_res, win, 'omitnan');
                stage2_mad = movmedian(abs(stage2_res - stage2_center), win, 'omitnan');
                stage2_sigma = max(1.4826 .* stage2_mad, floor_sigma);
                ci.stage2_sigma = stage2_sigma(:);
                ci.stage2_lower = stage2_fit - z .* stage2_sigma(:);
                ci.stage2_upper = stage2_fit + z .* stage2_sigma(:);
            end
        end

        function budget = compute_error_budget(result, ~)
            lb = result.loss_breakdown;
            budget = struct();
            budget.dominant_error_source = lb.dominant;
            budget.weighted_NRMSE = result.metrics.weighted_NRMSE;
            budget.highfreq_NRMSE = result.metrics.highfreq_NRMSE;
            budget.tail_overshoot_score = lb.tail_overshoot_score;
            budget.energy_dependent_ci_median = median(result.confidence_interval.half_width, 'omitnan');
            budget.energy_dependent_ci_p90 = EELS_PearsonSparse_FullFit.robust_quantile(result.confidence_interval.half_width, 0.90);
            budget.T_uncertainty_K = result.T_uncertainty;
            budget.reliability_score = result.reliability_score;
            budget.loss_fraction = lb.fraction;
            budget.loss_names = lb.names;
            if isfield(result.stage1, 'inner_solver') && isfield(result.stage1.inner_solver, 'solver')
                budget.inner_solver = result.stage1.inner_solver.solver;
            else
                budget.inner_solver = '';
            end
            if isfield(result, 'stage2') && isfield(result.stage2, 'temperature_policy')
                budget.temperature_policy = result.stage2.temperature_policy;
            else
                budget.temperature_policy = '';
            end
        end

        function regions = energy_region_masks(ene, settings)
            ene = ene(:);
            eabs = abs(ene);
            core = eabs < settings.core_downweight_width;
            shoulder = eabs >= settings.core_downweight_width & eabs < settings.em_min_meV;
            anti_signal = ene <= -settings.em_min_meV & eabs < settings.highfreq_start_meV;
            stokes_signal = ene >= settings.em_min_meV & eabs < settings.highfreq_start_meV;
            tail_start = settings.tail_overshoot_start_meV;
            if isfield(settings, 'stage2_far_tail_downweight_start_meV')
                tail_start = max(tail_start, settings.stage2_far_tail_downweight_start_meV);
            end
            highfreq = eabs >= settings.highfreq_start_meV & eabs < tail_start;
            tail = eabs >= tail_start;
            regions = struct('zlp_core', core, 'near_shoulder', shoulder, ...
                'anti_signal', anti_signal, 'stokes_signal', stokes_signal, ...
                'high_frequency', highfreq, 'far_tail', tail);
        end

        function be_sse = soft_be_penalty_value(result, settings)
            if ~isfield(result, 'PM_stokes') || ~isfield(result, 'PM_antistokes') || isempty(result.EM_final)
                be_sse = NaN;
                return;
            end
            kb = 0.08617333262;
            q = exp(-result.EM_final(:) ./ max(kb .* result.T_soft, eps));
            diffv = result.PM_antistokes(:) - result.PM_stokes(:) .* q;
            be_sse = settings.lambda_be_soft .* sum(diffv .^ 2, 'omitnan');
        end

        function z = normal_z_from_level(level)
            if level >= 0.99
                z = 2.576;
            elseif level >= 0.95
                z = 1.96;
            elseif level >= 0.90
                z = 1.645;
            else
                z = 1.0;
            end
        end

        function mask = stage2_transition_mask(stage1_result, ene, win, settings)
            ene = ene(:);
            if ~settings.stage2_transition_fit_enabled || numel(win) < 2
                mask = false(size(ene));
                return;
            end
            transition_start = settings.stage2_transition_min_meV;
            if isfield(stage1_result, 'zlp_region_policy') && isfield(stage1_result.zlp_region_policy, 'guard_start_meV')
                transition_start = max(transition_start, stage1_result.zlp_region_policy.guard_start_meV);
            elseif isfield(settings, 'zlp_signal_guard_start_meV')
                transition_start = max(transition_start, settings.zlp_signal_guard_start_meV);
            end
            transition_stop = win(1);
            if ~(isfinite(transition_start) && isfinite(transition_stop)) || transition_stop <= transition_start
                mask = false(size(ene));
                return;
            end
            eabs = abs(ene);
            mask = eabs >= transition_start & eabs < transition_stop & eabs <= win(2);
        end

        function attention = build_stage2_attention(stage1_result, settings)
            ene = stage1_result.ene_corrected(:);
            eabs = abs(ene);
            hard_win = sort(abs(settings.stage2_window_meV(:).'));
            if numel(hard_win) < 2, hard_win = [40, 150]; end
            att_win = sort(abs(settings.stage2_attention_window_meV(:).'));
            if numel(att_win) < 2, att_win = hard_win; end
            att_win(1) = min(att_win(1), hard_win(1));
            att_win(2) = max(att_win(2), hard_win(2));
            broad = eabs >= att_win(1) & eabs <= att_win(2);
            signal = abs(stage1_result.raw_data_norm(:) - stage1_result.stage1.zlp_norm(:));
            signal = movmedian(signal, 7, 'omitnan');
            grad_signal = [0; abs(diff(movmean(signal, 9, 'omitnan')))];
            resid = abs(stage1_result.residual_norm(:));

            pm_line = zeros(size(ene));
            EM = stage1_result.EM_final(:).';
            pm_s = stage1_result.PM_stokes(:).';
            pm_a = stage1_result.PM_antistokes(:).';
            sigma = max(settings.stage2_attention_peak_sigma_meV, eps);
            if ~isempty(EM)
                scale_pm = max([pm_s(:); pm_a(:); eps]);
                for ii = 1:numel(EM)
                    amp_s = pm_s(ii) ./ scale_pm;
                    amp_a = pm_a(ii) ./ scale_pm;
                    pm_line = pm_line + amp_s .* exp(-0.5 .* ((ene - EM(ii)) ./ sigma) .^ 2);
                    pm_line = pm_line + amp_a .* exp(-0.5 .* ((ene + EM(ii)) ./ sigma) .^ 2);
                end
            end

            feature = 0.48 .* EELS_PearsonSparse_FullFit.unit_scale(signal, broad) + ...
                0.20 .* EELS_PearsonSparse_FullFit.unit_scale(grad_signal, broad) + ...
                0.22 .* EELS_PearsonSparse_FullFit.unit_scale(pm_line, broad) + ...
                0.10 .* EELS_PearsonSparse_FullFit.unit_scale(resid, broad);
            anti_side = broad & ene < -hard_win(1);
            if any(anti_side)
                anti_signal = EELS_PearsonSparse_FullFit.unit_scale(signal, anti_side);
                anti_grad = EELS_PearsonSparse_FullFit.unit_scale(grad_signal, anti_side);
                anti_pm = EELS_PearsonSparse_FullFit.unit_scale(pm_line, anti_side);
                anti_resid = EELS_PearsonSparse_FullFit.unit_scale(resid, anti_side);
                anti_feature = 0.52 .* anti_signal + 0.20 .* anti_grad + 0.22 .* anti_pm + 0.06 .* anti_resid;
                anti_noise = EELS_PearsonSparse_FullFit.robust_noise(signal, anti_side);
                anti_pm_floor = 0.04 .* max(pm_line(anti_side), [], 'omitnan');
                if ~isfinite(anti_pm_floor), anti_pm_floor = 0; end
                anti_gate = anti_side & (signal >= settings.stage2_antistokes_noise_sigma .* anti_noise | pm_line >= anti_pm_floor);
                feature(anti_gate) = max(feature(anti_gate), ...
                    min(1, settings.stage2_antistokes_attention_boost .* anti_feature(anti_gate)));
            end
            mirrored = interp1(ene, feature, -ene, 'linear', 0);
            feature = max(feature, 0.65 .* mirrored);
            feature(~broad | ~isfinite(feature)) = 0;
            feature = min(max(movmean(feature, 5, 'omitnan'), 0), 1);
            active = broad & feature >= settings.stage2_attention_threshold;
            active_mirror = interp1(ene, double(active), -ene, 'nearest', 0) > 0;
            active = broad & (active | active_mirror);

            weight = ones(size(ene));
            weight(broad) = settings.stage2_attention_floor + settings.stage2_attention_gain .* feature(broad);
            weight = max(weight, settings.stage2_attention_floor);

            attention = struct();
            attention.energy = ene;
            attention.feature_score = feature;
            attention.weight = weight;
            attention.valid_mask = broad;
            attention.window_meV = att_win;
            attention.hard_be_window_meV = hard_win;
            attention.active_mask = active;
            attention.selected_fraction = nnz(active) ./ max(nnz(broad), 1);
            attention.bands_meV = EELS_PearsonSparse_FullFit.mask_to_bands(ene, active);
            attention = EELS_PearsonSparse_FullFit.update_stage2_attention_tail_masks(attention, ene, settings);
        end

        function attention = update_stage2_attention_tail_masks(attention, ene, settings)
            ene = ene(:);
            n = numel(ene);
            if ~(isstruct(attention) && n > 0)
                return;
            end
            score = zeros(n, 1);
            if isfield(attention, 'feature_score') && numel(attention.feature_score) == n
                score = max(score, min(max(attention.feature_score(:), 0), 1));
            end
            if isfield(attention, 'active_mask') && numel(attention.active_mask) == n
                score(logical(attention.active_mask(:))) = 1;
            end
            if isfield(attention, 'transition_mask') && numel(attention.transition_mask) == n
                score(logical(attention.transition_mask(:))) = max(score(logical(attention.transition_mask(:))), 0.5);
            end
            if isfield(attention, 'fit_candidate_alpha') && numel(attention.fit_candidate_alpha) == n
                score = max(score, min(max(attention.fit_candidate_alpha(:), 0), 1));
            elseif isfield(attention, 'fit_candidate_mask') && numel(attention.fit_candidate_mask) == n
                score(logical(attention.fit_candidate_mask(:))) = max(score(logical(attention.fit_candidate_mask(:))), ...
                    min(max(settings.stage2_attention_tail_signal_threshold, 0), 1));
            end
            valid = true(n, 1);
            if isfield(attention, 'valid_mask') && numel(attention.valid_mask) == n
                valid = logical(attention.valid_mask(:));
            end
            min_abs = settings.stage2_attention_tail_min_abs_meV;
            if isfield(settings, 'highfreq_start_meV')
                min_abs = max(min_abs, settings.highfreq_start_meV);
            end
            threshold = min(max(settings.stage2_attention_tail_signal_threshold, 0), 1);
            eabs = abs(ene);
            high = eabs >= min_abs;
            signal_mask = valid & high & score >= threshold;
            if isfield(attention, 'active_mask') && numel(attention.active_mask) == n
                signal_mask = signal_mask | (valid & high & logical(attention.active_mask(:)));
            end
            if isfield(attention, 'fit_candidate_mask') && numel(attention.fit_candidate_mask) == n
                signal_mask = signal_mask | (valid & high & logical(attention.fit_candidate_mask(:)));
            end
            tail_mask = high & ~signal_mask;
            attention.signal_score = score;
            attention.signal_mask = signal_mask;
            attention.tail_mask = tail_mask;
        end

        function score = stage2_attention_score_profile(attention)
            n = 0;
            if isstruct(attention) && isfield(attention, 'energy')
                n = numel(attention.energy);
            end
            score = zeros(n, 1);
            if n == 0
                return;
            end
            if isfield(attention, 'feature_score') && numel(attention.feature_score) == n
                score = max(score, attention.feature_score(:));
            end
            if isfield(attention, 'signal_score') && numel(attention.signal_score) == n
                score = max(score, attention.signal_score(:));
            end
            if isfield(attention, 'fit_candidate_alpha') && numel(attention.fit_candidate_alpha) == n
                score = max(score, attention.fit_candidate_alpha(:));
            end
            if isfield(attention, 'active_mask') && numel(attention.active_mask) == n
                score(logical(attention.active_mask(:))) = 1;
            end
            score(~isfinite(score)) = 0;
            score = min(max(score, 0), 1);
        end

        function basis = build_stage2_fixed_stage1_basis(stage1_result, attention, settings)
            base_EM = stage1_result.EM_final(:).';
            if isempty(base_EM) || any(~isfinite(base_EM)) || any(base_EM <= 0)
                basis = struct('energy_axis', stage1_result.ene_corrected(:), 'EM_all', zeros(1, 0), ...
                    'EM_hard', zeros(1, 0), 'EM_soft', zeros(1, 0), 'hard_idx', false(1, 0), ...
                    'soft_idx', false(1, 0), 'hard_fraction', 0, ...
                    'hard_selection_status', 'empty_stage1_em_grid');
                return;
            end
            pm_s_all = stage1_result.PM_stokes(:).';
            pm_a_all = stage1_result.PM_antistokes(:).';
            if numel(pm_s_all) ~= numel(base_EM)
                pm_s_all = EELS_PearsonSparse_FullFit.interpolate_pm(stage1_result.EM_final(:), stage1_result.PM_stokes(:), base_EM(:));
                pm_s_all = pm_s_all(:).';
            end
            if numel(pm_a_all) ~= numel(base_EM)
                pm_a_all = EELS_PearsonSparse_FullFit.interpolate_pm(stage1_result.EM_final(:), stage1_result.PM_antistokes(:), base_EM(:));
                pm_a_all = pm_a_all(:).';
            end
            T_stage1_seed = NaN;
            if isfield(stage1_result, 'T_soft') && isfinite(stage1_result.T_soft)
                T_stage1_seed = stage1_result.T_soft;
            elseif isfield(stage1_result, 'T_fit') && isfinite(stage1_result.T_fit)
                T_stage1_seed = stage1_result.T_fit;
            end
            if ~isfinite(T_stage1_seed)
                T_stage1_seed = settings.T_soft_init;
            end
            bounds = settings.T_soft_bounds(:).';
            if numel(bounds) >= 2 && all(isfinite(bounds(1:2)))
                T_stage1_seed = min(max(T_stage1_seed, min(bounds(1:2))), max(bounds(1:2)));
            end
            kb = 0.08617333262;
            q_seed = exp(-base_EM ./ max(kb .* T_stage1_seed, eps));
            attention_score = EELS_PearsonSparse_FullFit.stage2_attention_score_profile(attention);
            if isempty(attention_score) || ~isfield(attention, 'energy') || numel(attention.energy) ~= numel(attention_score)
                hard_feature_score = zeros(size(base_EM));
                hard_pair_score = zeros(size(base_EM));
            else
                att_pos = interp1(attention.energy(:), attention_score(:), base_EM(:), 'linear', 0);
                att_neg = interp1(attention.energy(:), attention_score(:), -base_EM(:), 'linear', 0);
                hard_feature_score = max(max(att_pos(:), 0), max(att_neg(:), 0)).';
                hard_pair_score = sqrt(max(att_pos(:), 0) .* max(att_neg(:), 0)).';
            end
            hard_idx = true(size(base_EM));
            soft_idx = false(size(base_EM));
            hard_selection_status = 'legacy_all_hard_be';
            threshold = 0.35;
            floor_threshold = 0.10;
            min_pairs = 0;
            if isfield(settings, 'stage2_hard_be_pair_threshold')
                threshold = min(max(settings.stage2_hard_be_pair_threshold, 0), 1);
            end
            if isfield(settings, 'stage2_hard_be_pair_floor')
                floor_threshold = min(max(settings.stage2_hard_be_pair_floor, 0), threshold);
            end
            if isfield(settings, 'stage2_ratio_min_pairs')
                min_pairs = max(0, round(settings.stage2_ratio_min_pairs));
            end
            if ~isfield(settings, 'stage2_attention_hardbe_only_enabled') || settings.stage2_attention_hardbe_only_enabled
                hard_idx = hard_pair_score >= threshold;
                hard_selection_status = 'threshold_ok';
                if nnz(hard_idx) < min_pairs
                    candidates = find(~hard_idx & hard_pair_score >= floor_threshold);
                    if ~isempty(candidates)
                        [~, ord] = sort(hard_pair_score(candidates), 'descend');
                        need = min(numel(candidates), min_pairs - nnz(hard_idx));
                        if need > 0
                            hard_idx(candidates(ord(1:need))) = true;
                            hard_selection_status = 'floor_backfill';
                        end
                    end
                end
                if nnz(hard_idx) < min_pairs
                    hard_selection_status = 'insufficient_hard_be_attention_pairs';
                end
                soft_idx = ~hard_idx;
            end
            hard_idx = logical(hard_idx(:).');
            soft_idx = logical(soft_idx(:).');
            stage1_soft_EM = base_EM(soft_idx);
            stage1_hard_EM = base_EM(hard_idx);
            hard_dense_EM = EELS_PearsonSparse_FullFit.stage2_hard_dense_em_grid( ...
                attention, settings, threshold);
            stage2_hard_EM = sort(unique([stage1_hard_EM(:); hard_dense_EM(:)]).');
            hard_dense_added_EM = setdiff(stage2_hard_EM(:).', base_EM(:).');
            soft_dense_EM = EELS_PearsonSparse_FullFit.stage2_soft_dense_em_grid( ...
                stage2_hard_EM, attention, settings, threshold);
            if ~isempty(soft_dense_EM)
                stage1_soft_EM = sort(unique([stage1_soft_EM(:); soft_dense_EM(:)]).');
            end
            stage2_EM_all = sort(unique([stage2_hard_EM(:); stage1_soft_EM(:)]).');
            if isempty(stage2_EM_all)
                stage2_EM_all = base_EM;
            end
            hard_idx_stage2 = EELS_PearsonSparse_FullFit.em_membership_mask(stage2_EM_all, stage2_hard_EM, 1e-8);
            soft_idx_stage2 = ~hard_idx_stage2;
            stage2_pair_score = zeros(size(stage2_EM_all));
            stage2_feature_score = zeros(size(stage2_EM_all));
            if ~isempty(attention_score) && isfield(attention, 'energy') && numel(attention.energy) == numel(attention_score)
                att_pos2 = interp1(attention.energy(:), attention_score(:), stage2_EM_all(:), 'linear', 0);
                att_neg2 = interp1(attention.energy(:), attention_score(:), -stage2_EM_all(:), 'linear', 0);
                stage2_feature_score = max(max(att_pos2(:), 0), max(att_neg2(:), 0)).';
                stage2_pair_score = sqrt(max(att_pos2(:), 0) .* max(att_neg2(:), 0)).';
            end
            pm_s_stage2 = EELS_PearsonSparse_FullFit.interpolate_pm(base_EM(:), pm_s_all(:), stage2_EM_all(:));
            pm_a_stage2 = EELS_PearsonSparse_FullFit.interpolate_pm(base_EM(:), pm_a_all(:), stage2_EM_all(:));
            q_seed_stage2 = exp(-stage2_EM_all ./ max(kb .* T_stage1_seed, eps));
            basis = struct();
            basis.energy_axis = stage1_result.ene_corrected(:);
            basis.EM_all = stage2_EM_all;
            basis.EM_stage1 = base_EM;
            basis.EM_hard = stage2_EM_all(hard_idx_stage2);
            basis.EM_soft = stage2_EM_all(soft_idx_stage2);
            basis.EM_hard_stage1 = base_EM(hard_idx);
            basis.EM_hard_dense_added = hard_dense_added_EM;
            basis.EM_soft_stage1 = base_EM(soft_idx);
            basis.EM_soft_dense_added = setdiff(stage1_soft_EM(:).', base_EM(:).');
            basis.hard_idx = hard_idx_stage2;
            basis.soft_idx = soft_idx_stage2;
            basis.stage1_hard_idx = hard_idx;
            basis.stage1_soft_idx = soft_idx;
            basis.hard_fraction = nnz(hard_idx_stage2) ./ max(numel(stage2_EM_all), 1);
            basis.stage1_PM_stokes_all = pm_s_stage2(:).';
            basis.stage1_PM_antistokes_all = pm_a_stage2(:).';
            basis.stage1_PM_stokes_hard = pm_s_stage2(hard_idx_stage2).';
            basis.stage1_PM_antistokes_hard = pm_a_stage2(hard_idx_stage2).';
            basis.stage1_PM_stokes_soft = pm_s_stage2(soft_idx_stage2).';
            basis.stage1_PM_antistokes_soft = pm_a_stage2(soft_idx_stage2).';
            basis.T_stage1_seed_K = T_stage1_seed;
            basis.PM_antistokes_initial_hardBE = pm_s_stage2(hard_idx_stage2).' .* q_seed_stage2(hard_idx_stage2);
            basis.hard_be_seed_ratio = q_seed_stage2(hard_idx_stage2);
            basis.q_seed_all = q_seed_stage2(:).';
            basis.hard_pair_score = stage2_pair_score(:).';
            basis.hard_pair_score_stage1 = hard_pair_score(:).';
            basis.hard_pair_score_hard = stage2_pair_score(hard_idx_stage2);
            basis.soft_pair_score = stage2_pair_score(soft_idx_stage2);
            basis.hard_feature_score = stage2_feature_score(:).';
            basis.hard_feature_score_stage1 = hard_feature_score(:).';
            basis.hard_selection_status = hard_selection_status;
            basis.hard_pair_threshold = threshold;
            basis.hard_pair_floor = floor_threshold;
            basis.min_hard_pairs = min_pairs;
            basis.num_hard_pairs = nnz(hard_idx_stage2);
            basis.num_soft_pairs = nnz(soft_idx_stage2);
            basis.num_stage1_em = numel(base_EM);
            basis.num_stage2_em = numel(stage2_EM_all);
            basis.num_hard_dense_added = numel(basis.EM_hard_dense_added);
            basis.num_soft_dense_added = numel(basis.EM_soft_dense_added);
        end

        function hard_EM = stage2_hard_dense_em_grid(attention, settings, hard_threshold)
            hard_EM = zeros(1, 0);
            if ~isfield(settings, 'stage2_hard_dense_em_enabled') || ~settings.stage2_hard_dense_em_enabled
                return;
            end
            lo_candidates = [];
            if isfield(settings, 'stage2_window_meV') && ~isempty(settings.stage2_window_meV)
                lo_candidates(end + 1) = min(abs(settings.stage2_window_meV(:)), [], 'omitnan'); %#ok<AGROW>
            end
            lo_candidates(end + 1) = max(0, settings.em_min_meV); %#ok<AGROW>
            lo = max(lo_candidates(isfinite(lo_candidates)), [], 'omitnan');
            if isfield(settings, 'stage2_hard_dense_min_meV') && ~isempty(settings.stage2_hard_dense_min_meV) && ...
                    isfinite(settings.stage2_hard_dense_min_meV)
                lo = max(0, settings.stage2_hard_dense_min_meV);
            end
            hi_candidates = [];
            if isfield(settings, 'stage2_attention_window_meV') && ~isempty(settings.stage2_attention_window_meV)
                hi_candidates(end + 1) = max(abs(settings.stage2_attention_window_meV(:)), [], 'omitnan'); %#ok<AGROW>
            end
            if isfield(settings, 'stage2_window_meV') && ~isempty(settings.stage2_window_meV)
                hi_candidates(end + 1) = max(abs(settings.stage2_window_meV(:)), [], 'omitnan'); %#ok<AGROW>
            end
            hi = max(hi_candidates(isfinite(hi_candidates)), [], 'omitnan');
            if isempty(hi) || ~isfinite(hi) || hi <= lo
                hi = 220;
            end
            if isfield(settings, 'stage2_hard_dense_extra_meV') && isfinite(settings.stage2_hard_dense_extra_meV)
                hi = hi + max(settings.stage2_hard_dense_extra_meV, 0);
            end
            if isfield(settings, 'stage2_hard_dense_max_meV') && ~isempty(settings.stage2_hard_dense_max_meV) && ...
                    isfinite(settings.stage2_hard_dense_max_meV)
                hi = min(hi, settings.stage2_hard_dense_max_meV);
            end
            if isstruct(attention) && isfield(attention, 'energy') && ~isempty(attention.energy)
                hi = min(hi, max(abs(attention.energy(:)), [], 'omitnan'));
            end
            step = 4.0;
            if isfield(settings, 'stage2_hard_dense_step_meV') && isfinite(settings.stage2_hard_dense_step_meV)
                step = max(0.5, settings.stage2_hard_dense_step_meV);
            end
            if ~(isfinite(lo) && isfinite(hi) && hi >= lo)
                return;
            end
            cand = lo:step:hi;
            cand = cand(isfinite(cand) & cand > 0);
            if isempty(cand)
                return;
            end
            attention_score = EELS_PearsonSparse_FullFit.stage2_attention_score_profile(attention);
            if isempty(attention_score) || ~isfield(attention, 'energy') || numel(attention.energy) ~= numel(attention_score)
                return;
            end
            att_pos = interp1(attention.energy(:), attention_score(:), cand(:), 'linear', 0);
            att_neg = interp1(attention.energy(:), attention_score(:), -cand(:), 'linear', 0);
            pair_score = sqrt(max(att_pos(:), 0) .* max(att_neg(:), 0)).';
            hard_EM = sort(unique(cand(pair_score >= hard_threshold)));
        end

        function soft_EM = stage2_soft_dense_em_grid(hard_EM, attention, settings, hard_threshold)
            soft_EM = zeros(1, 0);
            if ~isfield(settings, 'stage2_soft_dense_em_enabled') || ~settings.stage2_soft_dense_em_enabled
                return;
            end
            lo_candidates = [];
            if isfield(settings, 'stage2_window_meV') && ~isempty(settings.stage2_window_meV)
                lo_candidates(end + 1) = min(abs(settings.stage2_window_meV(:)), [], 'omitnan'); %#ok<AGROW>
            end
            lo_candidates(end + 1) = max(0, settings.em_min_meV); %#ok<AGROW>
            lo = max(lo_candidates(isfinite(lo_candidates)), [], 'omitnan');
            if isfield(settings, 'stage2_soft_dense_min_meV') && ~isempty(settings.stage2_soft_dense_min_meV) && ...
                    isfinite(settings.stage2_soft_dense_min_meV)
                lo = max(0, settings.stage2_soft_dense_min_meV);
            end
            hi_candidates = [];
            if isfield(settings, 'stage2_attention_window_meV') && ~isempty(settings.stage2_attention_window_meV)
                hi_candidates(end + 1) = max(abs(settings.stage2_attention_window_meV(:)), [], 'omitnan'); %#ok<AGROW>
            end
            if isfield(settings, 'stage2_window_meV') && ~isempty(settings.stage2_window_meV)
                hi_candidates(end + 1) = max(abs(settings.stage2_window_meV(:)), [], 'omitnan'); %#ok<AGROW>
            end
            hi = max(hi_candidates(isfinite(hi_candidates)), [], 'omitnan');
            if isempty(hi) || ~isfinite(hi) || hi <= lo
                hi = 220;
            end
            if isfield(settings, 'stage2_soft_dense_extra_meV') && isfinite(settings.stage2_soft_dense_extra_meV)
                hi = hi + max(settings.stage2_soft_dense_extra_meV, 0);
            end
            if isfield(settings, 'stage2_soft_dense_max_meV') && ~isempty(settings.stage2_soft_dense_max_meV) && ...
                    isfinite(settings.stage2_soft_dense_max_meV)
                hi = min(hi, settings.stage2_soft_dense_max_meV);
            end
            if isstruct(attention) && isfield(attention, 'energy') && ~isempty(attention.energy)
                hi = min(hi, max(abs(attention.energy(:)), [], 'omitnan'));
            end
            step = 4.0;
            if isfield(settings, 'stage2_soft_dense_step_meV') && isfinite(settings.stage2_soft_dense_step_meV)
                step = max(0.5, settings.stage2_soft_dense_step_meV);
            end
            if ~(isfinite(lo) && isfinite(hi) && hi >= lo)
                return;
            end
            cand = lo:step:hi;
            cand = cand(isfinite(cand) & cand > 0);
            if isempty(cand)
                return;
            end
            pair_score = zeros(size(cand));
            attention_score = EELS_PearsonSparse_FullFit.stage2_attention_score_profile(attention);
            if ~isempty(attention_score) && isfield(attention, 'energy') && numel(attention.energy) == numel(attention_score)
                att_pos = interp1(attention.energy(:), attention_score(:), cand(:), 'linear', 0);
                att_neg = interp1(attention.energy(:), attention_score(:), -cand(:), 'linear', 0);
                pair_score = sqrt(max(att_pos(:), 0) .* max(att_neg(:), 0)).';
            end
            max_pair = hard_threshold;
            if isfield(settings, 'stage2_soft_dense_pair_score_max') && isfinite(settings.stage2_soft_dense_pair_score_max)
                max_pair = min(max(settings.stage2_soft_dense_pair_score_max, 0), 1);
            end
            keep = pair_score < max_pair;
            hard_gap = 1.25;
            if isfield(settings, 'stage2_soft_dense_hard_exclusion_meV') && isfinite(settings.stage2_soft_dense_hard_exclusion_meV)
                hard_gap = max(0, settings.stage2_soft_dense_hard_exclusion_meV);
            end
            hard_EM = hard_EM(:).';
            if ~isempty(hard_EM) && hard_gap > 0
                for ii = 1:numel(cand)
                    if any(abs(cand(ii) - hard_EM) <= hard_gap)
                        keep(ii) = false;
                    end
                end
            end
            soft_EM = sort(unique(cand(keep)));
        end

        function mask = em_membership_mask(query_EM, ref_EM, tol)
            query_EM = query_EM(:).';
            ref_EM = ref_EM(:).';
            if nargin < 3 || isempty(tol) || ~isfinite(tol)
                tol = 1e-8;
            end
            mask = false(size(query_EM));
            if isempty(query_EM) || isempty(ref_EM)
                return;
            end
            for ii = 1:numel(query_EM)
                scale = max(1, abs(query_EM(ii)));
                mask(ii) = any(abs(query_EM(ii) - ref_EM) <= tol .* scale);
            end
        end

        function tf = stage2_fixed_grid_matches_stage1(stage2_EM, stage1_EM)
            stage2_EM = stage2_EM(:).';
            stage1_EM = stage1_EM(:).';
            tf = isequal(size(stage2_EM), size(stage1_EM));
            if tf && ~isempty(stage1_EM)
                tol = max(1e-10, 1e-10 .* max(abs(stage1_EM), [], 'omitnan'));
                tf = all(abs(stage2_EM - stage1_EM) <= tol);
            end
        end

        function [r2, pair_count] = stage2_hardbe_linearity_diagnostic(EM, pm_s, pm_a, T, weights)
            kb = 0.08617333262;
            EM = EM(:);
            pm_s = pm_s(:);
            pm_a = pm_a(:);
            if nargin < 5 || isempty(weights), weights = ones(size(EM)); end
            weights = weights(:);
            if numel(weights) ~= numel(EM), weights = ones(size(EM)); end
            log_ratio = log(max(pm_s, eps) ./ max(pm_a, eps));
            pred = EM ./ max(kb .* T, eps);
            good = isfinite(log_ratio) & isfinite(pred) & isfinite(weights) & weights > 0 & pm_s > 0 & pm_a > 0;
            pair_count = nnz(good);
            if pair_count < 2
                r2 = NaN;
                return;
            end
            w = weights(good) ./ max(sum(weights(good)), eps);
            y = log_ratio(good);
            ybar = sum(w .* y);
            ss_tot = sum(w .* (y - ybar) .^ 2);
            ss_res = sum(w .* (y - pred(good)) .^ 2);
            if ss_tot <= eps
                r2 = double(ss_res <= 1e-10);
            else
                r2 = 1 - ss_res ./ ss_tot;
            end
            r2 = min(max(r2, -Inf), 1);
        end

        function D = local_diff_matrix(x, order, max_gap)
            x = x(:);
            n = numel(x);
            D = [];
            if n < order + 1
                return;
            end
            if nargin < 3 || isempty(max_gap) || ~isfinite(max_gap)
                max_gap = Inf;
            end
            [x, sort_idx] = sort(x);
            if order == 1
                gaps = diff(x);
                keep = find(isfinite(gaps) & gaps > 0 & gaps <= max_gap);
                if isempty(keep), return; end
                Dsort = zeros(numel(keep), n);
                for kk = 1:numel(keep)
                    ii = keep(kk);
                    scale = 1 ./ max(gaps(ii), eps);
                    Dsort(kk, ii) = -scale;
                    Dsort(kk, ii + 1) = scale;
                end
            elseif order == 2
                gaps = diff(x);
                keep = find(isfinite(gaps(1:end-1)) & isfinite(gaps(2:end)) & ...
                    gaps(1:end-1) > 0 & gaps(2:end) > 0 & ...
                    gaps(1:end-1) <= max_gap & gaps(2:end) <= max_gap);
                if isempty(keep), return; end
                ref_gap = median(gaps(isfinite(gaps) & gaps > 0), 'omitnan');
                if ~isfinite(ref_gap) || ref_gap <= 0, ref_gap = 1; end
                Dsort = zeros(numel(keep), n);
                for kk = 1:numel(keep)
                    ii = keep(kk);
                    local_gap = 0.5 .* (gaps(ii) + gaps(ii + 1));
                    scale = (ref_gap ./ max(local_gap, eps)) .^ 2;
                    Dsort(kk, ii) = scale;
                    Dsort(kk, ii + 1) = -2 .* scale;
                    Dsort(kk, ii + 2) = scale;
                end
            else
                return;
            end
            D = zeros(size(Dsort));
            D(:, sort_idx) = Dsort;
        end

        function val = stage2_pm_second_diff_cost(EM, pm, settings)
            val = 0;
            if ~isfield(settings, 'stage2_pm_smooth_weight') || settings.stage2_pm_smooth_weight <= 0
                return;
            end
            EM = EM(:);
            pm = pm(:);
            if numel(EM) < 3 || numel(pm) ~= numel(EM)
                return;
            end
            D2 = EELS_PearsonSparse_FullFit.local_diff_matrix(EM, 2, settings.stage2_pm_smooth_max_gap_meV);
            if isempty(D2)
                return;
            end
            rr = D2 * pm;
            val = settings.stage2_pm_smooth_weight .* mean(rr .^ 2, 'omitnan');
        end



        function scan = stage2_scan_temperatures_adaptive(T_grid, basis, kernels, A, layout, fit_mask, y_norm, zlp_norm, sqrtW, stage1_result, settings)
            T_full = unique(T_grid(:).');
            T_full = T_full(isfinite(T_full) & T_full > 0);
            min_full_points = 25;
            if isfield(settings, 'stage2_adaptive_grid_min_full_points')
                min_full_points = max(3, round(settings.stage2_adaptive_grid_min_full_points));
            end
            if ~isfield(settings, 'stage2_adaptive_grid_scan_enabled') || ...
                    ~settings.stage2_adaptive_grid_scan_enabled || numel(T_full) <= min_full_points
                scan = EELS_PearsonSparse_FullFit.stage2_scan_temperatures(T_full, basis, kernels, A, layout, ...
                    fit_mask, y_norm, zlp_norm, sqrtW, stage1_result, settings);
                scan.scan_strategy = struct('mode', 'full_grid', 'requested_points', numel(T_full), ...
                    'coarse_points', 0, 'fine_points', 0, 'scanned_points', numel(scan.T_grid), ...
                    'coarse_best_T', NaN, 'fine_half_width_K', NaN);
                scan.temperature_profile.scan_strategy = scan.scan_strategy;
                return;
            end

            coarse_grid = EELS_PearsonSparse_FullFit.stage2_coarse_temperature_grid(T_full, settings);
            coarse_scan = EELS_PearsonSparse_FullFit.stage2_scan_temperatures(coarse_grid, basis, kernels, A, layout, ...
                fit_mask, y_norm, zlp_norm, sqrtW, stage1_result, settings);
            if ~any(isfinite(coarse_scan.cost_grid_smooth))
                scan = EELS_PearsonSparse_FullFit.stage2_scan_temperatures(T_full, basis, kernels, A, layout, ...
                    fit_mask, y_norm, zlp_norm, sqrtW, stage1_result, settings);
                scan.scan_strategy = struct('mode', 'full_grid_after_coarse_failure', ...
                    'requested_points', numel(T_full), 'coarse_points', numel(coarse_grid), ...
                    'fine_points', 0, 'scanned_points', numel(scan.T_grid), ...
                    'coarse_best_T', NaN, 'fine_half_width_K', NaN);
                scan.temperature_profile.scan_strategy = scan.scan_strategy;
                return;
            end

            coarse_best_T = coarse_scan.best.T_fit;
            if ~isfinite(coarse_best_T)
                [~, idx] = min(coarse_scan.cost_grid_smooth);
                coarse_best_T = coarse_scan.T_grid(idx);
            end
            half_width = 160;
            if isfield(settings, 'stage2_adaptive_grid_fine_half_width_K')
                half_width = max(0, settings.stage2_adaptive_grid_fine_half_width_K);
            end
            fine_grid = T_full(abs(T_full - coarse_best_T) <= half_width);
            if isempty(fine_grid)
                fine_grid = coarse_best_T;
            end
            fine_grid = setdiff(fine_grid(:).', coarse_scan.T_grid(:).');
            if isempty(fine_grid)
                scan = coarse_scan;
            else
                fine_scan = EELS_PearsonSparse_FullFit.stage2_scan_temperatures(fine_grid, basis, kernels, A, layout, ...
                    fit_mask, y_norm, zlp_norm, sqrtW, stage1_result, settings);
                scan = EELS_PearsonSparse_FullFit.merge_stage2_temperature_scans(coarse_scan, fine_scan, settings);
            end
            scan.scan_strategy = struct('mode', 'coarse_to_fine_grid', ...
                'requested_points', numel(T_full), 'coarse_points', numel(coarse_grid), ...
                'fine_points', numel(fine_grid), 'scanned_points', numel(scan.T_grid), ...
                'coarse_best_T', coarse_best_T, 'fine_half_width_K', half_width);
            scan.temperature_profile.scan_strategy = scan.scan_strategy;
        end

        function coarse_grid = stage2_coarse_temperature_grid(T_full, settings)
            T_full = unique(T_full(:).');
            if numel(T_full) <= 2
                coarse_grid = T_full;
                return;
            end
            diffs = diff(T_full);
            base_step = median(diffs(isfinite(diffs) & diffs > 0), 'omitnan');
            if ~isfinite(base_step) || base_step <= 0
                base_step = 20;
            end
            coarse_step = 100;
            if isfield(settings, 'stage2_adaptive_grid_coarse_step_K')
                coarse_step = max(settings.stage2_adaptive_grid_coarse_step_K, base_step);
            end
            stride = max(1, round(coarse_step ./ max(base_step, eps)));
            coarse_grid = T_full(1:stride:end);
            coarse_grid = unique([coarse_grid(:).', T_full(1), T_full(end)]);
        end

        function merged = merge_stage2_temperature_scans(scan_a, scan_b, settings)
            T = unique([scan_a.T_grid(:); scan_b.T_grid(:)]).';
            nT = numel(T);
            merged = scan_a;
            merged.T_grid = T;
            names = {'cost_grid', 'cost_grid_smooth'};
            for ii = 1:numel(names)
                merged.(names{ii}) = NaN(1, nT);
            end
            prof_names = fieldnames(scan_a.temperature_profile);
            prof = struct();
            for ii = 1:numel(prof_names)
                name = prof_names{ii};
                val = scan_a.temperature_profile.(name);
                if isnumeric(val) && isvector(val) && numel(val) == numel(scan_a.T_grid)
                    prof.(name) = NaN(1, nT);
                else
                    prof.(name) = val;
                end
            end
            prof.T_grid = T;
            candidates = cell(1, nT);
            [merged, prof, candidates] = EELS_PearsonSparse_FullFit.fill_stage2_scan_into_merge( ...
                merged, prof, candidates, T, scan_a);
            [merged, prof, candidates] = EELS_PearsonSparse_FullFit.fill_stage2_scan_into_merge( ...
                merged, prof, candidates, T, scan_b);
            merged.cost_grid_smooth = EELS_PearsonSparse_FullFit.smooth_stage2_cost_profile(T, merged.cost_grid, settings);
            prof.cost_smooth = merged.cost_grid_smooth;
            finite_cost = isfinite(merged.cost_grid_smooth);
            best = struct('cost', Inf, 'T_fit', NaN);
            if any(finite_cost)
                idx = find(finite_cost);
                [~, local] = min(merged.cost_grid_smooth(finite_cost));
                best_idx = idx(local);
                if best_idx <= numel(candidates) && ~isempty(candidates{best_idx})
                    best = candidates{best_idx};
                    best.cost_smoothed = merged.cost_grid_smooth(best_idx);
                    best.cost = merged.cost_grid_smooth(best_idx);
                end
            end
            merged.best = best;
            merged.temperature_profile = prof;
            merged.candidates = candidates;
        end

        function [merged, prof, candidates] = fill_stage2_scan_into_merge(merged, prof, candidates, T, src)
            [lia, loc] = ismember(src.T_grid(:).', T);
            loc = loc(lia);
            src_idx = find(lia);
            if isempty(loc)
                return;
            end
            merged.cost_grid(loc) = src.cost_grid(src_idx);
            merged.cost_grid_smooth(loc) = src.cost_grid_smooth(src_idx);
            if isfield(src, 'candidates') && numel(src.candidates) >= max(src_idx)
                for kk = 1:numel(loc)
                    candidates{loc(kk)} = src.candidates{src_idx(kk)};
                end
            end
            prof_names = fieldnames(src.temperature_profile);
            for ii = 1:numel(prof_names)
                name = prof_names{ii};
                if strcmp(name, 'T_grid') || strcmp(name, 'scan_strategy')
                    continue;
                end
                val = src.temperature_profile.(name);
                if isnumeric(val) && isvector(val) && numel(val) == numel(src.T_grid) && isfield(prof, name)
                    prof.(name)(loc) = val(src_idx);
                end
            end
        end

        function scan = stage2_refine_temperature_local(scan, T_grid, basis, kernels, A, layout, fit_mask, y_norm, zlp_norm, sqrtW, stage1_result, settings)
            if ~isfield(settings, 'stage2_temperature_refine_enabled') || ~settings.stage2_temperature_refine_enabled
                scan.temperature_refine = struct('enabled', false, 'status', 'disabled');
                scan.temperature_profile.temperature_refine = scan.temperature_refine;
                return;
            end
            if ~(isfield(scan, 'best') && isfield(scan.best, 'T_fit') && isfinite(scan.best.T_fit))
                scan.temperature_refine = struct('enabled', true, 'status', 'no_grid_best');
                scan.temperature_profile.temperature_refine = scan.temperature_refine;
                return;
            end
            T_full = unique(T_grid(:).');
            T_full = T_full(isfinite(T_full) & T_full > 0);
            if numel(T_full) < 2
                scan.temperature_refine = struct('enabled', true, 'status', 'insufficient_grid');
                scan.temperature_profile.temperature_refine = scan.temperature_refine;
                return;
            end
            grid_step = median(diff(T_full), 'omitnan');
            if ~isfinite(grid_step) || grid_step <= 0
                grid_step = 20;
            end
            half_width = min(max(settings.stage2_temperature_refine_half_width_K, 0), max(grid_step, eps));
            if half_width <= 0
                scan.temperature_refine = struct('enabled', true, 'status', 'zero_half_width');
                scan.temperature_profile.temperature_refine = scan.temperature_refine;
                return;
            end
            tol_K = max(settings.stage2_temperature_refine_tol_K, 0.05);
            step_K = min(max(settings.stage2_temperature_refine_initial_step_K, tol_K), half_width);
            max_evals = max(0, round(settings.stage2_temperature_refine_max_evals));
            best = scan.best;
            best_T = best.T_fit;
            grid_best_T = best_T;
            best_cost = EELS_PearsonSparse_FullFit.stage2_candidate_raw_cost(best);
            lower_T = max(min(T_full), best_T - half_width);
            upper_T = min(max(T_full), best_T + half_width);
            evals = 0;
            refined_scans = {};
            refined_T = zeros(1, 0);
            refined_cost = zeros(1, 0);
            while evals < max_evals && step_K >= tol_K
                trial_T = unique([best_T - step_K, best_T + step_K]);
                trial_T = trial_T(isfinite(trial_T) & trial_T >= lower_T & trial_T <= upper_T);
                improved = false;
                for ii = 1:numel(trial_T)
                    if evals >= max_evals
                        break;
                    end
                    T_try = trial_T(ii);
                    if any(abs([scan.T_grid(:).', refined_T] - T_try) <= 0.25 .* tol_K)
                        continue;
                    end
                    local_scan = EELS_PearsonSparse_FullFit.stage2_scan_temperatures(T_try, basis, kernels, A, layout, ...
                        fit_mask, y_norm, zlp_norm, sqrtW, stage1_result, settings);
                    evals = evals + 1;
                    refined_scans{end + 1} = local_scan; %#ok<AGROW>
                    refined_T(end + 1) = T_try; %#ok<AGROW>
                    if isfield(local_scan, 'best') && isfield(local_scan.best, 'T_fit') && isfinite(local_scan.best.T_fit)
                        c_try = EELS_PearsonSparse_FullFit.stage2_candidate_raw_cost(local_scan.best);
                    else
                        c_try = Inf;
                    end
                    refined_cost(end + 1) = c_try; %#ok<AGROW>
                    if isfinite(c_try) && c_try < best_cost
                        best = local_scan.best;
                        best_T = best.T_fit;
                        best_cost = c_try;
                        improved = true;
                    end
                end
                if ~improved
                    step_K = 0.5 .* step_K;
                end
            end
            for ii = 1:numel(refined_scans)
                scan = EELS_PearsonSparse_FullFit.merge_stage2_temperature_scans(scan, refined_scans{ii}, settings);
            end
            if isfinite(best_T)
                best.cost_raw = best_cost;
                best.cost = best_cost;
                best.cost_smoothed = best_cost;
                scan.best = best;
            end
            info = struct();
            info.enabled = true;
            info.status = 'ok';
            info.method = 'local_1d_descent_after_grid';
            info.start_T = grid_best_T;
            if isfield(scan, 'scan_strategy') && isfield(scan.scan_strategy, 'coarse_best_T')
                info.coarse_best_T = scan.scan_strategy.coarse_best_T;
            else
                info.coarse_best_T = NaN;
            end
            info.best_T = best_T;
            info.best_raw_cost = best_cost;
            info.evaluations = evals;
            info.trial_T = refined_T;
            info.trial_raw_cost = refined_cost;
            info.initial_step_K = settings.stage2_temperature_refine_initial_step_K;
            info.final_step_K = step_K;
            info.half_width_K = half_width;
            info.tol_K = tol_K;
            scan.temperature_refine = info;
            scan.temperature_profile.temperature_refine = info;
        end

        function c = stage2_candidate_raw_cost(cand)
            c = Inf;
            if isstruct(cand)
                if isfield(cand, 'cost_raw') && isfinite(cand.cost_raw)
                    c = cand.cost_raw;
                elseif isfield(cand, 'cost') && isfinite(cand.cost)
                    c = cand.cost;
                end
            end
        end

        function scan = stage2_scan_temperatures(T_grid, basis, kernels, A, layout, fit_mask, y_norm, zlp_norm, sqrtW, stage1_result, settings)
            %#ok<INUSD>
            n_all = numel(basis.EM_all);
            hard_idx = false(n_all, 1);
            soft_idx = false(n_all, 1);
            if isfield(basis, 'hard_idx') && numel(basis.hard_idx) == n_all
                hard_idx = logical(basis.hard_idx(:));
            else
                hard_idx(:) = true;
            end
            if isfield(basis, 'soft_idx') && numel(basis.soft_idx) == n_all
                soft_idx = logical(basis.soft_idx(:));
            else
                soft_idx = ~hard_idx;
            end
            hard_idx = hard_idx & ~soft_idx;
            soft_idx = soft_idx & ~hard_idx;
            all_stokes = kernels.stokes(:, 1:n_all);
            all_anti = kernels.anti(:, 1:n_all);
            hard_stokes = all_stokes(:, hard_idx);
            hard_anti = all_anti(:, hard_idx);
            soft_stokes = all_stokes(:, soft_idx);
            soft_anti = all_anti(:, soft_idx);
            n_hard = nnz(hard_idx);
            n_soft = nnz(soft_idx);
            sqrtW_stage2 = sqrtW(:);
            fit_mask = logical(fit_mask(:));
            data_mask = isfinite(y_norm(:)) & isfinite(zlp_norm(:)) & ...
                isfinite(sqrtW_stage2(:)) & sqrtW_stage2(:) > 0;
            y_stage2 = y_norm(:) - zlp_norm(:);
            transition_mask = false(size(fit_mask));
            zlp_refine = EELS_PearsonSparse_FullFit.stage2_zlp_refine_basis( ...
                stage1_result, zlp_norm, fit_mask, transition_mask, settings);
            n_zlp_refine = size(zlp_refine.B, 2);
            base_ncoef = n_hard + 2 .* n_soft;
            hard_cols = 1:n_hard;
            soft_stokes_cols = n_hard + (1:n_soft);
            soft_anti_cols = n_hard + n_soft + (1:n_soft);
            T_grid = unique(T_grid(:).');
            nT = numel(T_grid);
            cost_grid = NaN(1, nT);
            data_cost_grid = NaN(1, nT);
            attention_cost_grid = NaN(1, nT);
            joint_rmse_grid = NaN(1, nT);
            regularization_cost_grid = NaN(1, nT);
            pm_smooth_cost_grid = NaN(1, nT);
            zlp_delta_rms_grid = NaN(1, nT);
            linearity_R2_grid = NaN(1, nT);
            effective_pair_count_grid = zeros(1, nT);
            candidates = cell(1, max(nT, 1));
            kb = 0.08617333262;
            pair_weights = ones(n_hard, 1);
            if isfield(basis, 'hard_pair_score_hard') && numel(basis.hard_pair_score_hard) == n_hard
                pair_weights = basis.hard_pair_score_hard(:);
            elseif isfield(basis, 'hard_pair_score') && numel(basis.hard_pair_score) == n_all
                pair_weights = basis.hard_pair_score(:);
                pair_weights = pair_weights(hard_idx);
            end
            pair_floor = 0;
            if isfield(settings, 'stage2_effective_pair_weight_floor')
                pair_floor = min(max(settings.stage2_effective_pair_weight_floor, 0), 1);
            end
            if n_hard > 0
                pair_weights = min(max(pair_weights(:), pair_floor), 1);
            else
                pair_weights = zeros(0, 1);
            end
            for tt = 1:nT
                T = T_grid(tt);
                if ~isfinite(T) || T <= 0
                    continue;
                end
                q_hard = exp(-basis.EM_hard(:).' ./ max(kb .* T, eps));
                q_soft = exp(-basis.EM_soft(:).' ./ max(kb .* T, eps));
                hard_be_signal = hard_stokes + hard_anti .* q_hard;
                soft_signal = [soft_stokes, soft_anti];
                A2 = [hard_be_signal, soft_signal, zlp_refine.B];
                [Reg, reg_rhs] = EELS_PearsonSparse_FullFit.stage2_regularization_rows(basis, 0, settings, q_hard, q_soft);
                if n_zlp_refine > 0
                    Reg = [Reg, zeros(size(Reg, 1), n_zlp_refine)];
                    [zlp_Reg, zlp_rhs] = EELS_PearsonSparse_FullFit.stage2_zlp_refine_regularization( ...
                        zlp_refine, base_ncoef);
                    Reg = [Reg; zlp_Reg];
                    reg_rhs = [reg_rhs; zlp_rhs];
                end
                coef_bounds = struct('lb', zeros(base_ncoef, 1), 'ub', Inf(base_ncoef, 1), ...
                    'enabled', false, 'active_lower_count', 0, 'active_upper_count', 0);
                if n_zlp_refine > 0
                    coef_bounds = EELS_PearsonSparse_FullFit.stage2_append_zlp_refine_bounds(coef_bounds, zlp_refine);
                end
                l1 = zeros(size(A2, 2), 1);
                [coef, inner_info] = EELS_PearsonSparse_FullFit.solve_inner_linear(A2, y_stage2, ...
                    sqrtW_stage2, Reg, reg_rhs, basis.energy_axis(:), settings, l1, coef_bounds);
                if isempty(coef) || any(~isfinite(coef))
                    continue;
                end
                zlp_delta_norm = zeros(size(zlp_norm(:)));
                if n_zlp_refine > 0
                    zlp_delta_norm = zlp_refine.B * coef(base_ncoef + (1:n_zlp_refine));
                end
                model_norm = zlp_norm(:) + A2 * coef;
                pm_s = zeros(n_all, 1);
                pm_a = zeros(n_all, 1);
                hard_pm_s = zeros(n_hard, 1);
                hard_pm_a = zeros(n_hard, 1);
                soft_pm_s = zeros(n_soft, 1);
                soft_pm_a = zeros(n_soft, 1);
                vib_hard_norm = zeros(size(y_norm(:)));
                vib_soft_norm = zeros(size(y_norm(:)));
                if n_hard > 0
                    hard_pm_s = coef(hard_cols);
                    hard_pm_a = hard_pm_s(:) .* q_hard(:);
                    pm_s(hard_idx) = hard_pm_s;
                    pm_a(hard_idx) = hard_pm_a;
                    vib_hard_norm = hard_be_signal * hard_pm_s;
                end
                if n_soft > 0
                    soft_pm_s = coef(soft_stokes_cols);
                    soft_pm_a = coef(soft_anti_cols);
                    pm_s(soft_idx) = soft_pm_s;
                    pm_a(soft_idx) = soft_pm_a;
                    vib_soft_norm = soft_stokes * soft_pm_s + soft_anti * soft_pm_a;
                end
                vib_norm = vib_hard_norm + vib_soft_norm;
                full_residual = model_norm(:) - y_norm(:);
                rw = full_residual(data_mask) .* sqrtW_stage2(data_mask);
                wnorm = mean(sqrtW_stage2(data_mask) .^ 2, 'omitnan');
                data_cost = mean(rw .^ 2, 'omitnan') ./ max(wnorm, eps);
                joint_rmse = sqrt(mean(full_residual(data_mask) .^ 2, 'omitnan'));
                attention_cost = data_cost;
                attention_mask = fit_mask(:) & data_mask(:);
                if nnz(attention_mask) >= 3
                    arw = full_residual(attention_mask) .* sqrtW_stage2(attention_mask);
                    awnorm = mean(sqrtW_stage2(attention_mask) .^ 2, 'omitnan');
                    attention_cost = mean(arw .^ 2, 'omitnan') ./ max(awnorm, eps);
                end
                reg_res = Reg * coef - reg_rhs;
                reg_cost = mean(reg_res .^ 2, 'omitnan');
                if ~isfinite(reg_cost), reg_cost = 0; end
                pm_smooth_cost = EELS_PearsonSparse_FullFit.stage2_pm_second_diff_cost(basis.EM_hard, hard_pm_s, settings);
                if n_soft > 0
                    soft_settings = settings;
                    if isfield(settings, 'stage2_soft_pm_smooth_weight') && isfinite(settings.stage2_soft_pm_smooth_weight)
                        soft_settings.stage2_pm_smooth_weight = settings.stage2_soft_pm_smooth_weight;
                    end
                    pm_smooth_cost = pm_smooth_cost + ...
                        EELS_PearsonSparse_FullFit.stage2_pm_second_diff_cost(basis.EM_soft, soft_pm_s, soft_settings) + ...
                        EELS_PearsonSparse_FullFit.stage2_pm_second_diff_cost(basis.EM_soft, soft_pm_a, soft_settings);
                end
                zlp_delta_rms = sqrt(mean(zlp_delta_norm(fit_mask) .^ 2, 'omitnan'));
                [linearity_R2, pair_count] = EELS_PearsonSparse_FullFit.stage2_hardbe_linearity_diagnostic( ...
                    basis.EM_hard(:), hard_pm_s(:), hard_pm_a(:), T, pair_weights(:));
                cost = settings.stage2_profile_cost_weight .* ( ...
                    settings.stage2_joint_rmse_cost_weight .* data_cost + ...
                    settings.stage2_attention_cost_weight .* attention_cost) + ...
                    settings.stage2_regularization_cost_weight .* reg_cost;
                cost_grid(tt) = cost;
                data_cost_grid(tt) = data_cost;
                attention_cost_grid(tt) = attention_cost;
                joint_rmse_grid(tt) = joint_rmse;
                regularization_cost_grid(tt) = reg_cost;
                pm_smooth_cost_grid(tt) = pm_smooth_cost;
                zlp_delta_rms_grid(tt) = zlp_delta_rms;
                linearity_R2_grid(tt) = linearity_R2;
                effective_pair_count_grid(tt) = pair_count;
                cand = struct();
                cand.cost = cost;
                cand.cost_raw = cost;
                cand.data_cost = data_cost;
                cand.attention_cost = attention_cost;
                cand.joint_rmse = joint_rmse;
                cand.regularization_cost = reg_cost;
                cand.pm_smooth_cost = pm_smooth_cost;
                cand.T_fit = T;
                cand.coef = coef;
                cand.inner_solver = inner_info;
                cand.model_norm = model_norm;
                cand.zlp_norm = zlp_norm(:) + zlp_delta_norm(:);
                cand.zlp_delta_norm = zlp_delta_norm(:);
                cand.zlp_refine = EELS_PearsonSparse_FullFit.stage2_zlp_refine_result( ...
                    zlp_refine, coef(base_ncoef + (1:n_zlp_refine)), zlp_delta_norm, fit_mask);
                cand.vib_norm = vib_norm;
                cand.vib_hard_norm = vib_hard_norm;
                cand.vib_soft_norm = vib_soft_norm;
                cand.transition_hard_leakage_score = 0;
                cand.highfreq_overshoot_score = 0;
                cand.neighbor_prior_score = 0;
                cand.PM_stokes = pm_s(:).';
                cand.PM_antistokes = pm_a(:).';
                cand.PM_stokes_hard = hard_pm_s(:).';
                cand.PM_antistokes_hard = hard_pm_a(:).';
                cand.PM_stokes_soft = soft_pm_s(:).';
                cand.PM_antistokes_soft = soft_pm_a(:).';
                cand.bg_norm = zeros(size(y_norm(:)));
                cand.linearity_R2 = linearity_R2;
                cand.effective_pair_count = pair_count;
                candidates{tt} = cand;
            end
            cost_grid_smooth = EELS_PearsonSparse_FullFit.smooth_stage2_cost_profile(T_grid, cost_grid, settings);
            best = struct('cost', Inf, 'T_fit', NaN);
            finite_cost = isfinite(cost_grid_smooth);
            if any(finite_cost)
                idx = find(finite_cost);
                [~, local] = min(cost_grid_smooth(finite_cost));
                best_idx = idx(local);
                best = candidates{best_idx};
                best.cost_smoothed = cost_grid_smooth(best_idx);
                best.cost = cost_grid_smooth(best_idx);
            end
            scan = struct();
            scan.T_grid = T_grid;
            scan.cost_grid = cost_grid;
            scan.cost_grid_smooth = cost_grid_smooth;
            scan.best = best;
            scan.candidates = candidates;
            scan.temperature_profile = struct( ...
                'T_grid', T_grid, ...
                'cost_raw', cost_grid, ...
                'cost_smooth', cost_grid_smooth, ...
                'joint_data_mse', data_cost_grid, ...
                'attention_data_mse', attention_cost_grid, ...
                'joint_RMSE', joint_rmse_grid, ...
                'regularization_cost', regularization_cost_grid, ...
                'PM_smooth_penalty', pm_smooth_cost_grid, ...
                'ZLP_delta_rms', zlp_delta_rms_grid, ...
                'linearity_R2', linearity_R2_grid, ...
                'effective_pair_count', effective_pair_count_grid);
        end

        function mse = stage2_weighted_data_mse(model_norm, y_norm, sqrtW, fit_mask)
            model_norm = model_norm(:);
            y_norm = y_norm(:);
            sqrtW = sqrtW(:);
            fit_mask = logical(fit_mask(:));
            good = fit_mask & isfinite(model_norm) & isfinite(y_norm) & isfinite(sqrtW) & sqrtW > 0;
            if nnz(good) < 1
                mse = NaN;
                return;
            end
            r = (model_norm(good) - y_norm(good)) .* sqrtW(good);
            mse = mean(r .^ 2, 'omitnan');
        end

        function sqrtW = apply_stage2_far_tail_downweight(ene, sqrtW, settings, attention)
            if ~isfield(settings, 'stage2_far_tail_downweight_enabled') || ...
                    ~settings.stage2_far_tail_downweight_enabled
                return;
            end
            if ~isfield(settings, 'stage2_far_tail_weight_scale')
                return;
            end
            scale = min(max(settings.stage2_far_tail_weight_scale, 0.02), 1);
            if scale >= 0.999
                return;
            end
            ene = ene(:);
            sqrtW = sqrtW(:);
            if numel(ene) ~= numel(sqrtW)
                return;
            end
            if nargin >= 4 && isstruct(attention) && ...
                    isfield(settings, 'stage2_attention_tail_enabled') && settings.stage2_attention_tail_enabled
                n = numel(ene);
                score = zeros(n, 1);
                if isfield(attention, 'signal_score') && numel(attention.signal_score) == n
                    score = min(max(attention.signal_score(:), 0), 1);
                elseif isfield(attention, 'feature_score') && numel(attention.feature_score) == n
                    score = min(max(attention.feature_score(:), 0), 1);
                end
                threshold = min(max(settings.stage2_attention_tail_signal_threshold, 0), 1);
                soft_width = max(settings.stage2_attention_tail_soft_width, eps);
                release = 1 ./ (1 + exp(-(score - threshold) ./ soft_width));
                min_abs = settings.stage2_attention_tail_min_abs_meV;
                if isfield(settings, 'highfreq_start_meV')
                    min_abs = max(min_abs, settings.highfreq_start_meV);
                end
                eabs = abs(ene);
                tail_candidate = eabs >= min_abs;
                if isfield(attention, 'tail_mask') && numel(attention.tail_mask) == n
                    tail_candidate = tail_candidate | logical(attention.tail_mask(:));
                end
                profile = ones(n, 1);
                profile(tail_candidate) = scale + (1 - scale) .* release(tail_candidate);
                profile(~isfinite(profile)) = 1;
                sqrtW = sqrtW .* sqrt(min(max(profile, scale), 1));
                return;
            end
            start_meV = settings.stage2_far_tail_downweight_start_meV;
            width_meV = max(settings.stage2_far_tail_downweight_width_meV, eps);
            eabs = abs(ene);
            profile = scale + (1 - scale) ./ (1 + exp((eabs - start_meV) ./ width_meV));
            profile(~isfinite(profile)) = 1;
            sqrtW = sqrtW .* sqrt(min(max(profile, scale), 1));
        end

        function [Reg, rhs] = stage2_regularization_rows(basis, nbg, settings, q_hard, q_soft)
            %#ok<INUSD>
            n_hard = numel(basis.EM_hard);
            n_soft = numel(basis.EM_soft);
            ncoef = n_hard + 2 .* n_soft + nbg;
            Reg = zeros(0, ncoef);
            rhs = zeros(0, 1);
            hard_cols = 1:n_hard;
            soft_stokes_cols = n_hard + (1:n_soft);
            soft_anti_cols = n_hard + n_soft + (1:n_soft);
            if settings.stage2_pm_smooth_weight > 0 && n_hard >= 3
                D2 = EELS_PearsonSparse_FullFit.local_diff_matrix(basis.EM_hard(:), 2, settings.stage2_pm_smooth_max_gap_meV);
                if ~isempty(D2)
                    smooth = zeros(size(D2, 1), ncoef);
                    smooth(:, hard_cols) = D2;
                    Reg = [Reg; sqrt(settings.stage2_pm_smooth_weight) .* smooth]; %#ok<AGROW>
                    rhs = [rhs; zeros(size(D2, 1), 1)]; %#ok<AGROW>
                end
            end
            soft_smooth_weight = settings.stage2_pm_smooth_weight;
            if isfield(settings, 'stage2_soft_pm_smooth_weight') && isfinite(settings.stage2_soft_pm_smooth_weight)
                soft_smooth_weight = max(settings.stage2_soft_pm_smooth_weight, 0);
            end
            if soft_smooth_weight > 0 && n_soft >= 3
                D2 = EELS_PearsonSparse_FullFit.local_diff_matrix(basis.EM_soft(:), 2, settings.stage2_pm_smooth_max_gap_meV);
                if ~isempty(D2)
                    smooth_s = zeros(size(D2, 1), ncoef);
                    smooth_a = zeros(size(D2, 1), ncoef);
                    smooth_s(:, soft_stokes_cols) = D2;
                    smooth_a(:, soft_anti_cols) = D2;
                    Reg = [Reg; sqrt(soft_smooth_weight) .* smooth_s; sqrt(soft_smooth_weight) .* smooth_a]; %#ok<AGROW>
                    rhs = [rhs; zeros(2 .* size(D2, 1), 1)]; %#ok<AGROW>
                end
            end
            if n_soft <= 0
                return;
            end
            soft_l2_weight = 0;
            if isfield(settings, 'stage2_soft_pm_l2_weight')
                soft_l2_weight = max(settings.stage2_soft_pm_l2_weight, 0);
            end
            if soft_l2_weight > 0
                l2 = zeros(2 .* n_soft, ncoef);
                l2(1:n_soft, soft_stokes_cols) = eye(n_soft);
                l2(n_soft + (1:n_soft), soft_anti_cols) = eye(n_soft);
                Reg = [Reg; sqrt(soft_l2_weight) .* l2]; %#ok<AGROW>
                rhs = [rhs; zeros(2 .* n_soft, 1)]; %#ok<AGROW>
            end
            soft_prior_weight = 0;
            if isfield(settings, 'stage2_soft_pm_stage1_prior_weight')
                soft_prior_weight = max(settings.stage2_soft_pm_stage1_prior_weight, 0);
            end
            if soft_prior_weight > 0
                prior_s = zeros(n_soft, 1);
                prior_a = zeros(n_soft, 1);
                if isfield(basis, 'stage1_PM_stokes_soft') && numel(basis.stage1_PM_stokes_soft) == n_soft
                    prior_s = basis.stage1_PM_stokes_soft(:);
                end
                if isfield(basis, 'stage1_PM_antistokes_soft') && numel(basis.stage1_PM_antistokes_soft) == n_soft
                    prior_a = basis.stage1_PM_antistokes_soft(:);
                end
                prior_s(~isfinite(prior_s) | prior_s < 0) = 0;
                prior_a(~isfinite(prior_a) | prior_a < 0) = 0;
                P = zeros(2 .* n_soft, ncoef);
                P(1:n_soft, soft_stokes_cols) = eye(n_soft);
                P(n_soft + (1:n_soft), soft_anti_cols) = eye(n_soft);
                Reg = [Reg; sqrt(soft_prior_weight) .* P]; %#ok<AGROW>
                rhs = [rhs; sqrt(soft_prior_weight) .* [prior_s; prior_a]]; %#ok<AGROW>
            end
            soft_be_weight = 0;
            if isfield(settings, 'stage2_soft_be_weight')
                soft_be_weight = max(settings.stage2_soft_be_weight, 0);
            end
            if soft_be_weight > 0
                if nargin < 5 || isempty(q_soft) || numel(q_soft) ~= n_soft
                    q_soft = ones(1, n_soft);
                end
                q_soft = q_soft(:).';
                q_soft(~isfinite(q_soft) | q_soft < 0) = 0;
                be = zeros(n_soft, ncoef);
                for ii = 1:n_soft
                    be(ii, soft_stokes_cols(ii)) = -q_soft(ii);
                    be(ii, soft_anti_cols(ii)) = 1;
                end
                Reg = [Reg; sqrt(soft_be_weight) .* be]; %#ok<AGROW>
                rhs = [rhs; zeros(n_soft, 1)]; %#ok<AGROW>
            end
        end



        function unc = estimate_stage2_uncertainty(T_grid, cost_grid)
            valley = EELS_PearsonSparse_FullFit.stage2_cost_valley(T_grid, cost_grid);
            if ~valley.available
                unc = NaN;
                return;
            end
            unc = valley.uncertainty_K;
        end

        function cost_smooth = smooth_stage2_cost_profile(T_grid, cost_grid, settings)
            T_grid = T_grid(:).';
            cost_grid = cost_grid(:).';
            cost_smooth = cost_grid;
            good = isfinite(T_grid) & isfinite(cost_grid);
            if nnz(good) < 3
                return;
            end
            [~, order] = sort(T_grid(good));
            good_idx = find(good);
            good_idx = good_idx(order);
            c = cost_grid(good_idx);
            points = max(3, round(settings.stage2_cost_smooth_points));
            if mod(points, 2) == 0, points = points + 1; end
            points = min(points, numel(c));
            if points < 3
                return;
            end
            t = T_grid(good_idx).';
            c = c(:);
            cs = c;
            for ii = 1:numel(c)
                distance = abs(t - t(ii));
                [sorted_distance, nearest] = sort(distance, 'ascend');
                support_count = min(points, numel(c));
                bandwidth_idx = min(support_count + 1, numel(c));
                bandwidth = sorted_distance(bandwidth_idx);
                if ~isfinite(bandwidth) || bandwidth <= 0
                    continue;
                end
                use = nearest(1:support_count);
                x = t(use) - t(ii);
                u = min(abs(x) ./ bandwidth, 1);
                weight = (1 - u .^ 3) .^ 3;
                valid = isfinite(x) & isfinite(c(use)) & isfinite(weight) & weight > 0;
                if nnz(valid) < 3
                    continue;
                end
                x = x(valid);
                y = c(use(valid));
                weight = weight(valid);
                X = [ones(size(x)), x, x .^ 2];
                normal = X.' * (weight .* X);
                rhs = X.' * (weight .* y);
                if rcond(normal) < 1e-12
                    X = X(:, 1:2);
                    normal = X.' * (weight .* X);
                    rhs = X.' * (weight .* y);
                end
                if rcond(normal) >= 1e-12
                    beta = normal \ rhs;
                    local_value = beta(1);
                    if isfinite(local_value)
                        cs(ii) = min(max(local_value, min(y)), max(y));
                    end
                end
            end
            cs(~isfinite(cs)) = c(~isfinite(cs));
            cost_smooth(good_idx) = cs(:).';
        end

        function valley = stage2_cost_valley(T_grid, cost_grid)
            valley = struct();
            valley.available = false;
            valley.T_grid = [];
            valley.cost_grid = [];
            valley.relative_cost = [];
            valley.threshold = NaN;
            valley.relative_threshold = NaN;
            valley.in_valley = [];
            valley.best_T = NaN;
            valley.best_cost = NaN;
            valley.best_idx = NaN;
            valley.lower_K = NaN;
            valley.upper_K = NaN;
            valley.uncertainty_K = NaN;
            good = isfinite(T_grid) & isfinite(cost_grid);
            T_grid = T_grid(good);
            cost_grid = cost_grid(good);
            if numel(T_grid) < 3
                return;
            end
            [T_grid, order] = sort(T_grid(:).');
            cost_grid = cost_grid(order);
            [best_cost, best_idx] = min(cost_grid);
            threshold = best_cost + 0.10 .* max(abs(best_cost), eps);
            near = cost_grid <= threshold;
            scale = max(abs(best_cost), eps);
            valley.available = true;
            valley.T_grid = T_grid;
            valley.cost_grid = cost_grid;
            valley.relative_cost = (cost_grid - best_cost) ./ scale;
            valley.threshold = threshold;
            valley.relative_threshold = (threshold - best_cost) ./ scale;
            valley.in_valley = near;
            valley.best_T = T_grid(best_idx);
            valley.best_cost = best_cost;
            valley.best_idx = best_idx;
            if nnz(near) >= 2
                valley.lower_K = min(T_grid(near));
                valley.upper_K = max(T_grid(near));
                valley.uncertainty_K = 0.5 .* (valley.upper_K - valley.lower_K);
            else
                sorted_T = sort(T_grid);
                step = median(diff(sorted_T), 'omitnan');
                if ~isfinite(step) || step <= 0, step = 20; end
                valley.lower_K = max(min(T_grid), T_grid(best_idx) - step);
                valley.upper_K = min(max(T_grid), T_grid(best_idx) + step);
                valley.uncertainty_K = step;
            end
            if best_idx == 1 || best_idx == numel(T_grid)
                valley.uncertainty_K = max(valley.uncertainty_K, 0.25 .* (max(T_grid) - min(T_grid)));
            end
        end



        function line = best_ratio_origin_line(EM, log_stokes_over_anti, weight, settings)
            kb = 0.08617333262;
            EM = EM(:);
            log_stokes_over_anti = log_stokes_over_anti(:);
            weight = weight(:);
            n_input = numel(EM);
            line = struct('available', false, 'method', 'best_contiguous_origin_line', ...
                'T_K', NaN, 'slope_per_meV', NaN, 'R2_origin', NaN, ...
                'R2_origin_weighted', NaN, 'R2_origin_unweighted', NaN, ...
                'score', 0, 'num_points', 0, 'num_points_total', n_input, ...
                'weight_sum', 0, 'span_meV', NaN, 'EM_range_meV', [NaN, NaN], ...
                'used_mask', false(n_input, 1), 'residual_rms', NaN, ...
                'right_edge_residual_median', NaN, 'right_edge_residual_max', NaN, ...
                'excluded_points', n_input);
            if n_input < 2 || numel(log_stokes_over_anti) ~= n_input || numel(weight) ~= n_input
                return;
            end
            good = isfinite(EM) & isfinite(log_stokes_over_anti) & isfinite(weight) & ...
                EM > 0 & log_stokes_over_anti > 0 & weight > 0;
            if nnz(good) < 2
                return;
            end
            orig = find(good);
            x = EM(good);
            y = log_stokes_over_anti(good);
            w = weight(good);
            [x, order] = sort(x);
            y = y(order);
            w = w(order);
            orig = orig(order);
            w = w ./ max(max(w), eps);

            min_points = max(2, round(settings.stage2_ratio_best_line_min_points));
            min_span = max(0, settings.stage2_ratio_best_line_min_span_meV);
            ref_span = max(settings.stage2_ratio_best_line_ref_span_meV, eps);
            edge_points = max(1, round(settings.stage2_ratio_best_line_edge_points));
            right_edge_median_max = settings.stage2_ratio_best_line_right_edge_median_residual_max;
            right_edge_max_max = settings.stage2_ratio_best_line_right_edge_max_residual_max;
            if numel(x) < min_points
                min_points = numel(x);
                min_span = 0;
            end
            T_bounds = settings.T_soft_bounds(:).';
            if numel(T_bounds) < 2
                T_bounds = [0, Inf];
            end
            best_score = -Inf;
            best = struct();
            n = numel(x);
            for ii = 1:n
                for jj = (ii + min_points - 1):n
                    span = x(jj) - x(ii);
                    if span < min_span
                        continue;
                    end
                    idx = ii:jj;
                    fit = EELS_PearsonSparse_FullFit.ratio_origin_line_fit(x(idx), y(idx), w(idx), kb);
                    if ~fit.available || fit.T_K < T_bounds(1) || fit.T_K > T_bounds(2)
                        continue;
                    end
                    resid = y(idx) - fit.slope_per_meV .* x(idx);
                    edge_n = min(edge_points, numel(idx));
                    right_resid = resid((end - edge_n + 1):end);
                    right_resid_median = median(right_resid, 'omitnan');
                    right_resid_max = max(right_resid, [], 'omitnan');
                    if isfinite(right_edge_median_max) && isfinite(right_resid_median) && ...
                            right_resid_median > right_edge_median_max
                        continue;
                    end
                    if isfinite(right_edge_max_max) && isfinite(right_resid_max) && ...
                            right_resid_max > right_edge_max_max
                        continue;
                    end
                    n_score = min(1, numel(idx) ./ max(2 .* min_points, 1));
                    span_score = min(1, span ./ ref_span);
                    weight_ref = 0.35;
                    if isfield(settings, 'stage2_ratio_cost_weight_sum_ref')
                        weight_ref = settings.stage2_ratio_cost_weight_sum_ref;
                    end
                    weight_score = min(1, fit.weight_sum ./ max(weight_ref, eps));
                    r2_score = min(max(fit.R2_origin, 0), 1);
                    score = (r2_score .^ 4) .* (0.55 + 0.45 .* n_score) .* ...
                        (0.55 + 0.45 .* span_score) .* (0.65 + 0.35 .* weight_score);
                    if score > best_score
                        best_score = score;
                        best.fit = fit;
                        best.idx = idx;
                        best.span = span;
                        best.right_resid_median = right_resid_median;
                        best.right_resid_max = right_resid_max;
                    end
                end
            end
            if best_score < 0
                if ~settings.stage2_ratio_best_line_allow_all_fallback
                    return;
                end
                idx = 1:n;
                fit = EELS_PearsonSparse_FullFit.ratio_origin_line_fit(x, y, w, kb);
                if ~fit.available
                    return;
                end
                best_score = max(fit.R2_origin, 0);
                best.fit = fit;
                best.idx = idx;
                best.span = x(end) - x(1);
                resid = y(idx) - fit.slope_per_meV .* x(idx);
                edge_n = min(edge_points, numel(idx));
                right_resid = resid((end - edge_n + 1):end);
                best.right_resid_median = median(right_resid, 'omitnan');
                best.right_resid_max = max(right_resid, [], 'omitnan');
            end
            used_orig = orig(best.idx);
            line.available = true;
            line.T_K = best.fit.T_K;
            line.slope_per_meV = best.fit.slope_per_meV;
            line.R2_origin = best.fit.R2_origin;
            line.R2_origin_weighted = best.fit.R2_origin_weighted;
            line.R2_origin_unweighted = best.fit.R2_origin_unweighted;
            line.score = best_score;
            line.num_points = numel(best.idx);
            line.weight_sum = best.fit.weight_sum;
            line.span_meV = best.span;
            line.EM_range_meV = [min(x(best.idx)), max(x(best.idx))];
            line.used_mask(used_orig) = true;
            line.residual_rms = best.fit.residual_rms;
            line.right_edge_residual_median = best.right_resid_median;
            line.right_edge_residual_max = best.right_resid_max;
            line.excluded_points = n_input - nnz(line.used_mask);
        end

        function fit = ratio_origin_line_fit(x, y, w, kb)
            x = x(:);
            y = y(:);
            w = w(:);
            fit = struct('available', false, 'slope_per_meV', NaN, 'T_K', NaN, ...
                'R2_origin', NaN, 'R2_origin_weighted', NaN, 'R2_origin_unweighted', NaN, ...
                'weight_sum', 0, 'residual_rms', NaN);
            good = isfinite(x) & isfinite(y) & isfinite(w) & x > 0 & y > 0 & w > 0;
            x = x(good);
            y = y(good);
            w = w(good);
            if numel(x) < 2
                return;
            end
            slope = sum(w .* x .* y, 'omitnan') ./ max(sum(w .* x .^ 2, 'omitnan'), eps);
            if ~isfinite(slope) || slope <= 0
                return;
            end
            yhat = slope .* x;
            resid = y - yhat;
            sse = sum(w .* resid .^ 2, 'omitnan');
            sst0 = sum(w .* y .^ 2, 'omitnan');
            sse_unweighted = sum(resid .^ 2, 'omitnan');
            sst0_unweighted = sum(y .^ 2, 'omitnan');
            r2_weighted = 1 - sse ./ max(sst0, eps);
            r2_unweighted = 1 - sse_unweighted ./ max(sst0_unweighted, eps);
            fit.available = true;
            fit.slope_per_meV = slope;
            fit.T_K = 1 ./ max(kb .* slope, eps);
            fit.R2_origin_weighted = r2_weighted;
            fit.R2_origin_unweighted = r2_unweighted;
            fit.R2_origin = min(r2_weighted, r2_unweighted);
            fit.weight_sum = sum(w, 'omitnan');
            fit.residual_rms = sqrt(sum(w .* resid .^ 2, 'omitnan') ./ max(sum(w, 'omitnan'), eps));
        end

        function ratio_est = estimate_attention_ratio_temperature(stage1_result, attention, settings, basis)
            if nargin >= 4 && isfield(basis, 'EM_hard') && ~isempty(basis.EM_hard)
                EM = basis.EM_hard(:);
                pm_s = basis.stage1_PM_stokes_hard(:);
                pm_a = basis.stage1_PM_antistokes_hard(:);
            else
                EM = stage1_result.EM_final(:);
                pm_s = stage1_result.PM_stokes(:);
                pm_a = stage1_result.PM_antistokes(:);
            end
            kb = 0.08617333262;
            win = sort(abs(settings.stage2_window_meV(:).'));
            if numel(win) < 2, win = [40, 220]; end
            att_pos = interp1(attention.energy(:), attention.feature_score(:), EM, 'linear', 0);
            att_neg = interp1(attention.energy(:), attention.feature_score(:), -EM, 'linear', 0);
            pair_attention = sqrt(max(att_pos, 0) .* max(att_neg, 0));
            pm_strength = pm_s ./ max(max(pm_s), eps);
            weight = pair_attention(:) .* pm_strength(:);
            ratio = max(pm_a, eps) ./ max(pm_s, eps);
            log_stokes_over_anti = log(max(pm_s, eps) ./ max(pm_a, eps));
            T_pair = -EM ./ (kb .* log(ratio));
            valid = isfinite(T_pair) & T_pair > 0 & ratio > 0 & ratio < 0.98 & ...
                EM >= win(1) & EM <= win(2) & weight > 0.02;

            ratio_est = struct();
            ratio_est.reliable = false;
            ratio_est.T_ratio = NaN;
            ratio_est.spread_K = NaN;
            ratio_est.num_pairs = nnz(valid);
            ratio_est.weight_sum = sum(weight(valid), 'omitnan');
            ratio_est.reliability_score = 0;
            ratio_est.EM = zeros(0, 1);
            ratio_est.T_pair = zeros(0, 1);
            ratio_est.weight = zeros(0, 1);
            ratio_est.log_stokes_over_anti = zeros(0, 1);
            ratio_est.method = 'weighted_pair_median';
            ratio_est.best_line = struct();
            ratio_est.selected_mask = false(0, 1);
            ratio_est.R2_origin = NaN;
            ratio_est.used_EM_range_meV = [NaN, NaN];
            ratio_est.num_pairs_selected = 0;
            if nnz(valid) < settings.stage2_ratio_min_pairs || ratio_est.weight_sum <= 0
                return;
            end

            Tv = T_pair(valid);
            wv = weight(valid);
            EMv = EM(valid);
            logv = log_stokes_over_anti(valid);
            ratio_est.EM = EMv;
            ratio_est.T_pair = Tv(:);
            ratio_est.weight = wv(:);
            ratio_est.log_stokes_over_anti = logv;
            q25 = EELS_PearsonSparse_FullFit.weighted_quantile(Tv, wv, 0.25);
            q50 = EELS_PearsonSparse_FullFit.weighted_quantile(Tv, wv, 0.50);
            q75 = EELS_PearsonSparse_FullFit.weighted_quantile(Tv, wv, 0.75);
            ratio_est.T_ratio = q50;
            ratio_est.q25 = q25;
            ratio_est.q75 = q75;
            ratio_est.spread_K = q75 - q25;
            ratio_est.selected_mask = true(size(Tv(:)));
            if settings.stage2_ratio_best_line_enabled
                best_line = EELS_PearsonSparse_FullFit.best_ratio_origin_line(EMv, logv, wv, settings);
                if best_line.available
                    ratio_est.method = 'best_contiguous_origin_line';
                    ratio_est.best_line = best_line;
                    ratio_est.T_ratio = best_line.T_K;
                    ratio_est.R2_origin = best_line.R2_origin;
                    ratio_est.used_EM_range_meV = best_line.EM_range_meV;
                    ratio_est.num_pairs_selected = best_line.num_points;
                    ratio_est.selected_mask = best_line.used_mask(:);
                    Tv_sel = Tv(ratio_est.selected_mask);
                    wv_sel = wv(ratio_est.selected_mask);
                    if numel(Tv_sel) >= 2 && sum(wv_sel, 'omitnan') > 0
                        ratio_est.q25 = EELS_PearsonSparse_FullFit.weighted_quantile(Tv_sel, wv_sel, 0.25);
                        ratio_est.q75 = EELS_PearsonSparse_FullFit.weighted_quantile(Tv_sel, wv_sel, 0.75);
                        ratio_est.spread_K = ratio_est.q75 - ratio_est.q25;
                    end
                end
            end
            line_reliable = strcmp(ratio_est.method, 'best_contiguous_origin_line') && ...
                isfinite(ratio_est.R2_origin) && ratio_est.R2_origin >= settings.stage2_ratio_best_line_min_R2 && ...
                ratio_est.num_pairs_selected >= settings.stage2_ratio_best_line_min_points;
            if settings.stage2_ratio_best_line_enabled
                ratio_est.reliable = isfinite(ratio_est.T_ratio) && line_reliable;
            else
                ratio_est.reliable = isfinite(ratio_est.T_ratio) && ...
                    isfinite(ratio_est.spread_K) && ratio_est.spread_K <= settings.stage2_ratio_reliable_spread_K;
            end
            pair_score = min(1, ratio_est.num_pairs ./ max(2 .* settings.stage2_ratio_min_pairs, 1));
            weight_score = min(1, ratio_est.weight_sum ./ max(settings.stage2_ratio_cost_weight_sum_ref, eps));
            spread_score = 1 ./ (1 + (max(ratio_est.spread_K, 0) ./ max(settings.stage2_ratio_reliable_spread_K, eps)).^2);
            line_score = 0;
            if strcmp(ratio_est.method, 'best_contiguous_origin_line') && isfinite(ratio_est.R2_origin)
                selected_score = min(1, ratio_est.num_pairs_selected ./ max(2 .* settings.stage2_ratio_best_line_min_points, 1));
                line_score = min(max(ratio_est.R2_origin, 0), 1) .* (0.65 + 0.35 .* selected_score);
            end
            if settings.stage2_ratio_best_line_enabled
                if ~line_reliable
                    line_score = 0;
                end
                ratio_est.reliability_score = min(1, max(0, pair_score .* weight_score .* line_score));
            else
                ratio_est.reliability_score = min(1, max(0, pair_score .* weight_score .* spread_score));
            end
        end

        function analysis = stage1_pm_ratio_reliability(stage1_result, settings)
            if nargin < 2, settings = struct(); end
            settings = EELS_PearsonSparse_FullFit.defaults(settings);
            analysis = struct('source', 'stage1_pm_ratio', ...
                'standard', 'attention_weighted_stage1_PM_anti_over_stokes_ratio', ...
                'available', false, 'reliable', false, 'T_ratio', NaN, ...
                'uncertainty_K', NaN, 'uncertainty_method', 'stage1_pm_ratio_IQR_half_width_diagnostic', ...
                'spread_K', NaN, 'num_pairs', 0, 'weight_sum', 0, ...
                'method', 'unknown', 'ratio_fit_R2', NaN, 'num_pairs_selected', 0, ...
                'used_EM_range_meV', [NaN, NaN], ...
                'fit_quality_score', NaN, 'ratio_score', 0, 'reliability_score', 0, ...
                'reason', 'missing_stage1_PM_fields');
            needed = {'ene_corrected', 'EM_final', 'PM_stokes', 'PM_antistokes', 'stage1', 'raw_data_norm'};
            for ii = 1:numel(needed)
                if ~isfield(stage1_result, needed{ii})
                    return;
                end
            end
            try
                attention = EELS_PearsonSparse_FullFit.build_stage2_attention(stage1_result, settings);
                ratio_est = EELS_PearsonSparse_FullFit.estimate_attention_ratio_temperature( ...
                    stage1_result, attention, settings, struct());
            catch ME
                analysis.reason = ['analysis_failed_', char(string(ME.identifier))];
                return;
            end
            analysis.ratio_est = ratio_est;
            if isfield(ratio_est, 'T_ratio'), analysis.T_ratio = ratio_est.T_ratio; end
            if isfield(ratio_est, 'spread_K'), analysis.spread_K = ratio_est.spread_K; end
            if isfield(ratio_est, 'num_pairs'), analysis.num_pairs = ratio_est.num_pairs; end
            if isfield(ratio_est, 'weight_sum'), analysis.weight_sum = ratio_est.weight_sum; end
            if isfield(ratio_est, 'reliable'), analysis.reliable = ratio_est.reliable; end
            if isfield(ratio_est, 'reliability_score'), analysis.ratio_score = ratio_est.reliability_score; end
            if isfield(ratio_est, 'method'), analysis.method = ratio_est.method; end
            if isfield(ratio_est, 'R2_origin'), analysis.ratio_fit_R2 = ratio_est.R2_origin; end
            if isfield(ratio_est, 'num_pairs_selected'), analysis.num_pairs_selected = ratio_est.num_pairs_selected; end
            if isfield(ratio_est, 'used_EM_range_meV'), analysis.used_EM_range_meV = ratio_est.used_EM_range_meV; end
            if isfinite(analysis.spread_K)
                analysis.uncertainty_K = 0.5 .* analysis.spread_K;
            end
            analysis.available = isfinite(analysis.T_ratio) && analysis.num_pairs > 0 && analysis.weight_sum > 0;
            if ~analysis.available
                analysis.reason = 'insufficient_stage1_PM_ratio_pairs';
                return;
            end
            fit_score = NaN;
            if isfield(stage1_result, 'metrics') && isfield(stage1_result.metrics, 'R2') && isfinite(stage1_result.metrics.R2)
                fit_score = min(1, max(0, (stage1_result.metrics.R2 - 0.995) ./ 0.0045));
                if isfield(stage1_result.metrics, 'weighted_NRMSE') && isfinite(stage1_result.metrics.weighted_NRMSE)
                    wn_score = exp(-max(stage1_result.metrics.weighted_NRMSE, 0) ./ 0.010);
                    fit_score = 0.65 .* fit_score + 0.35 .* min(1, max(0, wn_score));
                end
            end
            analysis.fit_quality_score = fit_score;
            score = min(1, max(0, analysis.ratio_score));
            if isfinite(fit_score)
                score = score .* (0.75 + 0.25 .* fit_score);
            end
            analysis.reliability_score = min(1, max(0, score));
            analysis.reason = 'ok';
        end

        function source = reliability_source_from_result(res)
            source = 'unknown';
            if ~isstruct(res)
                return;
            end
            if isfield(res, 'reliability_source') && ~isempty(res.reliability_source)
                source = char(string(res.reliability_source));
                return;
            end
            if isfield(res, 'stage2') && isstruct(res.stage2) && ...
                    isfield(res.stage2, 'reliability_source') && ~isempty(res.stage2.reliability_source)
                source = char(string(res.stage2.reliability_source));
                return;
            end
            if isfield(res, 'stage1_pm_ratio_reliability') && isstruct(res.stage1_pm_ratio_reliability)
                source = 'stage1_pm_ratio';
            end
        end

        function [T_ratio, score, num_pairs, ratio_fit_R2, num_pairs_selected, selected_EM_lo, selected_EM_hi] = stage1_ratio_summary_from_result(res)
            T_ratio = NaN;
            score = NaN;
            num_pairs = NaN;
            ratio_fit_R2 = NaN;
            num_pairs_selected = NaN;
            selected_EM_lo = NaN;
            selected_EM_hi = NaN;
            ratio = [];
            if isstruct(res) && isfield(res, 'stage1_pm_ratio_reliability') && isstruct(res.stage1_pm_ratio_reliability)
                ratio = res.stage1_pm_ratio_reliability;
            elseif isstruct(res) && isfield(res, 'stage2') && isstruct(res.stage2) && ...
                    isfield(res.stage2, 'stage1_pm_ratio_reliability') && isstruct(res.stage2.stage1_pm_ratio_reliability)
                ratio = res.stage2.stage1_pm_ratio_reliability;
            end
            if isempty(ratio)
                return;
            end
            if isfield(ratio, 'T_ratio'), T_ratio = ratio.T_ratio; end
            if isfield(ratio, 'reliability_score'), score = ratio.reliability_score; end
            if isfield(ratio, 'num_pairs'), num_pairs = ratio.num_pairs; end
            if isfield(ratio, 'ratio_fit_R2'), ratio_fit_R2 = ratio.ratio_fit_R2; end
            if isfield(ratio, 'num_pairs_selected'), num_pairs_selected = ratio.num_pairs_selected; end
            if isfield(ratio, 'used_EM_range_meV') && numel(ratio.used_EM_range_meV) >= 2
                selected_EM_lo = ratio.used_EM_range_meV(1);
                selected_EM_hi = ratio.used_EM_range_meV(2);
            end
        end

        function [max_delta, fit_rms] = stage2_zlp_refine_summary_from_result(res)
            max_delta = NaN;
            fit_rms = NaN;
            if ~(isstruct(res) && isfield(res, 'stage2') && isstruct(res.stage2) && ...
                    isfield(res.stage2, 'zlp_refine') && isstruct(res.stage2.zlp_refine))
                return;
            end
            zr = res.stage2.zlp_refine;
            if isfield(zr, 'max_abs_delta_norm'), max_delta = zr.max_abs_delta_norm; end
            if isfield(zr, 'fit_rms_delta_norm'), fit_rms = zr.fit_rms_delta_norm; end
        end

        function score = stage2_reliability(stage2_result, stage1_result)
            if ~isfield(stage2_result, 'T_uncertainty') || ~isfinite(stage2_result.T_uncertainty)
                score = 0;
                return;
            end
            attention_fraction = 0;
            if isfield(stage2_result, 'attention') && isfield(stage2_result.attention, 'selected_fraction')
                attention_fraction = stage2_result.attention.selected_fraction;
            end
            attention_score = min(1, attention_fraction ./ 0.18);
            fit_score = 0.5;
            if isfield(stage1_result, 'metrics') && isfield(stage1_result.metrics, 'R2') && isfinite(stage1_result.metrics.R2)
                fit_score = min(1, max(0, stage1_result.metrics.R2 - 0.999) ./ 0.001);
            end
            uncertainty_score = exp(-stage2_result.T_uncertainty ./ 220);
            linearity_score = 0.5;
            if isfield(stage2_result, 'temperature_profile') && isfield(stage2_result.temperature_profile, 'linearity_R2') && ...
                    isfield(stage2_result, 'T_grid') && isfield(stage2_result, 'T_fit')
                [~, idx] = min(abs(stage2_result.T_grid(:) - stage2_result.T_fit));
                vals = stage2_result.temperature_profile.linearity_R2;
                if ~isempty(idx) && idx >= 1 && idx <= numel(vals) && isfinite(vals(idx))
                    linearity_score = min(max(vals(idx), 0), 1);
                end
            end
            bounds = [NaN, NaN];
            margin = 40;
            if isfield(stage2_result, 'settings')
                if isfield(stage2_result.settings, 'T_soft_bounds')
                    bounds = stage2_result.settings.T_soft_bounds(:).';
                end
                if isfield(stage2_result.settings, 'stage2_boundary_margin_K')
                    margin = stage2_result.settings.stage2_boundary_margin_K;
                end
            elseif isfield(stage1_result, 'settings') && isfield(stage1_result.settings, 'T_soft_bounds')
                bounds = stage1_result.settings.T_soft_bounds(:).';
            end
            boundary_score = 1;
            if numel(bounds) >= 2 && all(isfinite(bounds(1:2))) && isfield(stage2_result, 'T_fit') && isfinite(stage2_result.T_fit)
                lo = min(bounds(1:2));
                hi = max(bounds(1:2));
                if stage2_result.T_fit <= lo + margin || stage2_result.T_fit >= hi - margin
                    boundary_score = 0.25;
                end
            end
            score = max(0, min(1, 0.25 .* attention_score + 0.25 .* fit_score + ...
                0.25 .* uncertainty_score + 0.25 .* linearity_score));
            score = score .* boundary_score;
        end

        function [cube, ene, info] = load_datacube_mat(filename)
            [~, ~, ext] = fileparts(filename);
            if strcmpi(ext, '.dm3')
                imported = DM3Import(filename);
                cube = imported.asDatacube();
                ene = double(imported.ene(:));
                info = imported.info;
                info.source_file = filename;
                info.source_variable = 'DM3Import';
                info.source_class = class(imported);
                info.size = size(cube);
                info.energy_points = numel(ene);
                return;
            end
            S = load(filename);
            names = fieldnames(S);
            cube = [];
            ene = [];
            info = struct('source_file', filename, 'source_variable', '');
            for ii = 1:numel(names)
                val = S.(names{ii});
                if isobject(val) && isprop(val, 'data') && isprop(val, 'ene')
                    cube = double(val.data);
                    ene = double(val.ene(:));
                    info.source_variable = names{ii};
                    info.source_class = class(val);
                    break;
                elseif isnumeric(val) && ndims(val) == 3 && isempty(cube)
                    cube = double(val);
                    info.source_variable = names{ii};
                    info.source_class = class(val);
                elseif isnumeric(val) && isvector(val) && (contains(lower(names{ii}), 'ene') || contains(lower(names{ii}), 'energy'))
                    ene = double(val(:));
                end
            end
            if isempty(cube)
                error('No datacube-like variable found in %s.', filename);
            end
            if isempty(ene)
                E = load('eneref.mat', 'eneref');
                ene = double(E.eneref(:));
            end
            sz = size(cube);
            edim = find(sz == numel(ene), 1, 'last');
            if isempty(edim)
                error('Energy axis length does not match datacube dimensions in %s.', filename);
            end
            if edim ~= 3
                order = 1:ndims(cube);
                order([edim, 3]) = order([3, edim]);
                cube = permute(cube, order);
            end
            info.size = size(cube);
            info.energy_points = numel(ene);
        end

        function [cube_binned, info] = bin_datacube_spatial(cube, y_bin, x_bin)
            if nargin < 2 || isempty(y_bin), y_bin = 1; end
            if nargin < 3 || isempty(x_bin), x_bin = 1; end
            y_bin = max(1, round(double(y_bin)));
            x_bin = max(1, round(double(x_bin)));
            [ny, nx, ne] = size(cube);
            ny2 = ceil(ny ./ y_bin);
            nx2 = ceil(nx ./ x_bin);
            cube_binned = NaN(ny2, nx2, ne);
            for yy = 1:ny2
                y1 = (yy - 1) .* y_bin + 1;
                y2 = min(yy .* y_bin, ny);
                for xx = 1:nx2
                    x1 = (xx - 1) .* x_bin + 1;
                    x2 = min(xx .* x_bin, nx);
                    block = cube(y1:y2, x1:x2, :);
                    block = reshape(block, [], ne);
                    cube_binned(yy, xx, :) = mean(block, 1, 'omitnan');
                end
            end
            info = struct();
            info.y_bin = y_bin;
            info.x_bin = x_bin;
            info.input_size = [ny, nx, ne];
            info.output_size = size(cube_binned);
            info.method = 'mean';
        end

        function [cube_binned, ene_binned, info] = bin_datacube_energy(cube, ene, e_bin)
            if nargin < 3 || isempty(e_bin), e_bin = 1; end
            e_bin = max(1, round(double(e_bin)));
            ene = double(ene(:));
            [ny, nx, ne] = size(cube);
            if e_bin <= 1
                cube_binned = cube;
                ene_binned = ene;
                info = struct('energy_bin', 1, 'input_points', ne, ...
                    'output_points', ne, 'method', 'none');
                return;
            end
            ne2 = ceil(ne ./ e_bin);
            cube_binned = NaN(ny, nx, ne2);
            ene_binned = NaN(ne2, 1);
            for ee = 1:ne2
                e1 = (ee - 1) .* e_bin + 1;
                e2 = min(ee .* e_bin, ne);
                cube_binned(:, :, ee) = mean(cube(:, :, e1:e2), 3, 'omitnan');
                ene_binned(ee) = mean(ene(e1:e2), 'omitnan');
            end
            info = struct();
            info.energy_bin = e_bin;
            info.input_points = ne;
            info.output_points = ne2;
            info.method = 'mean consecutive energy channels';
        end

        function pixels = pick_datacube_pixels(ny, nx, settings)
            n_target = min(settings.datacube_sample_pixels, ny * nx);
            strategy = lower(string(settings.datacube_sample_strategy));
            if contains(strategy, "random")
                n_pick = min(ny * nx, max(0, round(n_target)));
                if n_pick <= 0
                    pixels = zeros(1, 0);
                    return;
                end
                if isfield(settings, 'datacube_random_seed') && ~isempty(settings.datacube_random_seed) && ...
                        all(isfinite(settings.datacube_random_seed))
                    old_rng = rng;
                    rng(double(settings.datacube_random_seed), 'twister');
                    pixels = randperm(ny * nx, n_pick);
                    rng(old_rng);
                else
                    pixels = randperm(ny * nx, n_pick);
                end
                return;
            end
            if contains(strategy, "neighbor") && n_target <= 12
                cy = round((ny + 1) ./ 2);
                cx = round((nx + 1) ./ 2);
                offsets = [0 0; 0 1; 1 0; 1 1; 0 -1; -1 0; -1 -1; 1 -1; -1 1; ...
                    -round(ny / 2) -round(nx / 2); round(ny / 2) -round(nx / 2); ...
                    -round(ny / 2) round(nx / 2); round(ny / 2) round(nx / 2)];
                yy = min(max(cy + offsets(:, 1), 1), ny);
                xx = min(max(cx + offsets(:, 2), 1), nx);
                pixels = unique(sub2ind([ny, nx], yy, xx), 'stable');
                pixels = pixels(1:min(numel(pixels), n_target));
                return;
            end
            if contains(strategy, "block")
                by = max(1, min(ny, round(sqrt(n_target))));
                bx = max(1, min(nx, ceil(n_target ./ by)));
                cy = round((ny + 1) ./ 2);
                cx = round((nx + 1) ./ 2);
                y1 = max(1, min(ny - by + 1, cy - floor(by ./ 2)));
                x1 = max(1, min(nx - bx + 1, cx - floor(bx ./ 2)));
                ys = y1:min(ny, y1 + by - 1);
                xs = x1:min(nx, x1 + bx - 1);
                [XX, YY] = meshgrid(xs, ys);
                pixels = sub2ind([ny, nx], YY(:), XX(:));
                pixels = pixels(1:min(numel(pixels), n_target));
                return;
            end
            switch strategy
                case "all"
                    pixels = 1:(ny * nx);
                    return;
                otherwise
                    gy = max(1, min(ny, round(sqrt(n_target))));
                    gx = max(1, min(nx, ceil(n_target / gy)));
                    ys = unique(round(linspace(1, ny, gy)));
                    xs = unique(round(linspace(1, nx, gx)));
                    [XX, YY] = meshgrid(xs, ys);
                    pixels = sub2ind([ny, nx], YY(:), XX(:));
                    center = sub2ind([ny, nx], round((ny + 1) ./ 2), round((nx + 1) ./ 2));
                    pixels = unique([center; pixels(:)], 'stable');
                    if contains(strategy, "neighbor")
                        [py, px] = ind2sub([ny, nx], pixels(:));
                        add_y = [];
                        add_x = [];
                        offsets = [0 0; 0 1; 1 0; 1 1; 0 -1; -1 0];
                        for ii = 1:numel(py)
                            for jj = 1:size(offsets, 1)
                                yy = min(max(py(ii) + offsets(jj, 1), 1), ny);
                                xx = min(max(px(ii) + offsets(jj, 2), 1), nx);
                                add_y(end + 1, 1) = yy; %#ok<AGROW>
                                add_x(end + 1, 1) = xx; %#ok<AGROW>
                            end
                        end
                        pixels = unique([pixels(:); sub2ind([ny, nx], add_y, add_x)], 'stable');
                    end
                    pixels = pixels(1:min(numel(pixels), n_target));
            end
        end

        function profile = bin_energy_profile(ene, vals, nbins)
            ene = ene(:);
            vals = vals(:);
            profile = NaN(1, nbins);
            if isempty(ene) || isempty(vals) || numel(ene) ~= numel(vals)
                return;
            end
            emax = max(abs(ene), [], 'omitnan');
            if ~isfinite(emax) || emax <= 0
                emax = 1;
            end
            edges = linspace(-emax, emax, nbins + 1);
            for bb = 1:nbins
                mask = ene >= edges(bb) & ene < edges(bb + 1) & isfinite(vals);
                if bb == nbins
                    mask = ene >= edges(bb) & ene <= edges(bb + 1) & isfinite(vals);
                end
                if nnz(mask) > 0
                    profile(bb) = mean(vals(mask), 'omitnan');
                end
            end
            if all(~isfinite(profile))
                profile(:) = 0;
            else
                fillv = median(profile(isfinite(profile)), 'omitnan');
                profile(~isfinite(profile)) = fillv;
                scale = max(abs(profile), [], 'omitnan');
                if isfinite(scale) && scale > 0
                    profile = profile ./ scale;
                end
            end
        end

        function delta = neighbor_profile_delta(profile_cube, yy, xx, y1, y2, x1, x2, keep)
            if isempty(profile_cube) || ~ismatrix(profile_cube(:, :, 1))
                delta = zeros(nnz(keep), 1);
                return;
            end
            center = squeeze(profile_cube(yy, xx, :));
            if isempty(center) || any(~isfinite(center))
                delta = zeros(nnz(keep), 1);
                return;
            end
            block = profile_cube(y1:y2, x1:x2, :);
            [rr, cc] = find(keep);
            delta = zeros(numel(rr), 1);
            for ii = 1:numel(rr)
                prof = squeeze(block(rr(ii), cc(ii), :));
                if isempty(prof) || any(~isfinite(prof))
                    delta(ii) = 0;
                else
                    delta(ii) = sqrt(mean((prof(:) - center(:)) .^ 2, 'omitnan'));
                end
            end
            delta(~isfinite(delta)) = 0;
        end

        function delta = neighbor_peak_shift_delta(profile_cube, yy, xx, y1, y2, x1, x2, keep)
            if isempty(profile_cube) || ndims(profile_cube) < 3
                delta = zeros(nnz(keep), 1);
                return;
            end
            center = squeeze(profile_cube(yy, xx, :));
            center_peak = EELS_PearsonSparse_FullFit.profile_peak_bin(center);
            if ~isfinite(center_peak)
                delta = zeros(nnz(keep), 1);
                return;
            end
            block = profile_cube(y1:y2, x1:x2, :);
            [rr, cc] = find(keep);
            nb = max(size(profile_cube, 3) - 1, 1);
            delta = zeros(numel(rr), 1);
            for ii = 1:numel(rr)
                prof = squeeze(block(rr(ii), cc(ii), :));
                pk = EELS_PearsonSparse_FullFit.profile_peak_bin(prof);
                if isfinite(pk)
                    delta(ii) = abs(pk - center_peak) ./ nb;
                end
            end
            delta(~isfinite(delta)) = 0;
        end

        function pk = profile_peak_bin(profile)
            profile = double(profile(:));
            good = isfinite(profile);
            if nnz(good) < 1
                pk = NaN;
                return;
            end
            profile(~good) = min(profile(good), [], 'omitnan');
            profile = profile - min(profile, [], 'omitnan');
            if max(profile, [], 'omitnan') <= 0
                pk = NaN;
                return;
            end
            [~, imax] = max(profile);
            lo = max(1, imax - 1);
            hi = min(numel(profile), imax + 1);
            w = profile(lo:hi);
            bins = (lo:hi).';
            if sum(w) > 0
                pk = sum(bins .* w) ./ sum(w);
            else
                pk = imax;
            end
        end

        function spatial = evaluate_spatial_consistency(T_map, R2_map, NRMSE_map, reliability_map, attention_map, status_map, settings, attention_profiles, signal_profiles)
            if nargin < 8, attention_profiles = []; end
            if nargin < 9, signal_profiles = []; end
            [ny, nx] = size(T_map);
            neighbor_loss = NaN(ny, nx);
            spectrum_neighbor_loss = NaN(ny, nx);
            local_median = NaN(ny, nx);
            local_scale = NaN(ny, nx);
            interface_score = zeros(ny, nx);
            temperature_neighbor_z = NaN(ny, nx);
            attention_peak_delta_map = NaN(ny, nx);
            signal_peak_delta_map = NaN(ny, nx);
            spatial_refit_candidate = false(ny, nx);
            low_R2_candidate = false(ny, nx);
            bad = false(ny, nx);
            warning_map = zeros(ny, nx);
            for yy = 1:ny
                for xx = 1:nx
                    t0 = T_map(yy, xx);
                    if ~isfinite(t0)
                        if status_map(yy, xx) ~= 0
                            bad(yy, xx) = true;
                            warning_map(yy, xx) = bitor(uint16(warning_map(yy, xx)), uint16(8));
                        end
                        continue;
                    end
                    y1 = max(1, yy - 1); y2 = min(ny, yy + 1);
                    x1 = max(1, xx - 1); x2 = min(nx, xx + 1);
                    nbr_T = T_map(y1:y2, x1:x2);
                    nbr_A = attention_map(y1:y2, x1:x2);
                    keep = isfinite(nbr_T);
                    keep(yy - y1 + 1, xx - x1 + 1) = false;
                    if nnz(keep) < 2
                        continue;
                    end
                    vals = nbr_T(keep);
                    att_vals = nbr_A(keep);
                    att0 = attention_map(yy, xx);
                    if ~isfinite(att0), att0 = median(att_vals(isfinite(att_vals)), 'omitnan'); end
                    if ~isfinite(att0), att0 = 0; end
                    att_vals(~isfinite(att_vals)) = att0;
                    profile_delta = EELS_PearsonSparse_FullFit.neighbor_profile_delta(attention_profiles, yy, xx, y1, y2, x1, x2, keep);
                    signal_delta = EELS_PearsonSparse_FullFit.neighbor_profile_delta(signal_profiles, yy, xx, y1, y2, x1, x2, keep);
                    attention_peak_delta = EELS_PearsonSparse_FullFit.neighbor_peak_shift_delta(attention_profiles, yy, xx, y1, y2, x1, x2, keep);
                    signal_peak_delta = EELS_PearsonSparse_FullFit.neighbor_peak_shift_delta(signal_profiles, yy, xx, y1, y2, x1, x2, keep);
                    combined_delta = max(abs(att_vals - att0), settings.spatial_interface_attention_weight .* profile_delta);
                    combined_delta = max(combined_delta, 0.5 .* signal_delta);
                    combined_delta = max(combined_delta, settings.spatial_peak_shift_weight .* attention_peak_delta);
                    combined_delta = max(combined_delta, settings.spatial_signal_peak_shift_weight .* signal_peak_delta);
                    edge_w = exp(-combined_delta ./ max(settings.spatial_edge_sigma, eps));
                    trusted = edge_w >= 0.25;
                    if nnz(trusted) >= 2
                        vals_use = vals(trusted);
                    else
                        vals_use = vals;
                    end
                    med = median(vals_use, 'omitnan');
                    madv = median(abs(vals_use - med), 'omitnan');
                    scale = max(1.4826 .* madv, settings.spatial_T_floor_K);
                    z = abs(t0 - med) ./ max(scale, eps);
                    edge_mean = mean(edge_w, 'omitnan');
                    z = z(1);
                    edge_mean = edge_mean(1);
                    temperature_neighbor_z(yy, xx) = z;
                    neighbor_loss(yy, xx) = settings.spatial_neighbor_weight * z * edge_mean;
                    spectrum_neighbor_loss(yy, xx) = settings.spatial_spectrum_smooth_weight .* mean(signal_delta .* edge_w, 'omitnan');
                    local_median(yy, xx) = med;
                    local_scale(yy, xx) = scale;
                    attention_peak_delta_map(yy, xx) = max(attention_peak_delta, [], 'omitnan');
                    signal_peak_delta_map(yy, xx) = max(signal_peak_delta, [], 'omitnan');
                    att_delta = combined_delta;
                    att_delta = att_delta(isfinite(att_delta));
                    if isempty(att_delta)
                        interface_score(yy, xx) = 0;
                    else
                        interface_score(yy, xx) = max(att_delta);
                    end
                    if interface_score(yy, xx) > settings.spatial_edge_sigma
                        warning_map(yy, xx) = bitor(uint16(warning_map(yy, xx)), uint16(2));
                    end
                    is_outlier = z > settings.spatial_bad_pixel_z && edge_mean > 0.35;
                    low_profile_change = max([mean(profile_delta, 'omitnan'), mean(signal_delta, 'omitnan')]) <= settings.spatial_refit_low_profile_delta_max;
                    low_peak_shift = max([attention_peak_delta_map(yy, xx), signal_peak_delta_map(yy, xx)]) <= settings.spatial_refit_peak_shift_max;
                    temp_jump = abs(t0 - med) >= settings.spatial_refit_min_temp_delta_K;
                    if temp_jump && low_profile_change && low_peak_shift && edge_mean > 0.30
                        spatial_refit_candidate(yy, xx) = true;
                        warning_map(yy, xx) = bitor(uint16(warning_map(yy, xx)), uint16(16));
                    end
                    low_R2_here = false;
                    if settings.spatial_bad_R2_enabled && isfinite(R2_map(yy, xx))
                        nbr_R2 = R2_map(y1:y2, x1:x2);
                        keep_R2 = isfinite(nbr_R2);
                        keep_R2(yy - y1 + 1, xx - x1 + 1) = false;
                        r2_med = median(nbr_R2(keep_R2), 'omitnan');
                        low_R2_here = R2_map(yy, xx) < settings.spatial_bad_R2_abs_threshold || ...
                            (isfinite(r2_med) && (r2_med - R2_map(yy, xx)) >= settings.spatial_bad_R2_local_drop);
                    end
                    if low_R2_here && edge_mean > 0.25
                        low_R2_candidate(yy, xx) = true;
                        spatial_refit_candidate(yy, xx) = true;
                        warning_map(yy, xx) = bitor(uint16(warning_map(yy, xx)), uint16(32));
                    end
                    weak_fit = (isfinite(reliability_map(yy, xx)) && reliability_map(yy, xx) < settings.spatial_reliability_floor) || ...
                        (isfinite(R2_map(yy, xx)) && R2_map(yy, xx) < 0.98) || ...
                        (isfinite(NRMSE_map(yy, xx)) && NRMSE_map(yy, xx) > 0.03);
                    if is_outlier
                        bad(yy, xx) = true;
                        warning_map(yy, xx) = bitor(uint16(warning_map(yy, xx)), uint16(1));
                    end
                    if weak_fit
                        warning_map(yy, xx) = bitor(uint16(warning_map(yy, xx)), uint16(4));
                    end
                    if status_map(yy, xx) < 0
                        bad(yy, xx) = true;
                        warning_map(yy, xx) = bitor(uint16(warning_map(yy, xx)), uint16(8));
                    end
                end
            end
            T_repaired = T_map;
            if settings.spatial_repair_isolated
                for yy = 1:ny
                    for xx = 1:nx
                        if ~bad(yy, xx) || ~isfinite(local_median(yy, xx))
                            continue;
                        end
                        y1 = max(1, yy - 1); y2 = min(ny, yy + 1);
                        x1 = max(1, xx - 1); x2 = min(nx, xx + 1);
                        nbr_bad = bad(y1:y2, x1:x2);
                        nbr_bad(yy - y1 + 1, xx - x1 + 1) = false;
                        if nnz(nbr_bad) <= 1
                            T_repaired(yy, xx) = local_median(yy, xx);
                        end
                    end
                end
            end
            T_spatial_refined = T_repaired;
            spatial_refine_mask = false(ny, nx);
            spatial_post_threshold = NaN;
            if settings.spatial_post_smooth_enabled
                finite_loss = neighbor_loss(isfinite(neighbor_loss));
                if ~isempty(finite_loss)
                    spatial_post_threshold = EELS_PearsonSparse_FullFit.robust_quantile( ...
                        finite_loss, settings.spatial_post_smooth_loss_quantile);
                end
                if isfinite(spatial_post_threshold)
                    n_iter = max(0, round(settings.spatial_post_smooth_iterations));
                    smooth_w = settings.spatial_post_smooth_weight;
                    if ~isscalar(smooth_w), smooth_w = smooth_w(1); end
                    smooth_w = min(max(smooth_w, 0), 1);
                    edge_sigma = settings.spatial_edge_sigma;
                    if ~isscalar(edge_sigma), edge_sigma = edge_sigma(1); end
                    edge_sigma = double(edge_sigma);
                    interface_release = settings.spatial_post_smooth_interface_release;
                    if ~isscalar(interface_release), interface_release = interface_release(1); end
                    interface_release = double(interface_release);
                    interface_limit = edge_sigma * max(interface_release, eps);
                    rel_floor = settings.spatial_post_smooth_reliability_floor;
                    if ~isscalar(rel_floor), rel_floor = rel_floor(1); end
                    rel_floor = double(rel_floor);
                    for iter = 1:n_iter
                        T_next = T_spatial_refined;
                        for yy = 1:ny
                            for xx = 1:nx
                                if ~isfinite(T_spatial_refined(yy, xx)) || status_map(yy, xx) < 0
                                    continue;
                                end
                                if ~(isfinite(neighbor_loss(yy, xx)) && neighbor_loss(yy, xx) >= spatial_post_threshold)
                                    continue;
                                end
                                if isfinite(reliability_map(yy, xx)) && reliability_map(yy, xx) < rel_floor
                                    continue;
                                end
                                if interface_score(yy, xx) > interface_limit
                                    continue;
                                end
                                y1 = max(1, yy - 1); y2 = min(ny, yy + 1);
                                x1 = max(1, xx - 1); x2 = min(nx, xx + 1);
                                nbr_T = T_spatial_refined(y1:y2, x1:x2);
                                nbr_A = attention_map(y1:y2, x1:x2);
                                keep = isfinite(nbr_T);
                                keep(yy - y1 + 1, xx - x1 + 1) = false;
                                if nnz(keep) < 2
                                    continue;
                                end
                                vals = nbr_T(keep);
                                att_vals = nbr_A(keep);
                                att0 = attention_map(yy, xx);
                                if ~isfinite(att0), att0 = median(att_vals(isfinite(att_vals)), 'omitnan'); end
                                if ~isfinite(att0), att0 = 0; end
                                att_vals(~isfinite(att_vals)) = att0;
                                profile_delta = EELS_PearsonSparse_FullFit.neighbor_profile_delta( ...
                                    attention_profiles, yy, xx, y1, y2, x1, x2, keep);
                                signal_delta = EELS_PearsonSparse_FullFit.neighbor_profile_delta( ...
                                    signal_profiles, yy, xx, y1, y2, x1, x2, keep);
                                attention_peak_delta = EELS_PearsonSparse_FullFit.neighbor_peak_shift_delta( ...
                                    attention_profiles, yy, xx, y1, y2, x1, x2, keep);
                                signal_peak_delta = EELS_PearsonSparse_FullFit.neighbor_peak_shift_delta( ...
                                    signal_profiles, yy, xx, y1, y2, x1, x2, keep);
                                combined_delta = max(abs(att_vals - att0), settings.spatial_interface_attention_weight .* profile_delta);
                                combined_delta = max(combined_delta, 0.5 .* signal_delta);
                                combined_delta = max(combined_delta, settings.spatial_peak_shift_weight .* attention_peak_delta);
                                combined_delta = max(combined_delta, settings.spatial_signal_peak_shift_weight .* signal_peak_delta);
                                edge_w = exp(-combined_delta ./ max(settings.spatial_edge_sigma, eps));
                                trusted = isfinite(vals) & isfinite(edge_w) & edge_w >= 0.20;
                                if nnz(trusted) < 2
                                    continue;
                                end
                                med = EELS_PearsonSparse_FullFit.weighted_quantile(vals(trusted), edge_w(trusted), 0.5);
                                if ~isfinite(med)
                                    continue;
                                end
                                T_next(yy, xx) = (1 - smooth_w) .* T_spatial_refined(yy, xx) + smooth_w .* med;
                                spatial_refine_mask(yy, xx) = true;
                            end
                        end
                        T_spatial_refined = T_next;
                    end
                end
            end
            interface_mask_final = interface_score > settings.spatial_edge_sigma;
            [region_label_map, region_median_map, region_delta_map, region_z_map, region_size_map] = ...
                EELS_PearsonSparse_FullFit.spatial_region_temperature_stats(T_map, interface_mask_final, status_map, settings);
            region_candidate = settings.spatial_region_refit_enabled & ...
                region_size_map >= settings.spatial_region_refit_min_pixels & ...
                region_delta_map >= settings.spatial_region_refit_delta_K & ...
                attention_peak_delta_map <= settings.spatial_refit_peak_shift_max & ...
                signal_peak_delta_map <= settings.spatial_refit_peak_shift_max;
            spatial_refit_candidate = spatial_refit_candidate | region_candidate;
            region_post_mask = false(size(T_map));
            if settings.spatial_region_post_smooth_enabled
                global_region_p90 = EELS_PearsonSparse_FullFit.robust_quantile(region_delta_map(:), 0.90);
                if isfinite(global_region_p90) && global_region_p90 >= settings.spatial_region_post_smooth_min_global_p90_K
                    [T_spatial_refined, region_post_mask] = EELS_PearsonSparse_FullFit.region_post_smooth_temperature( ...
                        T_spatial_refined, region_label_map, status_map, settings);
                    spatial_refine_mask = spatial_refine_mask | region_post_mask;
                end
            end

            spatial = struct();
            spatial.T_map_repaired = T_repaired;
            spatial.T_map_spatial_refined = T_spatial_refined;
            spatial.spatial_refine_delta_map = T_spatial_refined - T_map;
            spatial.spatial_refine_mask = spatial_refine_mask;
            spatial.spatial_post_threshold = spatial_post_threshold;
            spatial.bad_pixel_mask = bad;
            spatial.warning_map = double(warning_map);
            spatial.interface_score_map = interface_score;
            spatial.interface_mask = interface_mask_final;
            spatial.neighbor_loss_map = neighbor_loss;
            spatial.spectrum_neighbor_loss_map = spectrum_neighbor_loss;
            spatial.temperature_neighbor_z_map = temperature_neighbor_z;
            spatial.attention_peak_delta_map = attention_peak_delta_map;
            spatial.signal_peak_delta_map = signal_peak_delta_map;
            spatial.spatial_refit_candidate_mask = spatial_refit_candidate;
            spatial.low_R2_candidate_mask = low_R2_candidate;
            spatial.region_post_smooth_mask = region_post_mask;
            spatial.region_label_map = region_label_map;
            spatial.region_median_map = region_median_map;
            spatial.region_delta_map = region_delta_map;
            spatial.region_z_map = region_z_map;
            spatial.region_size_map = region_size_map;
            spatial.local_median_map = local_median;
            spatial.local_scale_map = local_scale;
            spatial.method = 'edge-aware local median plus optional high-neighbor-loss spectral refit and post smoothing; attention/signal peak shifts reduce coupling so interfaces are preserved';
        end

        function [label_map, median_map, delta_map, z_map, size_map] = spatial_region_temperature_stats(T_map, interface_mask, status_map, settings)
            [ny, nx] = size(T_map);
            label_map = zeros(ny, nx);
            median_map = NaN(ny, nx);
            delta_map = NaN(ny, nx);
            z_map = NaN(ny, nx);
            size_map = zeros(ny, nx);
            valid = isfinite(T_map) & status_map > 0 & ~interface_mask;
            label = 0;
            for yy0 = 1:ny
                for xx0 = 1:nx
                    if ~valid(yy0, xx0) || label_map(yy0, xx0) ~= 0
                        continue;
                    end
                    label = label + 1;
                    queue = zeros(ny * nx, 2);
                    head = 1;
                    tail = 1;
                    queue(tail, :) = [yy0, xx0];
                    label_map(yy0, xx0) = label;
                    members = zeros(ny * nx, 1);
                    nm = 0;
                    while head <= tail
                        yy = queue(head, 1);
                        xx = queue(head, 2);
                        head = head + 1;
                        nm = nm + 1;
                        members(nm) = sub2ind([ny, nx], yy, xx);
                        nbr = [yy - 1, xx; yy + 1, xx; yy, xx - 1; yy, xx + 1];
                        for kk = 1:size(nbr, 1)
                            y2 = nbr(kk, 1);
                            x2 = nbr(kk, 2);
                            if y2 < 1 || y2 > ny || x2 < 1 || x2 > nx
                                continue;
                            end
                            if valid(y2, x2) && label_map(y2, x2) == 0
                                tail = tail + 1;
                                queue(tail, :) = [y2, x2];
                                label_map(y2, x2) = label;
                            end
                        end
                    end
                    members = members(1:nm);
                    vals = T_map(members);
                    med = median(vals, 'omitnan');
                    madv = median(abs(vals - med), 'omitnan');
                    scale = max(1.4826 .* madv, settings.spatial_T_floor_K);
                    median_map(members) = med;
                    delta_map(members) = abs(vals - med);
                    z_map(members) = abs(vals - med) ./ max(scale, eps);
                    size_map(members) = nm;
                end
            end
        end

        function [T_out, smooth_mask] = region_post_smooth_temperature(T_in, label_map, status_map, settings)
            T_out = T_in;
            smooth_mask = false(size(T_in));
            n_iter = max(0, round(settings.spatial_region_post_smooth_iterations));
            if n_iter <= 0
                return;
            end
            labels = unique(label_map(:));
            labels = labels(labels > 0);
            if isempty(labels)
                return;
            end
            base_w = min(max(settings.spatial_region_post_smooth_weight, 0), 1);
            outlier_w = min(max(settings.spatial_region_post_smooth_outlier_weight, base_w), 1);
            z_thr = max(settings.spatial_region_post_smooth_z, 0);
            for iter = 1:n_iter
                T_next = T_out;
                for ll = 1:numel(labels)
                    idx = label_map == labels(ll) & status_map > 0 & isfinite(T_out);
                    if nnz(idx) < settings.spatial_region_refit_min_pixels
                        continue;
                    end
                    vals = T_out(idx);
                    med = median(vals, 'omitnan');
                    madv = median(abs(vals - med), 'omitnan');
                    scale = max(1.4826 .* madv, settings.spatial_T_floor_K);
                    region_p90 = EELS_PearsonSparse_FullFit.robust_quantile(abs(vals - med), 0.90);
                    if ~isfinite(region_p90) || region_p90 < settings.spatial_region_post_smooth_min_region_p90_K
                        continue;
                    end
                    z = abs(vals - med) ./ max(scale, eps);
                    w = base_w .* ones(size(vals));
                    outlier = z >= z_thr | abs(vals - med) >= settings.spatial_region_refit_delta_K;
                    w(outlier) = outlier_w;
                    old_vals = vals;
                    step = w .* (med - vals);
                    max_step = max(settings.spatial_region_post_smooth_max_step_K, 0);
                    if isfinite(max_step) && max_step > 0
                        step = min(max(step, -max_step), max_step);
                    end
                    vals = vals + step;
                    T_next(idx) = vals;
                    changed = abs(vals - old_vals) > 1e-6;
                    lin = find(idx);
                    smooth_mask(lin(changed)) = true;
                end
                T_out = T_next;
            end
        end

        function [results, maps, spatial, history] = spatial_refit_datacube(ene, cube, pixels, results, maps, settings, spatial)
            history = struct('iteration', {}, 'candidate_count', {}, 'attempted', {}, 'accepted', {});
            if ~settings.run_stage2
                return;
            end
            n_iter = max(0, round(settings.spatial_refit_iterations));
            if n_iter <= 0
                return;
            end
            [ny, nx, ~] = size(cube);
            nbins = settings.spatial_attention_profile_bins;
            for iter = 1:n_iter
                finite_loss = spatial.neighbor_loss_map(isfinite(spatial.neighbor_loss_map));
                if isempty(finite_loss)
                    break;
                end
                loss_threshold = EELS_PearsonSparse_FullFit.robust_quantile(finite_loss, settings.spatial_refit_loss_quantile);
                region_candidate = false(size(maps.T_map_raw));
                if settings.spatial_region_refit_enabled && isfield(spatial, 'region_delta_map')
                    region_candidate = spatial.region_size_map >= settings.spatial_region_refit_min_pixels & ...
                        spatial.region_delta_map >= settings.spatial_region_refit_delta_K & ...
                        spatial.attention_peak_delta_map <= settings.spatial_refit_peak_shift_max & ...
                        spatial.signal_peak_delta_map <= settings.spatial_refit_peak_shift_max;
                end
                candidate = maps.fit_status_map > 0 & isfinite(maps.T_map_raw) & ...
                    ~spatial.interface_mask & isfinite(spatial.local_median_map);
                low_R2_candidate = false(size(candidate));
                if isfield(spatial, 'low_R2_candidate_mask')
                    low_R2_candidate = spatial.low_R2_candidate_mask;
                end
                candidate = candidate & (spatial.spatial_refit_candidate_mask | ...
                    spatial.neighbor_loss_map >= loss_threshold | region_candidate | low_R2_candidate);
                candidate = candidate & (abs(maps.T_map_raw - spatial.local_median_map) >= settings.spatial_refit_min_temp_delta_K | ...
                    region_candidate | low_R2_candidate);
                cand_pix = find(candidate(:)).';
                if isempty(cand_pix)
                    history(end + 1) = struct('iteration', iter, 'candidate_count', 0, 'attempted', 0, 'accepted', 0); %#ok<AGROW>
                    break;
                end
                region_z = zeros(size(cand_pix));
                if isfield(spatial, 'region_z_map')
                    region_z = spatial.region_z_map(cand_pix);
                end
                score = spatial.neighbor_loss_map(cand_pix) + 0.08 .* spatial.temperature_neighbor_z_map(cand_pix) + ...
                    0.10 .* region_z;
                if any(low_R2_candidate(:))
                    low_R2_flat = low_R2_candidate(cand_pix);
                    r2_deficit = max(0, settings.spatial_bad_R2_abs_threshold - maps.R2_map(cand_pix));
                    score = score + settings.spatial_refit_low_R2_priority_bonus .* double(low_R2_flat) + ...
                        120 .* r2_deficit;
                end
                score(~isfinite(score)) = 0;
                [~, order] = sort(score, 'descend');
                cand_pix = cand_pix(order);
                max_attempts = min(numel(cand_pix), max(1, ceil(settings.spatial_refit_max_fraction .* nnz(candidate))));
                max_attempts = min(max_attempts, max(1, round(settings.spatial_refit_max_pixels_per_iter)));
                cand_pix = cand_pix(1:max_attempts);
                attempted = 0;
                accepted = 0;
                for jj = 1:numel(cand_pix)
                    pix = cand_pix(jj);
                    [yy, xx] = ind2sub([ny, nx], pix);
                    target_T = spatial.local_median_map(yy, xx);
                    if ~isfinite(target_T)
                        continue;
                    end
                    if settings.spatial_region_refit_enabled && isfield(spatial, 'region_median_map') && ...
                            isfinite(spatial.region_median_map(yy, xx)) && spatial.region_size_map(yy, xx) >= settings.spatial_region_refit_min_pixels
                        region_T = spatial.region_median_map(yy, xx);
                        blend = min(max(settings.spatial_region_refit_weight, 0), 1);
                        target_T = (1 - blend) .* target_T + blend .* region_T;
                    end
                    res_idx = find(pixels == pix, 1, 'first');
                    if isempty(res_idx) || res_idx > numel(results) || ~isstruct(results{res_idx}) || ~isfield(results{res_idx}, 'metrics')
                        continue;
                    end
                    old_res = results{res_idx};
                    if ~isfield(old_res, 'T_fit') || ~isfinite(old_res.T_fit)
                        continue;
                    end
                    attempted = attempted + 1;
                    pixel_settings = settings;
                    pixel_settings.initial_theta = old_res.theta;
                    spectrum = squeeze(cube(yy, xx, :));
                    try
                        new_res = EELS_PearsonSparse_FullFit.solve_stage1(ene, spectrum, pixel_settings);
                    catch
                        continue;
                    end
                    if ~isfield(new_res, 'metrics') || ~isfield(new_res, 'T_fit') || ~isfinite(new_res.T_fit)
                        continue;
                    end
                    old_delta = abs(old_res.T_fit - target_T);
                    new_delta = abs(new_res.T_fit - target_T);
                    old_R2 = old_res.metrics.R2;
                    new_R2 = new_res.metrics.R2;
                    old_NRMSE = old_res.metrics.NRMSE;
                    new_NRMSE = new_res.metrics.NRMSE;
                    improves_T = new_delta <= old_delta - settings.spatial_refit_accept_min_improvement_K;
                    improves_fit = new_R2 >= old_R2 + settings.spatial_refit_accept_min_R2_improvement || ...
                        new_NRMSE <= old_NRMSE - settings.spatial_refit_accept_min_NRMSE_improvement;
                    preserves_fit = new_R2 >= old_R2 - settings.spatial_refit_max_R2_drop && ...
                        new_NRMSE <= old_NRMSE + settings.spatial_refit_max_NRMSE_increase;
                    if (improves_T || improves_fit) && preserves_fit
                        new_res.spatial_refit = struct('accepted', true, 'iteration', iter, ...
                            'target_T_K', target_T, 'old_T_K', old_res.T_fit, ...
                            'old_delta_K', old_delta, 'new_delta_K', new_delta, ...
                            'old_R2', old_R2, 'new_R2', new_R2, ...
                            'old_NRMSE', old_NRMSE, 'new_NRMSE', new_NRMSE);
                        results{res_idx} = new_res;
                        maps = EELS_PearsonSparse_FullFit.update_datacube_maps_from_result( ...
                            maps, new_res, yy, xx, settings, nbins, 2);
                        accepted = accepted + 1;
                    end
                end
                spatial = EELS_PearsonSparse_FullFit.evaluate_spatial_consistency(maps.T_map_raw, maps.R2_map, ...
                    maps.NRMSE_map, maps.reliability_score_map, maps.attention_score_map, maps.fit_status_map, settings, ...
                    maps.attention_profile_cube, maps.signal_profile_cube);
                history(end + 1) = struct('iteration', iter, 'candidate_count', numel(cand_pix), ...
                    'attempted', attempted, 'accepted', accepted); %#ok<AGROW>
                if accepted == 0
                    break;
                end
            end
        end

        function maps = update_datacube_maps_from_result(maps, res, yy, xx, settings, nbins, status_code)
            maps.R2_map(yy, xx) = res.metrics.R2;
            maps.NRMSE_map(yy, xx) = res.metrics.NRMSE;
            maps.T_map_raw(yy, xx) = res.T_fit;
            maps.T_uncertainty_map(yy, xx) = res.T_uncertainty;
            maps.reliability_score_map(yy, xx) = res.reliability_score;
            if isfield(maps, 'reliability_source_map')
                maps.reliability_source_map{yy, xx} = EELS_PearsonSparse_FullFit.reliability_source_from_result(res);
            end
            if isfield(maps, 'stage1_PM_ratio_T_map')
                [maps.stage1_PM_ratio_T_map(yy, xx), maps.stage1_PM_ratio_reliability_map(yy, xx), ...
                    maps.stage1_PM_ratio_pairs_map(yy, xx), ratio_fit_R2, ratio_selected_pairs, ...
                    ratio_selected_EM_lo, ratio_selected_EM_hi] = EELS_PearsonSparse_FullFit.stage1_ratio_summary_from_result(res);
                if isfield(maps, 'stage1_PM_ratio_bestline_R2_map')
                    maps.stage1_PM_ratio_bestline_R2_map(yy, xx) = ratio_fit_R2;
                    maps.stage1_PM_ratio_selected_pairs_map(yy, xx) = ratio_selected_pairs;
                    maps.stage1_PM_ratio_selected_EM_lo_map(yy, xx) = ratio_selected_EM_lo;
                    maps.stage1_PM_ratio_selected_EM_hi_map(yy, xx) = ratio_selected_EM_hi;
                end
            end
            if isfield(maps, 'stage2_zlp_refine_max_delta_map')
                [maps.stage2_zlp_refine_max_delta_map(yy, xx), maps.stage2_zlp_refine_fit_rms_map(yy, xx)] = ...
                    EELS_PearsonSparse_FullFit.stage2_zlp_refine_summary_from_result(res);
            end
            if isfield(res, 'zlp_signal_guard')
                maps.zlp_guard_start_map(yy, xx) = res.zlp_signal_guard.energy_range_meV(1);
                maps.zlp_guard_fraction_map(yy, xx) = res.zlp_signal_guard.zlp_area_fraction;
            end
            if isfield(res, 'zlp_region_policy') && isfield(res.zlp_region_policy, 'core_downweight_width_meV')
                maps.zlp_core_width_map(yy, xx) = res.zlp_region_policy.core_downweight_width_meV;
            end
            maps.attention_profile_cube(yy, xx, :) = NaN;
            if isfield(res.stage2, 'attention') && isfield(res.stage2.attention, 'selected_fraction')
                maps.attention_score_map(yy, xx) = res.stage2.attention.selected_fraction;
                maps.attention_profile_cube(yy, xx, :) = EELS_PearsonSparse_FullFit.bin_energy_profile( ...
                    res.ene_corrected, res.stage2.attention.feature_score, nbins);
            elseif isfield(res, 'attention_summary')
                maps.attention_profile_cube(yy, xx, :) = EELS_PearsonSparse_FullFit.bin_energy_profile( ...
                    res.ene_corrected, res.attention_summary.attention_profile, nbins);
            end
            maps.signal_profile_cube(yy, xx, :) = EELS_PearsonSparse_FullFit.bin_energy_profile( ...
                res.ene_corrected, max(res.raw_data_norm(:) - res.stage1.zlp_norm(:), 0), nbins);
            maps.fit_status_map(yy, xx) = status_code;
        end

        function fig = plot_datacube_maps(cuberes)
            fig = figure('Visible', 'off', 'Color', 'w', 'Position', [60, 60, 1900, 1300]);
            if ~isfield(cuberes, 'zlp_guard_start_map')
                cuberes.zlp_guard_start_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'zlp_guard_fraction_map')
                cuberes.zlp_guard_fraction_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'zlp_core_width_map')
                cuberes.zlp_core_width_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'stage1_PM_ratio_bestline_R2_map')
                cuberes.stage1_PM_ratio_bestline_R2_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'stage1_PM_ratio_selected_pairs_map')
                cuberes.stage1_PM_ratio_selected_pairs_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'stage1_PM_ratio_selected_EM_lo_map')
                cuberes.stage1_PM_ratio_selected_EM_lo_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'stage1_PM_ratio_selected_EM_hi_map')
                cuberes.stage1_PM_ratio_selected_EM_hi_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'T_map_spatial_refined')
                cuberes.T_map_spatial_refined = cuberes.T_map_repaired;
            end
            if ~isfield(cuberes, 'spatial_refine_delta_map')
                cuberes.spatial_refine_delta_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'temperature_neighbor_z_map')
                cuberes.temperature_neighbor_z_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'attention_peak_delta_map')
                cuberes.attention_peak_delta_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'signal_peak_delta_map')
                cuberes.signal_peak_delta_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'spatial_refit_candidate_mask')
                cuberes.spatial_refit_candidate_mask = false(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'low_R2_candidate_mask')
                cuberes.low_R2_candidate_mask = false(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'region_delta_map')
                cuberes.region_delta_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'region_z_map')
                cuberes.region_z_map = NaN(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'region_label_map')
                cuberes.region_label_map = zeros(size(cuberes.R2_map));
            end
            if ~isfield(cuberes, 'region_post_smooth_mask')
                cuberes.region_post_smooth_mask = false(size(cuberes.R2_map));
            end
            maps = {cuberes.R2_map, cuberes.NRMSE_map, cuberes.T_map_raw, cuberes.T_map_repaired, ...
                cuberes.T_map_spatial_refined, cuberes.spatial_refine_delta_map, ...
                cuberes.T_uncertainty_map, cuberes.reliability_score_map, cuberes.attention_score_map, ...
                cuberes.zlp_guard_start_map, cuberes.zlp_guard_fraction_map, cuberes.zlp_core_width_map, ...
                cuberes.stage1_PM_ratio_bestline_R2_map, cuberes.stage1_PM_ratio_selected_pairs_map, ...
                cuberes.stage1_PM_ratio_selected_EM_lo_map, cuberes.stage1_PM_ratio_selected_EM_hi_map, ...
                cuberes.neighbor_loss_map, cuberes.temperature_neighbor_z_map, ...
                cuberes.spectrum_neighbor_loss_map, cuberes.attention_peak_delta_map, ...
                cuberes.signal_peak_delta_map, cuberes.warning_map, ...
                double(cuberes.interface_mask), double(cuberes.bad_pixel_mask), ...
                double(cuberes.spatial_refit_candidate_mask), double(cuberes.low_R2_candidate_mask), cuberes.region_delta_map, ...
                cuberes.region_z_map, cuberes.region_label_map, double(cuberes.region_post_smooth_mask)};
            titles = {'R2', 'NRMSE', 'T raw (K)', 'T repaired (K)', ...
                'T spatial refined (K)', 'spatial refine delta (K)', ...
                'T uncertainty (K)', 'reliability', 'attention fraction', ...
                'ZLP guard start (meV)', 'ZLP fraction in guard', 'ZLP core width (meV)', ...
                'Source-intensity ratio best-line R2', 'Source-intensity ratio selected pairs', ...
                'Source-intensity ratio selected mode-energy low (meV)', 'Source-intensity ratio selected mode-energy high (meV)', ...
                'T neighbor loss', 'temperature neighbor z', ...
                'spectrum neighbor loss', 'attention peak delta', ...
                'signal peak delta', 'warning map', ...
                'interface mask', 'bad pixel mask', 'spatial refit candidates', ...
                'low R2 candidates', 'region T delta (K)', 'region T z', 'region labels', 'region post smooth'};
            for ii = 1:numel(maps)
                subplot(8, 4, ii);
                imagesc(maps{ii});
                [map_ny, map_nx] = size(maps{ii});
                if map_nx / max(map_ny, 1) > 6
                    axis tight;
                    pbaspect([4, 1, 1]);
                else
                    axis image;
                end
                if contains(titles{ii}, 'mask')
                    clim([0, 1]);
                elseif strcmp(titles{ii}, 'warning map')
                    vmax = max(maps{ii}(:), [], 'omitnan');
                    clim([0, max(1, vmax)]);
                end
                colorbar;
                title(titles{ii}, 'Interpreter', 'none');
            end
            sgtitle(sprintf('%s | edge-aware spatial diagnostics', cuberes.filename), 'Interpreter', 'none');
        end

        function write_datacube_summary(filename, cuberes)
            fid = fopen(filename, 'w');
            if fid < 0, return; end
            cleaner = onCleanup(@() fclose(fid));
            fprintf(fid, 'Pearson sparse datacube summary\n');
            fprintf(fid, 'file: %s\n', cuberes.filename);
            fprintf(fid, 'sample_pixels: %d\n', numel(cuberes.sample_pixels));
            fprintf(fid, 'fit_success: %d\n', nnz(cuberes.fit_status_map > 0));
            fprintf(fid, 'fit_failed: %d\n', nnz(cuberes.fit_status_map < 0));
            fprintf(fid, 'median_R2: %.8g\n', median(cuberes.R2_map(:), 'omitnan'));
            fprintf(fid, 'median_NRMSE: %.8g\n', median(cuberes.NRMSE_map(:), 'omitnan'));
            fprintf(fid, 'median_T_raw_K: %.8g\n', median(cuberes.T_map_raw(:), 'omitnan'));
            fprintf(fid, 'median_reliability: %.8g\n', median(cuberes.reliability_score_map(:), 'omitnan'));
            if isfield(cuberes, 'stage2_zlp_refine_max_delta_map')
                fprintf(fid, 'median_stage2_zlp_refine_max_delta_norm: %.8g\n', median(cuberes.stage2_zlp_refine_max_delta_map(:), 'omitnan'));
                zvals = cuberes.stage2_zlp_refine_max_delta_map(:);
                zvals = zvals(isfinite(zvals));
                if isempty(zvals)
                    fprintf(fid, 'p90_stage2_zlp_refine_max_delta_norm: NaN\n');
                else
                    fprintf(fid, 'p90_stage2_zlp_refine_max_delta_norm: %.8g\n', prctile(zvals, 90));
                end
                fprintf(fid, 'median_stage2_zlp_refine_fit_rms_norm: %.8g\n', median(cuberes.stage2_zlp_refine_fit_rms_map(:), 'omitnan'));
            end
            if isfield(cuberes, 'reliability_source_map') && iscell(cuberes.reliability_source_map)
                sources = cuberes.reliability_source_map(:);
                sources = sources(~cellfun(@isempty, sources));
                if ~isempty(sources)
                    sources = cellfun(@(s) char(string(s)), sources, 'UniformOutput', false);
                    [unique_sources, ~, source_idx] = unique(sources);
                    for ss = 1:numel(unique_sources)
                        key = regexprep(unique_sources{ss}, '[^A-Za-z0-9_]', '_');
                        if isempty(key), key = 'empty'; end
                        fprintf(fid, 'reliability_source_count_%s: %d\n', key, nnz(source_idx == ss));
                    end
                end
            end
            if isfield(cuberes, 'zlp_guard_start_map')
                fprintf(fid, 'median_zlp_guard_start_meV: %.8g\n', median(cuberes.zlp_guard_start_map(:), 'omitnan'));
                fprintf(fid, 'median_zlp_guard_fraction: %.8g\n', median(cuberes.zlp_guard_fraction_map(:), 'omitnan'));
                fprintf(fid, 'median_zlp_core_width_meV: %.8g\n', median(cuberes.zlp_core_width_map(:), 'omitnan'));
            end
            if isfield(cuberes, 'cube_info') && isfield(cuberes.cube_info, 'spatial_bin')
                fprintf(fid, 'spatial_bin_yx: %d %d\n', cuberes.cube_info.spatial_bin.y_bin, cuberes.cube_info.spatial_bin.x_bin);
            end
            if isfield(cuberes, 'cube_info') && isfield(cuberes.cube_info, 'energy_bin')
                fprintf(fid, 'energy_bin: %d\n', cuberes.cube_info.energy_bin.energy_bin);
                fprintf(fid, 'energy_points_after_bin: %d\n', cuberes.cube_info.energy_bin.output_points);
            end
            if isfield(cuberes, 'stage1_sparse_pm_policy')
                p = cuberes.stage1_sparse_pm_policy;
                fprintf(fid, 'stage1_sparse_pm_status: %s\n', char(string(p.status)));
                if isfield(p, 'probe_pixel'), fprintf(fid, 'stage1_sparse_pm_probe_pixel: %.8g\n', p.probe_pixel); end
                if isfield(p, 'reference_count'), fprintf(fid, 'stage1_sparse_pm_reference_count: %.8g\n', p.reference_count); end
                if isfield(p, 'selected_count'), fprintf(fid, 'stage1_sparse_pm_selected_count: %.8g\n', p.selected_count); end
                if isfield(p, 'reference_R2'), fprintf(fid, 'stage1_sparse_pm_reference_R2: %.8g\n', p.reference_R2); end
                if isfield(p, 'selected_R2'), fprintf(fid, 'stage1_sparse_pm_selected_R2: %.8g\n', p.selected_R2); end
                if isfield(p, 'message'), fprintf(fid, 'stage1_sparse_pm_message: %s\n', char(string(p.message))); end
            end
            fprintf(fid, 'bad_pixel_fraction_sampled: %.8g\n', nnz(cuberes.bad_pixel_mask) ./ max(nnz(cuberes.fit_status_map ~= 0), 1));
            fprintf(fid, 'interface_fraction_sampled: %.8g\n', nnz(cuberes.interface_mask) ./ max(nnz(cuberes.fit_status_map ~= 0), 1));
            fprintf(fid, 'median_neighbor_loss: %.8g\n', median(cuberes.neighbor_loss_map(:), 'omitnan'));
            fprintf(fid, 'median_spectrum_neighbor_loss: %.8g\n', median(cuberes.spectrum_neighbor_loss_map(:), 'omitnan'));
            if isfield(cuberes, 'temperature_neighbor_z_map')
                fprintf(fid, 'median_temperature_neighbor_z: %.8g\n', median(cuberes.temperature_neighbor_z_map(:), 'omitnan'));
                fprintf(fid, 'median_attention_peak_delta: %.8g\n', median(cuberes.attention_peak_delta_map(:), 'omitnan'));
                fprintf(fid, 'median_signal_peak_delta: %.8g\n', median(cuberes.signal_peak_delta_map(:), 'omitnan'));
                fprintf(fid, 'spatial_refit_candidate_pixels: %d\n', nnz(cuberes.spatial_refit_candidate_mask));
                if isfield(cuberes, 'low_R2_candidate_mask')
                    fprintf(fid, 'low_R2_candidate_pixels: %d\n', nnz(cuberes.low_R2_candidate_mask));
                end
            end
            if isfield(cuberes, 'region_delta_map')
                fprintf(fid, 'median_region_T_delta_K: %.8g\n', median(cuberes.region_delta_map(:), 'omitnan'));
                fprintf(fid, 'p90_region_T_delta_K: %.8g\n', EELS_PearsonSparse_FullFit.robust_quantile(cuberes.region_delta_map(:), 0.90));
                fprintf(fid, 'median_region_T_z: %.8g\n', median(cuberes.region_z_map(:), 'omitnan'));
                fprintf(fid, 'num_spatial_regions: %d\n', max(cuberes.region_label_map(:), [], 'omitnan'));
                if isfield(cuberes, 'region_post_smooth_mask')
                    fprintf(fid, 'region_post_smooth_pixels: %d\n', nnz(cuberes.region_post_smooth_mask));
                end
            end
            if isfield(cuberes, 'spatial_refine_delta_map')
                fprintf(fid, 'spatial_refined_pixels: %d\n', nnz(cuberes.spatial_refine_mask));
                fprintf(fid, 'median_abs_spatial_refine_delta_K: %.8g\n', median(abs(cuberes.spatial_refine_delta_map(:)), 'omitnan'));
                fprintf(fid, 'p90_abs_spatial_refine_delta_K: %.8g\n', EELS_PearsonSparse_FullFit.robust_quantile(abs(cuberes.spatial_refine_delta_map(:)), 0.90));
            end
            if isfield(cuberes, 'spatial_refit_history') && ~isempty(cuberes.spatial_refit_history)
                for ii = 1:numel(cuberes.spatial_refit_history)
                    h = cuberes.spatial_refit_history(ii);
                    fprintf(fid, 'spatial_refit_iter_%d: candidates=%d attempted=%d accepted=%d\n', ...
                        h.iteration, h.candidate_count, h.attempted, h.accepted);
                end
            end
            if isfield(cuberes, 'benchmark')
                fprintf(fid, 'total_seconds: %.8g\n', cuberes.benchmark.total_seconds);
                fprintf(fid, 'median_seconds_per_pixel: %.8g\n', cuberes.benchmark.median_seconds_per_pixel);
                fprintf(fid, 'pixels_per_second: %.8g\n', cuberes.benchmark.pixels_per_second);
                if isfield(cuberes.benchmark, 'parallel_enabled')
                    fprintf(fid, 'parallel_enabled: %d\n', cuberes.benchmark.parallel_enabled);
                    fprintf(fid, 'parallel_workers: %d\n', cuberes.benchmark.parallel_workers);
                    fprintf(fid, 'parallel_message: %s\n', cuberes.benchmark.parallel_message);
                end
            end
            if isfield(cuberes, 'settings')
                fprintf(fid, 'datacube_fit_retries: %d\n', cuberes.settings.datacube_fit_retries);
                fprintf(fid, 'nonlinear_inner_solver: %s\n', string(cuberes.settings.nonlinear_inner_solver));
                fprintf(fid, 'stage2_scan_inner_solver: %s\n', string(cuberes.settings.stage2_scan_inner_solver));
            end
            fprintf(fid, 'spatial_method: %s\n', cuberes.spatial.method);
        end

        function y = unit_scale(x, mask)
            x = x(:);
            mask = mask(:) & isfinite(x);
            vals = x(mask);
            if isempty(vals)
                y = zeros(size(x));
                return;
            end
            lo = median(vals, 'omitnan');
            hi = EELS_PearsonSparse_FullFit.robust_quantile(vals, 0.98);
            if ~(isfinite(hi) && hi > lo)
                hi = max(vals);
            end
            span = max(hi - lo, eps);
            y = (x - lo) ./ span;
            y = min(max(y, 0), 1);
            y(~isfinite(y)) = 0;
        end

        function val = safe_mean(vals)
            vals = vals(:);
            vals = vals(isfinite(vals));
            if isempty(vals)
                val = NaN;
            else
                val = mean(vals, 'omitnan');
            end
        end

        function q = robust_quantile(vals, frac)
            vals = vals(:);
            vals = sort(vals(isfinite(vals)));
            if isempty(vals)
                q = NaN;
                return;
            end
            idx = max(1, min(numel(vals), round(1 + frac .* (numel(vals) - 1))));
            q = vals(idx);
        end

        function q = weighted_quantile(vals, weights, frac)
            vals = vals(:);
            weights = weights(:);
            good = isfinite(vals) & isfinite(weights) & weights > 0;
            vals = vals(good);
            weights = weights(good);
            if isempty(vals)
                q = NaN;
                return;
            end
            [vals, order] = sort(vals);
            weights = weights(order);
            c = cumsum(weights) ./ sum(weights);
            idx = find(c >= frac, 1, 'first');
            if isempty(idx), idx = numel(vals); end
            q = vals(idx);
        end

        function selected = select_spaced_em(EM, score, max_count, min_spacing)
            EM = double(EM(:).');
            score = double(score(:).');
            good = isfinite(EM) & isfinite(score);
            EM = EM(good);
            score = score(good);
            if isempty(EM)
                selected = zeros(1, 0);
                return;
            end
            max_count = max(1, round(double(max_count)));
            min_spacing = max(0, double(min_spacing));
            [~, order] = sort(score, 'descend');
            selected = zeros(1, 0);
            for kk = 1:numel(order)
                cand = EM(order(kk));
                if isempty(selected) || all(abs(cand - selected) >= min_spacing)
                    selected(end + 1) = cand; %#ok<AGROW>
                    if numel(selected) >= max_count
                        break;
                    end
                end
            end
            selected = sort(unique(selected));
        end

        function pm = interpolate_pm(EM_old, pm_old, EM_new)
            EM_old = EM_old(:);
            pm_old = pm_old(:);
            EM_new = EM_new(:);
            good = isfinite(EM_old) & isfinite(pm_old);
            EM_old = EM_old(good);
            pm_old = pm_old(good);
            if isempty(EM_old)
                pm = zeros(size(EM_new));
                return;
            end
            [EM_old, order] = sort(EM_old);
            pm_old = pm_old(order);
            pm = interp1(EM_old, pm_old, EM_new, 'linear', 'extrap');
            pm(~isfinite(pm) | pm < 0) = 0;
        end

        function alpha = smooth_binary_mask_by_energy(mask, ene, taper_meV)
            mask = logical(mask(:));
            ene = ene(:);
            alpha = double(mask);
            if numel(mask) < 3 || numel(ene) ~= numel(mask) || taper_meV <= 0
                return;
            end
            de = median(abs(diff(ene)), 'omitnan');
            if ~(isfinite(de) && de > 0)
                return;
            end
            bands = EELS_PearsonSparse_FullFit.mask_to_bands(ene, mask);
            if isempty(bands)
                alpha = zeros(size(mask));
                return;
            end
            tw = max([de, taper_meV ./ 4, eps]);
            alpha = zeros(size(ene));
            for ii = 1:size(bands, 1)
                lo = min(bands(ii, :));
                hi = max(bands(ii, :));
                left = 1 ./ (1 + exp(-(ene - lo) ./ tw));
                right = 1 ./ (1 + exp((ene - hi) ./ tw));
                alpha = max(alpha, left .* right);
            end
            alpha = min(max(alpha, 0), 1);
        end

        function bands = mask_to_bands(ene, mask)
            ene = ene(:);
            mask = mask(:);
            starts = find(mask & [true; ~mask(1:end-1)]);
            stops = find(mask & [~mask(2:end); true]);
            bands = zeros(numel(starts), 2);
            for ii = 1:numel(starts)
                bands(ii, :) = [ene(starts(ii)), ene(stops(ii))];
            end
        end

        function bands = stage2_temperature_bands(res)
            bands = struct('fit', [], 'attention', [], 'hard', [], 'transition', []);
            if ~isfield(res, 'stage2') || ~isfield(res.stage2, 'status') || ~strcmp(res.stage2.status, 'ok')
                return;
            end
            ene = res.ene_corrected(:);
            st2 = res.stage2;
            fit_mask = false(size(ene));
            if isfield(st2, 'fit_mask') && numel(st2.fit_mask) == numel(ene)
                fit_mask = logical(st2.fit_mask(:));
                bands.fit = EELS_PearsonSparse_FullFit.merge_energy_bands( ...
                    EELS_PearsonSparse_FullFit.mask_to_bands(ene, fit_mask), 0.5);
            end
            if isfield(st2, 'transition_mask') && numel(st2.transition_mask) == numel(ene)
                transition_mask = logical(st2.transition_mask(:));
                bands.transition = EELS_PearsonSparse_FullFit.merge_energy_bands( ...
                    EELS_PearsonSparse_FullFit.mask_to_bands(ene, transition_mask), 0.5);
            end
            if isfield(st2, 'attention') && isfield(st2.attention, 'active_mask') && ...
                    numel(st2.attention.active_mask) == numel(ene)
                active_mask = logical(st2.attention.active_mask(:));
                bands.attention = EELS_PearsonSparse_FullFit.merge_energy_bands( ...
                    EELS_PearsonSparse_FullFit.mask_to_bands(ene, active_mask), 3.0);
            elseif isfield(st2, 'attention') && isfield(st2.attention, 'bands_meV')
                bands.attention = EELS_PearsonSparse_FullFit.merge_energy_bands(st2.attention.bands_meV, 3.0);
            end
            if isfield(st2, 'EM_hard') && ~isempty(st2.EM_hard)
                half_width = 2.0;
                if isfield(res, 'settings')
                    candidates = half_width;
                    if isfield(res.settings, 'stage2_dense_em_step_meV')
                        candidates(end + 1) = 0.5 .* res.settings.stage2_dense_em_step_meV;
                    end
                    if isfield(res.settings, 'stage2_attention_peak_sigma_meV')
                        candidates(end + 1) = 0.35 .* res.settings.stage2_attention_peak_sigma_meV;
                    end
                    half_width = max(candidates);
                end
                centers = [st2.EM_hard(:); -st2.EM_hard(:)];
                centers = centers(isfinite(centers));
                hard_bands = [centers - half_width, centers + half_width];
                emin = min(ene, [], 'omitnan');
                emax = max(ene, [], 'omitnan');
                hard_bands(:, 1) = max(hard_bands(:, 1), emin);
                hard_bands(:, 2) = min(hard_bands(:, 2), emax);
                hard_bands = hard_bands(hard_bands(:, 2) > hard_bands(:, 1), :);
                bands.hard = EELS_PearsonSparse_FullFit.merge_energy_bands(hard_bands, 0.5);
            end
            if isempty(bands.attention) && ~isempty(bands.hard)
                bands.attention = bands.hard;
            end
        end

        function bands = merge_energy_bands(bands, max_gap)
            if nargin < 2, max_gap = 0; end
            if isempty(bands)
                bands = zeros(0, 2);
                return;
            end
            bands = double(bands);
            bands = bands(all(isfinite(bands), 2) & bands(:, 2) >= bands(:, 1), :);
            if isempty(bands)
                bands = zeros(0, 2);
                return;
            end
            bands = sortrows(bands, 1);
            merged = bands(1, :);
            for ii = 2:size(bands, 1)
                if bands(ii, 1) <= merged(end, 2) + max_gap
                    merged(end, 2) = max(merged(end, 2), bands(ii, 2));
                else
                    merged(end + 1, :) = bands(ii, :); %#ok<AGROW>
                end
            end
            bands = merged;
        end

        function shade_energy_bands(ax, bands, yl, color, alpha_value)
            if isempty(bands) || numel(yl) < 2 || ~(isfinite(yl(1)) && isfinite(yl(2))) || yl(2) <= yl(1)
                return;
            end
            hold(ax, 'on');
            for ii = 1:size(bands, 1)
                x1 = bands(ii, 1);
                x2 = bands(ii, 2);
                if ~(isfinite(x1) && isfinite(x2) && x2 > x1)
                    continue;
                end
                h = patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], color, ...
                    'FaceAlpha', alpha_value, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                try
                    uistack(h, 'bottom');
                catch
                end
            end
        end

        function txt = format_band_list(bands, max_count)
            if nargin < 2, max_count = 4; end
            if isempty(bands)
                txt = 'none';
                return;
            end
            n = size(bands, 1);
            nshow = min(n, max_count);
            parts = cell(1, nshow);
            for ii = 1:nshow
                parts{ii} = sprintf('[%.0f, %.0f]', bands(ii, 1), bands(ii, 2));
            end
            txt = strjoin(parts, ', ');
            if n > nshow
                txt = sprintf('%s +%d more', txt, n - nshow);
            end
        end

        function mask = values_in_energy_bands(vals, bands)
            vals = vals(:);
            mask = false(size(vals));
            if isempty(vals) || isempty(bands)
                return;
            end
            bands = double(bands);
            bands = bands(all(isfinite(bands), 2), :);
            if isempty(bands)
                return;
            end
            for ii = 1:size(bands, 1)
                lo = min(bands(ii, :));
                hi = max(bands(ii, :));
                mask = mask | (vals >= lo & vals <= hi);
            end
        end

        function pmview = pm_top_interpolated_spectrum(res, settings)
            if nargin < 2 || isempty(settings)
                if isstruct(res) && isfield(res, 'settings') && isstruct(res.settings)
                    settings = res.settings;
                else
                    settings = struct();
                end
            end
            settings = EELS_PearsonSparse_FullFit.defaults(settings);
            pmview = struct('available', false, 'reason', 'not_started', ...
                'source', 'none', 'ene', [], 'target_signal', [], ...
                'latent_signal', [], 'convolved_signal_unscaled', [], ...
                'convolved_signal', [], 'residual', [], 'zlp_norm', [], ...
                'bg_norm', [], 'reference_vib_norm', [], 'EM_meV', [], ...
                'PM_stokes_top', [], 'PM_antistokes_top', [], 'fit_mask', [], ...
                'sqrtW', [], 'R2_weighted', NaN, 'NRMSE_weighted', NaN, ...
                'scale_mode', char(settings.pm_interp_scale_mode), ...
                'scale_global', NaN, 'scale_loss', NaN, 'scale_gain', NaN, ...
                'peak_sigma_meV', settings.pm_interp_peak_sigma_meV, ...
                'interp_method', char(settings.pm_interp_method), ...
                'top_ratio_T_K', NaN, 'top_ratio_R2', NaN, ...
                'top_ratio_num_points', 0, 'num_centers', 0);
            if ~(isstruct(res) && isfield(res, 'ene_corrected'))
                pmview.reason = 'missing_result_energy_axis';
                return;
            end
            ene = res.ene_corrected(:);
            n = numel(ene);
            if isfield(res, 'raw_data_norm') && numel(res.raw_data_norm) >= n
                raw_norm = res.raw_data_norm(:);
                raw_norm = raw_norm(1:n);
            elseif isfield(res, 'raw_data') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                raw_norm = res.raw_data(:) ./ res.norm_val;
                raw_norm = raw_norm(1:n);
            else
                pmview.reason = 'missing_raw_data';
                return;
            end

            zlp_norm = zeros(n, 1);
            bg_norm = zeros(n, 1);
            reference_vib = zeros(n, 1);
            EM = [];
            pm_s = [];
            pm_a = [];
            source = 'stage1_PM';
            stage2_ok = isfield(res, 'stage2') && isstruct(res.stage2) && ...
                isfield(res.stage2, 'status') && strcmp(res.stage2.status, 'ok');
            if stage2_ok
                if isfield(res.stage2, 'zlp_norm') && numel(res.stage2.zlp_norm) >= n
                    tmp = res.stage2.zlp_norm(:); zlp_norm = tmp(1:n);
                end
                if isfield(res.stage2, 'bg_norm') && numel(res.stage2.bg_norm) >= n
                    tmp = res.stage2.bg_norm(:); bg_norm = tmp(1:n);
                end
                if isfield(res.stage2, 'vib_norm') && numel(res.stage2.vib_norm) >= n
                    tmp = res.stage2.vib_norm(:); reference_vib = tmp(1:n);
                end
                if isfield(res.stage2, 'EM_final') && isfield(res.stage2, 'PM_stokes') && ...
                        isfield(res.stage2, 'PM_antistokes')
                    EM = res.stage2.EM_final(:);
                    pm_s = res.stage2.PM_stokes(:);
                    pm_a = res.stage2.PM_antistokes(:);
                    source = 'stage2_PM';
                end
            end
            if isempty(EM) && isfield(res, 'EM_final') && isfield(res, 'PM_stokes') && isfield(res, 'PM_antistokes')
                EM = res.EM_final(:);
                pm_s = res.PM_stokes(:);
                pm_a = res.PM_antistokes(:);
            end
            if ~any(isfinite(zlp_norm) & abs(zlp_norm) > 0)
                if isfield(res, 'stage1') && isfield(res.stage1, 'zlp_norm') && numel(res.stage1.zlp_norm) >= n
                    tmp = res.stage1.zlp_norm(:); zlp_norm = tmp(1:n);
                elseif isfield(res, 'zlp_only') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                    tmp = res.zlp_only(:); zlp_norm = tmp(1:n) ./ res.norm_val;
                end
            end
            if ~any(isfinite(bg_norm) & abs(bg_norm) > 0)
                if isfield(res, 'stage1') && isfield(res.stage1, 'bg_norm') && numel(res.stage1.bg_norm) >= n
                    tmp = res.stage1.bg_norm(:); bg_norm = tmp(1:n);
                elseif isfield(res, 'bg_only') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                    tmp = res.bg_only(:); bg_norm = tmp(1:n) ./ res.norm_val;
                end
            end
            if ~any(isfinite(reference_vib) & abs(reference_vib) > 0)
                if isfield(res, 'stage1') && isfield(res.stage1, 'vib_norm') && numel(res.stage1.vib_norm) >= n
                    tmp = res.stage1.vib_norm(:); reference_vib = tmp(1:n);
                elseif isfield(res, 'vib_only') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                    tmp = res.vib_only(:); reference_vib = tmp(1:n) ./ res.norm_val;
                end
            end
            if isempty(EM) || numel(pm_s) ~= numel(EM) || numel(pm_a) ~= numel(EM)
                pmview.reason = 'missing_PM_dictionary';
                return;
            end

            target = raw_norm(:) - zlp_norm(:) - bg_norm(:);
            em_lo = max(0, settings.pm_interp_min_em_meV);
            em_hi = max(em_lo + eps, settings.pm_interp_max_em_meV);
            good = isfinite(EM) & isfinite(pm_s) & isfinite(pm_a) & EM(:) >= em_lo & EM(:) <= em_hi;
            EM = EM(good);
            pm_s = max(pm_s(good), 0);
            pm_a = max(pm_a(good), 0);
            if isempty(EM)
                pmview.reason = 'no_PM_tops_in_window';
                return;
            end
            pm_mass = max(pm_s(:), pm_a(:));
            max_mass = max(pm_mass, [], 'omitnan');
            if ~isfinite(max_mass) || max_mass <= 0
                pmview.reason = 'zero_PM_top_amplitude';
                return;
            end
            keep = pm_mass >= settings.pm_interp_min_pm_fraction .* max_mass;
            if nnz(keep) < 2
                [~, order] = sort(pm_mass, 'descend');
                keep(order(1:min(numel(order), 2))) = true;
            end
            cand_EM = EM(keep);
            cand_score = pm_mass(keep);
            max_centers = max(2, round(settings.pm_interp_max_centers));
            if numel(cand_EM) > max_centers
                selected_EM = EELS_PearsonSparse_FullFit.select_spaced_em(cand_EM(:).', cand_score(:).', ...
                    max_centers, settings.pm_interp_min_spacing_meV);
            else
                selected_EM = unique(cand_EM(:).');
            end
            selected_EM = selected_EM(isfinite(selected_EM) & selected_EM >= em_lo & selected_EM <= em_hi);
            selected_EM = sort(selected_EM(:).');
            if numel(selected_EM) < 2
                pmview.reason = 'insufficient_selected_PM_tops';
                return;
            end
            pm_s_sel = EELS_PearsonSparse_FullFit.interpolate_pm(EM(:), pm_s(:), selected_EM(:));
            pm_a_sel = EELS_PearsonSparse_FullFit.interpolate_pm(EM(:), pm_a(:), selected_EM(:));
            pm_s_sel = max(pm_s_sel(:), 0);
            pm_a_sel = max(pm_a_sel(:), 0);

            method = lower(string(settings.pm_interp_method));
            allowed = ["linear", "pchip", "makima", "nearest"];
            if ~any(method == allowed)
                method = "linear";
            end
            eabs = abs(ene(:));
            loss_mask = ene(:) >= 0 & eabs >= em_lo & eabs <= em_hi;
            gain_mask = ene(:) <= 0 & eabs >= em_lo & eabs <= em_hi;
            latent_loss = zeros(n, 1);
            latent_gain = zeros(n, 1);
            latent_loss(loss_mask) = interp1(selected_EM(:), pm_s_sel(:), eabs(loss_mask), char(method), 0);
            latent_gain(gain_mask) = interp1(selected_EM(:), pm_a_sel(:), eabs(gain_mask), char(method), 0);
            latent_loss(~isfinite(latent_loss) | latent_loss < 0) = 0;
            latent_gain(~isfinite(latent_gain) | latent_gain < 0) = 0;

            sigma = max(settings.pm_interp_peak_sigma_meV, eps);
            de = median(abs(diff(ene)), 'omitnan');
            if ~isfinite(de) || de <= 0
                de = 1;
            end
            half_n = max(1, ceil(5 .* sigma ./ max(de, eps)));
            grid = (-half_n:half_n).' .* de;
            g = exp(-0.5 .* (grid ./ max(sigma, eps)) .^ 2);
            g = g ./ max(sum(g, 'omitnan'), eps);
            conv_loss = conv(latent_loss, g, 'same');
            conv_gain = conv(latent_gain, g, 'same');
            conv_loss(~isfinite(conv_loss) | conv_loss < 0) = 0;
            conv_gain(~isfinite(conv_gain) | conv_gain < 0) = 0;

            fit_mask = abs(ene(:)) >= em_lo & abs(ene(:)) <= em_hi & isfinite(target(:));
            fit_mask = fit_mask & ((conv_loss + conv_gain) > 0);
            if nnz(fit_mask) < 8
                pmview.reason = 'insufficient_fit_samples';
                return;
            end
            scale_ref = prctile(abs(target(fit_mask)), 85);
            if ~isfinite(scale_ref) || scale_ref <= 0
                scale_ref = max(abs(target(fit_mask)), [], 'omitnan');
            end
            if ~isfinite(scale_ref) || scale_ref <= 0
                scale_ref = 1;
            end
            amp_weight = min(1, abs(target(:)) ./ max(scale_ref, eps));
            sqrtW = sqrt(settings.pm_interp_data_weight_floor + ...
                (1 - settings.pm_interp_data_weight_floor) .* amp_weight);
            if stage2_ok && isfield(res.stage2, 'attention') && isfield(res.stage2.attention, 'feature_score') && ...
                    numel(res.stage2.attention.feature_score) >= n
                att = res.stage2.attention.feature_score(:);
                sqrtW = sqrtW .* sqrt(0.45 + 0.55 .* max(min(att(1:n), 1), 0));
            end
            sqrtW(~isfinite(sqrtW)) = 0;

            mode = lower(string(settings.pm_interp_scale_mode));
            if mode == "two_sided_nnls"
                X = [conv_loss(:), conv_gain(:)];
                try
                    coef = lsqnonneg(bsxfun(@times, X(fit_mask, :), sqrtW(fit_mask)), ...
                        target(fit_mask) .* sqrtW(fit_mask));
                catch
                    coef = max((bsxfun(@times, X(fit_mask, :), sqrtW(fit_mask)) \ ...
                        (target(fit_mask) .* sqrtW(fit_mask))), 0);
                end
                scale_loss = coef(1);
                scale_gain = coef(2);
                scale_global = NaN;
                convolved = X * coef(:);
            else
                mode = "global_nnls";
                basis = conv_loss(:) + conv_gain(:);
                num = sum((sqrtW(fit_mask) .^ 2) .* basis(fit_mask) .* target(fit_mask), 'omitnan');
                den = sum((sqrtW(fit_mask) .^ 2) .* basis(fit_mask) .^ 2, 'omitnan');
                scale_global = max(num ./ max(den, eps), 0);
                scale_loss = scale_global;
                scale_gain = scale_global;
                convolved = scale_global .* basis(:);
            end
            residual = target(:) - convolved(:);
            w = sqrtW(fit_mask) .^ 2;
            sse = sum(w(:) .* residual(fit_mask) .^ 2, 'omitnan');
            sst = sum(w(:) .* target(fit_mask) .^ 2, 'omitnan');
            R2 = 1 - sse ./ max(sst, eps);
            nrmse = sqrt(sse ./ max(sum(w, 'omitnan'), eps)) ./ ...
                max(sqrt(sum(w(:) .* target(fit_mask) .^ 2, 'omitnan') ./ max(sum(w, 'omitnan'), eps)), eps);
            ratio_line = EELS_PearsonSparse_FullFit.best_ratio_origin_line( ...
                selected_EM(:), log(max(pm_s_sel(:), eps) ./ max(pm_a_sel(:), eps)), ...
                sqrt(max(pm_s_sel(:), 0) .* max(pm_a_sel(:), 0)), settings);

            pmview.available = true;
            pmview.reason = 'ok';
            pmview.source = source;
            pmview.ene = ene(:);
            pmview.target_signal = target(:);
            pmview.latent_signal = scale_loss .* latent_loss(:) + scale_gain .* latent_gain(:);
            pmview.convolved_signal_unscaled = conv_loss(:) + conv_gain(:);
            pmview.convolved_signal = convolved(:);
            pmview.residual = residual(:);
            pmview.zlp_norm = zlp_norm(:);
            pmview.bg_norm = bg_norm(:);
            pmview.reference_vib_norm = reference_vib(:);
            pmview.EM_meV = selected_EM(:);
            pmview.PM_stokes_top = pm_s_sel(:);
            pmview.PM_antistokes_top = pm_a_sel(:);
            pmview.fit_mask = fit_mask(:);
            pmview.sqrtW = sqrtW(:);
            pmview.R2_weighted = R2;
            pmview.NRMSE_weighted = nrmse;
            pmview.scale_mode = char(mode);
            pmview.scale_global = scale_global;
            pmview.scale_loss = scale_loss;
            pmview.scale_gain = scale_gain;
            pmview.peak_sigma_meV = sigma;
            pmview.interp_method = char(method);
            pmview.num_centers = numel(selected_EM);
            if ratio_line.available
                pmview.top_ratio_T_K = ratio_line.T_K;
                pmview.top_ratio_R2 = ratio_line.R2_origin;
                pmview.top_ratio_num_points = ratio_line.num_points;
                pmview.top_ratio_EM_range_meV = ratio_line.EM_range_meV;
            else
                pmview.top_ratio_EM_range_meV = [NaN, NaN];
            end
        end

        function rl = rl_deconvolve_background_subtracted(res, settings)
            if nargin < 2 || isempty(settings)
                if isstruct(res) && isfield(res, 'settings') && isstruct(res.settings)
                    settings = res.settings;
                else
                    settings = struct();
                end
            end
            settings = EELS_PearsonSparse_FullFit.defaults(settings);
            rl = struct('available', false, 'reason', 'not_started', ...
                'source', 'fitted_total_zlp_kernel', 'ene', [], 'target_signal', [], ...
                'positive_input', [], 'latent_signal', [], 'reconvolved_signal', [], ...
                'residual', [], 'zlp_norm', [], 'bg_norm', [], 'reference_vib_norm', [], ...
                'psf', [], 'psf_energy_meV', [], 'fit_mask', [], 'taper', [], ...
                'sqrtW', [], 'R2_weighted', NaN, 'NRMSE_weighted', NaN, ...
                'R2_positive_weighted', NaN, 'NRMSE_positive_weighted', NaN, ...
                'iterations', settings.rl_deconv_iterations, ...
                'kernel_width_scale', settings.rl_deconv_kernel_width_scale, ...
                'latent_smooth_sigma_meV', settings.rl_deconv_latent_smooth_sigma_meV, ...
                'final_scale', NaN, 'EM_meV', [], 'latent_stokes', [], ...
                'latent_antistokes', [], 'latent_ratio_T_K', NaN, ...
                'latent_ratio_R2', NaN, 'latent_ratio_num_points', 0);
            if ~(isstruct(res) && isfield(res, 'ene_corrected'))
                rl.reason = 'missing_result_energy_axis';
                return;
            end
            ene = res.ene_corrected(:);
            n = numel(ene);
            if n < 16
                rl.reason = 'too_few_energy_points';
                return;
            end
            if isfield(res, 'raw_data_norm') && numel(res.raw_data_norm) >= n
                raw_norm = res.raw_data_norm(:);
                raw_norm = raw_norm(1:n);
            elseif isfield(res, 'raw_data') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                raw_norm = res.raw_data(:) ./ res.norm_val;
                raw_norm = raw_norm(1:n);
            else
                rl.reason = 'missing_raw_data';
                return;
            end
            if isfield(res, 'stage1') && isstruct(res.stage1) && isfield(res.stage1, 'params')
                params = res.stage1.params;
            elseif isfield(res, 'zlp_params')
                params = res.zlp_params;
            else
                rl.reason = 'missing_zlp_kernel_params';
                return;
            end

            zlp_norm = zeros(n, 1);
            bg_norm = zeros(n, 1);
            reference_vib = zeros(n, 1);
            if isfield(res, 'stage2') && isstruct(res.stage2) && isfield(res.stage2, 'status') && ...
                    strcmp(res.stage2.status, 'ok')
                if isfield(res.stage2, 'zlp_norm') && numel(res.stage2.zlp_norm) >= n
                    tmp = res.stage2.zlp_norm(:); zlp_norm = tmp(1:n);
                end
                if isfield(res.stage2, 'bg_norm') && numel(res.stage2.bg_norm) >= n
                    tmp = res.stage2.bg_norm(:); bg_norm = tmp(1:n);
                end
                if isfield(res.stage2, 'vib_norm') && numel(res.stage2.vib_norm) >= n
                    tmp = res.stage2.vib_norm(:); reference_vib = tmp(1:n);
                end
            end
            if ~any(isfinite(zlp_norm) & abs(zlp_norm) > 0)
                if isfield(res, 'stage1') && isfield(res.stage1, 'zlp_norm') && numel(res.stage1.zlp_norm) >= n
                    tmp = res.stage1.zlp_norm(:); zlp_norm = tmp(1:n);
                elseif isfield(res, 'zlp_only') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                    tmp = res.zlp_only(:); zlp_norm = tmp(1:n) ./ res.norm_val;
                end
            end
            if ~any(isfinite(bg_norm) & abs(bg_norm) > 0)
                if isfield(res, 'stage1') && isfield(res.stage1, 'bg_norm') && numel(res.stage1.bg_norm) >= n
                    tmp = res.stage1.bg_norm(:); bg_norm = tmp(1:n);
                elseif isfield(res, 'bg_only') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                    tmp = res.bg_only(:); bg_norm = tmp(1:n) ./ res.norm_val;
                end
            end
            if ~any(isfinite(reference_vib) & abs(reference_vib) > 0)
                if isfield(res, 'stage1') && isfield(res.stage1, 'vib_norm') && numel(res.stage1.vib_norm) >= n
                    tmp = res.stage1.vib_norm(:); reference_vib = tmp(1:n);
                elseif isfield(res, 'vib_only') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                    tmp = res.vib_only(:); reference_vib = tmp(1:n) ./ res.norm_val;
                end
            end

            target = raw_norm(:) - zlp_norm(:) - bg_norm(:);
            eabs = abs(ene(:));
            e_lo = max(0, settings.rl_deconv_min_energy_meV);
            e_hi = max(e_lo + eps, settings.rl_deconv_max_energy_meV);
            fit_mask = eabs >= e_lo & eabs <= e_hi & isfinite(target(:));
            if nnz(fit_mask) < 16
                rl.reason = 'insufficient_fit_samples';
                return;
            end
            taper = EELS_PearsonSparse_FullFit.smooth_binary_mask_by_energy( ...
                fit_mask, ene(:), max(settings.rl_deconv_edge_taper_meV, 0));
            positive_input = max(target(:), 0) .* taper(:);
            peak = max(positive_input(fit_mask), [], 'omitnan');
            if ~isfinite(peak) || peak <= 0
                rl.reason = 'zero_positive_input_after_background_subtraction';
                return;
            end
            floor_val = max(settings.rl_deconv_floor_fraction, 0) .* peak;
            positive_input = max(positive_input, floor_val .* double(taper(:) > 0.01));

            de = median(abs(diff(ene)), 'omitnan');
            if ~isfinite(de) || de <= 0
                rl.reason = 'invalid_energy_spacing';
                return;
            end
            half_n = floor((n - 1) ./ 2);
            lag = (-half_n:half_n).' .* de;
            kernel_params = params;
            kscale = settings.rl_deconv_kernel_width_scale;
            if ~isscalar(kscale) || ~isfinite(kscale) || kscale <= 0
                kscale = 1;
            end
            if abs(kscale - 1) > 1e-12
                if isfield(kernel_params, 'sL')
                    kernel_params.sL = max(kernel_params.sL .* kscale, eps);
                end
                if isfield(kernel_params, 'sR')
                    kernel_params.sR = max(kernel_params.sR .* kscale, eps);
                end
            end
            % Stage 1/2 fitted-ZLP kernels use the amplitude-weighted total
            % fitted ZLP. Reuse that exact total here instead of rebuilding an
            % equal-weight sum from the nonlinear Pearson shape parameters.
            if any(isfinite(zlp_norm) & zlp_norm > 0)
                psf = interp1(ene(:), max(zlp_norm(:), 0), lag(:) ./ kscale, 'linear', 0);
            else
                psf = EELS_PearsonSparse_FullFit.signal_kernel(lag(:), kernel_params, settings);
            end
            psf = max(psf(:), 0);
            psf = psf ./ max(sum(psf, 'omitnan'), eps);
            psf_flip = flipud(psf(:));

            latent = max(positive_input(:), floor_val);
            latent = latent .* taper(:);
            smooth_sigma = max(settings.rl_deconv_latent_smooth_sigma_meV, 0);
            if smooth_sigma > 0
                shalf = max(1, ceil(4 .* smooth_sigma ./ max(de, eps)));
                sx = (-shalf:shalf).' .* de;
                sg = exp(-0.5 .* (sx ./ max(smooth_sigma, eps)) .^ 2);
                sg = sg ./ max(sum(sg, 'omitnan'), eps);
            else
                sg = [];
            end
            nit = max(1, round(settings.rl_deconv_iterations));
            eps_denom = max(floor_val, eps);
            for ii = 1:nit
                estimate = conv(latent, psf, 'same');
                estimate = max(estimate, eps_denom);
                ratio = positive_input(:) ./ estimate(:);
                ratio(~isfinite(ratio)) = 0;
                corr = conv(ratio, psf_flip, 'same');
                corr(~isfinite(corr) | corr < 0) = 0;
                latent = latent .* corr;
                latent(~isfinite(latent) | latent < 0) = 0;
                if ~isempty(sg)
                    latent = conv(latent, sg, 'same');
                end
                latent = latent .* taper(:);
            end
            recon = conv(latent, psf, 'same');
            final_scale = 1;
            if settings.rl_deconv_final_scale
                wtmp = double(fit_mask);
                num = sum(wtmp(:) .* recon(:) .* target(:), 'omitnan');
                den = sum(wtmp(:) .* recon(:) .^ 2, 'omitnan');
                final_scale = max(num ./ max(den, eps), 0);
                if isfinite(final_scale) && final_scale > 0
                    latent = latent .* final_scale;
                    recon = recon .* final_scale;
                else
                    final_scale = 1;
                end
            end
            residual = target(:) - recon(:);

            scale = prctile(abs(target(fit_mask)), 85);
            if ~isfinite(scale) || scale <= 0
                scale = max(abs(target(fit_mask)), [], 'omitnan');
            end
            if ~isfinite(scale) || scale <= 0
                scale = 1;
            end
            amp_weight = min(1, abs(target(:)) ./ max(scale, eps));
            sqrtW = sqrt(settings.rl_deconv_data_weight_floor + ...
                (1 - settings.rl_deconv_data_weight_floor) .* amp_weight);
            sqrtW(~isfinite(sqrtW)) = 0;
            w = sqrtW(fit_mask) .^ 2;
            sse = sum(w(:) .* residual(fit_mask) .^ 2, 'omitnan');
            sst = sum(w(:) .* target(fit_mask) .^ 2, 'omitnan');
            R2 = 1 - sse ./ max(sst, eps);
            nrmse = sqrt(sse ./ max(sum(w, 'omitnan'), eps)) ./ ...
                max(sqrt(sum(w(:) .* target(fit_mask) .^ 2, 'omitnan') ./ max(sum(w, 'omitnan'), eps)), eps);
            pres = positive_input(:) - recon(:);
            psse = sum(w(:) .* pres(fit_mask) .^ 2, 'omitnan');
            psst = sum(w(:) .* positive_input(fit_mask) .^ 2, 'omitnan');
            pR2 = 1 - psse ./ max(psst, eps);
            pnrmse = sqrt(psse ./ max(sum(w, 'omitnan'), eps)) ./ ...
                max(sqrt(sum(w(:) .* positive_input(fit_mask) .^ 2, 'omitnan') ./ max(sum(w, 'omitnan'), eps)), eps);

            EM = [];
            pm_s = [];
            pm_a = [];
            if isfield(res, 'stage2') && isstruct(res.stage2) && isfield(res.stage2, 'status') && ...
                    strcmp(res.stage2.status, 'ok') && isfield(res.stage2, 'EM_final') && ...
                    isfield(res.stage2, 'PM_stokes') && isfield(res.stage2, 'PM_antistokes')
                EM = res.stage2.EM_final(:);
                pm_s = res.stage2.PM_stokes(:);
                pm_a = res.stage2.PM_antistokes(:);
            elseif isfield(res, 'EM_final') && isfield(res, 'PM_stokes') && isfield(res, 'PM_antistokes')
                EM = res.EM_final(:);
                pm_s = res.PM_stokes(:);
                pm_a = res.PM_antistokes(:);
            end
            latent_s = [];
            latent_a = [];
            ratio_line = struct('available', false);
            selected_EM = [];
            if ~isempty(EM) && numel(pm_s) == numel(EM) && numel(pm_a) == numel(EM)
                good = isfinite(EM) & isfinite(pm_s) & isfinite(pm_a) & EM(:) >= e_lo & EM(:) <= e_hi;
                EMg = EM(good);
                pm_mass = max(max(pm_s(good), 0), max(pm_a(good), 0));
                if ~isempty(EMg) && any(pm_mass > 0)
                    keep = pm_mass >= 0.02 .* max(pm_mass, [], 'omitnan');
                    selected_EM = sort(unique(EMg(keep)));
                    if numel(selected_EM) > 60
                        selected_EM = EELS_PearsonSparse_FullFit.select_spaced_em( ...
                            selected_EM(:).', EELS_PearsonSparse_FullFit.interpolate_pm(EMg(:), pm_mass(:), selected_EM(:)).', ...
                            60, 0.8);
                    end
                    latent_s = interp1(ene(:), latent(:), selected_EM(:), 'linear', 0);
                    latent_a = interp1(ene(:), latent(:), -selected_EM(:), 'linear', 0);
                    latent_s = max(latent_s(:), 0);
                    latent_a = max(latent_a(:), 0);
                    ratio_line = EELS_PearsonSparse_FullFit.best_ratio_origin_line( ...
                        selected_EM(:), log(max(latent_s(:), eps) ./ max(latent_a(:), eps)), ...
                        sqrt(max(latent_s(:), 0) .* max(latent_a(:), 0)), settings);
                end
            end

            rl.available = true;
            rl.reason = 'ok';
            rl.ene = ene(:);
            rl.target_signal = target(:);
            rl.positive_input = positive_input(:);
            rl.latent_signal = latent(:);
            rl.reconvolved_signal = recon(:);
            rl.residual = residual(:);
            rl.zlp_norm = zlp_norm(:);
            rl.bg_norm = bg_norm(:);
            rl.reference_vib_norm = reference_vib(:);
            rl.psf = psf(:);
            rl.psf_energy_meV = lag(:);
            rl.fit_mask = fit_mask(:);
            rl.taper = taper(:);
            rl.sqrtW = sqrtW(:);
            rl.R2_weighted = R2;
            rl.NRMSE_weighted = nrmse;
            rl.R2_positive_weighted = pR2;
            rl.NRMSE_positive_weighted = pnrmse;
            rl.iterations = nit;
            rl.kernel_width_scale = kscale;
            rl.latent_smooth_sigma_meV = smooth_sigma;
            rl.final_scale = final_scale;
            rl.EM_meV = selected_EM(:);
            rl.latent_stokes = latent_s(:);
            rl.latent_antistokes = latent_a(:);
            if isfield(ratio_line, 'available') && ratio_line.available
                rl.latent_ratio_T_K = ratio_line.T_K;
                rl.latent_ratio_R2 = ratio_line.R2_origin;
                rl.latent_ratio_num_points = ratio_line.num_points;
                rl.latent_ratio_EM_range_meV = ratio_line.EM_range_meV;
            else
                rl.latent_ratio_EM_range_meV = [NaN, NaN];
            end
        end

        function deconv = pm_guided_deconvolution(res, settings)
            if nargin < 2 || isempty(settings)
                if isstruct(res) && isfield(res, 'settings') && isstruct(res.settings)
                    settings = res.settings;
                else
                    settings = struct();
                end
            end
            settings = EELS_PearsonSparse_FullFit.defaults(settings);
            deconv = struct('available', false, 'reason', 'not_started', ...
                'source', 'none', 'ene', [], 'target_signal', [], 'reconvolved_signal', [], ...
                'residual', [], 'zlp_norm', [], 'bg_norm', [], 'reference_vib_norm', [], ...
                'EM_meV', [], 'PM_stokes_prior', [], 'PM_antistokes_prior', [], ...
                'coef_stokes', [], 'coef_antistokes', [], 'confidence', [], ...
                'fit_mask', [], 'sqrtW', [], 'R2_weighted', NaN, 'NRMSE_weighted', NaN, ...
                'num_centers', 0, 'temperature_prior_K', NaN, 'coef_ratio_T_K', NaN, ...
                'coef_ratio_R2', NaN, 'coef_ratio_num_points', 0, 'solver', 'none', ...
                'kernel_width_scale', 1);
            if ~(isstruct(res) && isfield(res, 'ene_corrected'))
                deconv.reason = 'missing_result_energy_axis';
                return;
            end
            ene = res.ene_corrected(:);
            n = numel(ene);
            if isfield(res, 'raw_data_norm') && numel(res.raw_data_norm) >= n
                raw_norm = res.raw_data_norm(:);
                raw_norm = raw_norm(1:n);
            elseif isfield(res, 'raw_data') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                raw_norm = res.raw_data(:) ./ res.norm_val;
                raw_norm = raw_norm(1:n);
            else
                deconv.reason = 'missing_raw_data';
                return;
            end
            if isfield(res, 'stage1') && isstruct(res.stage1) && isfield(res.stage1, 'params')
                params = res.stage1.params;
            elseif isfield(res, 'zlp_params')
                params = res.zlp_params;
            else
                deconv.reason = 'missing_zlp_kernel_params';
                return;
            end

            zlp_norm = zeros(n, 1);
            bg_norm = zeros(n, 1);
            reference_vib = zeros(n, 1);
            EM = [];
            pm_s = [];
            pm_a = [];
            source = 'stage1_PM';
            stage2_ok = isfield(res, 'stage2') && isstruct(res.stage2) && ...
                isfield(res.stage2, 'status') && strcmp(res.stage2.status, 'ok');
            if stage2_ok
                if isfield(res.stage2, 'zlp_norm') && numel(res.stage2.zlp_norm) >= n
                    tmp = res.stage2.zlp_norm(:); zlp_norm = tmp(1:n);
                end
                if isfield(res.stage2, 'bg_norm') && numel(res.stage2.bg_norm) >= n
                    tmp = res.stage2.bg_norm(:); bg_norm = tmp(1:n);
                end
                if isfield(res.stage2, 'vib_norm') && numel(res.stage2.vib_norm) >= n
                    tmp = res.stage2.vib_norm(:); reference_vib = tmp(1:n);
                end
                if isfield(res.stage2, 'EM_final') && isfield(res.stage2, 'PM_stokes') && ...
                        isfield(res.stage2, 'PM_antistokes')
                    EM = res.stage2.EM_final(:);
                    pm_s = res.stage2.PM_stokes(:);
                    pm_a = res.stage2.PM_antistokes(:);
                    source = 'stage2_PM';
                end
            end
            if isempty(EM) && isfield(res, 'EM_final') && isfield(res, 'PM_stokes') && isfield(res, 'PM_antistokes')
                EM = res.EM_final(:);
                pm_s = res.PM_stokes(:);
                pm_a = res.PM_antistokes(:);
            end
            if ~any(isfinite(zlp_norm) & abs(zlp_norm) > 0)
                if isfield(res, 'stage1') && isfield(res.stage1, 'zlp_norm') && numel(res.stage1.zlp_norm) >= n
                    tmp = res.stage1.zlp_norm(:); zlp_norm = tmp(1:n);
                elseif isfield(res, 'zlp_only') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                    tmp = res.zlp_only(:); zlp_norm = tmp(1:n) ./ res.norm_val;
                end
            end
            if ~any(isfinite(bg_norm) & abs(bg_norm) > 0)
                if isfield(res, 'stage1') && isfield(res.stage1, 'bg_norm') && numel(res.stage1.bg_norm) >= n
                    tmp = res.stage1.bg_norm(:); bg_norm = tmp(1:n);
                elseif isfield(res, 'bg_only') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                    tmp = res.bg_only(:); bg_norm = tmp(1:n) ./ res.norm_val;
                end
            end
            if ~any(isfinite(reference_vib) & abs(reference_vib) > 0)
                if isfield(res, 'stage1') && isfield(res.stage1, 'vib_norm') && numel(res.stage1.vib_norm) >= n
                    tmp = res.stage1.vib_norm(:); reference_vib = tmp(1:n);
                elseif isfield(res, 'vib_only') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                    tmp = res.vib_only(:); reference_vib = tmp(1:n) ./ res.norm_val;
                end
            end
            if isempty(EM) || numel(pm_s) ~= numel(EM) || numel(pm_a) ~= numel(EM)
                deconv.reason = 'missing_PM_dictionary';
                return;
            end
            target = raw_norm(:) - zlp_norm(:) - bg_norm(:);
            em_lo = max(0, settings.pm_deconv_min_em_meV);
            em_hi = max(em_lo + eps, settings.pm_deconv_max_em_meV);
            good = isfinite(EM) & isfinite(pm_s) & isfinite(pm_a) & EM(:) >= em_lo & EM(:) <= em_hi;
            EM = EM(good);
            pm_s = max(pm_s(good), 0);
            pm_a = max(pm_a(good), 0);
            if isempty(EM)
                deconv.reason = 'no_PM_centers_in_window';
                return;
            end
            pm_mass = max(pm_s(:), pm_a(:));
            max_mass = max(pm_mass, [], 'omitnan');
            if ~isfinite(max_mass) || max_mass <= 0
                deconv.reason = 'zero_PM_prior';
                return;
            end
            keep = pm_mass >= settings.pm_deconv_min_pm_fraction .* max_mass;
            if nnz(keep) < 2
                [~, order] = sort(pm_mass, 'descend');
                keep(order(1:min(numel(order), 2))) = true;
            end
            cand_EM = EM(keep);
            cand_score = pm_mass(keep);
            max_centers = max(2, round(settings.pm_deconv_max_centers));
            if numel(cand_EM) > max_centers
                selected_EM = EELS_PearsonSparse_FullFit.select_spaced_em(cand_EM(:).', cand_score(:).', ...
                    max_centers, settings.pm_deconv_min_spacing_meV);
            else
                selected_EM = unique(cand_EM(:).');
            end
            selected_EM = selected_EM(isfinite(selected_EM) & selected_EM >= em_lo & selected_EM <= em_hi);
            selected_EM = sort(selected_EM(:).');
            if numel(selected_EM) < 2
                deconv.reason = 'insufficient_selected_PM_centers';
                return;
            end
            pm_s_sel = EELS_PearsonSparse_FullFit.interpolate_pm(EM(:), pm_s(:), selected_EM(:));
            pm_a_sel = EELS_PearsonSparse_FullFit.interpolate_pm(EM(:), pm_a(:), selected_EM(:));
            pm_s_sel = max(pm_s_sel(:), 0);
            pm_a_sel = max(pm_a_sel(:), 0);
            kernel_width_scale = settings.pm_deconv_kernel_width_scale;
            if ~isscalar(kernel_width_scale) || ~isfinite(kernel_width_scale) || kernel_width_scale <= 0
                kernel_width_scale = 1;
            end
            kernel_params = params;
            if abs(kernel_width_scale - 1) > 1e-12
                if isfield(kernel_params, 'sL')
                    kernel_params.sL = max(kernel_params.sL .* kernel_width_scale, eps);
                end
                if isfield(kernel_params, 'sR')
                    kernel_params.sR = max(kernel_params.sR .* kernel_width_scale, eps);
                end
            end
            [~, ~, ~, kernels] = EELS_PearsonSparse_FullFit.build_design_matrix(ene(:), selected_EM(:).', kernel_params, settings);
            C = [kernels.stokes, kernels.anti];
            if isempty(C) || size(C, 1) ~= n
                deconv.reason = 'empty_kernel_dictionary';
                return;
            end
            fit_mask = abs(ene(:)) >= em_lo & abs(ene(:)) <= em_hi & isfinite(target(:));
            fit_mask = fit_mask & any(isfinite(C) & abs(C) > 0, 2);
            if nnz(fit_mask) < max(8, numel(selected_EM))
                deconv.reason = 'insufficient_fit_samples';
                return;
            end
            scale = prctile(abs(target(fit_mask)), 85);
            if ~isfinite(scale) || scale <= 0
                scale = max(abs(target(fit_mask)), [], 'omitnan');
            end
            if ~isfinite(scale) || scale <= 0
                scale = 1;
            end
            amp_weight = min(1, abs(target(:)) ./ max(scale, eps));
            sqrtW = sqrt(settings.pm_deconv_data_weight_floor + ...
                (1 - settings.pm_deconv_data_weight_floor) .* amp_weight);
            if stage2_ok && isfield(res.stage2, 'attention') && isfield(res.stage2.attention, 'feature_score') && ...
                    numel(res.stage2.attention.feature_score) >= n
                att = res.stage2.attention.feature_score(:);
                sqrtW = sqrtW .* sqrt(0.45 + 0.55 .* max(min(att(1:n), 1), 0));
            end
            sqrtW(~isfinite(sqrtW)) = 0;
            Cw = bsxfun(@times, C(fit_mask, :), sqrtW(fit_mask));
            yw = target(fit_mask) .* sqrtW(fit_mask);
            nEM = numel(selected_EM);
            ncoef = 2 .* nEM;
            prior = [pm_s_sel(:); pm_a_sel(:)];
            conf = max(pm_s_sel(:), pm_a_sel(:)) ./ max(max(pm_s_sel(:), pm_a_sel(:)), [], 'omitnan');
            conf(~isfinite(conf)) = 0;
            coef_conf = [conf(:); conf(:)];
            Reg = zeros(0, ncoef);
            rhs = zeros(0, 1);
            if settings.pm_deconv_pm_prior_weight > 0
                wp = sqrt(settings.pm_deconv_pm_prior_weight .* (0.20 + 0.80 .* coef_conf(:)));
                Reg = [Reg; diag(wp)]; %#ok<AGROW>
                rhs = [rhs; wp(:) .* prior(:)]; %#ok<AGROW>
            end
            if settings.pm_deconv_l2_weight > 0
                wl2 = sqrt(settings.pm_deconv_l2_weight .* (1.25 - min(coef_conf(:), 1)));
                Reg = [Reg; diag(wl2)]; %#ok<AGROW>
                rhs = [rhs; zeros(ncoef, 1)]; %#ok<AGROW>
            end
            T_prior = NaN;
            if stage2_ok && isfield(res.stage2, 'T_fit') && isfinite(res.stage2.T_fit)
                T_prior = res.stage2.T_fit;
            elseif isfield(res, 'T_fit') && isfinite(res.T_fit)
                T_prior = res.T_fit;
            elseif isfield(res, 'T_soft') && isfinite(res.T_soft)
                T_prior = res.T_soft;
            end
            if settings.pm_deconv_be_prior_weight > 0 && isfinite(T_prior) && T_prior > 0
                kb = 0.08617333262;
                q = exp(-selected_EM(:) ./ max(kb .* T_prior, eps));
                B = zeros(nEM, ncoef);
                for ii = 1:nEM
                    B(ii, ii) = -q(ii);
                    B(ii, nEM + ii) = 1;
                end
                wb = sqrt(settings.pm_deconv_be_prior_weight .* (0.20 + 0.80 .* conf(:)));
                Reg = [Reg; bsxfun(@times, B, wb(:))]; %#ok<AGROW>
                rhs = [rhs; zeros(nEM, 1)]; %#ok<AGROW>
            end
            Aaug = [Cw; Reg];
            baug = [yw(:); rhs(:)];
            try
                coef = lsqnonneg(Aaug, baug);
                solver = 'lsqnonneg';
            catch
                coef = max(Aaug \ baug, 0);
                solver = 'backslash_nonnegative_clip';
            end
            coef(~isfinite(coef)) = 0;
            recon = C * coef(:);
            residual = target(:) - recon(:);
            w = sqrtW(fit_mask) .^ 2;
            sse = sum(w(:) .* residual(fit_mask) .^ 2, 'omitnan');
            sst = sum(w(:) .* target(fit_mask) .^ 2, 'omitnan');
            R2 = 1 - sse ./ max(sst, eps);
            nrmse = sqrt(sse ./ max(sum(w, 'omitnan'), eps)) ./ max(sqrt(sum(w(:) .* target(fit_mask) .^ 2, 'omitnan') ./ max(sum(w, 'omitnan'), eps)), eps);
            coef_s = coef(1:nEM);
            coef_a = coef(nEM + (1:nEM));
            ratio_line = EELS_PearsonSparse_FullFit.best_ratio_origin_line( ...
                selected_EM(:), log(max(coef_s(:), eps) ./ max(coef_a(:), eps)), ...
                sqrt(max(coef_s(:), 0) .* max(coef_a(:), 0)), settings);

            deconv.available = true;
            deconv.reason = 'ok';
            deconv.source = source;
            deconv.ene = ene(:);
            deconv.target_signal = target(:);
            deconv.reconvolved_signal = recon(:);
            deconv.residual = residual(:);
            deconv.zlp_norm = zlp_norm(:);
            deconv.bg_norm = bg_norm(:);
            deconv.reference_vib_norm = reference_vib(:);
            deconv.EM_meV = selected_EM(:);
            deconv.PM_stokes_prior = pm_s_sel(:);
            deconv.PM_antistokes_prior = pm_a_sel(:);
            deconv.coef_stokes = coef_s(:);
            deconv.coef_antistokes = coef_a(:);
            deconv.confidence = conf(:);
            deconv.fit_mask = fit_mask(:);
            deconv.sqrtW = sqrtW(:);
            deconv.R2_weighted = R2;
            deconv.NRMSE_weighted = nrmse;
            deconv.num_centers = nEM;
            deconv.temperature_prior_K = T_prior;
            deconv.solver = solver;
            deconv.kernel_width_scale = kernel_width_scale;
            if ratio_line.available
                deconv.coef_ratio_T_K = ratio_line.T_K;
                deconv.coef_ratio_R2 = ratio_line.R2_origin;
                deconv.coef_ratio_num_points = ratio_line.num_points;
                deconv.coef_ratio_EM_range_meV = ratio_line.EM_range_meV;
            else
                deconv.coef_ratio_EM_range_meV = [NaN, NaN];
            end
        end

        function fig = plot_rl_deconvolution_diagnostic(res, title_text, settings)
            if nargin < 2 || isempty(title_text), title_text = 'RL deconvolution of background-subtracted spectrum'; end
            if nargin < 3, settings = []; end
            rl = EELS_PearsonSparse_FullFit.rl_deconvolve_background_subtracted(res, settings);
            fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 60, 1450, 1050]);
            tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            if ~rl.available
                text(0.05, 0.55, sprintf('RL deconvolution unavailable: %s', rl.reason), ...
                    'Units', 'normalized', 'Interpreter', 'none');
                axis off;
                sgtitle(title_text, 'Interpreter', 'none');
                return;
            end
            ene = rl.ene(:);
            nexttile;
            plot(ene, rl.target_signal, 'k-', 'LineWidth', 0.9); hold on;
            plot(ene, rl.positive_input, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.9);
            plot(ene, rl.reconvolved_signal, 'Color', [0.55 0.10 0.78], 'LineWidth', 1.35);
            plot(ene, rl.reference_vib_norm, '--', 'Color', [0.25 0.45 0.90], 'LineWidth', 0.9);
            yline(0, ':', 'Color', [0.5 0.5 0.5]);
            title(sprintf('Target vs RL reconvolved | R2 %.3f, positive R2 %.3f', ...
                rl.R2_weighted, rl.R2_positive_weighted), 'Interpreter', 'none');
            xlabel('Energy loss (meV)'); ylabel('norm intensity'); grid on;
            legend({'raw - ZLP - bg', 'RL nonnegative input', 'RL reconvolved', 'current vib'}, 'Location', 'best');

            nexttile;
            plot(ene, rl.residual, 'Color', [0.25 0.25 0.25], 'LineWidth', 0.85); hold on;
            yline(0, ':', 'Color', [0.5 0.5 0.5]);
            if any(rl.fit_mask)
                yl = ylim;
                patch([min(ene(rl.fit_mask)) max(ene(rl.fit_mask)) max(ene(rl.fit_mask)) min(ene(rl.fit_mask))], ...
                    [yl(1) yl(1) yl(2) yl(2)], [0.95 0.90 1.0], 'FaceAlpha', 0.16, 'EdgeColor', 'none');
                plot(ene, rl.residual, 'Color', [0.25 0.25 0.25], 'LineWidth', 0.85);
            end
            title(sprintf('Residual after RL reconvolution | NRMSE %.3g', rl.NRMSE_weighted), 'Interpreter', 'none');
            xlabel('Energy loss (meV)'); ylabel('target - recon'); grid on;

            nexttile;
            plot(ene, rl.latent_signal, 'Color', [0.55 0.10 0.78], 'LineWidth', 1.1); hold on;
            plot(ene, rl.positive_input, ':', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.85);
            yline(0, ':', 'Color', [0.5 0.5 0.5]);
            title(sprintf('RL latent spectrum | iter %d | smooth %.2g meV', ...
                rl.iterations, rl.latent_smooth_sigma_meV), 'Interpreter', 'none');
            xlabel('Energy loss (meV)'); ylabel('latent amplitude'); grid on;
            legend({'RL latent', 'positive input'}, 'Location', 'best');

            nexttile;
            plot(rl.psf_energy_meV, rl.psf, 'Color', [0.35 0.05 0.55], 'LineWidth', 1.15);
            xlim([-80, 80]);
            title(sprintf('Fitted ZLP PSF for RL | width scale %.2g', rl.kernel_width_scale), 'Interpreter', 'none');
            xlabel('Energy lag (meV)'); ylabel('normalized PSF'); grid on;

            nexttile;
            if numel(rl.EM_meV) >= 2 && numel(rl.latent_stokes) == numel(rl.EM_meV)
                stem(rl.EM_meV, rl.latent_stokes, 'Color', [0.85 0.25 0.05], ...
                    'Marker', 'none', 'LineWidth', 0.95); hold on;
                stem(-rl.EM_meV, rl.latent_antistokes, 'Color', [0.05 0.25 0.85], ...
                    'Marker', 'none', 'LineWidth', 0.95);
                title(sprintf('RL latent sampled at sparse mode energies | BE T %.1f K, R2 %.3f', ...
                    rl.latent_ratio_T_K, rl.latent_ratio_R2), 'Interpreter', 'none');
                xlabel('EM (meV), gain side negative'); ylabel('latent amplitude'); grid on;
                legend({'loss sampled', 'gain sampled'}, 'Location', 'best');
            else
                text(0.05, 0.55, 'No paired sparse-mode samples for RL BE ratio', ...
                    'Units', 'normalized', 'Interpreter', 'none');
                axis off;
            end

            nexttile;
            lines = {
                sprintf('method: Richardson-Lucy')
                sprintf('input: max(raw - ZLP - bg, 0)')
                sprintf('iterations: %d', rl.iterations)
                sprintf('kernel width scale: %.3g', rl.kernel_width_scale)
                sprintf('latent smooth sigma: %.3g meV', rl.latent_smooth_sigma_meV)
                sprintf('final scale: %.5g', rl.final_scale)
                sprintf('R2 weighted: %.5g', rl.R2_weighted)
                sprintf('NRMSE weighted: %.5g', rl.NRMSE_weighted)
                sprintf('positive-input R2: %.5g', rl.R2_positive_weighted)
                sprintf('positive-input NRMSE: %.5g', rl.NRMSE_positive_weighted)
                sprintf('latent BE T: %.1f K', rl.latent_ratio_T_K)
                sprintf('latent BE R2: %.4g', rl.latent_ratio_R2)
                };
            text(0.02, 0.92, strjoin(lines, newline), 'Units', 'normalized', ...
                'Interpreter', 'none', 'VerticalAlignment', 'top', 'FontName', 'Consolas');
            axis off;
            sgtitle(title_text, 'Interpreter', 'none');
        end

        function fig = plot_rl_deconvolution_datacube_summary(cuberes, title_text, settings)
            if nargin < 2 || isempty(title_text), title_text = 'RL deconvolution datacube smoke'; end
            if nargin < 3, settings = []; end
            if ~(isstruct(cuberes) && isfield(cuberes, 'results') && isfield(cuberes, 'sample_pixels'))
                error('plot_rl_deconvolution_datacube_summary expects a datacube result struct.');
            end
            valid = find(cellfun(@(r) isstruct(r) && ~isfield(r, 'error'), cuberes.results));
            valid = valid(1:min(numel(valid), 6));
            nrow = max(1, numel(valid));
            fig = figure('Visible', 'off', 'Color', 'w', 'Position', [70, 50, 1500, 260 .* nrow + 120]);
            tiledlayout(nrow, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            for rr = 1:nrow
                if rr > numel(valid)
                    nexttile; axis off; nexttile; axis off; continue;
                end
                ii = valid(rr);
                res = cuberes.results{ii};
                rl = EELS_PearsonSparse_FullFit.rl_deconvolve_background_subtracted(res, settings);
                pix = cuberes.sample_pixels(ii);
                if isfield(cuberes, 'R2_map')
                    [yy, xx] = ind2sub(size(cuberes.R2_map), pix);
                    pix_label = sprintf('pixel %d (%d,%d)', pix, yy, xx);
                else
                    pix_label = sprintf('pixel %d', pix);
                end
                nexttile;
                if rl.available
                    ene = rl.ene(:);
                    mask = abs(ene) <= max(abs(ene(rl.fit_mask)), [], 'omitnan') + 35;
                    if nnz(mask) < 10, mask = true(size(ene)); end
                    plot(ene(mask), rl.target_signal(mask), 'k-', 'LineWidth', 0.8); hold on;
                    plot(ene(mask), rl.reconvolved_signal(mask), 'Color', [0.55 0.10 0.78], 'LineWidth', 1.2);
                    yline(0, ':', 'Color', [0.5 0.5 0.5]);
                    title(sprintf('%s target/RL reconv | R2 %.3f, NRMSE %.3g', ...
                        pix_label, rl.R2_weighted, rl.NRMSE_weighted), 'Interpreter', 'none');
                    xlabel('Energy loss (meV)'); ylabel('norm intensity'); grid on;
                else
                    text(0.05, 0.55, sprintf('%s unavailable: %s', pix_label, rl.reason), ...
                        'Units', 'normalized', 'Interpreter', 'none');
                    axis off;
                end
                nexttile;
                if rl.available
                    ene = rl.ene(:);
                    plot(ene, rl.latent_signal, 'Color', [0.55 0.10 0.78], 'LineWidth', 1.0); hold on;
                    if numel(rl.EM_meV) >= 2
                        stem(rl.EM_meV, rl.latent_stokes, 'Color', [0.85 0.25 0.05], ...
                            'Marker', 'none', 'LineWidth', 0.75);
                        stem(-rl.EM_meV, rl.latent_antistokes, 'Color', [0.05 0.25 0.85], ...
                            'Marker', 'none', 'LineWidth', 0.75);
                    end
                    title(sprintf('RL latent | iter %d | BE T %.1f K', ...
                        rl.iterations, rl.latent_ratio_T_K), 'Interpreter', 'none');
                    xlabel('Energy loss / EM (meV)'); ylabel('latent amplitude'); grid on;
                else
                    axis off;
                end
            end
            sgtitle(title_text, 'Interpreter', 'none');
        end

        function fig = plot_pm_top_interpolation_diagnostic(res, title_text, settings)
            if nargin < 2 || isempty(title_text), title_text = 'Sparse source-intensity interpolation spectrum'; end
            if nargin < 3, settings = []; end
            pmview = EELS_PearsonSparse_FullFit.pm_top_interpolated_spectrum(res, settings);
            fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 60, 1450, 1050]);
            tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            if ~pmview.available
                text(0.05, 0.55, sprintf('Sparse source-intensity interpolation unavailable: %s', pmview.reason), ...
                    'Units', 'normalized', 'Interpreter', 'none');
                axis off;
                sgtitle(title_text, 'Interpreter', 'none');
                return;
            end

            ene = pmview.ene(:);
            nexttile;
            plot(ene, pmview.target_signal, 'k-', 'LineWidth', 0.9); hold on;
            plot(ene, pmview.convolved_signal, 'Color', [0.05 0.55 0.20], 'LineWidth', 1.35);
            plot(ene, pmview.reference_vib_norm, '--', 'Color', [0.25 0.45 0.90], 'LineWidth', 0.95);
            yline(0, ':', 'Color', [0.5 0.5 0.5]);
            title(sprintf('Target vs sparse source-intensity small-kernel spectrum | R2 %.3f', pmview.R2_weighted), 'Interpreter', 'none');
            xlabel('Energy loss (meV)'); ylabel('norm intensity'); grid on;
            legend({'raw - ZLP - bg', 'source-intensity convolution', 'current vib'}, 'Location', 'best');

            nexttile;
            plot(ene, pmview.residual, 'Color', [0.25 0.25 0.25], 'LineWidth', 0.85); hold on;
            yline(0, ':', 'Color', [0.5 0.5 0.5]);
            if any(pmview.fit_mask)
                yl = ylim;
                patch([min(ene(pmview.fit_mask)) max(ene(pmview.fit_mask)) max(ene(pmview.fit_mask)) min(ene(pmview.fit_mask))], ...
                    [yl(1) yl(1) yl(2) yl(2)], [0.9 1.0 0.9], 'FaceAlpha', 0.16, 'EdgeColor', 'none');
                plot(ene, pmview.residual, 'Color', [0.25 0.25 0.25], 'LineWidth', 0.85);
            end
            title(sprintf('Residual after source-intensity interpolation | NRMSE %.3g', pmview.NRMSE_weighted), 'Interpreter', 'none');
            xlabel('Energy loss (meV)'); ylabel('target - PMtop'); grid on;

            nexttile;
            plot(ene, pmview.latent_signal, 'Color', [0.15 0.65 0.20], 'LineWidth', 1.10); hold on;
            plot(ene, pmview.convolved_signal, 'Color', [0.05 0.35 0.80], 'LineWidth', 0.95);
            yline(0, ':', 'Color', [0.5 0.5 0.5]);
            title(sprintf('Latent source intensity and %.2g meV Gaussian convolution', pmview.peak_sigma_meV), ...
                'Interpreter', 'none');
            xlabel('Energy loss (meV)'); ylabel('scaled amplitude'); grid on;
            legend({'latent interp', 'small-kernel conv'}, 'Location', 'best');

            nexttile;
            stem(pmview.EM_meV, pmview.PM_stokes_top, 'Color', [0.85 0.25 0.05], ...
                'Marker', 'none', 'LineWidth', 1.0); hold on;
            stem(-pmview.EM_meV, pmview.PM_antistokes_top, 'Color', [0.05 0.25 0.85], ...
                'Marker', 'none', 'LineWidth', 1.0);
            title(sprintf('Selected sparse source intensities, n=%d', pmview.num_centers), 'Interpreter', 'none');
            xlabel('Mode energy (meV), gain side negative'); ylabel('Source intensity'); grid on;
            legend({'loss/Stokes intensity', 'gain/anti-Stokes intensity'}, 'Location', 'best');

            nexttile;
            good_ratio = pmview.PM_stokes_top > 0 & pmview.PM_antistokes_top > 0;
            if nnz(good_ratio) >= 2
                y = log(max(pmview.PM_stokes_top(good_ratio), eps) ./ max(pmview.PM_antistokes_top(good_ratio), eps));
                w = sqrt(max(pmview.PM_stokes_top(good_ratio), 0) .* max(pmview.PM_antistokes_top(good_ratio), 0));
                scatter(pmview.EM_meV(good_ratio), y, 28, w, 'filled'); hold on;
                if isfinite(pmview.top_ratio_T_K)
                    kb = 0.08617333262;
                    xx = linspace(0, max(pmview.EM_meV(good_ratio), [], 'omitnan'), 120);
                    plot(xx, xx ./ max(kb .* pmview.top_ratio_T_K, eps), 'k-', 'LineWidth', 1.2);
                end
                colorbar;
                title(sprintf('Source-intensity BE ratio: T %.1f K, R2 %.3f, n=%d', ...
                    pmview.top_ratio_T_K, pmview.top_ratio_R2, pmview.top_ratio_num_points), 'Interpreter', 'none');
                xlabel('EM (meV)'); ylabel('ln(loss/gain)'); grid on;
            else
                text(0.05, 0.55, 'Too few paired source intensities for BE ratio', 'Units', 'normalized', 'Interpreter', 'none');
                axis off;
            end

            nexttile;
            lines = {
                sprintf('source: %s', pmview.source)
                sprintf('interp method: %s', pmview.interp_method)
                sprintf('peak sigma: %.3g meV', pmview.peak_sigma_meV)
                sprintf('scale mode: %s', pmview.scale_mode)
                sprintf('scale global: %.5g', pmview.scale_global)
                sprintf('scale loss/gain: %.5g / %.5g', pmview.scale_loss, pmview.scale_gain)
                sprintf('centers: %d', pmview.num_centers)
                sprintf('R2 weighted: %.5g', pmview.R2_weighted)
                sprintf('NRMSE weighted: %.5g', pmview.NRMSE_weighted)
                sprintf('Source-intensity BE T: %.1f K', pmview.top_ratio_T_K)
                sprintf('Source-intensity BE R2: %.4g', pmview.top_ratio_R2)
                };
            text(0.02, 0.92, strjoin(lines, newline), 'Units', 'normalized', ...
                'Interpreter', 'none', 'VerticalAlignment', 'top', 'FontName', 'Consolas');
            axis off;
            sgtitle(title_text, 'Interpreter', 'none');
        end

        function fig = plot_pm_top_interpolation_datacube_summary(cuberes, title_text, settings)
            if nargin < 2 || isempty(title_text), title_text = 'Sparse source-intensity interpolation datacube smoke'; end
            if nargin < 3, settings = []; end
            if ~(isstruct(cuberes) && isfield(cuberes, 'results') && isfield(cuberes, 'sample_pixels'))
                error('plot_pm_top_interpolation_datacube_summary expects a datacube result struct.');
            end
            valid = find(cellfun(@(r) isstruct(r) && ~isfield(r, 'error'), cuberes.results));
            valid = valid(1:min(numel(valid), 6));
            nrow = max(1, numel(valid));
            fig = figure('Visible', 'off', 'Color', 'w', 'Position', [70, 50, 1500, 260 .* nrow + 120]);
            tiledlayout(nrow, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            for rr = 1:nrow
                if rr > numel(valid)
                    nexttile; axis off; nexttile; axis off; continue;
                end
                ii = valid(rr);
                res = cuberes.results{ii};
                pmview = EELS_PearsonSparse_FullFit.pm_top_interpolated_spectrum(res, settings);
                pix = cuberes.sample_pixels(ii);
                if isfield(cuberes, 'R2_map')
                    [yy, xx] = ind2sub(size(cuberes.R2_map), pix);
                    pix_label = sprintf('pixel %d (%d,%d)', pix, yy, xx);
                else
                    pix_label = sprintf('pixel %d', pix);
                end
                nexttile;
                if pmview.available
                    ene = pmview.ene(:);
                    mask = abs(ene) <= max(abs(pmview.EM_meV), [], 'omitnan') + 35;
                    if nnz(mask) < 10, mask = true(size(ene)); end
                    plot(ene(mask), pmview.target_signal(mask), 'k-', 'LineWidth', 0.8); hold on;
                    plot(ene(mask), pmview.convolved_signal(mask), 'Color', [0.05 0.55 0.20], 'LineWidth', 1.2);
                    yline(0, ':', 'Color', [0.5 0.5 0.5]);
                    title(sprintf('%s target/source intensity | R2 %.3f, NRMSE %.3g', ...
                        pix_label, pmview.R2_weighted, pmview.NRMSE_weighted), 'Interpreter', 'none');
                    xlabel('Energy loss (meV)'); ylabel('norm intensity'); grid on;
                else
                    text(0.05, 0.55, sprintf('%s unavailable: %s', pix_label, pmview.reason), ...
                        'Units', 'normalized', 'Interpreter', 'none');
                    axis off;
                end
                nexttile;
                if pmview.available
                    ene = pmview.ene(:);
                    plot(ene, pmview.latent_signal, 'Color', [0.15 0.65 0.20], 'LineWidth', 1.0); hold on;
                    stem(pmview.EM_meV, pmview.scale_loss .* pmview.PM_stokes_top, ...
                        'Color', [0.85 0.25 0.05], 'Marker', 'none', 'LineWidth', 0.75);
                    stem(-pmview.EM_meV, pmview.scale_gain .* pmview.PM_antistokes_top, ...
                        'Color', [0.05 0.25 0.85], 'Marker', 'none', 'LineWidth', 0.75);
                    title(sprintf('Sparse source intensity | sigma %.2g meV | BE T %.1f K', ...
                        pmview.peak_sigma_meV, pmview.top_ratio_T_K), 'Interpreter', 'none');
                    xlabel('Energy loss / EM (meV)'); ylabel('scaled amplitude'); grid on;
                else
                    axis off;
                end
            end
            sgtitle(title_text, 'Interpreter', 'none');
        end

        function fig = plot_pm_deconvolution_diagnostic(res, title_text, settings)
            if nargin < 2 || isempty(title_text), title_text = 'Source-intensity-guided deconvolution'; end
            if nargin < 3, settings = []; end
            deconv = EELS_PearsonSparse_FullFit.pm_guided_deconvolution(res, settings);
            fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 60, 1450, 1050]);
            tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            if ~deconv.available
                text(0.05, 0.55, sprintf('Source-intensity-guided deconvolution unavailable: %s', deconv.reason), ...
                    'Units', 'normalized', 'Interpreter', 'none');
                axis off;
                sgtitle(title_text, 'Interpreter', 'none');
                return;
            end
            ene = deconv.ene(:);
            nexttile;
            plot(ene, deconv.target_signal, 'k-', 'LineWidth', 0.9); hold on;
            plot(ene, deconv.reconvolved_signal, 'Color', [0.9 0.25 0.05], 'LineWidth', 1.35);
            plot(ene, deconv.reference_vib_norm, '--', 'Color', [0.25 0.45 0.90], 'LineWidth', 0.95);
            yline(0, ':', 'Color', [0.5 0.5 0.5]);
            xlim([min(ene), max(ene)]);
            title(sprintf('Target residual vs reconvolved source intensity | kernel %.2gx | R2 %.3f', ...
                deconv.kernel_width_scale, deconv.R2_weighted), 'Interpreter', 'none');
            xlabel('Energy loss (meV)'); ylabel('norm intensity'); grid on;
            legend({'raw - ZLP - bg', 'source intensity reconvolved', 'current vib'}, 'Location', 'best');

            nexttile;
            plot(ene, deconv.residual, 'Color', [0.25 0.25 0.25], 'LineWidth', 0.85); hold on;
            yline(0, ':', 'Color', [0.5 0.5 0.5]);
            if any(deconv.fit_mask)
                yl = ylim;
                patch([min(ene(deconv.fit_mask)) max(ene(deconv.fit_mask)) max(ene(deconv.fit_mask)) min(ene(deconv.fit_mask))], ...
                    [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 1.0], 'FaceAlpha', 0.16, 'EdgeColor', 'none');
                plot(ene, deconv.residual, 'Color', [0.25 0.25 0.25], 'LineWidth', 0.85);
            end
            title(sprintf('Residual after forward deconvolution | NRMSE %.3g', deconv.NRMSE_weighted), 'Interpreter', 'none');
            xlabel('Energy loss (meV)'); ylabel('target - recon'); grid on;

            nexttile;
            stem(deconv.EM_meV, deconv.coef_stokes, 'Color', [0.85 0.25 0.05], 'Marker', 'none', 'LineWidth', 1.0); hold on;
            stem(-deconv.EM_meV, deconv.coef_antistokes, 'Color', [0.05 0.25 0.85], 'Marker', 'none', 'LineWidth', 1.0);
            stem(deconv.EM_meV, deconv.PM_stokes_prior, ':', 'Color', [0.95 0.55 0.25], 'Marker', 'none', 'LineWidth', 0.8);
            stem(-deconv.EM_meV, deconv.PM_antistokes_prior, ':', 'Color', [0.25 0.55 0.95], 'Marker', 'none', 'LineWidth', 0.8);
            title(sprintf('Excitation sticks from sparse mode dictionary, n=%d', deconv.num_centers), 'Interpreter', 'none');
            xlabel('EM (meV), gain side negative'); ylabel('amplitude'); grid on;
            legend({'deconv loss', 'deconv gain', 'loss source prior', 'gain source prior'}, 'Location', 'best');

            nexttile;
            scatter(deconv.PM_stokes_prior, deconv.coef_stokes, 28, deconv.EM_meV, 'filled'); hold on;
            scatter(deconv.PM_antistokes_prior, deconv.coef_antistokes, 28, deconv.EM_meV, 'd', 'filled');
            limv = max([deconv.PM_stokes_prior(:); deconv.PM_antistokes_prior(:); deconv.coef_stokes(:); deconv.coef_antistokes(:)], [], 'omitnan');
            if ~isfinite(limv) || limv <= 0, limv = 1; end
            plot([0 limv], [0 limv], 'k--', 'LineWidth', 0.8);
            xlim([0 limv]); ylim([0 limv]);
            cb = colorbar; cb.Label.String = 'EM (meV)';
            title('Source-intensity prior vs deconvolved coefficients', 'Interpreter', 'none');
            xlabel('Source-intensity prior'); ylabel('deconv amplitude'); grid on;

            nexttile;
            good_ratio = deconv.coef_stokes > 0 & deconv.coef_antistokes > 0;
            if nnz(good_ratio) >= 2
                y = log(max(deconv.coef_stokes(good_ratio), eps) ./ max(deconv.coef_antistokes(good_ratio), eps));
                w = sqrt(max(deconv.coef_stokes(good_ratio), 0) .* max(deconv.coef_antistokes(good_ratio), 0));
                scatter(deconv.EM_meV(good_ratio), y, 28, w, 'filled'); hold on;
                if isfinite(deconv.coef_ratio_T_K)
                    kb = 0.08617333262;
                    xx = linspace(0, max(deconv.EM_meV(good_ratio), [], 'omitnan'), 120);
                    plot(xx, xx ./ max(kb .* deconv.coef_ratio_T_K, eps), 'k-', 'LineWidth', 1.2);
                end
                colorbar;
                title(sprintf('Deconv BE ratio: T %.1f K, R2 %.3f, n=%d', ...
                    deconv.coef_ratio_T_K, deconv.coef_ratio_R2, deconv.coef_ratio_num_points), 'Interpreter', 'none');
                xlabel('EM (meV)'); ylabel('ln(loss/gain)'); grid on;
            else
                text(0.05, 0.55, 'Too few paired coefficients for BE ratio', 'Units', 'normalized', 'Interpreter', 'none');
                axis off;
            end

            nexttile;
            lines = {
                sprintf('source: %s', deconv.source)
                sprintf('solver: %s', deconv.solver)
                sprintf('kernel width scale: %.3g', deconv.kernel_width_scale)
                sprintf('centers: %d', deconv.num_centers)
                sprintf('R2 weighted: %.5g', deconv.R2_weighted)
                sprintf('NRMSE weighted: %.5g', deconv.NRMSE_weighted)
                sprintf('T prior: %.1f K', deconv.temperature_prior_K)
                sprintf('coef BE T: %.1f K', deconv.coef_ratio_T_K)
                sprintf('coef BE R2: %.4g', deconv.coef_ratio_R2)
                };
            text(0.02, 0.92, strjoin(lines, newline), 'Units', 'normalized', ...
                'Interpreter', 'none', 'VerticalAlignment', 'top', 'FontName', 'Consolas');
            axis off;
            sgtitle(title_text, 'Interpreter', 'none');
        end

        function fig = plot_pm_deconvolution_datacube_summary(cuberes, title_text, settings)
            if nargin < 2 || isempty(title_text), title_text = 'Source-intensity-guided deconvolution datacube smoke'; end
            if nargin < 3, settings = []; end
            if ~(isstruct(cuberes) && isfield(cuberes, 'results') && isfield(cuberes, 'sample_pixels'))
                error('plot_pm_deconvolution_datacube_summary expects a datacube result struct.');
            end
            valid = find(cellfun(@(r) isstruct(r) && ~isfield(r, 'error'), cuberes.results));
            valid = valid(1:min(numel(valid), 6));
            nrow = max(1, numel(valid));
            fig = figure('Visible', 'off', 'Color', 'w', 'Position', [70, 50, 1500, 260 .* nrow + 120]);
            tiledlayout(nrow, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            for rr = 1:nrow
                if rr > numel(valid)
                    nexttile; axis off; nexttile; axis off; continue;
                end
                ii = valid(rr);
                res = cuberes.results{ii};
                deconv = EELS_PearsonSparse_FullFit.pm_guided_deconvolution(res, settings);
                pix = cuberes.sample_pixels(ii);
                if isfield(cuberes, 'R2_map')
                    [yy, xx] = ind2sub(size(cuberes.R2_map), pix);
                    pix_label = sprintf('pixel %d (%d,%d)', pix, yy, xx);
                else
                    pix_label = sprintf('pixel %d', pix);
                end
                nexttile;
                if deconv.available
                    ene = deconv.ene(:);
                    mask = abs(ene) <= max(abs(deconv.EM_meV), [], 'omitnan') + 35;
                    if nnz(mask) < 10, mask = true(size(ene)); end
                    plot(ene(mask), deconv.target_signal(mask), 'k-', 'LineWidth', 0.8); hold on;
                    plot(ene(mask), deconv.reconvolved_signal(mask), 'Color', [0.9 0.25 0.05], 'LineWidth', 1.2);
                    yline(0, ':', 'Color', [0.5 0.5 0.5]);
                    title(sprintf('%s target/reconv | kernel %.2gx | R2 %.3f, NRMSE %.3g', ...
                        pix_label, deconv.kernel_width_scale, deconv.R2_weighted, deconv.NRMSE_weighted), ...
                        'Interpreter', 'none');
                    xlabel('Energy loss (meV)'); ylabel('norm intensity'); grid on;
                else
                    text(0.05, 0.55, sprintf('%s unavailable: %s', pix_label, deconv.reason), ...
                        'Units', 'normalized', 'Interpreter', 'none');
                    axis off;
                end
                nexttile;
                if deconv.available
                    stem(deconv.EM_meV, deconv.coef_stokes, 'Color', [0.85 0.25 0.05], 'Marker', 'none', 'LineWidth', 0.9); hold on;
                    stem(-deconv.EM_meV, deconv.coef_antistokes, 'Color', [0.05 0.25 0.85], 'Marker', 'none', 'LineWidth', 0.9);
                    stem(deconv.EM_meV, deconv.PM_stokes_prior, ':', 'Color', [0.95 0.55 0.25], 'Marker', 'none', 'LineWidth', 0.75);
                    stem(-deconv.EM_meV, deconv.PM_antistokes_prior, ':', 'Color', [0.25 0.55 0.95], 'Marker', 'none', 'LineWidth', 0.75);
                    title(sprintf('sticks | centers %d | coef BE T %.1f K', deconv.num_centers, deconv.coef_ratio_T_K), ...
                        'Interpreter', 'none');
                    xlabel('EM (meV), gain negative'); ylabel('amplitude'); grid on;
                else
                    axis off;
                end
            end
            sgtitle(title_text, 'Interpreter', 'none');
        end

        function lgfit = spectrum_loss_gain_ln_ratio_fit(res)
            lgfit = struct('available', false, 'fit_available', false, ...
                'E_meV', [], 'ln_loss_gain', [], 'I_loss', [], 'I_gain', [], ...
                'weight', [], 'attention_pair', [], 'slope_per_meV', NaN, ...
                'T_K', NaN, 'R2_origin', NaN, 'num_points', 0, ...
                'num_points_selected', 0, 'selected_mask', [], ...
                'used_E_range_meV', [NaN, NaN], 'method', 'best_contiguous_origin_line', ...
                'source', 'spectrum_minus_zlp_bg', 'window_meV', [NaN, NaN]);
            if ~(isstruct(res) && isfield(res, 'ene_corrected'))
                return;
            end
            ene = res.ene_corrected(:);
            if isfield(res, 'raw_data_norm') && ~isempty(res.raw_data_norm)
                raw = res.raw_data_norm(:);
            elseif isfield(res, 'raw_data') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                raw = res.raw_data(:) ./ res.norm_val;
            elseif isfield(res, 'raw_data')
                raw = res.raw_data(:);
            else
                return;
            end
            n = min(numel(ene), numel(raw));
            ene = ene(1:n);
            raw = raw(1:n);
            zlp = zeros(n, 1);
            bg = zeros(n, 1);
            if isfield(res, 'stage2') && isstruct(res.stage2) && ...
                    isfield(res.stage2, 'status') && strcmp(res.stage2.status, 'ok')
                if isfield(res.stage2, 'zlp_norm') && numel(res.stage2.zlp_norm) >= n
                    tmp = res.stage2.zlp_norm(:);
                    zlp = tmp(1:n);
                end
                if isfield(res.stage2, 'bg_norm') && numel(res.stage2.bg_norm) >= n
                    tmp = res.stage2.bg_norm(:);
                    bg = tmp(1:n);
                end
            end
            if ~any(isfinite(zlp) & abs(zlp) > 0)
                if isfield(res, 'stage1') && isfield(res.stage1, 'zlp_norm') && numel(res.stage1.zlp_norm) >= n
                    tmp = res.stage1.zlp_norm(:);
                    zlp = tmp(1:n);
                elseif isfield(res, 'zlp_only') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                    tmp = res.zlp_only(:);
                    zlp = tmp(1:n) ./ res.norm_val;
                end
            end
            if ~any(isfinite(bg) & abs(bg) > 0)
                if isfield(res, 'stage1') && isfield(res.stage1, 'bg_norm') && numel(res.stage1.bg_norm) >= n
                    tmp = res.stage1.bg_norm(:);
                    bg = tmp(1:n);
                elseif isfield(res, 'bg_only') && isfield(res, 'norm_val') && isfinite(res.norm_val) && res.norm_val > 0
                    tmp = res.bg_only(:);
                    bg = tmp(1:n) ./ res.norm_val;
                end
            end
            signal = raw(:) - zlp(:) - bg(:);
            if numel(signal) >= 5
                signal = movmedian(signal, 3, 'omitnan');
            end
            pos = ene(ene > 0);
            neg_abs = -ene(ene < 0);
            if isempty(pos) || isempty(neg_abs)
                return;
            end
            full_lo = max(min(pos, [], 'omitnan'), min(neg_abs, [], 'omitnan'));
            full_hi = min(max(pos, [], 'omitnan'), max(neg_abs, [], 'omitnan'));
            min_em = 8;
            if isfield(res, 'settings') && isfield(res.settings, 'em_min_meV') && isfinite(res.settings.em_min_meV)
                min_em = max(4, res.settings.em_min_meV);
            end
            full_lo = max(full_lo, min_em);
            if ~(isfinite(full_lo) && isfinite(full_hi) && full_hi > full_lo)
                return;
            end
            lo = full_lo;
            hi = min(full_hi, max(full_lo + 20, 220));
            if isfield(res, 'stage2') && isstruct(res.stage2) && isfield(res.stage2, 'window_meV')
                win = sort(abs(res.stage2.window_meV(:).'));
                if numel(win) >= 2 && all(isfinite(win(1:2)))
                    lo = max(full_lo, win(1));
                    hi = min(full_hi, win(2));
                end
            elseif isfield(res, 'settings') && isfield(res.settings, 'stage2_attention_window_meV')
                win = sort(abs(res.settings.stage2_attention_window_meV(:).'));
                if numel(win) >= 2 && all(isfinite(win(1:2)))
                    lo = max(full_lo, win(1));
                    hi = min(full_hi, win(2));
                end
            end
            E = unique(pos(pos >= lo & pos <= hi));
            if numel(E) < 3
                E = unique(pos(pos >= full_lo & pos <= full_hi));
                lo = full_lo;
                hi = full_hi;
            end
            if numel(E) < 2
                return;
            end
            Iloss = interp1(ene, signal, E, 'linear', NaN);
            Igain = interp1(ene, signal, -E, 'linear', NaN);
            noise_mask = abs(ene) >= full_lo & abs(ene) <= full_hi;
            noise = EELS_PearsonSparse_FullFit.robust_noise(signal, noise_mask);
            raw_scale = max(abs(raw), [], 'omitnan');
            if ~isfinite(raw_scale) || raw_scale <= 0, raw_scale = 1; end
            floor_val = max([eps, 0.25 .* max(noise, 0), 1e-8 .* raw_scale]);
            good = isfinite(E) & isfinite(Iloss) & isfinite(Igain) & Iloss > floor_val & Igain > floor_val;
            if nnz(good) < 2
                good = isfinite(E) & isfinite(Iloss) & isfinite(Igain) & Iloss > eps & Igain > eps;
            end
            if nnz(good) < 2
                return;
            end
            x = E(good);
            loss = Iloss(good);
            gain = Igain(good);
            y = log(loss ./ gain);
            w = sqrt(max(loss, 0) .* max(gain, 0));
            att_pair = ones(size(x));
            if isfield(res, 'stage2') && isstruct(res.stage2) && isfield(res.stage2, 'attention') && ...
                    isfield(res.stage2.attention, 'feature_score')
                att = res.stage2.attention.feature_score(:);
                if numel(att) >= n
                    att_pos = interp1(ene, att(1:n), x, 'linear', 0);
                    att_neg = interp1(ene, att(1:n), -x, 'linear', 0);
                    att_pair = sqrt(max(att_pos, 0) .* max(att_neg, 0));
                    w = w .* max(att_pair, 0.15);
                end
            end
            good2 = isfinite(x) & isfinite(y) & isfinite(w) & x > 0 & w > 0;
            x = x(good2);
            y = y(good2);
            loss = loss(good2);
            gain = gain(good2);
            w = w(good2);
            att_pair = att_pair(good2);
            if numel(x) < 2
                return;
            end
            w = w ./ max(max(w), eps);
            lgfit.available = true;
            lgfit.E_meV = x(:);
            lgfit.ln_loss_gain = y(:);
            lgfit.I_loss = loss(:);
            lgfit.I_gain = gain(:);
            lgfit.weight = w(:);
            lgfit.attention_pair = att_pair(:);
            lgfit.selected_mask = false(size(x(:)));
            best_line = EELS_PearsonSparse_FullFit.best_ratio_origin_line(x, y, w, EELS_PearsonSparse_FullFit.defaults(struct()));
            if best_line.available
                lgfit.fit_available = true;
                lgfit.slope_per_meV = best_line.slope_per_meV;
                lgfit.T_K = best_line.T_K;
                lgfit.R2_origin = best_line.R2_origin;
                lgfit.num_points_selected = best_line.num_points;
                lgfit.selected_mask = best_line.used_mask(:);
                lgfit.used_E_range_meV = best_line.EM_range_meV;
            end
            lgfit.num_points = numel(x);
            lgfit.window_meV = [lo, hi];
        end

        function befit = hard_be_ln_ratio_fit(res)
            befit = struct('available', false, 'EM_meV', [], 'ln_loss_gain', [], ...
                'weight', [], 'slope_per_meV', NaN, 'T_K', NaN, 'R2_origin', NaN, ...
                'num_points', 0, 'source', 'none', 'label', 'Attention source intensity', ...
                'num_attention_points', 0, 'num_both_side_points', 0, ...
                'method', 'weighted_origin_line', 'selected_mask', [], ...
                'num_points_selected', 0, 'used_EM_range_meV', [NaN, NaN], ...
                'num_points_raw', 0, 'num_outliers_hidden', 0, ...
                'outlier_residual_threshold', NaN);
            kb = 0.08617333262;
            EM = [];
            PMs = [];
            PMa = [];
            hard_EM = [];
            attention_bands = [];
            if isfield(res, 'stage2') && isstruct(res.stage2) && ...
                    isfield(res.stage2, 'EM_final') && isfield(res.stage2, 'PM_stokes') && ...
                    isfield(res.stage2, 'PM_antistokes') && ~isempty(res.stage2.EM_final)
                EM = res.stage2.EM_final(:);
                PMs = res.stage2.PM_stokes(:);
                PMa = res.stage2.PM_antistokes(:);
                if isfield(res.stage2, 'EM_hard')
                    hard_EM = res.stage2.EM_hard(:);
                end
                bands = EELS_PearsonSparse_FullFit.stage2_temperature_bands(res);
                if isfield(bands, 'attention') && ~isempty(bands.attention)
                    attention_bands = bands.attention;
                elseif isfield(bands, 'fit') && ~isempty(bands.fit)
                    attention_bands = bands.fit;
                end
                befit.source = 'stage2_attention_PM';
            elseif isfield(res, 'EM_final') && isfield(res, 'PM_stokes') && isfield(res, 'PM_antistokes')
                EM = res.EM_final(:);
                PMs = res.PM_stokes(:);
                PMa = res.PM_antistokes(:);
                befit.source = 'stage1_PM_fallback';
            end
            if isempty(EM) || numel(PMs) ~= numel(EM) || numel(PMa) ~= numel(EM)
                return;
            end
            keep = true(size(EM));
            both_side = false(size(EM));
            if ~isempty(attention_bands)
                pos_in = EELS_PearsonSparse_FullFit.values_in_energy_bands(EM, attention_bands);
                neg_in = EELS_PearsonSparse_FullFit.values_in_energy_bands(-EM, attention_bands);
                keep = pos_in | neg_in;
                both_side = pos_in & neg_in;
                befit.label = 'Attention-window source intensity';
                befit.num_attention_points = nnz(keep);
                befit.num_both_side_points = nnz(both_side);
            elseif ~isempty(hard_EM)
                keep = false(size(EM));
                tol = max(0.75, 0.45 .* median(diff(unique(sort(EM))), 'omitnan'));
                if ~isfinite(tol) || tol <= 0, tol = 1.5; end
                for ii = 1:numel(hard_EM)
                    [dmin, idx] = min(abs(EM - hard_EM(ii)));
                    if isfinite(dmin) && dmin <= tol
                        keep(idx) = true;
                    end
                end
                both_side = keep;
                befit.label = 'Hard-pair source intensity';
                befit.num_attention_points = nnz(keep);
                befit.num_both_side_points = nnz(both_side);
            else
                befit.label = 'All source intensities';
                befit.num_attention_points = numel(EM);
                befit.num_both_side_points = numel(EM);
            end
            good = keep & isfinite(EM) & EM > 0 & isfinite(PMs) & isfinite(PMa) & PMs > 0 & PMa > 0;
            if nnz(good) < 2
                return;
            end
            x = EM(good);
            y = log(PMs(good) ./ PMa(good));
            w = sqrt(max(PMs(good), 0) .* max(PMa(good), 0));
            both_good = both_side(good);
            good2 = isfinite(x) & isfinite(y) & isfinite(w) & x > 0 & w > 0;
            x = x(good2);
            y = y(good2);
            w = w(good2);
            both_good = both_good(good2);
            if numel(x) < 2
                return;
            end
            w = w ./ max(max(w), eps);
            if any(both_good)
                w(both_good) = min(1.0, 1.5 .* w(both_good));
            end
            plot_settings = EELS_PearsonSparse_FullFit.defaults(struct());
            best_line = EELS_PearsonSparse_FullFit.best_ratio_origin_line(x, y, w, plot_settings);
            if best_line.available
                slope = best_line.slope_per_meV;
                befit.method = best_line.method;
                befit.num_points_selected = best_line.num_points;
                befit.used_EM_range_meV = best_line.EM_range_meV;
                befit.R2_origin = best_line.R2_origin;
                befit.T_K = best_line.T_K;
                befit.available = true;
            else
                return;
            end
            befit.slope_per_meV = slope;
            befit.num_points_raw = numel(x);
            selected_mask_raw = best_line.used_mask(:);
            residual = y(:) - slope .* x(:);
            ref_residual = residual(selected_mask_raw);
            sigma = 1.4826 .* median(abs(ref_residual - median(ref_residual, 'omitnan')), 'omitnan');
            if ~isfinite(sigma) || sigma <= 0
                sigma = sqrt(mean(ref_residual .^ 2, 'omitnan'));
            end
            flyer_threshold = max(1.25, 2.8 .* max(sigma, eps));
            befit.outlier_residual_threshold = flyer_threshold;
            hide_flyers = false(size(x(:)));
            hide_flyers_enabled = true;
            if isfield(plot_settings, 'stage2_pm_ratio_plot_hide_flyers')
                hide_flyers_enabled = plot_settings.stage2_pm_ratio_plot_hide_flyers;
            end
            flyer_weight_fraction = 0.30;
            if isfield(plot_settings, 'stage2_pm_ratio_plot_flyer_weight_fraction')
                flyer_weight_fraction = plot_settings.stage2_pm_ratio_plot_flyer_weight_fraction;
            end
            if hide_flyers_enabled
                hide_flyers = ~selected_mask_raw & abs(residual) > flyer_threshold & ...
                    w(:) <= flyer_weight_fraction;
            end
            plot_mask = ~hide_flyers;
            if nnz(plot_mask) < max(2, nnz(selected_mask_raw))
                plot_mask = true(size(plot_mask));
                hide_flyers = false(size(hide_flyers));
            end
            befit.num_outliers_hidden = nnz(hide_flyers);
            befit.EM_meV = x(plot_mask);
            befit.ln_loss_gain = y(plot_mask);
            befit.weight = w(plot_mask);
            befit.selected_mask = selected_mask_raw(plot_mask);
            befit.num_points = nnz(plot_mask);
        end

        function fig = plot_result(res, title_text)
            if nargin < 2, title_text = 'Pearson sparse full-spectrum fit'; end
            fig = figure('Visible', 'off', 'Color', 'w', 'Position', [80, 35, 1550, 1450]);
            tiledlayout(6, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            ene = res.ene_corrected(:);
            ci_color = [0.52 0.66 0.96];
            st2_color = [0.62 0.1 0.72];
            has_ci = isfield(res, 'confidence_interval') && isfield(res.confidence_interval, 'fit_lower') && ...
                isfield(res.confidence_interval, 'fit_upper');
            has_stage2 = isfield(res, 'stage2') && isfield(res.stage2, 'status') && ...
                strcmp(res.stage2.status, 'ok') && isfield(res.stage2, 'model_norm');
            stage2_bands = struct('fit', [], 'attention', [], 'hard', [], 'transition', []);
            if has_stage2
                stage2_bands = EELS_PearsonSparse_FullFit.stage2_temperature_bands(res);
            end
            if has_stage2
                plot_fit = res.stage2.model_norm(:) .* res.norm_val;
                if isfield(res.stage2, 'zlp_norm')
                    plot_zlp = res.stage2.zlp_norm(:) .* res.norm_val;
                else
                    plot_zlp = res.zlp_only(:);
                end
                if isfield(res.stage2, 'vib_norm')
                    plot_vib = res.stage2.vib_norm(:) .* res.norm_val;
                else
                    plot_vib = res.vib_only(:);
                end
                if isfield(res.stage2, 'vib_hard_norm')
                    plot_vib_hard = res.stage2.vib_hard_norm(:) .* res.norm_val;
                else
                    plot_vib_hard = [];
                end
                if isfield(res.stage2, 'vib_soft_norm')
                    plot_vib_soft = res.stage2.vib_soft_norm(:) .* res.norm_val;
                else
                    plot_vib_soft = [];
                end
                if isfield(res.stage2, 'bg_norm')
                    plot_bg = res.stage2.bg_norm(:) .* res.norm_val;
                else
                    plot_bg = res.bg_only(:);
                end
                plot_total_label = 'Stage 2 hard-BE';
                plot_zlp_label = 'Stage 2 ZLP';
                plot_vib_label = 'Stage 2 vib';
                plot_bg_label = 'Stage 2 bg';
                plot_ci_label = 'Stage 2 95% CI';
                has_plot_ci = has_ci && isfield(res.confidence_interval, 'stage2_lower') && ...
                    isfield(res.confidence_interval, 'stage2_upper');
                if has_plot_ci
                    plot_ci_lower = res.confidence_interval.stage2_lower(:);
                    plot_ci_upper = res.confidence_interval.stage2_upper(:);
                else
                    plot_ci_lower = [];
                    plot_ci_upper = [];
                end
                if has_ci && isfield(res.confidence_interval, 'stage2_sigma')
                    plot_ci_half_width = EELS_PearsonSparse_FullFit.normal_z_from_level(res.confidence_interval.level) .* ...
                        res.confidence_interval.stage2_sigma(:);
                elseif has_ci && isfield(res.confidence_interval, 'half_width')
                    plot_ci_half_width = res.confidence_interval.half_width(:);
                else
                    plot_ci_half_width = [];
                end
                title_tail = sprintf('Tfit %.1f +/- %.1f K, rel %.2f', ...
                    res.stage2.T_fit, res.stage2.T_uncertainty, res.stage2.reliability_score);
            else
                plot_fit = res.fitted_spectrum(:);
                plot_zlp = res.zlp_only(:);
                plot_vib = res.vib_only(:);
                plot_vib_hard = [];
                plot_vib_soft = [];
                plot_bg = res.bg_only(:);
                plot_total_label = 'Stage 1 total';
                plot_zlp_label = 'ZLP';
                plot_vib_label = 'Stage 1 sparse conv';
                plot_bg_label = 'Background';
                plot_ci_label = 'Stage 1 95% CI';
                has_plot_ci = has_ci;
                if has_plot_ci
                    plot_ci_lower = res.confidence_interval.fit_lower(:);
                    plot_ci_upper = res.confidence_interval.fit_upper(:);
                else
                    plot_ci_lower = [];
                    plot_ci_upper = [];
                end
                if has_ci && isfield(res.confidence_interval, 'half_width')
                    plot_ci_half_width = res.confidence_interval.half_width(:);
                else
                    plot_ci_half_width = [];
                end
                title_tail = sprintf('Tsoft %.1f K', res.T_soft);
            end
            plot_residual = res.raw_data(:) - plot_fit(:);
            plot_has_bg = any(isfinite(plot_bg(:)) & abs(plot_bg(:)) > 1e-10 .* max(max(abs(plot_fit(:))), 1));
            zlp_component_label = 'ZLP components';
            if isfield(res, 'settings') && isfield(res.settings, 'zlp_model_type')
                zlp_component_label = sprintf('%s ZLP components', strrep(res.settings.zlp_model_type, '_', '-'));
            end
            nexttile([1 2]);
            hmain = gobjects(0);
            lmain = {};
            if has_plot_ci
                h = fill([ene; flipud(ene)], [plot_ci_lower(:); flipud(plot_ci_upper(:))], ...
                    ci_color, 'FaceAlpha', 0.22, 'EdgeColor', 'none');
                hmain(end + 1) = h;
                lmain{end + 1} = plot_ci_label;
                hold on;
                plot(ene, plot_ci_lower(:), '--', 'Color', ci_color .* 0.75, 'LineWidth', 0.8, 'HandleVisibility', 'off');
                plot(ene, plot_ci_upper(:), '--', 'Color', ci_color .* 0.75, 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
            h = plot(ene, res.raw_data, 'k.', 'MarkerSize', 4); hold on;
            hmain(end + 1) = h; lmain{end + 1} = 'Raw';
            if has_stage2
                h = plot(ene, plot_fit, '-', 'Color', st2_color, 'LineWidth', 1.8);
            else
                h = plot(ene, plot_fit, 'k-', 'LineWidth', 1.8);
            end
            hmain(end + 1) = h; lmain{end + 1} = plot_total_label;
            h = plot(ene, plot_zlp, 'Color', [0.05 0.45 0.18], 'LineWidth', 1.2);
            hmain(end + 1) = h; lmain{end + 1} = plot_zlp_label;
            if isfield(res, 'zlp_components') && ~isempty(res.zlp_components)
                zc = res.zlp_components;
                if size(zc, 1) ~= numel(ene) && size(zc, 2) == numel(ene)
                    zc = zc.';
                end
                if size(zc, 1) == numel(ene)
                    h = plot(ene, zc(:, 1), '-', 'Color', [0.55 0.78 0.28], 'LineWidth', 0.7);
                    hmain(end + 1) = h; lmain{end + 1} = zlp_component_label;
                    for kk = 2:size(zc, 2)
                        plot(ene, zc(:, kk), '-', 'Color', [0.55 0.78 0.28], ...
                            'LineWidth', 0.7, 'HandleVisibility', 'off');
                    end
                end
            end
            h = plot(ene, plot_vib, 'Color', [0.9 0.25 0.05], 'LineWidth', 1.2);
            hmain(end + 1) = h; lmain{end + 1} = plot_vib_label;
            if has_stage2 && ~isempty(plot_vib_hard)
                h = plot(ene, plot_vib_hard, ':', 'Color', [0.78 0.18 0.02], 'LineWidth', 0.9);
                hmain(end + 1) = h; lmain{end + 1} = 'Stage 2 hard vib part';
            end
            if has_stage2 && ~isempty(plot_vib_soft)
                h = plot(ene, plot_vib_soft, ':', 'Color', [0.15 0.35 0.95], 'LineWidth', 0.9);
                hmain(end + 1) = h; lmain{end + 1} = 'Stage 2 soft vib part';
            end
            if plot_has_bg
                h = plot(ene, plot_bg, '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 1.1);
                hmain(end + 1) = h; lmain{end + 1} = plot_bg_label;
            end
            annotation_lines = {};
            if has_stage2
                yl = ylim;
                EELS_PearsonSparse_FullFit.shade_energy_bands(gca, stage2_bands.fit, yl, [0.72 0.72 0.72], 0.055);
                EELS_PearsonSparse_FullFit.shade_energy_bands(gca, stage2_bands.transition, yl, [0.35 0.70 1.00], 0.075);
                EELS_PearsonSparse_FullFit.shade_energy_bands(gca, stage2_bands.attention, yl, [1.0 0.72 0.08], 0.11);
                ylim(yl);
                soft_count = 0;
                if isfield(res.stage2, 'EM_soft')
                    soft_count = numel(res.stage2.EM_soft);
                end
                annotation_lines{end + 1} = sprintf('Stage2 feature attention: %s | hard/soft EM %d/%d | policy %s', ...
                    EELS_PearsonSparse_FullFit.format_band_list(stage2_bands.attention, 4), ...
                    numel(res.stage2.EM_hard), soft_count, ...
                    strrep(res.stage2.temperature_policy, '_', ' '));
                if ~isempty(stage2_bands.transition)
                    leak_score = NaN;
                    if isfield(res.stage2, 'transition_hard_leakage_score')
                        leak_score = res.stage2.transition_hard_leakage_score;
                    end
                    annotation_lines{end + 1} = sprintf('Stage2 soft transition fit: %s | hard-leak %.3g', ...
                        EELS_PearsonSparse_FullFit.format_band_list(stage2_bands.transition, 3), ...
                        leak_score);
                end
            end
            if isfield(res, 'zlp_signal_guard') && isfield(res.zlp_signal_guard, 'enabled') && res.zlp_signal_guard.enabled
                method = 'static';
                if isfield(res, 'zlp_region_policy') && isfield(res.zlp_region_policy, 'method')
                    method = res.zlp_region_policy.method;
                end
                annotation_lines{end + 1} = sprintf('ZLP guard: %s [%.1f, %.1f] meV | ZLP fraction %.2f', ...
                    method, res.zlp_signal_guard.energy_range_meV(1), ...
                    res.zlp_signal_guard.energy_range_meV(2), res.zlp_signal_guard.zlp_area_fraction);
                if isfield(res.zlp_signal_guard, 'lowfreq_guard_enabled') && res.zlp_signal_guard.lowfreq_guard_enabled
                    annotation_lines{end + 1} = sprintf('Low-freq ZLP guard: [%.1f, %.1f] meV | weight %.2g', ...
                        res.zlp_signal_guard.lowfreq_guard_range_meV(1), ...
                        res.zlp_signal_guard.lowfreq_guard_range_meV(2), ...
                        res.zlp_signal_guard.lowfreq_guard_weight);
                end
            end
            if isfield(res, 'stage1_zlp_alignment') && isstruct(res.stage1_zlp_alignment)
                za = res.stage1_zlp_alignment;
                if isfield(za, 'applied_shift_meV') && isfield(za, 'status')
                    residual_shift = NaN;
                    if isfield(res, 'zlp_params') && isfield(res.zlp_params, 'E_shift')
                        residual_shift = res.zlp_params.E_shift;
                    end
                    annotation_lines{end + 1} = sprintf('Energy align: pre-shift %.3g meV (%s), residual Eshift %.3g meV', ...
                        za.applied_shift_meV, za.status, residual_shift);
                end
            end
            if isfield(res, 'zlp_params') && isfield(res.zlp_params, 'sL') && isfield(res.zlp_params, 'sR')
                widths = 0.5 .* (res.zlp_params.sL(:).' + res.zlp_params.sR(:).');
                if ~isempty(widths)
                    annotation_lines{end + 1} = sprintf('ZLP component mean widths: %s meV', ...
                        strjoin(arrayfun(@(v) sprintf('%.1f', v), widths, 'UniformOutput', false), ', '));
                end
            end
            if ~isempty(annotation_lines)
                text(0.012, 0.965, strjoin(annotation_lines, newline), 'Units', 'normalized', 'Interpreter', 'none', ...
                    'FontSize', 8.5, 'VerticalAlignment', 'top', 'BackgroundColor', [1 1 1], ...
                    'EdgeColor', [0.70 0.70 0.70], 'Margin', 3);
            end
            legend(hmain, lmain, 'Location', 'northeastoutside');
            title(sprintf('%s | R2 %.6g, NRMSE %.4g, wNRMSE %.4g, %s', ...
                title_text, res.metrics.R2, res.metrics.NRMSE, res.metrics.weighted_NRMSE, title_tail), 'Interpreter', 'none');
            ylabel('Intensity'); xlabel('Energy loss (meV)'); grid on;

            nexttile;
            zoom_mask = abs(ene) >= max(10, res.settings.em_min_meV);
            core_clip_meV = max([25, res.settings.optimization_keep_core_meV, 2 .* res.settings.core_downweight_width]);
            zoom_scale_mask = zoom_mask & abs(ene) >= core_clip_meV;
            if nnz(zoom_scale_mask) < 10
                zoom_scale_mask = zoom_mask;
            end
            if has_plot_ci
                zz = zoom_mask(:);
                fill([ene(zz); flipud(ene(zz))], [plot_ci_lower(zz); flipud(plot_ci_upper(zz))], ...
                    ci_color, 'FaceAlpha', 0.28, 'EdgeColor', 'none', 'HandleVisibility', 'off'); hold on;
                plot(ene(zz), plot_ci_lower(zz), '--', 'Color', ci_color .* 0.75, 'LineWidth', 0.8, 'HandleVisibility', 'off');
                plot(ene(zz), plot_ci_upper(zz), '--', 'Color', ci_color .* 0.75, 'LineWidth', 0.8, 'HandleVisibility', 'off');
            end
            plot(ene(zoom_mask), res.raw_data(zoom_mask), 'k.', 'MarkerSize', 4); hold on;
            if has_stage2
                plot(ene(zoom_mask), plot_fit(zoom_mask), '-', 'Color', st2_color, 'LineWidth', 1.5);
            else
                plot(ene(zoom_mask), plot_fit(zoom_mask), 'k-', 'LineWidth', 1.5);
            end
            plot(ene(zoom_mask), plot_zlp(zoom_mask), 'Color', [0.05 0.45 0.18], 'LineWidth', 1.1);
            if isfield(res, 'zlp_components') && ~isempty(res.zlp_components)
                zc = res.zlp_components;
                if size(zc, 1) ~= numel(ene) && size(zc, 2) == numel(ene)
                    zc = zc.';
                end
                if size(zc, 1) == numel(ene)
                    for kk = 1:size(zc, 2)
                        plot(ene(zoom_mask), zc(zoom_mask, kk), '-', ...
                            'Color', [0.55 0.78 0.28], 'LineWidth', 0.55, 'HandleVisibility', 'off');
                    end
                end
            end
            plot(ene(zoom_mask), plot_vib(zoom_mask), 'Color', [0.9 0.25 0.05], 'LineWidth', 1.1);
            if has_stage2 && ~isempty(plot_vib_hard)
                plot(ene(zoom_mask), plot_vib_hard(zoom_mask), ':', 'Color', [0.78 0.18 0.02], 'LineWidth', 0.85);
            end
            if has_stage2 && ~isempty(plot_vib_soft)
                plot(ene(zoom_mask), plot_vib_soft(zoom_mask), ':', 'Color', [0.15 0.35 0.95], 'LineWidth', 0.85);
            end
            if plot_has_bg
                plot(ene(zoom_mask), plot_bg(zoom_mask), '--', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.9);
            end
            zoom_vals = [res.raw_data(zoom_scale_mask); plot_fit(zoom_scale_mask); plot_vib(zoom_scale_mask)];
            if has_plot_ci
                zoom_vals = [zoom_vals; plot_ci_lower(zoom_scale_mask); plot_ci_upper(zoom_scale_mask)];
            end
            zoom_lims = EELS_PearsonSparse_FullFit.robust_limits(zoom_vals, 0.01, 0.995);
            if all(isfinite(zoom_lims)) && zoom_lims(2) > zoom_lims(1)
                ylim(zoom_lims);
            end
            if has_stage2
                yl = ylim;
                EELS_PearsonSparse_FullFit.shade_energy_bands(gca, stage2_bands.fit, yl, [0.72 0.72 0.72], 0.065);
                EELS_PearsonSparse_FullFit.shade_energy_bands(gca, stage2_bands.transition, yl, [0.35 0.70 1.00], 0.10);
                EELS_PearsonSparse_FullFit.shade_energy_bands(gca, stage2_bands.attention, yl, [1.0 0.72 0.08], 0.14);
                ylim(yl);
                text(0.012, 0.94, sprintf('Feature attention bands: %s', ...
                    EELS_PearsonSparse_FullFit.format_band_list(stage2_bands.attention, 5)), ...
                    'Units', 'normalized', 'Interpreter', 'none', 'FontSize', 8.5, ...
                    'VerticalAlignment', 'top', 'BackgroundColor', [1 1 1], ...
                    'EdgeColor', [0.70 0.70 0.70], 'Margin', 3);
            end
            if has_stage2 && isfield(res.stage2, 'attention')
                yyaxis right;
                plot(ene(zoom_mask), res.stage2.attention.feature_score(zoom_mask), '-', 'Color', [0.35 0.35 0.9], 'LineWidth', 0.9);
                ylabel('Attention');
                ylim([0, 1]);
                yyaxis left;
            end
            title(sprintf('Signal and tail zoom, y-scale outside |E| < %.0f meV', core_clip_meV));
            ylabel('Intensity'); xlabel('Energy loss (meV)'); grid on;

            nexttile;
            if plot_has_bg
                raw_log = log10(max(res.raw_data - plot_bg, eps));
            else
                raw_log = log10(max(res.raw_data, eps));
            end
            zlp_log = log10(max(plot_zlp, eps));
            vib_log = log10(max(plot_vib, eps));
            plot(ene, raw_log, 'k.', 'MarkerSize', 4); hold on;
            plot(ene, zlp_log, 'Color', [0.05 0.45 0.18], 'LineWidth', 1.2);
            plot(ene, vib_log, 'Color', [0.9 0.25 0.05], 'LineWidth', 1.2);
            log_lims = EELS_PearsonSparse_FullFit.robust_limits([raw_log; zlp_log; vib_log], 0.02, 0.995);
            if all(isfinite(log_lims)) && log_lims(2) > log_lims(1)
                ylim(log_lims);
            end
            if plot_has_bg
                raw_log_label = 'log raw-bg';
            else
                raw_log_label = 'log raw';
            end
            if has_stage2
                legend({raw_log_label, 'log Stage 2 ZLP', 'log Stage 2 vib'}, 'Location', 'best');
            else
                legend({raw_log_label, 'log ZLP', 'log sparse conv'}, 'Location', 'best');
            end
            title('Log-scale signal vs ZLP'); ylabel('log10'); xlabel('Energy loss (meV)'); grid on;

            nexttile;
            hold on;
            pm_handles = gobjects(0);
            pm_labels = {};
            stage2_guarded_display = false;
            pm_ymax = 0;
            if isfield(res, 'EM_final') && isfield(res, 'PM_stokes') && isfield(res, 'PM_antistokes') && ...
                    ~isempty(res.EM_final)
                h = stem(res.EM_final, res.PM_stokes, 'Color', [0.95 0.45 0.12], ...
                    'Marker', 'none', 'LineWidth', 0.8);
                pm_handles(end + 1) = h; pm_labels{end + 1} = 'Stage 1 Stokes intensity';
                h = stem(-res.EM_final, res.PM_antistokes, 'Color', [0.1 0.50 0.95], ...
                    'Marker', 'none', 'LineWidth', 0.8);
                pm_handles(end + 1) = h; pm_labels{end + 1} = 'Stage 1 anti-Stokes intensity';
                pm_ymax = max([pm_ymax; res.PM_stokes(:); res.PM_antistokes(:)], [], 'omitnan');
            end
            if has_stage2 && ~stage2_guarded_display && isfield(res.stage2, 'PM_stokes') && ...
                    isfield(res.stage2, 'PM_antistokes') && ~isempty(res.stage2.PM_stokes)
                if isfield(res.stage2, 'EM_final') && ~isempty(res.stage2.EM_final)
                    em2 = res.stage2.EM_final;
                else
                    em2 = res.EM_final;
                end
                h = stem(em2, res.stage2.PM_stokes, 'Color', [0.75 0.02 0.02], ...
                    'Marker', 'none', 'LineWidth', 1.15);
                pm_handles(end + 1) = h; pm_labels{end + 1} = 'Stage 2 Stokes intensity';
                h = stem(-em2, res.stage2.PM_antistokes, 'Color', [0.0 0.20 0.85], ...
                    'Marker', 'none', 'LineWidth', 1.15);
                pm_handles(end + 1) = h; pm_labels{end + 1} = 'Stage 2 anti-Stokes intensity';
                pm_ymax = max([pm_ymax; res.stage2.PM_stokes(:); res.stage2.PM_antistokes(:)], [], 'omitnan');
            end
            if ~isfinite(pm_ymax) || pm_ymax <= 0
                pm_ymax = 1;
            end
            ylim([0, pm_ymax .* 1.20]);
            if has_stage2 && isfield(res.stage2, 'EM_hard') && ~isempty(res.stage2.EM_hard)
                plot(res.stage2.EM_hard, 0 .* res.stage2.EM_hard + pm_ymax .* 1.10, ...
                    'm.', 'MarkerSize', 7, 'HandleVisibility', 'off');
                plot(-res.stage2.EM_hard, 0 .* res.stage2.EM_hard + pm_ymax .* 1.04, ...
                    'c.', 'MarkerSize', 7, 'HandleVisibility', 'off');
            end
            if ~isempty(pm_handles)
                legend(pm_handles, pm_labels, 'Location', 'best');
            end
            if has_stage2
                title(sprintf('Sparse source intensity and hard-BE mode markers, feature windows %s', ...
                    EELS_PearsonSparse_FullFit.format_band_list(stage2_bands.attention, 3)), 'Interpreter', 'none');
            else
                title('Stage 1 independent source intensity');
            end
            ylabel('Source intensity'); xlabel('Mode energy (meV)'); grid on;

            nexttile;
            lgfit = EELS_PearsonSparse_FullFit.spectrum_loss_gain_ln_ratio_fit(res);
            if lgfit.available
                scatter(lgfit.E_meV, lgfit.ln_loss_gain, 24, lgfit.weight, 'filled'); hold on;
                if isfield(lgfit, 'selected_mask') && numel(lgfit.selected_mask) == numel(lgfit.E_meV) && any(lgfit.selected_mask)
                    scatter(lgfit.E_meV(lgfit.selected_mask), lgfit.ln_loss_gain(lgfit.selected_mask), ...
                        42, 'o', 'MarkerEdgeColor', [0.02 0.02 0.02], 'LineWidth', 0.95);
                end
                if lgfit.fit_available
                    xmax = max(lgfit.E_meV, [], 'omitnan');
                    xx = linspace(0, max(xmax, eps), 120);
                    plot(xx, lgfit.slope_per_meV .* xx, 'k-', 'LineWidth', 1.25);
                end
                yline(0, ':', 'Color', [0.45 0.45 0.45], 'HandleVisibility', 'off');
                cb = colorbar;
                cb.Label.String = 'spectrum pair weight';
                xlabel('|E| (meV)');
                ylabel('ln(I_{loss}/I_{gain})');
                if lgfit.fit_available
                    diag_Tfit = lgfit.T_K;
                    if has_stage2 && isfield(res.stage2, 'T_fit') && isfinite(res.stage2.T_fit)
                        diag_Tfit = res.stage2.T_fit;
                    end
                    title(sprintf('Spectrum BE consistency at Tfit %.1f K: R2 %.3f, n=%d/%d, [%.0f %.0f] meV', ...
                        diag_Tfit, lgfit.R2_origin, lgfit.num_points_selected, lgfit.num_points, ...
                        lgfit.used_E_range_meV(1), lgfit.used_E_range_meV(2)), ...
                        'Interpreter', 'none');
                else
                    title(sprintf('Spectrum ln(Iloss/Igain): n=%d, no accepted linear segment', lgfit.num_points), ...
                        'Interpreter', 'none');
                end
                grid on;
            else
                text(0.05, 0.55, 'No valid mirrored ZLP/bg-subtracted Iloss/Igain points', ...
                    'Units', 'normalized', 'Interpreter', 'none');
                axis off;
            end

            nexttile;
            befit = EELS_PearsonSparse_FullFit.hard_be_ln_ratio_fit(res);
            if befit.available
                scatter(befit.EM_meV, befit.ln_loss_gain, 28, befit.weight, 'filled'); hold on;
                if isfield(befit, 'selected_mask') && numel(befit.selected_mask) == numel(befit.EM_meV) && any(befit.selected_mask)
                    scatter(befit.EM_meV(befit.selected_mask), befit.ln_loss_gain(befit.selected_mask), ...
                        48, 'o', 'MarkerEdgeColor', [0.02 0.02 0.02], 'LineWidth', 1.0);
                end
                xmax = max(befit.EM_meV, [], 'omitnan');
                xx = linspace(0, max(xmax, eps), 120);
                plot(xx, befit.slope_per_meV .* xx, 'k-', 'LineWidth', 1.35);
                cb = colorbar;
                cb.Label.String = 'Source-intensity pair weight';
                xlabel('EM (meV)');
                ylabel('ln(I_{source,loss}/I_{source,gain})');
                hidden_suffix = '';
                if isfield(befit, 'num_outliers_hidden') && befit.num_outliers_hidden > 0
                    hidden_suffix = sprintf(', hidden %d flyers', befit.num_outliers_hidden);
                end
                diag_Tfit = befit.T_K;
                if has_stage2 && isfield(res.stage2, 'T_fit') && isfinite(res.stage2.T_fit)
                    diag_Tfit = res.stage2.T_fit;
                end
                title(sprintf('%s hard-BE consistency at Tfit %.1f K: R2 %.3f, n=%d/%d%s, [%.0f %.0f] meV', ...
                    befit.label, diag_Tfit, befit.R2_origin, befit.num_points_selected, befit.num_points, ...
                    hidden_suffix, befit.used_EM_range_meV(1), befit.used_EM_range_meV(2)), ...
                    'Interpreter', 'none');
                grid on;
            else
                text(0.05, 0.55, 'No valid attention-window Stokes/anti-Stokes source-intensity pairs for BE ratio fit', ...
                    'Units', 'normalized', 'Interpreter', 'none');
                axis off;
            end

            nexttile;
            if has_stage2 && isfield(res.stage2, 'T_grid') && isfield(res.stage2, 'cost_grid')
                if isfield(res.stage2, 'cost_grid_smooth')
                    valley_cost = res.stage2.cost_grid_smooth;
                else
                    valley_cost = res.stage2.cost_grid;
                end
                valley = EELS_PearsonSparse_FullFit.stage2_cost_valley(res.stage2.T_grid, valley_cost);
                if valley.available
                    Tplot = valley.T_grid(:);
                    Cplot = valley.relative_cost(:);
                    in_valley = valley.in_valley(:);
                    plot(Tplot, Cplot, '-', 'Color', [0.55 0.55 0.55], 'LineWidth', 0.8, ...
                        'HandleVisibility', 'off'); hold on;
                    scatter(Tplot(~in_valley), Cplot(~in_valley), 18, [0.62 0.62 0.62], ...
                        'filled', 'MarkerFaceAlpha', 0.55, 'DisplayName', 'sampled T');
                    scatter(Tplot(in_valley), Cplot(in_valley), 30, [0.62 0.10 0.72], ...
                        'filled', 'MarkerFaceAlpha', 0.85, 'DisplayName', 'inside valley');
                    plot(valley.best_T, 0, 'p', 'Color', [0.05 0.15 0.05], ...
                        'MarkerFaceColor', [1.0 0.85 0.10], 'MarkerSize', 10, 'DisplayName', 'cost-min T');
                    final_T = NaN;
                    if isfield(res.stage2, 'T_fit')
                        final_T = res.stage2.T_fit;
                    end
                    final_C = NaN;
                    if isfinite(final_T)
                        final_C = interp1(Tplot, Cplot, final_T, 'linear', 'extrap');
                        plot(final_T, final_C, 'd', 'Color', [0.0 0.35 0.95], ...
                            'MarkerFaceColor', [0.35 0.75 1.0], 'MarkerSize', 8, 'DisplayName', 'final Tfit');
                    end
                    yline(valley.relative_threshold, '--', 'Color', [0.85 0.25 0.10], ...
                        'LineWidth', 0.9, 'DisplayName', 'valley threshold');
                    ylims = ylim;
                    if isfinite(valley.lower_K)
                        plot([valley.lower_K valley.lower_K], ylims, '--', 'Color', [0.62 0.10 0.72], ...
                            'LineWidth', 1.0, 'HandleVisibility', 'off');
                    end
                    if isfinite(valley.upper_K)
                        plot([valley.upper_K valley.upper_K], ylims, '--', 'Color', [0.62 0.10 0.72], ...
                            'LineWidth', 1.0, 'HandleVisibility', 'off');
                    end
                    if isfinite(final_T)
                        plot([final_T final_T], ylims, '-', 'Color', [0.0 0.35 0.95], ...
                            'LineWidth', 1.0, 'HandleVisibility', 'off');
                    end
                    ylim(ylims);
                    xlabel('Temperature (K)');
                    ylabel('\Delta cost / cost_{min}');
                    title(sprintf('Stage 2 cost valley [%.0f, %.0f] K, half-width %.1f K', ...
                        valley.lower_K, valley.upper_K, valley.uncertainty_K), 'Interpreter', 'none');
                    legend('Location', 'best');
                    grid on;
                    if isfinite(final_T)
                        label_text = sprintf('cost-min %.1f K | final %.1f K | offset %.1f K | valley points %d/%d', ...
                            valley.best_T, final_T, final_T - valley.best_T, nnz(in_valley), numel(in_valley));
                    else
                        label_text = sprintf('cost-min %.1f K | valley points %d/%d', ...
                            valley.best_T, nnz(in_valley), numel(in_valley));
                    end
                    text(0.012, 0.94, label_text, ...
                        'Units', 'normalized', 'Interpreter', 'none', 'FontSize', 8.5, ...
                        'VerticalAlignment', 'top', 'BackgroundColor', [1 1 1], ...
                        'EdgeColor', [0.70 0.70 0.70], 'Margin', 3);
                else
                    text(0.05, 0.55, 'No valid Stage 2 cost valley points', ...
                        'Units', 'normalized', 'Interpreter', 'none');
                    axis off;
                end
            else
                text(0.05, 0.55, 'Stage 2 cost valley is not available', ...
                    'Units', 'normalized', 'Interpreter', 'none');
                axis off;
            end

            nexttile;
            if isfield(res, 'attention_summary')
                vals = res.attention_summary.mean_attention;
                labels = strrep(res.attention_summary.region_names, '_', ' ');
                barh(vals, 'FaceColor', [0.25 0.50 0.95]);
                xlim([0, max(1, 1.05 .* max(vals, [], 'omitnan'))]);
                set(gca, 'YTick', 1:numel(labels), 'YTickLabel', labels);
                xlabel('Mean attention'); title('Attention by spectral region'); grid on;
            else
                text(0.05, 0.5, 'No attention summary available'); axis off;
            end

            nexttile;
            if isfield(res, 'loss_breakdown')
                vals = 100 .* res.loss_breakdown.fraction;
                labels = strrep(res.loss_breakdown.names, '_', ' ');
                barh(vals, 'FaceColor', [0.95 0.48 0.20]); 
                if any(isfinite(vals))
                    xlim([0, max(1, 1.15 .* max(vals, [], 'omitnan'))]);
                end
                set(gca, 'YTick', 1:numel(labels), 'YTickLabel', labels);
                xlabel('Loss (%)'); title(sprintf('Loss budget, dominant = %s', strrep(res.loss_breakdown.dominant, '_', ' '))); grid on;
            else
                text(0.05, 0.5, 'No loss breakdown available'); axis off;
            end

            nexttile([1 2]);
            area(ene, plot_residual, 'FaceColor', [0.45 0.62 0.95], 'EdgeColor', [0.1 0.35 0.95]); hold on;
            if ~isempty(plot_ci_half_width)
                yyaxis right;
                plot(ene, plot_ci_half_width(:), '-', 'Color', [0.1 0.4 0.95], 'LineWidth', 0.9);
                ylabel('CI half-width');
                yyaxis left;
            end
            yline(0, 'k-');
            xlabel('Energy loss (meV)'); ylabel('Data - fit'); title('Residuals and energy-dependent CI width'); grid on;
        end

        function idx = pick_test_indices(T, n_spec, n_pick)
            n_pick = min(n_pick, n_spec);
            if n_pick <= 1
                idx = 1;
                return;
            end
            [~, order] = sort(T(:));
            pos = round(linspace(1, numel(order), n_pick));
            idx = order(pos).';
        end

        function write_summary(filename, testres)
            fid = fopen(filename, 'w');
            if fid < 0, return; end
            cleaner = onCleanup(@() fclose(fid));
            fprintf(fid, 'Pearson sparse full-spectrum Stage 1 summary\n');
            fprintf(fid, 'num_spectra: %d\n', numel(testres.indices));
            fprintf(fid, 'median_R2: %.8g\n', median(testres.R2_all, 'omitnan'));
            fprintf(fid, 'mean_R2: %.8g\n', mean(testres.R2_all, 'omitnan'));
            fprintf(fid, 'median_NRMSE: %.8g\n', median(testres.NRMSE_all, 'omitnan'));
            fprintf(fid, 'mean_NRMSE: %.8g\n', mean(testres.NRMSE_all, 'omitnan'));
            fprintf(fid, 'median_weighted_NRMSE: %.8g\n', median(testres.weighted_NRMSE_all, 'omitnan'));
            if isfield(testres, 'T_fit_all') && any(isfinite(testres.T_fit_all))
                fprintf(fid, 'stage2_median_abs_T_error_K: %.8g\n', median(abs(testres.T_error_all), 'omitnan'));
                fprintf(fid, 'stage2_mean_abs_T_error_K: %.8g\n', mean(abs(testres.T_error_all), 'omitnan'));
            end
            fprintf(fid, 'failed: %d\n\n', nnz(testres.fit_status < 0));
            for ii = 1:numel(testres.indices)
                if isfield(testres, 'T_fit_all') && isfinite(testres.T_fit_all(ii))
                    fprintf(fid, 'index %d | T true %.2f | T fit %.2f | T err %.2f | R2 %.8g | NRMSE %.8g | wNRMSE %.8g | status %d\n', ...
                        testres.indices(ii), testres.T_true(ii), testres.T_fit_all(ii), testres.T_error_all(ii), ...
                        testres.R2_all(ii), testres.NRMSE_all(ii), testres.weighted_NRMSE_all(ii), testres.fit_status(ii));
                else
                    fprintf(fid, 'index %d | R2 %.8g | NRMSE %.8g | wNRMSE %.8g | status %d\n', ...
                        testres.indices(ii), testres.R2_all(ii), testres.NRMSE_all(ii), ...
                        testres.weighted_NRMSE_all(ii), testres.fit_status(ii));
                end
                if testres.fit_status(ii) > 0
                    hist = testres.results{ii}.iter_history;
                    for rr = 1:numel(hist)
                        fprintf(fid, '  %s: R2 %.8g, NRMSE %.8g, high %.8g, EM %d\n', ...
                            hist(rr).round_name, hist(rr).R2, hist(rr).NRMSE, hist(rr).highfreq_NRMSE, hist(rr).num_EM);
                    end
                    if isfield(testres.results{ii}, 'stage2') && isfield(testres.results{ii}.stage2, 'attention') && ...
                            strcmp(testres.results{ii}.stage2.status, 'ok')
                        st2 = testres.results{ii}.stage2;
                        best_linearity = NaN;
                        effective_pairs = NaN;
                        if isfield(st2, 'temperature_profile') && isfield(st2.temperature_profile, 'linearity_R2') && ...
                                isfield(st2, 'T_grid')
                            [~, best_idx] = min(abs(st2.T_grid(:) - st2.T_fit));
                            if best_idx >= 1 && best_idx <= numel(st2.temperature_profile.linearity_R2)
                                best_linearity = st2.temperature_profile.linearity_R2(best_idx);
                            end
                            if isfield(st2.temperature_profile, 'effective_pair_count') && ...
                                    best_idx <= numel(st2.temperature_profile.effective_pair_count)
                                effective_pairs = st2.temperature_profile.effective_pair_count(best_idx);
                            end
                        end
                        soft_count = 0;
                        if isfield(st2, 'EM_soft')
                            soft_count = numel(st2.EM_soft);
                        end
                        hard_status = '';
                        if isfield(st2, 'hard_be_attention_status')
                            hard_status = st2.hard_be_attention_status;
                        end
                        fprintf(fid, '  stage2: status %s, T %.2f, unc %.2f, rel %.3f, attention_fraction %.3f, hard_EM %d, soft_EM %d, hard_status %s, EM_matches_stage1 %d, linearity_R2 %.4f, effective_pairs %.0f, policy %s\n', ...
                            st2.status, st2.T_fit, ...
                            st2.T_uncertainty, st2.reliability_score, st2.attention.selected_fraction, ...
                            numel(st2.EM_hard), soft_count, hard_status, logical(st2.stage2_EM_matches_stage1), ...
                            best_linearity, effective_pairs, st2.temperature_policy);
                    end
                end
            end
        end

        function noise = robust_noise(y, mask)
            yy = y(mask);
            yy = yy(:);
            yy = yy(isfinite(yy));
            if numel(yy) < 8
                yy = y(isfinite(y));
                yy = yy(:);
            end
            if numel(yy) < 8
                noise = eps;
                return;
            end
            med = median(yy, 'omitnan');
            noise = 1.4826 .* median(abs(yy - med), 'omitnan');
            noise = noise(1);
            if ~(isfinite(noise) && noise > 0)
                noise = std(yy, 'omitnan');
            end
            noise = noise(1);
            if ~(isfinite(noise) && noise > 0)
                noise = eps;
            end
        end

        function lims = robust_limits(vals, lo_frac, hi_frac)
            vals = vals(:);
            vals = vals(isfinite(vals));
            if isempty(vals)
                lims = [NaN, NaN];
                return;
            end
            vals = sort(vals);
            n = numel(vals);
            ilo = max(1, min(n, round(1 + lo_frac .* (n - 1))));
            ihi = max(1, min(n, round(1 + hi_frac .* (n - 1))));
            lo = vals(ilo);
            hi = vals(ihi);
            if hi <= lo
                span = max(abs([lo, hi]));
                if span <= 0, span = 1; end
                lo = lo - 0.05 .* span;
                hi = hi + 0.05 .* span;
            else
                pad = 0.08 .* (hi - lo);
                lo = lo - pad;
                hi = hi + pad;
            end
            lims = [lo, hi];
        end
    end
end
