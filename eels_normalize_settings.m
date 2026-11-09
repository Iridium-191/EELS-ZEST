function settings = eels_normalize_settings(settings)
% Keep the current configuration schema and supported runtime context.
names=fieldnames(settings);
if isempty(names),return;end
if isfield(settings,'use_bg')&&logical(settings.use_bg)
    error('EELS:UnsupportedBackground','No independent background term is supported; use_bg must be false.');
end
persistent supported
if isempty(supported)
    % Empty input exits above before looking up the schema.
    supported=fieldnames(EELS_PearsonSparse_FullFit.defaults(struct()));
    optional={'name','output_dir','spatial_bin','datacube_spatial_bin','energy_bin', ...
        'datacube_input_file','datacube_parallel_effective_workers', ...
        'native_objective_gradient','stage1_sparse_pm_policy', ...
        'stage2_attention_class_prune_enabled','stage2_dense_em_step_meV', ...
        'stage2_fine_half_width_K','stage2_fine_step_K','stage2_force_raw_display', ...
        'stage2_hard_em_max_centers','stage2_shape_guard_enabled', ...
        'stage2_single_fast_T_grid','stage2_single_fast_enabled', ...
        'stage2_single_fast_fine_half_width_K','stage2_single_fast_fine_step_K', ...
        'stage2_single_fast_hard_em_max_centers'};
    supported=union(supported,optional(:));
end
unknown=setdiff(names,supported);
if ~isempty(unknown),settings=rmfield(settings,unknown);end
if isfield(settings,'stage2_temperature_algorithm')
    settings.stage2_temperature_algorithm='grid_scan';
end
end
