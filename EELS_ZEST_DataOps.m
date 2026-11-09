classdef EELS_ZEST_DataOps
    %EELS_ZEST_DATAOPS Reproducible post-fit data transformations.

    methods (Static)
        function [cutoff_meV, meta] = pdos_cutoff_from_json(filename, sigma_meV, relative_floor)
            if nargin < 2 || isempty(sigma_meV), sigma_meV = 2; end
            if nargin < 3 || isempty(relative_floor), relative_floor = 0.005; end
            assert(isfile(filename), 'EELS:PDOSFileMissing', ...
                'Materials Project JSON file does not exist: %s', filename);
            data = jsondecode(fileread(filename));
            assert(isfield(data, 'phonon') && isfield(data.phonon, 'dos_frequencies') && ...
                isfield(data.phonon, 'ph_dos'), 'EELS:InvalidPDOSJSON', ...
                'JSON must contain phonon.dos_frequencies and phonon.ph_dos.');
            energy = double(data.phonon.dos_frequencies(:)) ./ 8.065543937;
            dos = max(double(data.phonon.ph_dos(:)), 0);
            keep = isfinite(energy) & isfinite(dos) & energy >= 0;
            energy = energy(keep); dos = dos(keep);
            [energy, order] = sort(energy); dos = dos(order);
            [energy, unique_idx] = unique(energy, 'stable'); dos = dos(unique_idx);
            assert(numel(energy) >= 3 && any(dos > 0), 'EELS:EmptyPDOS', ...
                'JSON does not contain usable positive-frequency PDOS data.');

            step = median(diff(energy), 'omitnan');
            if ~isfinite(step) || step <= 0, step = 0.25; end
            grid = (min(energy):step:max(energy)).';
            broadened = interp1(energy, dos, grid, 'linear', 0);
            if sigma_meV > 0
                half_width = max(1, ceil(5 .* sigma_meV ./ step));
                xk = (-half_width:half_width).' .* step;
                kernel = exp(-0.5 .* (xk ./ sigma_meV) .^ 2);
                kernel = kernel ./ max(sum(kernel), eps);
                broadened = conv(broadened, kernel, 'same');
            end
            peak = max(broadened, [], 'omitnan');
            support = isfinite(broadened) & broadened > relative_floor .* peak;
            assert(any(support), 'EELS:EmptyPDOSSupport', ...
                'No PDOS support remains above the selected relative floor.');
            cutoff_meV = ceil(max(grid(support)) .* 2) ./ 2;

            meta = struct('source_file', char(filename), 'sigma_meV', sigma_meV, ...
                'relative_floor', relative_floor, 'frequency_unit', 'cm^-1', ...
                'energy_unit', 'meV', 'conversion', 'meV = cm^-1 / 8.065543937', ...
                'cutoff_meV', cutoff_meV);
            metadata = struct();
            if isfield(data, 'metadata')
                metadata = data.metadata;
            elseif isfield(data.phonon, 'metadata')
                metadata = data.phonon.metadata;
            end
            if ~isempty(fieldnames(metadata))
                fields = {'material_id', 'formula', 'point_group'};
                for ii = 1:numel(fields)
                    if isfield(metadata, fields{ii})
                        meta.(fields{ii}) = metadata.(fields{ii});
                    end
                end
            end
        end

        function cutoff_meV = attention_cutoff(res, tail_meV)
            if nargin < 2 || isempty(tail_meV), tail_meV = 20; end
            cutoff_meV = NaN;
            if ~(isstruct(res) && isfield(res, 'stage2') && isstruct(res.stage2))
                return;
            end
            st2 = res.stage2;
            support = [];
            if isfield(st2, 'EM_hard')
                support = double(st2.EM_hard(:));
            end
            support = support(isfinite(support) & support > 0);
            if isempty(support) && isfield(st2, 'attention') && ...
                    isfield(st2.attention, 'bands_meV')
                bands = double(st2.attention.bands_meV);
                if size(bands, 2) >= 2
                    support = abs(bands(:));
                    support = support(isfinite(support) & support > 0);
                end
            end
            if ~isempty(support)
                cutoff_meV = max(support) + max(0, double(tail_meV));
            end
        end

        function [signal, detail] = reconstruct_sparse_signal(res, stage, cutoff_meV)
            if nargin < 2 || isempty(stage), stage = 'stage2'; end
            if nargin < 3 || isempty(cutoff_meV), cutoff_meV = Inf; end
            [ene, em, pm_s, pm_a, kernel, norm_val, source] = ...
                EELS_ZEST_DataOps.sparse_components(res, stage);
            assert(isfinite(cutoff_meV) && cutoff_meV > 0 || isinf(cutoff_meV), ...
                'EELS:InvalidCutoff', 'Cutoff must be a positive finite value or Inf.');
            retained = em <= cutoff_meV;
            pm_s_cut = pm_s; pm_a_cut = pm_a;
            pm_s_cut(~retained) = 0; pm_a_cut(~retained) = 0;
            stokes_kernel = EELS_PearsonSparse_FullFit.shift_kernel_matrix(ene, kernel, em.');
            anti_kernel = EELS_PearsonSparse_FullFit.shift_kernel_matrix(ene, kernel, -em.');
            signal_norm = stokes_kernel * pm_s_cut + anti_kernel * pm_a_cut;
            signal = signal_norm(:) .* norm_val;
            signal(~isfinite(signal)) = NaN;
            detail = struct('stage', source, 'energy_meV', ene(:), 'EM_meV', em(:), ...
                'PM_stokes_original', pm_s(:), 'PM_antistokes_original', pm_a(:), ...
                'PM_stokes_cut', pm_s_cut(:), 'PM_antistokes_cut', pm_a_cut(:), ...
                'retained_mask', retained(:), 'cutoff_meV', cutoff_meV, ...
                'kernel_norm', kernel(:), 'signal_norm', signal_norm(:), ...
                'norm_val', norm_val, 'vibrational_component', signal(:));
        end

        function [ene, em, pm_s, pm_a, kernel, norm_val, source] = sparse_components(res, stage)
            assert(isstruct(res) && ~isfield(res, 'error'), 'EELS:InvalidFitResult', ...
                'A valid fitted spectrum result is required.');
            stage = lower(string(stage));
            ene = EELS_ZEST_DataOps.field_or(res, 'ene_corrected', ...
                EELS_ZEST_DataOps.field_or(res, 'ene_raw', []));
            norm_val = EELS_ZEST_DataOps.field_or(res, 'norm_val', 1);
            assert(~isempty(ene) && isfinite(norm_val) && norm_val > 0, ...
                'EELS:InvalidFitResult', 'Result energy axis or normalization is invalid.');
            if stage == "stage2"
                assert(isfield(res, 'stage2') && isstruct(res.stage2) && ...
                    isfield(res.stage2, 'EM_final') && isfield(res.stage2, 'PM_stokes') && ...
                    isfield(res.stage2, 'PM_antistokes'), 'EELS:MissingStage2Source', ...
                    'Result does not contain Stage 2 sparse phonon source intensities.');
                em = double(res.stage2.EM_final(:));
                pm_s = double(res.stage2.PM_stokes(:));
                pm_a = double(res.stage2.PM_antistokes(:));
                source = 'stage2';
            else
                assert(isfield(res, 'EM_final') && isfield(res, 'PM_stokes') && ...
                    isfield(res, 'PM_antistokes'), 'EELS:MissingStage1Source', ...
                    'Result does not contain Stage 1 sparse phonon source intensities.');
                em = double(res.EM_final(:));
                pm_s = double(res.PM_stokes(:));
                pm_a = double(res.PM_antistokes(:));
                source = 'stage1';
            end
            assert(numel(em) == numel(pm_s) && numel(em) == numel(pm_a), ...
                'EELS:InvalidSourceGrid', 'Mode-energy and source-intensity sizes differ.');
            kernel = [];
            if stage == "stage2" && isfield(res, 'stage2') && isstruct(res.stage2) && ...
                    isfield(res.stage2, 'zlp_norm')
                kernel = double(res.stage2.zlp_norm(:));
            elseif isfield(res, 'stage1') && isstruct(res.stage1) && isfield(res.stage1, 'zlp_norm')
                kernel = double(res.stage1.zlp_norm(:));
            elseif isfield(res, 'zlp_only')
                kernel = double(res.zlp_only(:)) ./ norm_val;
            end
            assert(numel(kernel) == numel(ene), 'EELS:MissingSignalKernel', ...
                'Result does not contain a fitted-ZLP kernel on the result energy axis.');
            kernel(~isfinite(kernel) | kernel < 0) = 0;
            kernel = kernel ./ max(sum(kernel), eps);
            ene = double(ene(:));
        end

        function multiplier = sta_multiplier(energy_meV, temperature_K)
            energy_meV = double(energy_meV(:));
            temperature_K = double(temperature_K);
            assert(all(isfinite(temperature_K(:)) | isnan(temperature_K(:))) && ...
                all(temperature_K(isfinite(temperature_K)) > 0), ...
                'EELS:InvalidTemperature', 'STA temperatures must be positive.');
            kb = 0.08617333262;
            multiplier = 0.5 .* energy_meV .* ...
                (1 - exp(-energy_meV ./ max(kb .* temperature_K, eps)));
            multiplier(~isfinite(multiplier)) = NaN;
        end

        function [temperature_map, fallback_mask, invalid_mask] = ...
                resolve_sta_temperature(source, map_size, manual_T, T_stage1, T_stage2)
            source = lower(char(string(source)));
            assert(numel(map_size) == 2 && all(map_size > 0), 'EELS:InvalidMapSize', ...
                'STA requires two positive spatial dimensions.');
            map_size = double(map_size(:).');
            fallback_mask = false(map_size);
            switch source
                case 'manual'
                    assert(isscalar(manual_T) && isfinite(manual_T) && manual_T > 0, ...
                        'EELS:InvalidTemperature', 'Enter a positive manual temperature.');
                    temperature_map = repmat(double(manual_T), map_size);
                case 'stage1'
                    assert(isequal(size(T_stage1), map_size), 'EELS:TemperatureMapSizeMismatch', ...
                        'Stage 1 temperature map does not match the Data page.');
                    temperature_map = double(T_stage1);
                case 'stage2_fallback_stage1'
                    assert(isequal(size(T_stage2), map_size), 'EELS:TemperatureMapSizeMismatch', ...
                        'Stage 2 temperature map does not match the Data page.');
                    temperature_map = double(T_stage2);
                    bad = ~isfinite(temperature_map) | temperature_map <= 0;
                    if isequal(size(T_stage1), map_size)
                        use = bad & isfinite(T_stage1) & T_stage1 > 0;
                        temperature_map(use) = T_stage1(use);
                        fallback_mask = use;
                    end
                otherwise
                    error('EELS:InvalidTemperatureSource', 'Unknown STA temperature source: %s', source);
            end
            invalid_mask = ~isfinite(temperature_map) | temperature_map <= 0;
        end

        function [cube_out, ene_out, info] = apply_sta_cube(cube, ene, temperature_map)
            ene = double(ene(:));
            assert(ndims(cube) == 3 && size(cube, 3) == numel(ene), ...
                'EELS:InvalidDatacube', 'Datacube energy dimension does not match the energy axis.');
            [ny, nx, ~] = size(cube);
            if isscalar(temperature_map)
                temperature_map = repmat(double(temperature_map), ny, nx);
            else
                temperature_map = double(temperature_map);
            end
            assert(isequal(size(temperature_map), [ny, nx]), 'EELS:TemperatureMapSizeMismatch', ...
                'Temperature map must match the datacube spatial dimensions.');
            loss = ene > 0;
            assert(any(loss), 'EELS:NoLossChannels', 'No positive-energy loss channels are available.');
            ene_out = ene(loss);
            invalid = ~isfinite(temperature_map) | temperature_map <= 0;
            cube_out = NaN(ny, nx, numel(ene_out));
            bytes_per_row = max(1, nx .* numel(ene_out) .* 8 .* 3);
            rows_per_block = max(1, floor((128 .* 1024 .^ 2) ./ bytes_per_row));
            energy3 = reshape(ene_out, 1, 1, []);
            for y1 = 1:rows_per_block:ny
                y2 = min(ny, y1 + rows_per_block - 1);
                Tblock = temperature_map(y1:y2, :);
                factor = 0.5 .* energy3 .* ...
                    (1 - exp(-energy3 ./ max(0.08617333262 .* Tblock, eps)));
                block = double(cube(y1:y2, :, loss)) .* factor;
                bad = invalid(y1:y2, :);
                if any(bad(:))
                    tmp = reshape(block, numel(bad), []);
                    tmp(bad(:), :) = NaN;
                    block = reshape(tmp, size(block));
                end
                cube_out(y1:y2, :, :) = block;
            end
            info = struct('temperature_map_K', temperature_map, ...
                'valid_temperature_mask', ~invalid, 'loss_only', true, ...
                'formula', 'I_STA(E,T)=I(E)*E/2*(1-exp(-E/(kB*T)))', ...
                'kB_meV_per_K', 0.08617333262, 'rows_per_block', rows_per_block);
        end

        function [T_stage1, T_stage2, valid_stage2] = temperature_maps(stage1_cube, stage2_cube, map_size)
            if nargin < 3 || isempty(map_size)
                map_size = EELS_ZEST_DataOps.result_map_size(stage2_cube);
                if isempty(map_size), map_size = EELS_ZEST_DataOps.result_map_size(stage1_cube); end
            end
            assert(numel(map_size) == 2 && all(map_size > 0), 'EELS:MissingResultMapSize', ...
                'Could not determine fitted-result map dimensions.');
            T_stage1 = NaN(map_size); T_stage2 = NaN(map_size);
            T_stage1 = EELS_ZEST_DataOps.fill_temperature_map(T_stage1, stage1_cube, 'T_soft');
            T_stage2 = EELS_ZEST_DataOps.fill_temperature_map(T_stage2, stage2_cube, 'T_fit');
            valid_stage2 = isfinite(T_stage2) & T_stage2 > 0;
        end

        function context = fit_context_from_info(info, map_size)
            context = struct();
            if ~(isstruct(info) && isfield(info, 'eels_zest_fit_context') && ...
                    isstruct(info.eels_zest_fit_context))
                return;
            end
            saved = info.eels_zest_fit_context;
            map_size = double(map_size(:).');
            if numel(map_size) ~= 2 || any(map_size < 1)
                return;
            end
            context.map_size = map_size;
            if isfield(saved, 'result_file') && (ischar(saved.result_file) || isstring(saved.result_file))
                context.result_file = char(saved.result_file);
            else
                context.result_file = '';
            end
            map_fields = {'T_stage1_soft', 'T_stage2_fit', 'stage2_valid_mask'};
            for ii = 1:numel(map_fields)
                name = map_fields{ii};
                if isfield(saved, name) && isequal(size(saved.(name)), map_size)
                    context.(name) = saved.(name);
                end
            end
            if isfield(saved, 'source_stage')
                context.source_stage = saved.source_stage;
            end
        end

        function info = append_history(info, operation, parameters)
            if nargin < 1 || ~isstruct(info), info = struct(); end
            if nargin < 3 || isempty(parameters), parameters = struct(); end
            entry = struct('operation', char(operation), ...
                'timestamp', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
                'parameters', parameters);
            if ~isfield(info, 'processing_history') || isempty(info.processing_history)
                info.processing_history = entry;
            else
                history = info.processing_history;
                if iscell(history), history = [history{:}]; end
                info.processing_history = history;
                info.processing_history(end + 1) = entry;
            end
        end

        function obj = make_eels(cube, ene, info)
            if nargin < 3 || isempty(info), info = struct(); end
            obj = EELS(cube, ene, info);
        end

        function save_eels(filename, cube, ene, info)
            eels_data = EELS_ZEST_DataOps.make_eels(cube, ene, info); %#ok<NASGU>
            save(filename, 'eels_data', '-v7.3');
        end
    end

    methods (Static, Access = private)
        function value = field_or(s, name, fallback)
            value = fallback;
            if isstruct(s) && isfield(s, name) && ~isempty(s.(name)), value = s.(name); end
        end

        function map_size = result_map_size(cuberes)
            map_size = [];
            if ~isstruct(cuberes), return; end
            names = {'R2_map', 'fit_status_map', 'T_map_raw'};
            for ii = 1:numel(names)
                if isfield(cuberes, names{ii}) && ismatrix(cuberes.(names{ii}))
                    map_size = size(cuberes.(names{ii})); return;
                end
            end
        end

        function map = fill_temperature_map(map, cuberes, field_name)
            if ~(isstruct(cuberes) && isfield(cuberes, 'results')), return; end
            results = cuberes.results;
            for ii = 1:numel(results)
                if iscell(results), res = results{ii}; else, res = results(ii); end
                if ~(isstruct(res) && ~isfield(res, 'error')), continue; end
                value = NaN;
                if strcmp(field_name, 'T_fit') && isfield(res, 'stage2') && ...
                        isstruct(res.stage2) && isfield(res.stage2, 'T_fit')
                    value = res.stage2.T_fit;
                elseif isfield(res, field_name)
                    value = res.(field_name);
                end
                if isfield(cuberes, 'sample_pixels') && ii <= numel(cuberes.sample_pixels) && ...
                        isfinite(cuberes.sample_pixels(ii)) && cuberes.sample_pixels(ii) >= 1 && ...
                        cuberes.sample_pixels(ii) <= numel(map)
                    [yy, xx] = ind2sub(size(map), round(cuberes.sample_pixels(ii)));
                elseif isfield(res, 'pixel_y') && isfield(res, 'pixel_x')
                    yy = round(res.pixel_y); xx = round(res.pixel_x);
                elseif numel(results) == numel(map)
                    [yy, xx] = ind2sub(size(map), ii);
                else
                    continue;
                end
                if yy >= 1 && xx >= 1 && yy <= size(map, 1) && xx <= size(map, 2) && ...
                        isfinite(value) && value > 0
                    map(yy, xx) = value;
                end
            end
        end

    end
end
