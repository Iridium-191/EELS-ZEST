classdef EELS_ZEST_Export
    %EELS_ZEST_EXPORT File writers shared by GUI export actions.

    methods (Static)
        function save_payload(filename, payload, format, variable_name)
            if nargin < 3 || isempty(format)
                [~, ~, ext] = fileparts(filename);
                format = erase(lower(ext), '.');
            end
            if nargin < 4 || isempty(variable_name), variable_name = 'eels_zest_export'; end
            format = lower(char(string(format)));
            switch format
                case 'mat'
                    S = struct();
                    S.(matlab.lang.makeValidName(variable_name)) = payload;
                    save(filename, '-struct', 'S', '-v7.3');
                case 'csv'
                    table_out = EELS_ZEST_Export.payload_table(payload);
                    writetable(table_out, filename);
                otherwise
                    error('EELS:UnsupportedExportFormat', 'Unsupported payload format: %s', format);
            end
        end

        function save_figure(filename, fig, format)
            assert(isgraphics(fig, 'figure'), 'EELS:InvalidFigure', 'A valid figure is required.');
            if nargin < 3 || isempty(format)
                [~, ~, ext] = fileparts(filename); format = erase(lower(ext), '.');
            end
            switch lower(char(string(format)))
                case 'fig'
                    savefig(fig, filename);
                case 'svg'
                    try
                        exportgraphics(fig, filename, 'ContentType', 'vector');
                    catch
                        print(fig, filename, '-dsvg');
                    end
                case 'png'
                    exportgraphics(fig, filename, 'Resolution', 360);
                case {'tif', 'tiff'}
                    exportgraphics(fig, filename, 'Resolution', 600);
                otherwise
                    error('EELS:UnsupportedExportFormat', 'Unsupported figure format: %s', format);
            end
        end

        function T = payload_table(payload)
            if istable(payload), T = payload; return; end
            if isnumeric(payload) && ismatrix(payload)
                T = array2table(payload); return;
            end
            if isstruct(payload) && isfield(payload, 'ene') && isfield(payload, 'spectrum')
                ene = double(payload.ene(:));
                values = payload.spectrum;
                if iscell(values), values = cat(2, values{:}); end
                values = double(values);
                if size(values, 1) ~= numel(ene), values = values.'; end
                assert(size(values, 1) == numel(ene), 'EELS:CSVShapeMismatch', ...
                    'Spectrum size does not match its energy axis.');
                T = array2table([ene, values], 'VariableNames', ...
                    [{'energy_meV'}, compose('spectrum_%d', 1:size(values, 2))]);
                return;
            end
            error('EELS:CSVPayloadUnsupported', ...
                'This payload is not naturally representable as a CSV table.');
        end
    end
end
