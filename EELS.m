classdef EELS
    %EELS Lightweight data container for EELS spectra and datacubes.

    properties
        data {mustBeNumeric}
        ene {mustBeNumeric, mustBeReal}
        info = struct()
    end

    methods
        function obj = EELS(data, ene, info)
            if nargin >= 1
                obj.data = data;
            end
            if nargin >= 2
                obj.ene = ene(:);
            end
            if nargin >= 3 && ~isempty(info)
                obj.info = info;
            end
            obj.info = EELS.fillInfoDefaults(obj.info, obj.ene);
        end

        function dim = xdim(obj)
            dim = size(obj.data, 2);
        end

        function dim = ydim(obj)
            dim = size(obj.data, 1);
        end

        function dim = edim(obj)
            dim = size(obj.data, 3);
        end

        function scale = xscale(obj)
            scale = EELS.infoValue(obj.info, 'xscale', 1);
        end

        function scale = yscale(obj)
            scale = EELS.infoValue(obj.info, 'yscale', 1);
        end

        function scale = escale(obj)
            scale = EELS.infoValue(obj.info, 'escale', EELS.defaultEscale(obj.ene));
        end

        function tf = isaligned(obj)
            tf = EELS.infoValue(obj.info, 'isaligned', false);
        end

        function tf = isnormalized(obj)
            tf = EELS.infoValue(obj.info, 'isnormalized', false);
        end

        function tf = issaturated(obj)
            tf = EELS.infoValue(obj.info, 'issaturated', false);
        end

        function sumspec = sumall(obj)
            data3 = obj.asDatacube();
            sumspec = squeeze(sum(sum(data3, 1, 'omitnan'), 2, 'omitnan'));
        end

        function idx = zlp(obj)
            [~, idx] = min(abs(obj.ene(:)));
        end

        function cube = asDatacube(obj)
            data0 = double(obj.data);
            if isvector(data0)
                cube = reshape(data0(:), 1, 1, []);
            elseif ismatrix(data0)
                if numel(obj.ene) == size(data0, 2)
                    cube = reshape(data0, size(data0, 1), 1, size(data0, 2));
                elseif numel(obj.ene) == size(data0, 1)
                    cube = reshape(data0.', 1, size(data0, 2), size(data0, 1));
                else
                    cube = reshape(data0, 1, size(data0, 1), size(data0, 2));
                end
            else
                cube = data0;
            end
        end

        function plot(obj)
            plot(obj.ene(:), obj.sumall());
            grid on;
            xlabel('Energy loss (meV)');
            ylabel('Intensity (a.u.)');
        end

        function semilogy(obj)
            semilogy(obj.ene(:), max(obj.sumall(), eps));
            grid on;
            xlabel('Energy loss (meV)');
            ylabel('Intensity (a.u.)');
        end
    end

    methods (Static)
        function obj = fromDM3(filename)
            obj = DM3Import(filename);
        end

        function info = fillInfoDefaults(info, ene)
            if nargin < 1 || isempty(info)
                info = struct();
            end
            if ~isfield(info, 'xscale'), info.xscale = 1; end
            if ~isfield(info, 'yscale'), info.yscale = 1; end
            if ~isfield(info, 'escale'), info.escale = EELS.defaultEscale(ene); end
            if ~isfield(info, 'isaligned'), info.isaligned = false; end
            if ~isfield(info, 'isnormalized'), info.isnormalized = false; end
            if ~isfield(info, 'issaturated'), info.issaturated = false; end
        end

        function val = infoValue(info, field, fallback)
            if isstruct(info) && isfield(info, field) && ~isempty(info.(field))
                val = info.(field);
            else
                val = fallback;
            end
        end

        function escale = defaultEscale(ene)
            if nargin < 1 || numel(ene) < 2
                escale = 1;
            else
                escale = median(diff(double(ene(:))), 'omitnan');
                if ~isfinite(escale) || escale == 0
                    escale = 1;
                end
            end
        end
    end
end
