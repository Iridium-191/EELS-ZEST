function imported = DM3Import(dm3filename)
%DM3IMPORT Import a DM3 file into the lightweight EELS container.
%
%   obj = DM3Import(filename) returns an EELS object with fields:
%       obj.data  Ny x Nx x Ne datacube or 1 x 1 x Ne spectrum
%       obj.ene   Ne x 1 energy axis in meV when calibration is available
%       obj.info  scale/status metadata
%
if nargin < 1 || isempty(dm3filename)
    error('DM3Import:MissingFile', 'A DM3 filename is required.');
end

raw = DM3Import0(dm3filename);
if isempty(raw)
    error('DM3Import:ReadFailed', 'DM3Import0 returned an empty result for %s.', dm3filename);
end

data = [];
if isfield(raw, 'data') && ~isempty(raw.data)
    data = raw.data;
elseif isfield(raw, 'image_data') && ~isempty(raw.image_data)
    data = raw.image_data;
elseif isfield(raw, 'spectra_data') && ~isempty(raw.spectra_data)
    data = raw.spectra_data{1};
end
if isempty(data)
    error('DM3Import:NoData', 'No image or spectrum data were found in %s.', dm3filename);
end

if isfield(raw, 'ene') && ~isempty(raw.ene)
    ene = raw.ene(:);
elseif isfield(raw, 'zaxis') && isstruct(raw.zaxis) && isfield(raw.zaxis, 'scale')
    ne = size(data, ndims(data));
    origin = 0;
    if isfield(raw.zaxis, 'origin'), origin = raw.zaxis.origin; end
    ene = ((1:ne).' - origin) .* raw.zaxis.scale;
elseif isfield(raw, 'xaxis') && isstruct(raw.xaxis) && isfield(raw.xaxis, 'scale') && isvector(data)
    ne = numel(data);
    origin = 0;
    if isfield(raw.xaxis, 'origin'), origin = raw.xaxis.origin; end
    ene = ((1:ne).' - origin) .* raw.xaxis.scale;
else
    ne = size(data, ndims(data));
    ene = (0:(ne - 1)).';
end

info = struct();
info.source_file = char(dm3filename);
info.importer = 'DM3Import0 compatibility parser';
info.xscale = localField(raw, 'xscale', 1);
info.yscale = localField(raw, 'yscale', 1);
info.escale = localField(raw, 'escale', EELS.defaultEscale(ene));
info.isaligned = false;
info.isnormalized = false;
info.issaturated = false;

imported = EELS(double(data), double(ene(:)), info);
end

function val = localField(S, name, fallback)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    val = S.(name);
else
    val = fallback;
end
end
