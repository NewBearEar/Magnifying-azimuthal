function varargout = variable_ortho(varargin)
mproj.default = @orthoDefault;
mproj.forward = @orthoFwd;
mproj.inverse = @orthoInv;
mproj.auxiliaryLatitudeType = 'geodetic';
mproj.classCode = 'Azim';

varargout = applyMagnifyingAzimuthalProjection(mproj, varargin{:});

%--------------------------------------------------------------------------

function mstruct = orthoDefault(mstruct)

[mstruct.trimlat, mstruct.trimlon] ...
          = fromDegrees(mstruct.angleunits, [-Inf 89], [-180 180]);
mstruct.mapparallels = [];
mstruct.nparallels   = 0;
mstruct.fixedorient  = [];

%--------------------------------------------------------------------------

function [x, y] = orthoFwd(mstruct, rng, az)

a = ellipsoidprops(mstruct);
r = a * sin(rng);

x = r .* sin(az);
y = r .* cos(az);

%--------------------------------------------------------------------------

function [rng, az] = orthoInv(mstruct, x, y)

a = ellipsoidprops(mstruct);

rng = asin(hypot(x,y) / a);
az = atan2(x,y);
