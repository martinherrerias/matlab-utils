function s = sph2cartV(az,el)
% S = SPH2CARTV(AZ,EL) - Transform spherical to Cartesian coordinates on a unit circle,
%   returning an N·3 array of coordinates S = [X,Y,Z] for N-vectors of azimuth and elevation
%   angles (degrees!) AZ, EL.
%
%   See also SPH2CART

    if size(az,2)~=1, az = az(:); end
    if size(el,2)~=1, el = el(:); end
    az = double(az);
    el = double(el);
    
    z = sind(el);
    if isscalar(z) && numel(az) > 1, z = repelem(z,numel(az),1); end
    
    s = [cosd(az).*cosd(el),sind(az).*cosd(el),z];
end