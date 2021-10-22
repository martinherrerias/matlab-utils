function S = meteofilename(Loc,t_start,t_end,dt,varargin)
% S = METEOFILENAME(LOC,A,B,DT) - Generate a string of the form: 
%
%   PREFIX_LAT_LON_PDT_INTERVAL[SUFFIX], e.g. 'name_21N08_104W61_PT60M_2018'[...]
%
% Where PREFIX (if any) is taken from LOC.name, or overriden by 'prefix',X; LAT, LON are 
% abbreviated decimal degrees e.g. 21.08°N, -104.61°E -> 21N08_104W61; PDT is ISODURATION(DT), 
% and INTERVAL is PRETTYINTERVAL(A,B,..).
%
% 'decimals',N - changes the number of decimals for LAT,LON representations
% 'suffix',Y - concatenates a free string Y at the end of S, e.g. a file extension
% 'tol',X - pass optional tolerance to PRETTYINTERVAL

    parselocation(Loc);
    if ~isfield(Loc,'name'), Loc.name = ''; end

    opt.decimals = 2;
    opt.suffix = '';
    opt.prefix = Loc.name;
    opt.tol = [];
    opt = getpairedoptions(varargin,opt);
    
    if isfield(Loc,'TimeZone'), TZ = Loc.TimeZone; else, TZ = 'keep'; end

    t_start = parsetime(t_start,'TimeZone',TZ);
    t_end = parsetime(t_end,'TimeZone',TZ);
    
    interval = prettyinterval(t_start,t_end,opt.tol);
    dtlbl = isoduration(dt);
    
    if ~isempty(opt.prefix)
        opt.prefix = [matlab.lang.makeValidName(opt.prefix) '_'];
    end
    
    if opt.decimals > 0
        fmt = ['%d%c%0' num2str(opt.decimals) 'd'];
        K = 10^opt.decimals;
        deg2dd = @(x,C) sprintf(fmt,floor(abs(x)),(x<0)*(C(2)-C(1))+C(1),round(K*rem(abs(x),1)));  
    else
        deg2dd = @(x,C) sprintf('%d%c',round(abs(x)),(x<0)*(C(2)-C(1))+C(1)); 
    end
    loclbl = [deg2dd(Loc.latitude,'NS') '_' deg2dd(Loc.longitude,'EW')];

    S = [opt.prefix loclbl '_' dtlbl '_' interval opt.suffix];
end
