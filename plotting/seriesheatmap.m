function varargout = seriesheatmap(t,z,varargin)
% SERIESHEATMAP(t,z) - Plot time series (t,Z) as a color-coded image (IMAGESC), with pixel coords.
%   given by day and hour, i.e. x = floor(t), y = 24·mod(t,1), color ~ Z.
%
% SERIESHEATMAP(..,'clim',[LO,HI]) - Specify color scaling of IMAGESC.
%
% SERIESHEATMAP(..,'alphadata',A) - provide size(Z) alpha channel. Unless 'alphadata' is provided
%   explicitly, pixels with no data will be rendered transparent.
%
% [Z,x,y,F,ih] = SERIESHEATMAP(..,'name',VAL) - Provide additional arguments to IMAGESC or 
%   SERIES2DAYMATRIX, and/or get the output of both functions.
%
% INPUT:
%   t - N-vector of DATENUM timestamps 
%   z - N-vector of values
%   'name',VAL - optional arguments to IMAGESC or SERIES2DAYMATRIX.
%
% OUTPUT: 
%   Z, x, y, F - output of SERIES2DAYMATRIX
%   ih - image handle (not the same as axis/figure handle), returned by IMAGESC
%
% See also: SERIES2DAYMATRIX, IMAGESC

    narginchk(2,Inf);
    
    opt.alphadata = [];
    opt.clim = [];
    
    % options for PARSETIME
    opt.gridunits = [];
    opt.tolerance = 0.5;
    opt.interval = 'c';
    opt.step = [];
    opt.TimeZone = '';
    opt.InputFormat = '';
    
    if ~isempty(varargin) && isscalar(varargin{end}) && ishandle(varargin{end}) && ...
       isa(varargin{end},'matlab.graphics.axis.Axes')
        ax = varargin{end};
        varargin(end) = [];
    else
        ax = gca(); 
    end
    [opt,varargin] = getpairedoptions(varargin,opt);

    clim = opt.clim;
    if nargin < 3 || isempty(clim)
        clim = [min(z,[],'omitnan'),max(z,[],'omitnan')];
    end
    
    alpha = opt.alphadata;
    opt = rmfield(opt,{'clim','alphadata'});
    opt = [fieldnames(opt),struct2cell(opt)]';
    
    [Z,x,y,F] = series2daymatrix(t,z,opt{:});
    
    if isempty(alpha)
        alpha = ~isnan(Z)*1;
    else
        assert(numel(Z) == numel(alpha),'Inconsistent alphadata');
    end
    
    ih = imagesc(ax,x,y*24,Z,varargin{:},'AlphaData',alpha,clim);
    set(ax,'YDir','normal');
    datetick(ax,'x','keeplimits');

    if nargout > 0, varargout = {Z,x,y*24,F,ih}; end
end
