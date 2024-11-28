function varargout = multiscatter(X,Y,N,M,varargin)
% MULTISCATTER(X,Y,N,M,..,NAME,VAL,...) - Generate N x M tiled scatterplots of the data in X, Y,
%   using common axis limits for all plots, and ommitting x/y-ticks except on the left- and lower-
%   edges.
%
%   Everything will be plot on the current figure handle, GCF/CLF beforehand, if required.
%
% H = MULTISCATTER([],[],N,M,...) - create and label empty subplots.
%
% INPUT:
%   N,M - number of rows and columns of SUBPLOT(N,M)
%   X,Y - Equal-sized numerical arrays. Each column X(:,j),Y(:,j) will be plotted on SUBPLOT(N,M,j)
%       If the number of columns SIZE(X,2) > M*N, the rest will be ignored.
%   
%   ..,'limits',[XLO,XHI,YLO,YHI] - Axis limits. Default is min/max for all data.
%   ..,'titles',C - provide a size(X,2) cellstring of titles, one for for each SUBPLOT
%   ..,'xlabel',A,'ylabel',B,'grid','on','size',S,'filled',F,.. - have their usual meaning. The
%       only difference is they will be set equal for all subplots.
%   ..,'oversize',[1.2 1.2] - Multiplies axes dimensions [W,H] to reduce margins between tiles.
%   ..,'symmetric',true - assume X is meant to be equal to Y plot identity line.
%   ..,'errors',true - (default if symmetric) write error metrics NSTD(y-x) and NMBE(y-x)
%   ..,'densityplot',true,NAME,VAL - Use DENSITYPLOT instead of SCATTER, and optionally provide
%       custom options (resolution, bandwith,..).
%
% EXAMPLE:
%
%     X = rand(1000,6);               % 6 test cases of an independent variable X
%     Y = X.^2 + randn(1000,6)*0.05;  % Measurements results for the processes Y(X)
%     Z = X.^2;                       % A model that is meant to describe the data
% 
%     % Plot independent X vs measured Y on a 2x3 tile, using DENSITYPLOT
%     figure(1); clf();
%     multiscatter(X,Y,2,3,'densityplot',true,'xlabel','X','ylabel','Y');
% 
%     % Plot model Z vs measurements Y on a 2x3 tile
%     figure(2); clf();
%     multiscatter(Y,Z,2,3,'symmetric',true,'xlabel','Measured','ylabel','Model');
% 
%     % Plot only the first column of Y vs Z
%     figure(3); clf();
%     multiscatter(Y(:,1),Z(:,1),1,1,'densityplot',true,'symmetric',true,'size',5,...
%        'limits',[0 1 0 1],'xlabel','Measured','ylabel','Model','titles',{'Test Case 1'});
%
% See also: DENSITYPLOT, SCATTER, SUBPLOT

assert(M > 0 && N > 0 && all(mod([M,N],1) == 0),'Expecting nonzero integers M, N');

[X,Y] = compatiblesize(X,Y);
dummy = isempty(X);
if dummy
    MN = M*N;
else
    MN = size(X,2);
end
assert(MN <= M*N);

% Set defaults
OPT.errors = true;
OPT.units = '%';
OPT.symmetric = true;
OPT.xlabel = 'X';
OPT.ylabel = 'Y';
OPT.limits = [];
OPT.titles = arrayfun(@num2str,1:MN,'unif',0);
OPT.grid = 'on';

OPT.size = 2;
OPT.color = [0.1,0.5,1];
OPT.filled = true;
OPT.oversize = [1.2 1.2];
if size(X,2) <= 1, OPT.oversize(:) = 1; end

OPT.densityplot = false;
OPT.resolution = 100;
OPT.maxpts = 10000;
OPT.bandwidth = [];
OPT.cmap = [];
    
% Parse name,value pairs. Anything not recognized will be passed directly to SCATTER
[OPT,varargin,isdef] = getpairedoptions(varargin,OPT);

if ~dummy && isempty(OPT.limits)
   OPT.limits = [min(X(:)),max(X(:)),min(Y(:)),max(Y(:))];
   if OPT.symmetric
       OPT.limits = repmat([min(OPT.limits([1,3])),max(OPT.limits([2,4]))],1,2);
   end
end
OPT.limits = double(OPT.limits);
OPT.symmetric = OPT.symmetric && ~isempty(OPT.limits);

% Collect Density-Plot args
if OPT.densityplot
    OPT.args = {'resolution','maxpts','bandwidth','cmap','limits'};
    OPT.args(2,:) = cellfun(@(f) OPT.(f),OPT.args,'unif',0);
    OPT.args = OPT.args(:)';
    if OPT.filled, OPT.args{end+1} = 'filled'; end
end

if OPT.filled, OPT.filled = {'filled'}; end

OPT.errors = ~dummy && OPT.symmetric && OPT.errors;
if OPT.errors
    f = isfinite(X) & isfinite(Y);
    X(~f) = NaN;
    Y(~f) = NaN;

    MBE = mean(Y-X,1,'omitnan');
    STD = std(Y-X,1,'omitnan');
    OPT.units = strtrim(OPT.units);
    if OPT.units == '%'
        AVG = mean(X,1,'omitnan');
        MBE = MBE./AVG*100;
        STD = STD./AVG*100;
        tags = {'nMBE','nSTD'};
    else
        OPT.units = [' ' OPT.units];
        tags = {'MBE','STD'};
    end
end

if ~iscell(OPT.xlabel) || (isscalar(OPT.xlabel) && M > 1), OPT.xlabel = repmat({OPT.xlabel},1,M); end
if ~iscell(OPT.ylabel) || (isscalar(OPT.ylabel) && N > 1), OPT.ylabel = repmat({OPT.ylabel},1,N); end

parsestruct(OPT,{'oversize','size'},'-n','-Z','-p','-r');
if isscalar(OPT.oversize), OPT.oversize(2) = OPT.oversize; end
assert(numel(OPT.oversize) == 2,'Expecting 1-2 valued OPT.oversize');
K = OPT.oversize-1;
K = [eye(2),zeros(2);-diag(K)/2,diag(1+K)]; % centered scaling matrix

if ~isdef.size, compatiblesize(OPT.size,X); end

if ~isdef.color
    if OPT.densityplot
        warning('Ignoring color information');
    else
        try
            switch size(OPT.color,1)
            case {1,size(X,1)}, assert(any(size(OPT.color,2) == [1,3,MN]));
            case {numel(X),MN}, assert(any(size(OPT.color,2) == [1,3]));
            end
        catch
            error('Incompatible color size'); 
        end
    end
end


gcf();
H = arrayfun(@(j) subplot(N,M,j),1:MN);

for i = 1:MN
    [c,r] = ind2sub([M,N],i);
    % j = sub2ind([M,N],r,c);
    % H(j) = subplot(N,M,i); hold on; 
    j = i;
    hold(H(j),'on');
    if ~isempty(OPT.limits), axis(H(j),OPT.limits); end
    
    if ~isempty(OPT.titles), title(H(j),OPT.titles{i}); end
    if c == 1
        ylabel(H(j),OPT.ylabel{r});
    else
        yticklabels(H(j),{});
    end
    
    if r == N
        xlabel(H(j),OPT.xlabel{c});
    else
        xticklabels(H(j),{});
    end
    set(H(j),'position',get(H(j),'position')*K);
    
    if ~dummy
        if OPT.densityplot
            try
                densityplot(X(:,j),Y(:,j),OPT.size,'ax',H(j),OPT.args{:});
            catch
                warning('densityplot(X(:,j),Y(:,j) error at j = %d',j);
                scatter(H(j),X(:,j),Y(:,j),OPT.size,OPT.filled{:});
            end
        else
            if size(OPT.size,2) == MN, s = OPT.size(:,j); else, s = OPT.size; end
            if size(OPT.color,2) == MN
                c = OPT.color(:,j); 
            elseif size(OPT.color,2) == 3 && size(OPT.color,1) == numel(X)
                c = OPT.color((1:size(X,1)) + (j-1)*size(X,1),:);
            else
                c = compatiblesize(OPT.color,X(:,j)); 
            end
            scatter(H(j),X(:,j),Y(:,j),s,c,OPT.filled{:});
        end

        if OPT.errors
            txt = sprintf('\n %s: %0.2f%s\n %s: %0.2f%s\n',tags{2},STD(j),OPT.units,tags{1},MBE(j),OPT.units);
            text(H(j),0,OPT.limits(2),txt,'fontsize',10,'verticalalignment','top');
        end
    end
    
    if OPT.symmetric
        line(H(j),OPT.limits(1:2),OPT.limits(3:4),'color',[1 1 1]*0.8);
    end
    grid(H(j),OPT.grid);
    
    if ~isempty(varargin), set(H(j),varargin{:}); end
end
if nargout > 0, varargout{1} = H; end
