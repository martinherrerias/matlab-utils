function [L,c,n,lb,P,ub] = findlag(x,y,lags,varargin)
% [L,C,N,LB,P] = FINDLAG(X,Y,LAGS) - Estimate the lag between two shifted- but otherwise similar
%   vectors X, Y, by calculating the cross-correlations C = corr( shift(X,LAGS(j)), Y ) (*)
%   and identifying L = LAGS(k) such that LB(k) > LB(j) for every other j. Here LB is a lower-
%   bound estimate for C, accounting for uncertainty in the calculated C. 
%
% (*) by default, shift(x,+n) = [ NaN,..,NaN ,x(1),x(2),..,x(end-n) ] and
%                 shift(x,-n) = [ x(n+1),x(n+2),..,x(end),NaN,..,NaN ]
% [..] = FINDLAG(..'-periodic') makes shift(x,n) equivalent to circshift(x,n)
%
%   The use of max(LB) instead of max(C) provides robustness when working with vectors X, Y that 
%   contain many NaNs (e.g. when using large shifts). Spurious correlation values can be calculated 
%   for lags at which very few valid samples of Y and shift(X,LAGS(j)) happen to align. Introducing
%   a variance proportional to 1/(n-3) penalizes these lags, and favors those with a good data 
%   overlap. Use ..,'P', 0, .. to avoid this behavior (see below).
%
%   P(j) is the probability that C(k) > C(j) for k ~= j, i.e. that the correlation found for the
%   'best lag' L = LAGS(k) is not statistically different than that of each LAGS(j). 
%   NOTE: for k == j, the trivial P( C(k) > C(k) ) is replaced by P(k) = min( P(j) | j~=k) ),
%   that is, the probability that the correlation for the 'best' lag is not statistically
%   different from that of the 'next-best' lag.
%
% [..] = FINDLAG(..'P',P) - adjust the probability for confidence interval LB. The default is
%   P = 0.95, i.e. LB is a lower 95% confidence bound for C.
%   In general,  LB = tanh(atanh(C) - norminv((1+P)/2) / sqrt(n-3) )    (see [1])
%
% INPUT 
%   X, Y : equal-size vectors
%   LAGS: n-vector of integers
% OUTPUT
%   L : scalar integer, argmax(k) LB(k). L > 0 means Y is ahead of X by L steps.
%   C : n-vector, correlation coefficients corresponding to each lag
%   N : n-vector, number of finite, overlaping samples for each lag
%   LB: n-vector, lower-bound estimate for C
%   P : n-vector, probability that C(k) > C(j) | j ~= k, and P(k) = min( P(j) | j~=k) )
%
% TODO: search method when lags = [], maybe on frequency space? (see PLOMB)
%
% REF: [1] Wasserman, "All of Statistics" (2004) p.242
%
% See also: CHECKUTCOFFSET, CIRCSHIFT, CORR

    if nargin == 0, test(); return; end
    
    narginchk(3,9);
    [opt,varargin] = getflagoptions(varargin,{'-periodic','-plot'});
    opt.figh = [];
    opt.P = 0.95;
    opt = getpairedoptions(varargin,opt,'restchk');
    
    validateattributes(x,{'numeric'},{'vector','real'});
    validateattributes(x,{'numeric'},{'vector','real','numel',numel(y)});
    validateattributes(lags,{'numeric'},{'vector','real','integer'});
    validateattributes(opt.P,{'numeric'},{'vector','real','scalar','>=',0,'<=',1});
    assert(isequal(lags,unique(lags,'stable')),'Lags must be unique');
    
    if size(x,2) > 1, x = x'; y = y'; end
    x = double(x); y = double(y);
    fx = isfinite(x);
    fy = isfinite(y);
    
    [c,n] = arrayfun(@lagged_corr,lags);
    assert(any(n > 3),'Not enough samples');
    
    if opt.P ~= 0
    % Confidence interval using Fisher-transform [1]
        z = atanh(c);
        s = norminv((1+opt.P)/2)./sqrt(max(0,n-3));
        lb = tanh(z-s);
        ub = tanh(z+s);
    else
        lb = c;
        ub = c;
    end
       
    [~,best] = max(lb); 
    L = lags(best);
    
    if nargout > 1 || opt.plot
    % Probability that c(idx) > c(j) for every other j, using Fischer Z test
        % s = sqrt(1./(n(best)-3)+1./(n-3));
        s = sqrt(1./max(0,n(best)-3))+sqrt(1./max(0,n-3)); % c ~ 1
        z = (atanh(c(best))-atanh(c))./s;
        % P = 2*normcdf(z)-1;
        P = normcdf(z);
        P(best) = min(P(lags ~= L));
    end
    
    if opt.plot, findlag_plot(); end

    function [c,n] = lagged_corr(lag)
        
        xs = circshift(x,lag);
        useful = circshift(fx,lag) & fy;
        if ~opt.periodic
            if lag > 0, useful(1:lag) = 0; elseif lag < 0, useful(end+lag+1:end) = 0; end
        end
        n = nnz(useful);
        if n < 2, c = NaN; return; end
        c = corr(xs(useful),y(useful));
    end

    function findlag_plot()
        P0 = 1 - P(lags == 0);
        if isempty(P0), P0 = NaN; end

        if isempty(opt.figh) || ~ishandle(opt.figh) || ~isa(opt.figh,'matlab.ui.Figure')
            opt.figh = GUIfigure('findlag','findlag-plot','3:1'); 
        end
        figure(opt.figh); clf();
        ax = subplot(1,3,1:2); hold on;
        ax2 = subplot(1,3,3); hold on;
        
        title(ax,sprintf('detected lag = %d, P_{0} = %0.0f%%, P_B = %0.0f%%',...
            L,P0*100, P(best)*100));
        plot(ax,x);
        plot(ax,y);
        plot(ax,circshift(x,L));
        legend(ax,'x','y',sprintf('shift(x,%d)',L),'box','off')
        
        grid(ax2,'on');
        plot(ax2,lags,c,'r');
        plot(ax2,lags,lb,'r:');
        plot(ax2,lags,ub,'r:');
        ylabel('correlation');
        xlabel('lag');
        
        yyaxis(ax2,'right');
        set(ax2,'yscale','log')
        plot(ax2,lags,P,'o-');
        ylabel(sprintf('P( corr(%d) > corr(x) )',L));
    end
end

function test()
    N = 3;  
    m = 100;
    offset = randi(10)-5;
    noise = 0.1;
    missing = 0.2;
    
    x = sin(linspace(0,2*N*pi,m*N));
    y = circshift(x,offset) + randn(size(x))*noise;
    
    x(rand(size(x)) < missing) = NaN;
    y(rand(size(y)) < missing) = NaN;

    lags = max([5,abs(offset)*2]);
    lags = -lags:lags;
    findlag(x,y,lags,'-plot');
    h = findobj('type','axes');
    h = h(arrayfun(@(h) contains(h.Title.String,'detected lag ='),h));
    h.Title.String = [sprintf('true lag = %d, ',offset),h.Title.String];
end
