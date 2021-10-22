function [x,exitflag] = quickfzero(f,df,x0,tolx,MaxIter)
% X = QUICKFZERO(F,dF,X0,TOLX,MAXITER) - A minimalistic solver for element-wise functions of the 
%   form F(x) = 0, where X in an array of any size. With known function derivatives dF = F'(x), 
%   using Newton-Raphson method with seed X0, and ensuring tolerance TOLX(x) (note that TOLX is
%   expected to be a function, e.g. a result of PARSETOLERANCE).
%
%   QUICKFZERO allways iterates with size(X0) arrays, which makes it less efficient, but allows
%   function handles that use external size(X0) element-wise operations.
%
%   QUICKFZERO doesn't check dF for coherence. However, If the calculated step dx = -F(x)/dF(x) doesn't
%   improve results, 
%   For most cases, consider using FZERO / FSOLVE / LSQNONLIN.

    narginchk(3,5);
    if nargin < 4, [~,tolx] = parsetolerance(0); end
    if nargin < 5, MaxIter = getSimOption('MaxIter'); end
    
    exitflag = zeros(size(x0)); % 0 = Too many function evaluations or iterations
    
    xp = x0;
    fxp = f(xp);
    dydx = df(xp);
    dx = -fxp./dydx; % tentative correction
    
    iter = 0;
    
    % Termination criteria: fx = x (solution), x=xp / fx=fxp (no further changes in x/f(x))
    % while any(abs(cat(d,fx,fx-fxp,x-xp)) > tolx(x),'all')
    while iter <= MaxIter
        
        x = xp + dx;
        fx = f(x);
        dydx = df(x);
        tx = tolx(x);
        ty = max(eps(fx),abs(tx.*dydx));
        
        ok = abs(fx) < abs(fxp); % any improvement?
        toofar = ~ok;

        exitflag(abs(fx-fxp) <= ty) = 3;   % 3 = Change in residual small enough
        exitflag(abs(fx) <= ty) = 1;       % 1 = Error small enough
        exitflag(abs(dx) <= tx) = 2;       % 2 = Change in X small enough
        
        if all(exitflag > 0,'all')
            x(toofar) = xp(toofar);
            break; 
        end

        % back 1/2 a step, don't update last step
        dx(toofar) = dx(toofar)/2;

        % update last step, and take 1 step forward 
        xp(ok) = x(ok);
        fxp(ok) = fx(ok);
        dx(ok) = -fx(ok)./dydx(ok);

        iter = iter+1;
    end
    if iter > MaxIter
        error('quickfzero:maxiter','Could not converge to a solution'); 
    end
end