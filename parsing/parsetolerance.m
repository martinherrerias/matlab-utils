function [t,tx,ty] = parsetolerance(t,varargin)
% TOL = PARSETOLERANCE(T) - Return a unified tolerance vector [ABSTOLX,ABSTOLY,RELTOL] from
%   incomplete vector T, by completing/bounding with default values:
%
%  INPUT: T - incomplete tolerance vector:
%
%       [] - (empty), returns TOL = [mAT,mAT,dRT]
%       tR - scalar relative tolerance, returns TOL = [mAT,mAT,max(R,mRT)]
%       [tx,ty] - absolute tolerance in x, y. Returns [max(tx,mAT),max(ty,mAT),mRT]
%
%   Where default minimum absolute & relative tolerances (mAT,mRT) and relative tolerance dRT are 
%   set (persistently!) from global SimOptions.{'minabstol','minreltol','RelTol'}, or can be set
%   individually using name-value pairs, i.e.
%
%       TOL = PARSETOLERANCE(..,'reltol',dRT)
%       TOL = PARSETOLERANCE(..,'minabstol',mAT)
%       TOL = PARSETOLERANCE(..,'minreltol',mRT)
%
% OUTPUT: TOL = [ABSTOLX,ABSTOLY,RELTOL] (3-vector). The vector is meant to be used to determine
%   convergence based on RELTOL i.e. TOL(3) except when x~0 or y~0, where it should be bound by 
%   minimum absolute values [ABSTOLX,ABSTOLY].
%
% [TOL @TX @TY] = PARSETOLERANCE(T) - Additionally returns tolerance-function handles:
%
%       TX = @(x) max(|x|·RELTOL,ABSTOLX,NEPS·eps(x))
%       TY = @(y) max(|y|·RELTOL,ABSTOLY,NEPS·eps(y))

    persistent Def
    if isempty(Def)
        Def.RelTol = getSimOption('RelTol');
        Def.minabstol = getSimOption('minabstol');
        Def.minreltol = getSimOption('minreltol');
        Def.NEPS = getSimOption('NEPS');
    end
    if nargin > 1
        opt = getpairedoptions(varargin,{'RelTol','minabstol','minreltol','NEPS'},...
            {Def.RelTol,Def.minabstol,Def.minreltol,Def.NEPS},'restchk');
    else, opt = Def;
    end

    % Resolve/bound tolerances
    if nargin < 1 || isempty(t), t = opt.RelTol; end
    switch numel(t)
        case 1, t = [0,0,t];        % reltol ~= 0
        case 2, t(3) = 0;           % [maxtolx, maxtoly]
        case 3  % all set           % [maxtolx, maxtoly, reltol]
        otherwise, error('parsetolerance:ntol','Expecting 1-3 vector tolerance');
    end
    t(1:2) = max(t(1:2),opt.minabstol,'includenan');
    t(3) = max(t(3),opt.minreltol,'includenan');
    assert(all(t>=0),'parsetolerance:tolrange','Tolerances must be real-positive');
    
    if nargout > 1
        tx = @(x) min(max(max(abs(x)*t(3),opt.NEPS*eps(x)),t(1)),realmax);
        ty = @(y) min(max(max(abs(y)*t(3),opt.NEPS*eps(y)),t(2)),realmax);
    end
end