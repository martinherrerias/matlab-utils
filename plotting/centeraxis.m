function centeraxis(ax,X0)
% centeraxis(ax) - set all ruler CrossoverValue properties to 0

    if nargin < 1, ax = gca(); end
    if nargin < 2, X0 = [0,0,0]; end
    
    ax.XRuler.FirstCrossoverValue  = X0(2);
    ax.XRuler.SecondCrossoverValue = X0(3);
    ax.YRuler.FirstCrossoverValue  = X0(1);
    ax.YRuler.SecondCrossoverValue = X0(3);
    ax.ZRuler.FirstCrossoverValue  = X0(1);
    ax.ZRuler.SecondCrossoverValue = X0(2);
end