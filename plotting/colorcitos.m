function map = colorcitos(m)
% MAP = COLORCITOS(M) - a warmer, dirtier PARULA.

    C = [0.50 0.20 0.65;
         0.35 0.40 0.80;
         0.25 0.60 0.60;
         0.50 0.65 0.30;
         0.75 0.60 0.10;
         0.80 0.60 0.20;
         1.00 0.80 0.30];

    x = [0 1/4 1/2 2/3 3/4 5/6 1];

    map = interp1(x,C,linspace(0,1,m),'pchip');
end