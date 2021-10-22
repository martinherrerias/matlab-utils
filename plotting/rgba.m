function [cdata,alpha] = rgba(colorinfo)
% RGBA Converts some color format provided by the user in colorinfo into a n·3 matrix cdata and an
% n-vector alpha.
%
% colorinfo: can be a vertical vector of scalar indices for the current ax.ColorOrder, character
%   keys for basic colors {'y','m','c','r','g','b','w','k'} or n·3/n·4 matrixes for (R,G,B,[alpha]).
%   Special string 'none' is recognized as [0,0,0,0], the same as integer 0.
%   cell-arrays of the above will be resolved individually and vertically-concatenated.
%
% See also COLORMAP, GCA

    if iscell(colorinfo)
        for j = numel(colorinfo):-1:1
            [cdata{j},alpha{j}] = rgba(colorinfo{j});
        end
        cdata = vertcat(cdata{:});
        alpha = vertcat(alpha{:});
        return;
    end
        
    if ischar(colorinfo)
        if isscalar(colorinfo) 
            cdata = findcolor(colorinfo); % 'y','m','c','r','g'...
            alpha = 1;
        elseif strcmpi(colorinfo,'none'), cdata = [0 0 0]; alpha = 0; % transparent 
        else, cdata = NaN; alpha = NaN;
        end
    elseif isnumeric(colorinfo)
        switch size(colorinfo,2)
            case 1 
                cdata = zeros(size(colorinfo,1),3);
                alpha = zeros(size(colorinfo,1),1);
                col = get(gca,'ColorOrder'); % Ordered color, from ax.ColorOrder
                nc = size(col,1);
                nz = colorinfo~=0; % leave zeros as they are
                colorinfo = mod(colorinfo(nz)-1,nc)+1; % recycle colors
                cdata(nz,:) = col(colorinfo,:); 
                alpha(nz) = 1;
            case 3, cdata = colorinfo; alpha = 1; % [R,G,B]
            case 4, cdata = colorinfo(:,1:3); alpha = colorinfo(:,4); % [R,G,B,alpha]
            otherwise, cdata = NaN; alpha = NaN;
        end
    else
        cdata = NaN; alpha = NaN;
    end
end

function c = findcolor(a)
    CD = [1 1 0 1 0 0 1 0;
         1 0 1 0 1 0 1 0;
         0 1 1 0 0 1 1 0]';
    c = CD(a=='ymcrgbwk',:); 
end
