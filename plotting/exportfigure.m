function exportfigure(sufix,figh)
% EXPORTFIGURE(sufix,figh) - save a copy of figure FIGH as PATH/fig/PRJ_SUFIX.png

    narginchk(1,2);
    if nargin < 2, figh = gcf(); end

    [path,prj] = fileparts(getSimOption('prjname'));
    path = fullfile(path,'fig');
    if isempty(dir(path)), mkdir(path); end
    
    figname = fullfile(path,[prj '_' sufix '.png']);
    resolution = getSimOption('exportplotres');
    if ~isempty(dir(figname)), delete(figname); end
    print(figh,'-dpng',sprintf('-r%d',resolution),figname);
end
