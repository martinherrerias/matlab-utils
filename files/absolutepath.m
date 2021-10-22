function abspaths = absolutepath(relpaths,rootpath)
% RELPATHS = ABSOLUTEPATH(FILEPATHS) - Return absolute paths for file(s)/directory(ies) FILEPATHS.
% See also: RELATIVEPATH

    if nargin > 1 && ~isempty(rootpath)
        basedir = pwd();
        lastwill = onCleanup(@() cd(basedir));
        cd(rootpath);
    end

    cellinput = iscell(relpaths);
    if ~cellinput, relpaths = {relpaths}; end
    
    abspaths = cell(size(relpaths));
    for j = 1:numel(relpaths)
        d = dir(relpaths{j});
        assert(~isempty(d),'file/directory not found');
        if d(1).isdir
            abspaths{j} = d(1).folder;
        else
            abspaths{j} = fullfile(d(1).folder,d(1).name);
        end
    end
    
    if ~cellinput, abspaths = abspaths{:};
end
