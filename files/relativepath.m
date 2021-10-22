function filepaths = relativepath(filepaths,rootpath,maxdepth)
% RELPATHS = RELATIVEPATH(ABSPATHS,ROOTPATH,MAXDEPTH)
% Try to write each element of FILEPATHS (cell-array of strings), as a relative path from ROOTPATH.
% See also: ABSOLUTEPATH

    if nargin < 2, rootpath = pwd(); end
    if nargin < 3, maxdepth = Inf; end
  
    roottree = strsplit(rootpath,{'\','/'});
    if isempty(roottree{end}), roottree = roottree(1:end-1); end
    if isempty(roottree), return; end
    nr = numel(roottree);
    
    cellinput = iscell(filepaths);
    if ~cellinput, filepaths = {filepaths}; end
    
    for j = 1:numel(filepaths)
        filetree = strsplit(filepaths{j},{'\','/'});
        endslash = ~isempty(filetree) && isempty(filetree{end});
        if endslash, filetree(end) = []; end 
        if isempty(filetree), filepaths{j} = ''; continue; end
        nf = numel(filetree);
        
        % Get index of last matching folder name in root- and file-tree
        k = 0; 
        while strcmpi(filetree{k+1},roottree{k+1})
            k = k + 1;
            if k == nr || k == nf, break; end
        end
        if k == 0, continue; end    % nothing in common
     
        if nr - k > maxdepth
        % keep using absolute path
            filepaths{j} =  fullfile(filetree{:});
        else
        % go back nr-k folders, then forward nf-k-1 into the divergent file path
            if nr > k, backpath = repmat({'..'},1,nr-k); else, backpath = {'.'}; end
            filepaths{j} = fullfile(backpath{:},filetree{k+1:nf});
        end
        
        if endslash, filepaths{j} = [filepaths{j},filesep()]; end % put back ending /
    end
    
    if ~cellinput, filepaths = filepaths{:};
end
