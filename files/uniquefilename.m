function filename = uniquefilename(filename)
% UFILE = UNIQUEFILENAME(FILE) - check if FILE already exists, and return a unique name with the 
%   same path and extension.
%
%  EXAMPLE:
%   !touch ./foo.bar
%   uniquefilename('./foo.bar') % returns './foo_1.bar'
%
% See also: BACKUPDELETE, RIGHT2OVERWRITE

    if ~isfile(filename), return; end
    
    [path,filename,ext] = fileparts(filename);

    ots = pickfile(fullfile(path,['*' ext]),Inf);
    [~,ots] = cellfun(@fileparts,ots,'unif',0);
    
    filename = matlab.lang.makeUniqueStrings(filename,ots);
    filename = fullfile(path,[filename ext]);
end