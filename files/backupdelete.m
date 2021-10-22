function varargout = backupdelete(filename,varargin)
% BACKUP = BACKUPDELETE(FILE) - delete FILE, making sure a BACKUP copy is kept somewhere prior to
%   deletion. Rules for BACKUP naming are the following:
%
%   a) PATH/FILE.EXT --> TEMPDIR()/MATLAB_Files_../FILE.EXT     if  RECYCLE() == 'on'
%   b) PATH/FILE.EXT --> PATH/FILE_EXT.bak                      if  RECYCLE() == 'off'
%   c) .../FILE.EXT --> .../FILE_EXT_CPUTIME.bak  if either of the above matches an existing file
%
% See also: DELETE, RECYCLE, TEMPDIR

    [opt,varargin] = getflagoptions(varargin,{'-warning'});
    opt.suffixes{1} = '.bak';
    opt.suffixes{2} = [num2str(cputime*1000) '.bak'];
    opt = getpairedoptions(varargin,opt);

    if isempty(dir(filename))
       warning('File (%s) not found',filename);
       if nargout > 0, varargout{1} = ''; end
       return
    end

    recycling = strcmp(recycle,'on') && isfolder(tempdir());
    if recycling
        d = dir(fullfile(tempdir,'MATLAB_Files_*'));
        [~,idx] = max([d.datenum]);
        backupname = fullfile(tempdir(),d(idx).name,filename);
    else
        [path,backupname,ext] = fileparts(filename);
        backupname = fullfile(path,[backupname '_' ext(2:end) opt.suffixes{1}]);
    end
    if ~isempty(dir(backupname))
        [path,backupname,ext] = fileparts(backupname);
        if ~recycling,ext = ''; end
        backupname = fullfile(path,[backupname '_' ext(2:end) opt.suffixes{2}]);
        if recycling
            lastwill = onCleanup(@() recycle('on'));
            recycle('off');
            recycling = false;
        end
    end
    
    if ~recycling
        [success,msg] = copyfile(filename,backupname);
        assert(success,msg);
    end
    
    delete(filename);
    if nargout > 0, varargout{1} = backupname; end
    
    if opt.warning
        warning('Conflicting file (%s) moved to: %s',relativepath(filename),...
                                                     relativepath(backupname));
    end
end