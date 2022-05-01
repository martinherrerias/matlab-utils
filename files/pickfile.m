function filenames = pickfile(pattern,varargin)
% PICKFILE(PATTERN,...) - Search the current directory for files whose name matches the provided 
%   PATTERN (e.g. '*.ext'), and return A string or cell-array of strings with the complete path 
%   and file name(s). PATTERN can be anything that works with DIR(), e.g. '../**/*.txt'
%
% F = PICKFILE(PATTERN) expects a single file, if none/too-many files match PATTERN, the user is 
%   asked to pick manually (or an error is thrown when RUNNINGFROMUI sets a non-interactive mode).
%   The result is a string with the name (with relative path) of the file.
%
% F = PICKFILE(PATTERN,PROMPT) ... use a custom prompt when/if the user has to pick.
%
% C = PICKFILE(PATTERN,N) ask-for/expect N files, use Inf for 'any number' of files.
%   The result is a cell-array of strings, except when N = 1. Note that with N = Inf, a cell-array
%   will be returned even if only one file is found.
% C = PICKFILE(PATTERN,[N1,N2]) expect a range [N1,N2] number of files. The syntax above is 
%   equivalent to [N,N], and [0,Inf], respectively.
%
% C = PICKFILE(PATTERN,N,PROMPT) ask for N files with custom prompt.
%
% X = PICKFILE(..,'ui',1) force user-interactive pick, even if number of candidates matches n
% X = PICKFILE(..,'ui',0) default, user-interactive pick only if numel(candidates ) ~= n
% X = PICKFILE(..,'ui',-1) blocks user-interactive pick, throws error if numel(candidates) ~= n
%   NOTE: if the code is running in non-interactive mode (see RUNNINGFROMUI), PICKFILE will work 
%   the same as with 'ui' = -1, throwing an error when numel(candidates) ~= n.
%
% X = PICKFILE(..,'fullpath',true) attach full-path to each returned file's name. If set to false
%   (default), local (or relative up to ../) file-names will be returned when possible.
% 
% Examples: 
%   F = PICKFILE('*.txt') - look for a text file, ask the user to pick if there are none/too many
%   C = PICKFILE('*.txt',Inf,'ui',1) - ask the user to select any number of text files
%   C = PICKFILE('*.txt',Inf,'fullpath',1) - return full paths for any and all text files in dir()
%
%   TODO: use predefined global seach PATH from SimOptions? would allow projects with files spread
%   over multiple directories.
%
% See also: DIR, UIGETFILE, RUNNINGFROMUI, WHICH

    narginchk(0,7);
    [opt,remargs] = getflagoptions(varargin,{'-fullpath','-ui','-soft'});
    opt.ui = double(opt.ui);
    [opt,remargs] = getpairedoptions(remargs,opt);
    prompt = 'Select a file';
    n = 1;
    switch numel(remargs) + 1*(nargin > 0)
        case 0, pattern = '*';
        case 1 % other defaults set above
        case 2, if ischar(remargs{1}), prompt = remargs{1}; else, n = remargs{1}; end
        case 3, n = remargs{1}; prompt = remargs{2};
        otherwise, error('pickfile:nargs','Expecting up to 3 arguments + 2 property-value pairs');
    end
    if isscalar(n)
        if isinf(n), n = [0,Inf]; else, n = [n,n]; end
    end
    if n(2) <= 0, filenames = {}; return; end
 
    % 0. Look for candidates
    if iscell(pattern)
        if ~opt.ui
        	for j = numel(pattern):-1:1, candidates{j} = dir(pattern{j}); end
            candidates = cat(1,candidates{:});
        else
            candidates = struct('name',{},'folder',{},'isdir',{});
        end
    else
        assert(isempty(regexp(pattern,'[?"<>|]|[\\/]{2}','once')),'Invalid file name or PATTERN');
        candidates = dir(pattern);
        candidates = candidates(~[candidates.isdir]); % remove directories
    end
    
    if ~(opt.ui > 0) && ( numel(candidates) >= n(1) && numel(candidates) <= n(2))
    % 1. not-user-interactive pick
        if isempty(candidates), filenames = {}; return; end % n = 0, or n = Inf, make no fuzz
        
        filenames = {candidates.name};
        filepaths = {candidates.folder}; 
    else        
        % 2. No match, but not user-interactive: throw error
        if ~runningfromUI() || opt.ui < 0
            if numel(candidates) > n(2), msg = 'Too many candidates'; 
            elseif numel(candidates) < n(1), msg = 'File(s) not found';
            else, msg = 'Forced UI-pick'; 
            end
            error('pickfile:uipick','%s, Cannot use UI-pick in not-interactive mode!', msg);
        end
        
        % 3. User-interactive pick
        if n(2) == 1, [filenames, filepaths] = uigetfile(pattern,prompt);
        else, [filenames, filepaths] = uigetfile(pattern,prompt,'MultiSelect','on');
        end
    end
    
    if ~iscell(filenames)
    % Make sure filename is a cell-array of strings (just for now)
       if filenames == 0, filenames = {};  else, filenames = {filenames}; end
    end 
    if ~opt.soft && (numel(filenames) < n(1) || numel(filenames) > n(2))
    % Throw error if the result is not what we expected
        if n(2) > n(1), msg = sprintf('%d-%d files.',n(1),n(2)); 
        else, msg = nthings(n(1),'file'); 
        end
        error('pickfile:n','pickfile: you were supposed to pick %s',msg);
    end
    if isempty(filenames), return; end
    
    % Make sure filepaths is a cell-array of (absolute/relative) paths, based on opt.fullpath
    if ~opt.fullpath, filepaths = relativepath(filepaths); end
    if ischar(filepaths), filepaths = {filepaths}; end
    
    if numel(filenames) > numel(filepaths), filepaths(1:numel(filenames)) = filepaths; end

    % Put filenames and filepaths together
    for j = 1:numel(filenames)
        filenames{j} = fullfile(filepaths{j}, filenames{j}); 
    end
    if isscalar(filenames) && n(2) == 1, filenames = filenames{1}; end
end
