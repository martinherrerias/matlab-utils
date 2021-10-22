function [interactive,version,id] = runningfromUI()
% [ISINTERACTIVE, VERSION, ID] = RUNNINGFROMUI Analyzes the current directory/workspace in order 
% to guess what version of the code to use, mainly to choose between graphical/log-based display
% when running in interactive vs unattended modes, but also e.g. to provide temporal backwards 
% compatibility when deploying new versions.
%
% A specific result can be 'forced' by creating a RUNFROM = VERSION; variable in base workspace.
% With secondary priority, a global SimOptions structure is scanned for SimOptions.version. 
% Current versions include:
%
%  ID: VERSION [ISINTERACTIVE]
%  -1: 'Debug' [?]
%   0: 'Unknown' [?]
%   1: 'PythonGUI' [true] - there is a *.json file in the current directory
%   2: 'splitGUI' [true] - basic interactive/debugging mode
%   3: 'splitscript' [false] - log-based, unattended mode (-nodisplay Matlab)
%   4: 'parallel' [false] - meant to run in parallel-clusters   
%
% For versions in which ISINTERACTIVE is not defined, RUNNINGFROMUI scans for jDesktop clients.
% In any case, ISINTERACTIVE is overriden if ~feature('ShowFigureWindows')
%
% See also: SPLITGUI, SPLITSCRIPT, GETSIMOPTION

    VERSIONS = {'Debug','Unknown','PythonGUI','splitGUI','splitscript','parallel'};
    VERSION_IDS = [-1,0,1,2,3,4];
    ISINTERACTIVE = [1,-1,1,1,0,0]; % -1 for who-knows

    % 1. DEBUGGING back-door
    if ~evalin('base','exist(''RUNFROM'')'), version = '';
    else, version = evalin('base','RUNFROM');
    end

    % 2. Get version from SimOptions
    if isempty(version) && any(strcmp(who('global'),'SimOptions')), version = getSimOption('version'); end

    % 3. Missing SimOptions structure and a *.json file somewhere point to old pythonGUI
    if isempty(version) && ~isempty(dir('*.json')), version = 'pythonGUI'; end

    % Match version to known VERSIONS
    if ischar(version)
       id = find(strcmpi(version,VERSIONS),1);
       if ~isempty(id), id = VERSION_IDS(id); else, id = 0; end
    else
       id = version;
    end
    if id == 0 || ~any(id==VERSION_IDS)
        if ~isempty(version), warning('Unknown version ID: %s',version); end
        id = 0; 
    end
    version = VERSIONS{id == VERSION_IDS};

    % Define if version is interactive
    interactive = ISINTERACTIVE(id == VERSION_IDS);
    if interactive < 0, interactive = usejava('desktop'); end
    
    interactive = interactive && feature('ShowFigureWindows');
end
