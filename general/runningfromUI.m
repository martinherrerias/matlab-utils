function [interactive,graphic,version] = runningfromUI()
% [ISINTERACTIVE, ISGRAPHIC, VERSION] = RUNNINGFROMUI is meant to provide support for a program
%   to have several 'modes': Graphical User Interface (GUI), Command Line Interface (CLI), or
%   unattended 'batch' mode.
%
% A specific version can be forced by setting the 'version' field of global SimOptions, e.g.
% with setSimOption('version','GUI'). If no version is defined, RUNNINGFROMUI scans for jDesktop 
% clients and feature('ShowFigureWindows') to choose between GUI/batch mode.
%
% See also: SPLITGUI, SPLITSCRIPT, GETSIMOPTION

    VERSIONS = {'GUI', {'GUI','UI'};
                'CLI', {'CLI'};
                'batch', {'batch','parallel'}}; 
                % 'API'?
    ISINTERACTIVE = [1,1,0];
    ISGRAPHIC = [1,0,0];

    % Get version from SimOptions
    if any(strcmp(who('global'),'SimOptions'))
        version = getSimOption('version');
    else
        version = '';
    end

    id = 0;
    if ~isempty(version)
        [version,id] = parselist(version,VERSIONS,'-soft');
        if id == 0
            warning('Unknown version ID: %s',version); 
        end
    end
    if id == 0
        version = 'Unknown';
        graphic = usejava('desktop') && feature('ShowFigureWindows');
        interactive = graphic;
    else
        interactive = ISINTERACTIVE(id);
        graphic = ISGRAPHIC(id);
    end
end
