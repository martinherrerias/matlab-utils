function D = DefaultsFrom(varargin)
% D = DEFAULTSFROM(PATH) - call a specific function of DEFAULTFUNCTIONS, found in PATH. 
%     Meant to allow DEFAULTOPTIONS to be overriden/extended.
% D = DEFAULTSFROM(BASE,REL,..) - Use a series of relative paths.
    
    here = pwd(); 
    foo = onCleanup(@() cd(here));
    cellfun(@cd,varargin);
    assert(isfile('DefaultOptions.m'),'Failed to find DefaultOptions.m in PATH');
    D = DefaultOptions();
end