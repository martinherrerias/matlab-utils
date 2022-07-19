function varargout = completeoptions(options,varargin)
% OPTIONS = COMPLETEOPTIONS(OPT) - Complete the options structure OPT with default 
%   values from existing global SimOptions, and from DEFAULTOPTIONS. This call will NOT change 
%   the existing global SimOptions.
%
% [OPTIONS,FLAGS,MSG] = COMPLETEOPTIONS(OPT) - return comparison with current defaults.
%   FLAGS is boolean structure with the same fields as OPTIONS, where FLAGS.(FLD) = FALSE signals
%   non-default option FLD. MSG is a comparison string of these custom settings.
%
% COMPLETEOPTIONS(OPT) - Set global SimOptions and OptionFlags to the resulting OPTIONS,FLAGS.
%
% See also: COMPLETESTRUCT, GETSIMOPTION, SETSIMOPTION

    global SimOptions
    global OptionFlags
    if isempty(SimOptions), SimOptions = struct(); end
    
    if nargin < 1 || isempty(options), options = struct(); end

    options = completestruct(options,SimOptions,varargin{:});
    
    DEF = DefaultOptions();
        
    % ... Remove obsolete options
    [val,names] = nestedstruct2cell(options);
    obsolete = ~isnestedfield(DEF,names);
    if any(obsolete)
        warning('Removed obsolete %s',shortliststr(names(obsolete),'option'));
        options = cell2nestedstruct(val(~obsolete),names(~obsolete));
    end

    % ... Use defaults for any new options
    options = completestruct(options,DEF,varargin{:});
    
    [custom,~,txt] = comparestruct(options,DEF);
    
    names = nestedfieldnames(options);
    flags = cell2nestedstruct(repmat({true},numel(names),1),names);
    nondef = nestedfieldnames(custom);
    for j = 1:numel(nondef)
        flags = setnestedfield(flags,nondef{j},false);
    end
    
    if nargout == 0
        SimOptions = options; 
        OptionFlags = flags;
    else
        varargout = {options,flags,txt};
    end
end