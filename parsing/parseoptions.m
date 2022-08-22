function [par,rem,isdef] = parseoptions(args,flags,varargin)
% [PAR,REM,ISDEF] = PARSEOPTIONS(ARGS,FLAGS,NAMES,[DEFVAL,...]) - Equivalent to GETFLAGOPTIONS
%   followed by GETPAIREDOPTIONS. Collects all flags, and name-value pairs in structure PAR,
%   any remaining arguments in REM.
%
% NOTES:
%   (§) Flags take precedence over name-value pairs, i.e. for arguments {..,'-flag','flag',0,..},
%       the result will be PAR.flag = 1.
%   ($) Positional argument order (for use with 'dealrest' option of GETPAIREDOPTIONS) takes 
%       name-value pairs first, then flags.
%
% EXAMPLES:
%   parseoptions({'-foo','bam',3},{'-foo','-bar'},{'bam'})
%   parseoptions({3,1,0},{'-foo','-bar'},{'bam'},'dealrest') % see note ($)
% 
% See also: GETPAIREDOPTIONS, GETFLAGOPTIONS

    if isempty(flags)
        [par,rem,isdef] = getpairedoptions(args,varargin{:});
        return;
    end

    [opt,args] = getflagoptions(args,flags);
    fld = fieldnames(opt);
    if ~isempty(varargin) && isstruct(varargin{1})
        varargin{1} = completestruct(varargin{1},cell2struct(repmat({false},numel(fld),1),fld));
    else
        if ~isempty(varargin) && iscell(varargin{1})
            names = varargin{1};
            n = numel(names);
            names = unique([names(:);fieldnames(opt)],'stable');
        else
            names = fieldnames(opt);
        end
        varargin{1} = names;
        
        if numel(varargin) > 1 && iscell(varargin{2})
            defval = varargin{2};
            if isscalar(defval) && n > 1, defval = repmat(defval,n,1); end
            [defval{end+1:numel(names)}] = deal(false);
            varargin{2} = defval;
        end
    end
    [par,rem,isdef] = getpairedoptions(args,varargin{:});
    
    % merge parameter-pair and flag options
    for j = 1:numel(fld)
       if ~isfield(opt,fld{j}), continue; end
       if ~isfield(par,fld{j}) || isempty(par.(fld{j})) || isequal(par.(fld{j}),false)
           par.(fld{j}) = opt.(fld{j});
           isdef.(fld{j}) = ~opt.(fld{j});
       end
    end
end