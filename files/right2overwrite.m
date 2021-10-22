function [doit,msginfo,reallyasked] = right2overwrite(filename,forceoverwrite)
% Checks the current directory to see if filename exists, if so it asks the user whether or not to
% replace it.
%   filename: can be a string, or a cell-array of strings.
%   if forceoverwrite == true, doit is set to true without even checking.
%   doit: boolean result, true means user agreed to overwrite file (or forceoverwrite is true)
%   reallyasked: boolean flag, used to avoid asking repeatedly for several files.
%   msginfo: user feedbak string.
 
    if nargin < 2, forceoverwrite = false; end
    
    doit = true;
    reallyasked = false;
    
    if forceoverwrite
        msginfo = 'Force-overwrite is on';
    else
        if iscell(filename) % for several files, ask the user only once, if at all.
            alreadyasked = false;
            for j = 1:numel(filename)
                if ~alreadyasked 
                    [doit,msginfo,alreadyasked] = right2overwrite(filename{j},false);
                    if ~doit, break; end 
                end
            end
        else
            if ~isempty(dir(filename))
                reallyasked = true;
                switch optquestdlg('Do you want to overwrite the existing file(s)?',...
                                                        'File(s) already exists','Yes','No','No')
                    case 'Yes'
                        msginfo = 'Overwriting existing file';
                        doit = true;
                    case 'No'
                        msginfo = 'Existing file, you asked me not to touch it';
                        doit = false;
                    otherwise
                        error('Failed to get answer');
                end
            else
                msginfo = 'Creating new file';
            end
        end
    end
end
