function varargout = crashdump()
% CRASHDUMP() - Save the workspace of calling function into an auto-named file, for later 
%   debugging. File names follow the pattern:
%
%       $HOME/CALLER_yyyymmddHHMMSS_XXXX.mat'
%
%   where $HOME is the user's home directory, CALLER is the name of the calling function
%   yyyymmddHHMMSS the formatted date, and XXXX a random alphanumeric key.
%
%   The result of dbstack() is also appended to the file, as a variable named DBSTACK. 
%   Finally a warning is issued: 'Crash-dump saved to ...'
%
% EXAMPLE:
%     if something_terrible
%         if debugging(), keyboard(); 
%         else
%             crashdump();
%             warning('Attempting to continue...');
%         end
%     end
%
% See also: DEBUGGING

    DBSTACK = dbstack(1);
    if isempty(DBSTACK), source = 'cmd'; 
    else
        source = matlab.lang.makeValidName(DBSTACK(1).name); 
    end
    date = datestr(now(),'yyyymmddHHMMSS');
    uid = dec2base(randi(36^4-1),36,4);
    path = getuserdir();
    file = fullfile(path,[strjoin({source,date,uid},'_') '.mat']);
    
    evalin('caller',['save(''' file ''');']);
    save(file,'DBSTACK','-append');
    warning('Crash-dump saved to %s',file);
    
    if nargout > 0, varargout{1} = file; end
end

function userDir = getuserdir
%GETUSERDIR   return the user home directory.
%   USERDIR = GETUSERDIR returns the user home directory using the registry
%   on windows systems and using Java on non windows systems as a string
%
% Source: Sven Probst (2020). Get user home directory 
%   (https://www.mathworks.com/matlabcentral/fileexchange/15885-get-user-home-directory), 
%   MATLAB Central File Exchange. Retrieved August 19, 2020. 

    if ispc
        userDir = winqueryreg('HKEY_CURRENT_USER',...
            ['Software\Microsoft\Windows\CurrentVersion\' ...
             'Explorer\Shell Folders'],'Personal');
    else
        userDir = char(java.lang.System.getProperty('user.home'));
    end
end