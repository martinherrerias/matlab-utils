function shelled(code)
% SHELLED(CODE) - eval('base',CODE) inside a try-catch block, then exit(f) MATLAB with the 
%   appropiate flag. Meant to simplify command-line MATLAB calls, providing the shell/console
%   an exit flag, and avoiding that MATLAB stays open after crash.
%
% NOTE: make darn sure SHELLED is in your current PATH! otherwise MATLAB will still open, throw an 
%   error complaining that SHELLED wasn't found, and stay open waiting for input.
%
% EXAMPLE:
%     #!/bin/bash
%     matlab -nodisplay -r "shelled('dostuff()')"
%     [ $? == "0" ] && echo "Hey, it actually worked!" || echo "Nope, it crashed"

    if nargin < 1 || isempty(code), code = "disp('shelled: Nothing to do')"; end
    
    try
        evalin('base',code);
        exit(0);
    catch ERR
        disp(getReport(ERR));
        exit(1);
    end
end