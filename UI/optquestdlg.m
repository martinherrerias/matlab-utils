function ButtonName = optquestdlg(varargin)
% OPTQUESTDLG opens a question dialog box ONLY if running MATLAB in interactive mode, logging both
%   the question and the answer into the console/diary; otherwise it directly returns the Default 
%   option and issues a warning. 
% See also: QUESTDLG

    if nargin < 1, error(message('MATLAB:questdlg:TooFewArguments')); end
    
    if runningfromUI()
        ButtonName = questdlg(varargin{:});
        fprintf('Questdlg: %s - %s\n',varargin{1},ButtonName);
    else
        if nargin < 3, ButtonName = 'Yes';
        else, ButtonName = varargin{end}; end
        if isstruct(ButtonName), ButtonName = ButtonName.Default; end
        
        warning('optquestdlg:default','Using default answer: %s - ''%s''',varargin{1},ButtonName);
    end
end
