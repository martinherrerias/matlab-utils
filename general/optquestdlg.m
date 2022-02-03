function ButtonName = optquestdlg(varargin)
% OPTQUESTDLG opens a question dialog box ONLY if running MATLAB in interactive mode, logging both
%   the question and the answer into the console/diary; otherwise it directly returns the Default 
%   option and issues a warning. 
% See also: QUESTDLG

    LBL = 'Questdlg:';

    if nargin < 1, error(message('MATLAB:questdlg:TooFewArguments')); end
    
    if nargin < 2 || isempty(varargin{2}), question = varargin{1};
    else, question = ['(' varargin{2} ') ' varargin{1}];
    end
    
    [interactive,hasdisplay] = runningfromUI();
    if hasdisplay
        ButtonName = questdlg(varargin{:});
        fprintf('%s %s - %s\n',LBL,question,ButtonName);
    elseif interactive
        ButtonName = CLIquestion([LBL ' ' question],varargin{:});       
    else
        if nargin < 3, ButtonName = 'Yes';
        else, ButtonName = varargin{end}; end
        if isstruct(ButtonName), ButtonName = ButtonName.Default; end
        
        warning('optquestdlg:default','Using default answer: %s - ''%s''',question,ButtonName);
    end
end

function answer = CLIquestion(question,varargin)

    switch nargin
        case {1,2}, options = {'Yes','No','Cancel'}; default = 'Yes';
        case 3, options = {'Yes','No','Cancel'}; default = varargin{3};
        otherwise, options = varargin(3:end-1); default = varargin{end}; 
    end
    if ~ismember(default,options)
        warning('optquestdlg:StringMismatch',...
            'Default character vector does not match any button character vector name.');
    end
    
    prompt = [question ' ' strjoin(options,'/') ' (' default '): '];
    match = 0;
    
    while nnz(match) ~= 1
        answer = input(prompt,'s');
        if isempty(answer)
            answer = default; return;
        else
            match = contains(options,answer);
            if ~any(match), match = contains(options,answer,'IgnoreCase',true); end
            if nnz(match) == 1
                answer = options{match}; return;
            else
                fprintf('Invalid/ambiguous answer\n'); 
            end
        end
    end
end

