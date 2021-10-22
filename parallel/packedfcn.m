function varargout = packedfcn(fcnfile,varargin)
% PACKEDFCN(FILE) - run function fcn(args{:}) where both the function F.fcn and its arguments 
%   F.args are included in a -mat file FILE. Results will be appended to the same FILE as F.res.
% 	[NOTE F.x is used here as a shortcut for FILE variables, i.e. F = load(FILE)]
%
% FLAG = PACKEDFCN(FILE,'-exit') - will exit(FLAG) the current MATLAB session, regardless of the
%   result. It's designed to be used on shell/console calls of the sort:
%
%     #!/bin/bash
%     matlab -nodisplay -r "packedfcn('foo.mat','-exit')"
%     [ $? == "0" ] && echo "Hey, it actually worked!" || echo "Nope, it crashed"
%
% PACKEDFCN(FILE,'-overwrite') - re-run res = fcn(args{:}) and save the result to F.res even if
%   res already exists in the provided FILE.
%
% PACKEDFCN(FILE,'-backup') - include autosave-crash-resume functionality. For this to work, fcn
%    must recognize the name-value pair ..,'backup',X. Specifically the calls:
%
%       fcn(args{:},'backup',FILE~) - where FILE~ is the name of a backup file to use
%       fcn('backup',S) - where S = load(FILE~,'-mat') contains the contents of an interrupted run
%
%   On a fresh run (first syntax), fcn is informed that is should start dumping any provisional 
%   results into a MAT file named FILE~. On a recovery run (i.e. if PACKEDFCN finds an existing
%   FILE~), fcn is expected to recover and resume its task from the contents of S.
%
% NOTE: the '-backup' flag will be set automatically if the function arguments F.args already
%   contain a name,value pair 'backup',X for ischar(X). To override this behavior (e.g. if fcn
%   uses 'backup' for something else) use PACKEDFCN(..,'backup',false).
%
% INPUT:
%   FILE - file name (usually *.mat), containing the following variables:
%
%     args - cell-array of arguments to be passed to function fcn.
%
%     fcn - function to be executed as R = fcn(args{:}), and optionally as R = fcn(S) where S is a 
%       structure with contents of a backup file: S = load(FILE~). fcn can be a string recongizable 
%       by STR2FUNC or directly a function-handle (internal-references seem to work). A single
%       output is required. Storage of backup files 
%
%     opt - If included, COMPLETEOPTIONS(opt) will be set as global SimOptions prior to execution.
%
%     res - Meant as the output variable for the results of fcn(args{:}), if found on the input
%       (and unless -overwrite flag is used), calculation will be aborted. This functionality is
%       avoids unnecessary calculations when resuming interrupted/crashed RUNPARALLEL calls.
%
% OUTPUT: results of the function evaluation are saved to FILE.res. 
%   FLAG = PACKEDFCN(FILE,..) returns a simple binary flag for success (0) failure (1).
%
% EXAMPLE:
%     args = {6,7};
%     fcn = @times;
%     save('foo.mat','args','fcn');
%     !matlab -nodisplay -r "packedfcn('foo.mat','-exit')"  % run packedfcn from shell
%     F = load('foo.mat');
%     fprintf('The answer to life, the universe, and everything else: %d\n',F.res)
%
% See also: RUNPARALLEL, PARALLELCONFIG

    global SimOptions
    % VARNAMES = args, fcn, opt, res - just check against RUNPARALLEL
    ERR_FLAG = 1;
    OK_FLAG = 0;
        
    if debugging() && nargin == 0, test(); return; end
    
    try
        narginchk(1,4);

        feature('DefaultCharacterSet','UTF-8');
        
        warning_resetter = naptime(); %#ok<NASGU>
        warning on all; warning on backtrace; warning on verbose
        warning('error','MATLAB:load:variableNotFound'); %#ok<CTPCT>

        [opt,varargin] = getflagoptions(varargin,{'-exit','-backup','-overwrite'});
        if ~opt.backup, opt.backup = []; end
        opt = getpairedoptions(varargin,opt,'restchk');

        [path,fname,~] = fileparts(fcnfile);
        diary(fullfile(path,[fname '.log']));
        
        if ~isempty(SimOptions)
        % Backup current SimOptions, restore once packedfcn is done
            option_resetter = restoreSimOptions(); %#ok<NASGU>
        end
        
        try
        % Load/complete SimOptions first (for use in loadobj methods)
        % NOTE: using vars = whos('-file',fcnfile); if ismember('opt',vars)...
        %   loads all variables (calling loaaobj methods) in the process, defeating the purpose
        %   the same goes for using F = matfile(fcnfile); ...
            
            F = load(fcnfile,'-mat','opt');
            SimOptions = [];
            SimOptions = completeoptions(F.opt);
        catch
        % Quietly delete SimOptions. Let the first call to getSimOption (possibly by some loadobj)
        % restore it to defaults with a warning, if required.
            SimOptions = [];
        end
     
        % Load everything else
        F = load(fcnfile,'-mat');
        assert(isfield(F,'args') && iscell(F.args),'Missing/bad "args" variable');
        assert(isfield(F,'fcn'),'Missing "fcn" variable');
        
        if isfield(F,'res') && ~opt.overwrite
        % Don't re-run function if result is already available
            warning('Skipping completed function-file: %s',fcnfile);
            if nargout > 0, varargout{1} = OK_FLAG; end
            if opt.exit, exit(OK_FLAG); end
            return;
        end

        if ischar(F.fcn) || isstring(F.fcn), F.fcn = str2func(F.fcn); end
        assert(isa(F.fcn,'function_handle'),'Failed to get function handle from "fcn"');

        if isempty(opt.backup) || opt.backup
            try
                [x,y] = getpairedoptions(F.args,{'backup'},{''});
                backupfile = x.backup;
                assert(ischar(backupfile));
                F.args = y;
                clear x y
            catch
                backupfile = '';
            end
            
            if isempty(opt.backup), opt.backup = ~isempty(backupfile);
            elseif isempty(backupfile)
                DEF = parallelconfig();
                backupfile = [fcnfile DEF.backup];
            end
        end

        if opt.backup
            F.args = [F.args(:)',{'backup',backupfile}];
        end
            
        if ~opt.backup || ~isfile(backupfile)
        % Fresh function call to fcn(args{:},['backup',file])
            fprintf('Running %s on: %s\n',func2str(F.fcn),fcnfile);
            F.res = F.fcn(F.args{:});
        else
        % Attempt recovery from backup

            fprintf('Attempting to resume interrupted calculation from:\n\t%s\n',backupfile);
            try
            % Attempt to read SUBFILE~ (backup) and resume from contents
                B = load(backupfile,'-mat');  
                F.res = F.fcn('backup',B);
            catch ERR
            % If this fails...
                fprintf('... resume failed: %s\n', getReport(ERR)); 
                oldbackup = backupdelete(backupfile);
                fprintf('Backup file moved to %s\n',oldbackup);
                F.res = F.fcn(F.args{:},'backup',backupfile);
            end
        end

        warning_resetter = naptime(); %#ok<NASGU>  
        try
            warning('error','MATLAB:onCleanup:DoNotSave'); %#ok<CTPCT>
            save(fcnfile,'-struct','F');
        catch ERR
            dbstack()
            getReport(ERR);
            rethrow(ERR);
        end
        flag = OK_FLAG; % all went right
        
        if opt.backup && isfile(backupfile), delete(backupfile); end
    catch ERR
        disp(getReport(ERR));
        flag = ERR_FLAG;
    end
    
    if nargout > 0, varargout{1} = flag; end
    if opt.exit, exit(flag); end
end

function test()
    NAME = 'foo.mat';
    DEF = parallelconfig();
    BACK = [NAME DEF.backup];

    if isfile(NAME), delete(NAME); end
    if isfile(BACK), delete(BACK); end

    fprintf('\n\nPACKEDFCN test\n');
    args = {6,7};
    fcn = @times;
    save(NAME,'args','fcn');
    system(['matlab -nodisplay -r "packedfcn(''' NAME ''',''-exit'')"']);  % run from shell
    F = load(NAME);
    assert(isequal(F.res,42));
    fprintf('Shell call works\n');
    
    args = {};
    opt = completeoptions();
    opt.prjname = 'packedfcn_test';
    opt.RelTol = 3.1416;
    fcn = @testoptions;
    save(NAME,'args','fcn','opt');
    packedfcn(NAME);
    F = load(NAME);
    assert(isequal(F.res,true));
    fprintf('Loading of SimOptions works\n');
    
    fcn = @testbackup;
    save(NAME,'args','fcn');
    packedfcn(NAME,'-backup'); % call to dummy('backup',BACK)
    packedfcn(NAME,'-backup'); % call to dummy('backup',load(BACK))
    F = load(NAME);
    assert(isequal(F.res,'Hmm...'))
    fprintf('Backup works\n');
    
    args = {'backup',42};
    save(NAME,'args','fcn');
    packedfcn(NAME);
    F = load(NAME);
    assert(isequal(F.res,42));
    
    args = {'backup','42!'};
    save(NAME,'args','fcn');
    packedfcn(NAME,'backup',false);
    F = load(NAME);
    assert(isequal(F.res,'42!'));
    fprintf('Confusing calls work\n');

    delete(NAME);
    fprintf('\nPACKEDFCN tests passed with flying colors\n\n');

    function r = testoptions(varargin)
       r = isequal(getSimOption('prjname'),'packedfcn_test') && ...
          isequal(getSimOption('RelTol'),3.1416);
    end

    function r = testbackup(varargin)
        
       S = getpairedoptions(varargin,{'backup'},{[]},'restchk');
       if ~isempty(S.backup)
          if ischar(S.backup) && S.backup(end) == '~'
              very_important_provisional_results = 'Hmm...';
              save(S.backup,'very_important_provisional_results');
              error('Oh no, a terrible crash! (don''t worry, it is part of the test)');
          elseif isstruct(S.backup)
              r = S.backup.very_important_provisional_results;
              return;
          else
              r = S.backup;
          end
       end
    end
end