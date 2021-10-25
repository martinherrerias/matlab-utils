function Defaults = DefaultOptions()
% Returns a structure with default simulation options
% TODO: should be replaced by a language-neutral configuration file (XML, YAML, ...)

    Defaults.prjname = '';
    Defaults.version = 'GUI'; % (see RUNNINGFROMUI)
    Defaults.verbose = true;    
    
    Defaults.outliers.P = 1e-6;          % Probability to consider a point an outlier
    Defaults.outliers.warning = 0.25;    % fraction of outlier data to cause warning
    Defaults.outliers.error = 0.95;      % fraction of outlier data to cause error

    % Precision and Tolerances (see PARSETOLERANCE)
    Defaults.RelTol = 1e-3;             % Relative tolerance
    Defaults.MaxIter = 100;             % Used for iterative solutions
    Defaults.minabstol = 1e-12;         % Used to resolve machine precision errors (§)
    Defaults.minreltol = 1e-6;          % Used to resolve machine precision errors (§)
    Defaults.NEPS = 8;                  % Use NEPS·eps(x) as minimum relative tolerance
    
    % (§) NOTE: most operations are performed on double-precision floats, which in principle have
    % a minimum absolute tolerance of eps(0) = 2^(-1074), and a minimum relative tolerance of 
    % eps(1) = 2^(-52). In practice numerical errors stack up, and setting relative and absolute
    % tolerances of 1e-8 and 1e-12 is already ambitious.
    
    % Parallel execution
    Defaults.runparallel = false;       % see RUNPARALLEL
    Defaults.maxthreads = 24;           % limit the max. number of instances (see getnworkers)
    Defaults.chunksize = 1000;
                                        
    % Plotting
    Defaults.plotting = false;           % Plot every shading/solver step (editable on runtime)
    Defaults.exportplots = false;        % Save a copy of each plot as PNG image
    Defaults.exportplotres = 300;        % resolution
    
    % File Importing
    Defaults.DelimiterPriority = [char(9),';, ']; % delimiter priority: (tab) ; , (space)
end
