function [M,Np] = getnworkers(M0,varargin)
% [M,Np] = GETNWORKERS(M0) - determine a maximum recommended number M of parallel workers/threads
%       given a preferred M0, available memmory, number of processors, and SimOptions.
%
% [M,Np] = GETNWORKERS(..,'workersize',S0,'buffer',B) - the max memmory-bound number of workers is
%   set as floor[ free-memmory·(1-B)/S0 ]. Defaults are 120e6 (Bytes) for S0, and 0.2 for B.
%
% [M,Np] = GETNWORKERS(..,'vars',S,'factor',F) - where S is the output of WHOS('a','b',..) adds 
%   sum([S.bytes])·F to the memmory requirements S0 of each worker.
%
% [M,Np] = GETNWORKERS(..,'create') - create/update parallel pool, if Parallel Toolbox is available
    
    global SimOptions;

    if nargin < 1 || isempty(M0), M0 = NaN; end
    assert(isscalar(M0) && isnumeric(M0),'Bad M');
    
    [opt,varargin] = getflagoptions(varargin,{'create'});
    opt.vars = [];
    opt.factor = 3.0;
    opt.buffer = 0.2;
    opt.workersize = 1.4e9; 
    opt.max = SimOptions.maxthreads;
    opt = getpairedoptions(varargin,opt,'restchk');

    % Get number of processors from the system
    if ispc, [status,Np] = system('%NUMBER_OF_PROCESSORS%');  % tested on wine cmd only
    elseif ismac, [status,Np] = system('sysctl -n hw.ncpu');  % not tested at all
    else, [status,Np] = system('nproc --all');
    end
    if status ~= 0, Np = ''; end
    Np = str2double(Np);
    % opt.max = round(Np^0.72); % M = {1 2 3 4 6 10 16 .. } for P = {1 2 4 8 12 24 48 ..}
    
    if ~SimOptions.runparallel
        if M0 > 1
            warning('Single-thread override due to SimOptions.runparallel = false'); 
        end 
        M = 1;
        return;
    end
    if ~isfinite(M0), M0 = Inf; end

    if ~isempty(opt.vars) && isstruct(opt.vars)
        opt.workersize = sum([opt.vars.bytes])*opt.factor + opt.workersize;
    end
    available = (maxarraysize('double')*64/8)*(1-opt.buffer);
    mbM = floor(available/opt.workersize) + 1;
    M = min([M0,opt.max,Np/2,mbM]);

    % Parse M, check that it's reasonable
    assert(isfinite(M) && mod(M,1) == 0, '%f is not a reasonable number of threads',M);
    
    % Create/update parallel pool of workers, if required 
    if opt.create && M > 1
        lic = ver();
        if ~any(strcmpi({lic.Name},'Parallel Computing Toolbox'))
            warning('Parallel Computing Toolbox is not available - using single-thread');
            M = 1;
            return
        end
        P = gcp('nocreate');
        if isempty(P) || P.NumWorkers < M
           delete(P);
           parpool(M,'EnvironmentVariables','SimOptions');
        end
    end
end