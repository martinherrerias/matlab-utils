function OPT = parallelconfig(name)
% OPT = PARALLELCONFIG(PRJNAME) - returns an options structure with unifying naming criteria for 
% RUNPARALLEL, PACKEDFCN, and MERGERESULTS. Check documentation for those functions for details.
% They are centralized here to resolve platform-dependent issues easily.
%
% In general terms: RUNPARALLEL will generate a set of files as input/output for PACKEDFCN, the
% exact names of which are defined by OPT = PARALLELCONFIG(PRJNAME):
%
%   PATH/.par/PREFIX.idx - Index file: containing Sub-file names, and a full index. Full path
%                          name is given by: OPT.path/OPT.subdir/[OPT.prefix,OPT.index]
%
%   PATH/.par/PREFIX_01.part - Sub files: each containing input-output for PACKEDFCN, as well
%   PATH/.par/PREFIX_02.part   as the row-indices for that given file. File names are given by                
%   ...                        OPT.path/OPT.subdir/[OPT.prefix_XX,OPT.part] where XX is a part
%   PATH/.par/PREFIX_NN.part   number.
%
% Variable names inside .part files {args,fcn,opt,res,idx} are now fixed (see PACKEDFCN), just
% like variable names inside .idx file {index,parts}
%
% The functions wrapped in PACKEDFCN are expected to write provisional backup files named:
%
%   PATH/.par/PREFIX_XX.part~ - Where the suffix ~ is defined by OPT.backup
% 
% See also: RUNPARALLEL, PACKEDFCN, MERGERESULTS

    if nargin < 1, name = datestr(now,'job_yymmddHHMMSS'); end
    
    [OPT.path,OPT.prefix] = fileparts(name);
    OPT.prefix = matlab.lang.makeValidName(OPT.prefix);
    if isempty(OPT.path), OPT.path = pwd(); end
    
    OPT.subdir = '.par';
    OPT.parts = '.part';
    OPT.index = '.idx';
    OPT.backup = '~';
    
    OPT.option_overrides = struct('solverlog',false,'resultsxls',false,'version','parallel');
end