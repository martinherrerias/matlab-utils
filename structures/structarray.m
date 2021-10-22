function S = structarray(C,d,FILL)
% S = STRUCTARRAY(C,DIM,FILL) - Concatenate cell-array of structures C along dimension DIM,
%   Filling missing fields of any given C{j} with FILL, to avoid a disimilar-structure error
%   (MATLAB:catenate:structFieldBad).
%
%   NOTE: Fields from empty structures will be copied (and possibly filled with FILL) in all
%   other structures, even when the empty structure doesn't show in the concatenated array.
%
%   Defaults are DIM = 2, and FILL = [].

    if isempty(C), S = struct.empty; return; end
    if nargin < 2 || isempty(d), d = 2; end
    if nargin < 3, FILL = []; end
    
    ff = cellfun(@fieldnames,C,'unif',0);
    ff = uniquecell(cat(1,ff{:}))';
    [ff{2,:}] = deal(FILL);
    
    empties = cellfun(@isempty,C);
    C = C(~empties);
    
    C = cellfun(@(S) arrayfun(@(s) completestruct(s,struct(ff{:})),S),C,'unif',0);
    S = cat(d,C{:});
end