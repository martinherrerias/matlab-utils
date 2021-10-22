function varargout = renamefields(T,rules,varargin)
% R = RENAMEFIELDS(S,RULES,...) - Rename fields of structure/table S according to given RULES.
% [R,I] = RENAMEFIELDS(S,RULES,...) - Return also an index I of the actual substitutions.
% [M,I,U] = RENAMEFIELDS(S,RULES,...) - Return a partial copy M and index I with fields that have
%   been matched/renamed, and a residual U with all other fields. 
%
% INPUT:
%   S - (nested) structure
%   RULES - 2-column cell of the form {alias, target; ...}, or [nested] structure with fields 
%       arranged as RL.target = alias. Each alias can be in a cellstring, allowing for multiple 
%       fields to be parsed into the same target. e.g.
%
%       RENAMEFIELDS(S,{ 'w','width'; {'L','len'},'length'} }) % two replacement rules
%
%   ..,'-regexp' - read rule aliases as regex expressions, e.g. {'^len.*','length';...}
%
%   ..,'-regexrep' - use both alias and target as regex expressions, e.g. {'^len(.*)','L$1';...}
%
%   ..,'-nested' - allow nested field names as both alias & targets
%
%   ..,'-ignorecase' - ignore case when comparing aliases
%
%   ..,'-overwrite' - by default, whenever a the target field-name for a matching rule already
%       exists in the original structure, said replacement rule will be skipped.
%
%   ..,'cat',N - Allow non-injective rules, i.e. multiple original fields (sources) parsed into  
%       the same target field. The sources will then be concatenated along the Nth dimension, 
%       i.e. they must be compatible so that R.target = cat(N,S.A,S.B,...) for sources {A, B,..}.
%       The index I in such a case will take the form I.target = {'A.1','A.2',..,'B','C.1',..}, 
%       where A, B, C are the original fields, with suffixes .1, .2, .. added to any fields x 
%       with size(x,N) > 1, so that size(R.target,N) == numel(I.target).
%
%   ..,'scale',X - where X is an n-vector (n = number of rules), applies factor X(j) to any
%       matching fields of S that are copied onto M. 'scale' does not work with '-regexprep'.
% 
% OUTPUT:
%   R - (nested) structure, with matching fields replaced
%   I - (nested) structure, index with fields I.target = alias, for each replacement. If 'dim',N
%       with N > 0 is used, I.target = {alias1,alias2} indicates concatenation of several fields
%   M - (nested) structure, with replaced fields only
%   U - (nested) structue, with non-matching fields only
%
%  See RENAMEFIELDS.TEST for examples.
%
% See also: PARSELIST, NESTEDSTRUCT, SETNESTEDFIELD, NESTEDFIELDNAMES, REORDERSTRUCTURE

%     f = strsplit(field,'.');
%     v = getfield(S,{1},f{:});
    if nargin == 0, test(); return; end

    [opt,varargin] = getflagoptions(varargin,{'-regexp','-regexprep','-ignorecase','-nested','-overwrite'});
    opt.cat = 0;
    opt.scale = [];
    opt = getpairedoptions(varargin,opt,'restchk');
    
    if isstruct(rules)
        if opt.nested
           [src,tgt] = nestedstruct2cell(rules);
        else
           tgt = fieldnames(rules);
           src = struct2cell(rules);
        end
    elseif iscell(rules)
        src = rules(:,1);
        if iscellstr(rules) && size(rules,2) == 1
            tgt = src;
        else
            assert(size(rules,2) == 2,'Expecting structure or Nx2 cell array of rules')
            tgt = rules(:,2);
        end
    end

    try
        tgt = cellstr(tgt);
        try
            src = cellstr(src);
            ir = (1:numel(src))';
        catch
            src = cellfun(@(x) reshape(cellstr(x),1,[]),src,'unif',0);
            ir = repelem(1:numel(src),cellfun(@numel,src))';
            src = cat(2,src{:})';
        end
        [tgt,~,ir] = unique(tgt(ir)); % different rules might still assign the same target
    catch
        error('Invalid replacement rules');
    end
    n = numel(src);
    
    if ~isempty(opt.scale)
       assert(isnumeric(opt.scale) && isvector(opt.scale) && numel(opt.scale) == n,...
           'Expecting %d numeric vector as scale',n);
       assert(~opt.regexprep,'''scale'' does not work with ''-regexprep''');
       opt.scale = opt.scale(:);
    else
        opt.scale = ones(n,1);
    end
    
    % in case of early return (no match)
    if nargout > 2, varargout = {struct(),struct(),T}; else, varargout = {T,struct()}; end
    if isempty(src) || isempty(T), return; end
    
    if isa(T,'tabular')
        oldnames = T.Properties.VariableNames';
        
        warning_resetter = naptime('MATLAB:table:ModifiedVarnames'); %#ok<NASGU>
        S = struct2cell(table2struct(T,'ToScalar',true));
        
        if istimetable(T)
            S = [{T.Properties.RowTimes};S];
            oldnames = [T.Properties.DimensionNames{1};oldnames];
        end
        if nargout > 2, varargout{1} = T(:,[]); end
    elseif isstruct(T)
        S = T;
        if opt.nested
            [S,oldnames] = nestedstruct2cell(S);
        else
            oldnames = fieldnames(S);
            S = struct2cell(S);
        end
    else
       error('Expecting table or struct'); 
    end
    if isempty(oldnames), return; end

    % Parse rules to generate a nxm matrix idx, where idx(j,k) = 1 means oldnames(j) -> tgt(k)
    if opt.regexprep
        if opt.ignorecase, tgt = regexprep(oldnames(:),src,tgt(ir),'ignorecase');
        else, tgt = regexprepi(oldnames(:),src,tgt(ir));
        end
        n = numel(tgt);
        tochange = ~strcmp(tgt,oldnames);
        idx = eye(n);
        idx(~tochange,:) = 0;
        opt.scale = ones(1,n);
    else
        [fld,src] = meshgrid(oldnames,src);
        if opt.regexp  
            if opt.ignorecase, idx = ~cellfun(@isempty,regexpi(fld,src,'start'));
            else, idx = ~cellfun(@isempty,regexp(fld,src,'start'));
            end
        else
            if opt.ignorecase, idx = strcmpi(fld,src);
            else, idx = strcmp(fld,src);
            end
        end
        idx = groupsummary(idx,ir,'sum')'; % group columns that point to same target
        tochange = any(idx,2);
    end
    if ~any(tochange), return; end

    used = any(idx,1);
    tgt = tgt(used);
    idx = idx(:,used);
    opt.scale = opt.scale(used);
    src = oldnames;

    overwritten = ismember(tgt,oldnames(~tochange));
    
    if any(overwritten)
        if opt.overwrite
            warning('renamefields:overwrite','Replacing %s',...
                shortliststr(tgt(overwritten),'existing field'))
        else
            warning('renamefields:overwrite','Skipping %s, use ''-overwrite'' to replace',...
                shortliststr(tgt(overwritten),'existing field'))
            tgt(overwritten) = [];
            opt.scale(overwritten) = [];
            idx(:,overwritten) = [];
        end
        tochange = any(idx,2);
        if ~any(tochange), return; end
    end

    ambiguous = sum(idx,2) > 1;
    if any(ambiguous)
        error('%s map to more than one target', shortliststr(oldnames(ambiguous),'Field'));
    end

    notunique = sum(idx,1) > 1;
    merged = false(size(idx,1),1);
    if any(notunique)
    % non injective map: if a-> b and c -> b, then b = cat(opt.cat,a,b)
    
        if ~(opt.cat > 0)
            error(['Multiple matches for %s, use ''dim'',N to concatenate along ',...
            'a given dimension'],shortliststr(tgt(notunique),'target'));
        end
        
        for j = find(notunique)
            ia = find(idx(:,j));
            try
                v = cat(opt.cat,S{ia});
                src{ia(1)} = mergedsrc(oldnames(ia),cellfun(@(x) size(x,opt.cat),S(ia)));
                S{ia(1)} = v;
                merged(ia(2:end)) = true;
                [S{ia(2:end)}] = deal([]);
            catch
                warning('Failed to concatenate %s along dimension %d',...
                    shortliststr(oldnames(ia),'field'),opt.cat);
                idx(ia,j) = 0;
                notunique(j) = false;
            end
            idx = idx & ~merged;
        end
        tochange = any(idx,2);
        if ~any(tochange), return; end
    end
    assert(all(sum(idx,1) <= 1) && all(sum(idx,2) <= 1)); % debug

    [ia,ib] = find(idx);
    s = opt.scale(ib);
    toscale = (s ~= 1);
    if any(toscale)
        S(ia(toscale)) = cellfun(@times,num2cell(s(toscale)),S(ia(toscale)),'unif',0);
    end
    tgt = tgt(ib);
    tgt_idx = 1:numel(tgt);
    tgt_idx = tgt_idx(ib);
    
    tokeep = ~tochange & ~merged;
    if nargout < 3
        ic = find(tokeep);
        ia = [ia;ic];
        tgt = [tgt;oldnames(ic)];
        tgt_idx(end+(1:numel(ic))) = 0;
        [ia,ib] = sort(ia);
        tgt = tgt(ib);
        tgt_idx = tgt_idx(ib);
        src = src(ia);
        tochange(:) = true;
        tokeep(:) = false;
    else
        src = src(tochange & ~merged);
    end
    index = cell2nestedstruct(src,tgt);

    if isstruct(T)
        if opt.nested
            M = cell2nestedstruct(S(ia),tgt);
            R = cell2nestedstruct(S(tokeep),oldnames(tokeep));
        else
            M = cell2struct(S(ia),tgt);
            R = cell2struct(S(tokeep),oldnames(tokeep)); 
        end
        if isempty(R), R = struct(); end
        if isempty(M), M = struct(); end
    else
        if istimetable(T) && tochange(1)
            M = T(:,ia(2:end)-1);
            M.Properties.VariableNames = tgt(2:end);
            M.Properties.DimensionNames{1} = tgt{1};
        else
            M = T(:,ia);
            M.Properties.VariableNames = tgt;
        end
        R = T(:,tokeep);
        
        if any(notunique)
            for j = find(notunique)
                ib = find(idx(:,j),1);
                M.(tgt{tgt_idx == j}) = S{ib};
            end
        end
    end
    varargout = {M,index,R};
end

function m = mergedsrc(c,n)
    if all(n == 1), m = c(:)'; return; end
    for j = 1:numel(c)
       switch n(j)
       case 0, c{j} = {};
       case 1, c{j} = c{j};
       otherwise, c{j} = arrayfun(@(k) sprintf('%s.%d',c{j},k),1:n,'unif',0);    
       end
    end
    m = cat(2,c{:});
end

function test()
    
    x = struct('width',1:3,'len',4:6,'height',7:9,'color','purple');
    M = struct('W',1:3,'L',4:6,'H',7:9);
    R = struct('color','purple');

    [m,idx,r] = renamefields(x,{'Width','W';'len','L';'Height','H'});
    assert(isequal(m,struct('L',4:6)) && isequal(r,rmfield(x,'len')))
    assert(isequal(m,renamefields(rmfield(x,fieldnames(r)),idx)));

    [m,idx,r] = renamefields(x,{'Width','W';'len','L';'Height','H'},'-ignorecase');
    assert(isequal(m,M) && isequal(R,r))
    assert(isequal(m,renamefields(rmfield(x,fieldnames(r)),idx)));
    
    [m,idx,r] = renamefields(x,{'^([wlh]).*','${upper($1)}'},'-regexprep','-ignorecase');
    assert(isequal(m,M) && isequal(R,r))
    assert(isequal(m,renamefields(rmfield(x,fieldnames(r)),idx)));
    
    [m,idx,r] = renamefields(x,{'width','W';'len','L';'height','H'},'scale',[1,10,100]);
    assert(isequal(m,struct('W',1:3,'L',(4:6)*10,'H',(7:9)*100)) && isequal(R,r));
    assert(isequal(m,renamefields(rmfield(x,fieldnames(r)),idx,'scale',[1,10,100])));
    
    [m,idx,r] = renamefields(x,{'^[wlh].*','dims'},'-regexp','-ignorecase','cat',1);
    assert(isequal(m,struct('dims',[1:3;4:6;7:9])) && isequal(R,r));
    assert(isequal(m,renamefields(rmfield(x,fieldnames(r)),idx,'cat',1)));
    
    repl.dims.W = {'width','w'};
    repl.dims.L = {'len','length','L'};
    repl.dims.H = {'height','h'};
    [m,idx,r] = renamefields(x,repl,'-nested','-ignorecase');
    assert(isequal(m.dims,M) && isequal(R,r));
    assert(isequal(m,renamefields(rmfield(x,fieldnames(r)),idx,'-nested')));

    warning('off','renamefields:overwrite');
    x.clr = 'green';
    [m,idx,r] = renamefields(x,{'clr','color'});
    assert(isequal(r,x) && isequal(m,struct()));
    assert(isequal(m,renamefields(rmfield(x,fieldnames(r)),idx)));
    
    [m,idx,r] = renamefields(x,{'clr','color'},'-overwrite');
    assert(isequal(m,struct('color','green')) && isequal(r,rmfield(x,'clr')));
    assert(isequal(m,renamefields(rmfield(x,fieldnames(r)),idx)));
end