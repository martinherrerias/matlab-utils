function [A,fromB] = completestruct(A,B,varargin)
% [A,fromB] = COMPLETESTRUCT(A,B) - Complete (nested) structure A with any missing fields in B.
% fromB has the same structure as A, but boolean field values indicating if A.(j) comes from B.(j)
%
% [..] = COMPLETESTRUCT(..,'valid',@f) - Treat fields from A/B that don't comply with @f as missing. 
%
% [..] = COMPLETESTRUCT(..,'warning','A>B') - Issue a warning if any field of A is not to B.
% [..] = COMPLETESTRUCT(..,'warning','B>A') - Issue a warning if any field in B is not in A.
% [..] = COMPLETESTRUCT(..,'warning','B~A') - Both of the above. 
%
% See also: GETSIMOPTIONS, GETNESTEDFIELD, SETNESTEDFIELD

    opt = getpairedoptions(varargin,{'warning','valid'},'restchk');
    if ~isfield(opt,'warning'), opt.warning = ''; end
    if isfield(opt,'valid')
        validateattributes(opt.valid,{'function_handle'},{'nonempty'});
        validfcn = opt.valid;
    else
        validfcn = [];
    end

    switch lower(opt.warning)
        case {'a>b','b<a','anotinb'}, warnif.AnotinB = true; warnif.BnotinA = false;
        case {'a<b','b>a','bnotina'}, warnif.AnotinB = false; warnif.BnotinA = true;
        case {'b~a','b~=a','a~b','a~=b','all'}, warnif.AnotinB = true; warnif.BnotinA = true;
        case {'',false,'none','off'}, warnif.AnotinB = false; warnif.BnotinA = false;
        otherwise
            error('Unknown warning type: %s',opt.warning);
    end
    
    [Bval,Bfields,implicitB] = nestedstruct2cell(B);
    if ~isempty(validfcn)
        validB = cellfun(validfcn,Bval);
    else
        validB = true(size(Bval));
    end

    if isempty(A)
        A = B;
        % if warnif.BnotinA
            warning('completestruct:BnotinA','Empty structure, set to non-empty defaults');
        % end
        if nargout > 1, fromB = cell2nestedstruct(num2cell(true(numel(Bfields),1)),Bfields); end
    else
        [Aval,Afields,implicitA] = nestedstruct2cell(A);
        if ~isempty(validfcn)
            validA = cellfun(validfcn,Aval);
            A = cell2nestedstruct(Aval(validA),Afields(validA));
        else
            validA = true(size(Aval));
        end
        if nargout > 1, fromB = cell2nestedstruct(num2cell(false(nnz(validA),1)),Afields(validA,1)); end

        AnotB = setdiff(Afields(validA),[Bfields(validB);implicitB]);
        
        if warnif.AnotinB && ~isempty(AnotB)
            warning('completestruct:AnotinB','%s not found in defaults',shortliststr(AnotB,'Field'));
        end

        BnotA = setdiff(Bfields(validB),[Afields(validA);implicitA]);
        [~,ib] = ismember(BnotA,Bfields);
        
        warning_resetter = naptime('setnestedfield:nostruct','error'); %#ok<NASGU>
        conflicts = false(size(BnotA));
        for j = 1:numel(BnotA)
            try
            % NOTE: if A.b is not a structure, it should not be removed to make place for B.b.c
                A = setnestedfield(A,BnotA{j},Bval{ib(j)});
                if nargout > 1, fromB = setnestedfield(fromB,BnotA{j},1); end
            catch ERR
                if ~strcmp(ERR.identifier,'setnestedfield:nostruct'), rethrow(ERR); end
                conflicts(j) = true;
            end
        end
        if any(conflicts)
            warning('completestruct:nostruct',...
                'Ignoring %s from defaults, as parent field(s) in base are not structures',...
                shortliststr(BnotA,'subfield'));
        end
        
        BnotA(conflicts) = [];
        if ~isempty(BnotA) && warnif.BnotinA
            warning('completestruct:BnotinA','Copying %s from defaults',shortliststr(BnotA,'field'));
        end
        
        if any(~validA)
            [ABfields,implicitAB] = nestedfieldnames(A);
            AnotAB = setdiff(Afields,[ABfields;implicitAB]);
            [~,ia] = ismember(AnotAB,Afields);
            
            for j = 1:numel(AnotAB)
                try
                    A = setnestedfield(A,AnotAB{j},Aval{ia(j)});
                    if nargout > 1, fromB = setnestedfield(fromB,AnotAB{j},1); end
                catch ERR
                    error('You should not be here, call for help: %s',getReport(ERR));
                end
            end
        end
    end
end