function varargout = comparestruct(A,B,op,missing,tags)
% [C,D] = COMPARESTRUCT(A,B,OP,FILL) - Return partial nested sructures C, D that contain only
%   the fields and subfields of A,B that are different.
%
%   OP = 'xor' (default) - C,D will also contain any fields of A or B that are not found or the
%              the other, filled with FILL (default []);
%   OP = 'and' - C,D will only include (sub)fields that both A and B have in common
%   OP = 'diff' - C,D will only include fields of A that are not in B.
%
% [C,D,TXT] = COMPARESTRUCT(A,B,OP,FILL,TAGS) - where TAGS is a 2-cell-string, returns a side-by-
%   side summary of C, D with the headers in TAGS, as cellstring TXT.
%
% See also NESTEDFIELDNAMES, DISPNESTED

    if nargin == 0, error(['Usage: ' evalc('test()')]); end
    if nargin < 3, op = 'xor'; end
    if nargin < 4, missing = []; end
    if nargin < 5, tags = {}; end
    
    [A,names_A] = nestedstruct2cell(A);
    [B,names_B] = nestedstruct2cell(B);

    extra = ~ismember(names_A,names_B);
    if strcmpi(op,'and')
        A(extra) = [];
        names_A(extra) = [];
    else
        names_B = cat(1,names_B,names_A(extra));
        B = cat(1,B,repmat({missing},nnz(extra),1));
    end
    extra = ~ismember(names_B,names_A);
    switch lower(op)
    case {'diff','and'}
        B(extra) = [];
        names_B(extra) = [];
    case 'xor'
        extra = ~ismember(names_B,names_A);
        names_A = cat(1,names_A,names_B(extra));
        A = cat(1,A,repmat({missing},nnz(extra),1));
    otherwise
        error('Expecting diff/xor as OP keywords');
    end
    [~,idx] = ismember(names_A,names_B);
    B = B(idx);

    C = struct(); D = struct(); TXT={};
    
    matching = cellfun(@isequal,A,B);
    if all(matching)
        if nargout > 0, varargout = {C,D,TXT}; end
        return;
    end
    
    C = cell2nestedstruct(A(~matching),names_A(~matching));

    if nargout > 1 || nargout == 0
        D = cell2nestedstruct(B(~matching),names_A(~matching));
    end
    if nargout > 2 || nargout == 0
        TXT = side2side(C,D,tags);
    end
    
    if nargout == 0, fprintf('%s\n',TXT{:});
    else, varargout = {C,D,TXT};
    end
end

function varargout = side2side(A,B,tags)
    if nargin < 3, tags = {}; end

    txtA = strsplit(dispnested(A),newline())';
    txtB = strsplit(dispnested(B),newline())';
    if ~isempty(tags)
        txtA = [tags(1);txtA];
        txtB = [tags(2);txtB];
    end
    txtA = char(txtA);
    txtB = char(txtB);
    n = [size(txtA,1),size(txtB,1)];
    N = max(n);
    
    if n(1) < N, txtA(end+1:N,:) = ' '; end
    if n(2) < N, txtB(end+1:N,:) = ' '; end
    txt = cellstr([txtA,repmat(' ',N,4),txtB]);
    if nargout == 0
        fprintf('%s\n',txt{:});
    else
        varargout{1} = txt; 
    end
end

function test(varargin) %#ok<DEFNU>
   A = struct('a',1,'b',struct('b_1',2.1,'b_2',2.2),'c',3,'d',4);
   B = struct('a',1,'b',struct('b_1',2.1,'b_2',22,'b_3',3.1),'c',30);
   side2side(A,B,{'A =','B = '});
   [C,D] = comparestruct(A,B,'xor');
   fprintf('\ncomparestruct(A,B,''xor'') = \n'); side2side(C,D);
   [C,D] = comparestruct(B,A,'xor');
   fprintf('\ncomparestruct(B,A,''xor'') = \n'); side2side(C,D);
   [C,D] = comparestruct(A,B,'diff');
   fprintf('\ncomparestruct(A,B,''diff'') = \n'); side2side(C,D);
   [C,D] = comparestruct(B,A,'diff');
   fprintf('\ncomparestruct(B,A,''diff'') = \n'); side2side(C,D);
   [C,D] = comparestruct(A,B,'and');
   fprintf('\ncomparestruct(A,B,''and'') = \n'); side2side(C,D);
   [C,D] = comparestruct(B,A,'and');
   fprintf('\ncomparestruct(B,A,''and'') = \n'); side2side(C,D);
end