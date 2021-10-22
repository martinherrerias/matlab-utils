function varargout = uniquecell(A,varargin)
% [U,ia,ic] = UNIQUECELL(C,...) - Overload of MATLAB's UNIQUE so that it works with cell elements 
%   of arbitrary class. For numeric, logical, and cell-strings, UNIQUECELL simply calls built-in
%   UNIQUE (with the subtlelty that 'rows' works with cell-strings -- by default it doesn't).
%   For other/mixed class elements, an exhaustive (slow!) search generates an index of elements
%   first, and then behaves as [~,ic,..] = UNIQUE(idx,...), U = A(u).

    if iscellstr(A) %#ok<ISCLSTR>
        rowsarg = find(strcmpi('rows',varargin),1); % check if 'rows' is an argument
        if iscell(A) && ~isempty(rowsarg)
            % Call unique() with all arguments except 'rows'
            [C,~,ic] = unique(A(:),varargin{setdiff(1:numel(varargin),rowsarg)});

            % Replace A with an array of unique indices, call unique(...,'rows') upon this array...
            ic = reshape(ic,size(A));
            [varargout{1:nargout}] = unique(ic,varargin{:});

            % ... and replace the indices back with the original cell-elements
            varargout{1} = C(varargout{1});
        else
            [varargout{1:nargout}] = unique(A,varargin{:}); 
        end
    else
    % Exhaustive search, when there's no other choice
    % TODO: for 'rows' option, row-by-row comparison (not element by element) is sure faster
    
        [opt,varargin] = getflagoptions(varargin,{'sorted','stable'});
        if opt.sorted
           warning('uniquecell:sorted','sorting not defined for cell-arrays, will be ignored'); 
        end
   
        ic = zeros(size(A));
        nu = 0;
        ia = zeros(size(A));
        for j = 1:numel(A)
            match = false;
            for k = 1:nu
                match = isequal(A(j),A(ic(k)));
                if match, break; end
            end
            if match
                ia(j) = k;
            else
                ic(nu+1) = j;
                ia(j) = nu+1;
                nu = nu+1;
            end
        end
        
        if isempty(varargin)
            ic = ic(1:nu);
            U = A(ic);
        else
        % handle 'legacy' and 'rows' options
            idx = ic(1:nu);
            [u,ic,ia] = unique(ia,varargin{:});
            U = A(idx(u));
        end
        [varargout{1:3}] = deal(U,ic,ia);
    end
end
