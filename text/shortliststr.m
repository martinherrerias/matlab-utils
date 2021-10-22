function thelist = shortliststr(x,prefix,varargin)
% L =  SHORTLISTSTR(X,THING,[MAXN])
% Return a string with the (possibly truncated) list X (cell-array of strings, or numeric array);
% joined by commas, and possibly 'and', 'to', or 'and others'; and (optionally) preceded by THING
% or THINGs (according to the number of X). 
% By default, repeated elements will appear just once (see 'unique' option), and integer-ranges 
% will be shortened, i.e. 'rows 1 to 4' instead of 'rows 1, 2, 3, and 4' (see 'mingroup' option).
%
%   THING: the name of the stuff represented on the list (e.g. 'row', 'element', 'input', ...)
%       the plural form (think irregular plurals) can be passed as the second element of a cell-
%       array, eg. {'radius','radii'}, and an uncountable THING can be specified as a single 
%       cell-array, e.g. {'stuff'}.
%       For char-arrays, 's' will be added to words not already ending in 's'. i.e. using
%       'row' is equivalent to {'row','rows'}. Default is '', i.e. no prefix before the list.
%   MAXN: will truncate the list if it exceeds the provided number of elements. Note that series of
%       integers will try to be synthesized using 'to' (e.g. 1:10 -> 1 to 10), and MAXN will then
%       limit the number of these 'lists'.
%
% L =  SHORTLISTSTR(X) - use no prefix, list up to 5 elements.
%
% L =  SHORTLISTSTR(..,'empty',M) - By default, an empty X with an empty THING will result on an
%   empty string L, but an empty list X with a non empty THING (e.g. 'row') will result in 
%   a negated plural, e.g. 'no rows'. Use 'empty' to change this behavior (see examples).
% L =  SHORTLISTSTR(..,'unique',false) - don't filter repeated elements
% L =  SHORTLISTSTR(..,'mingroup',N) - don't shorten integer-ranges smaller than N elements
% L =  SHORTLISTSTR(..,'colon',S) - insert string S (e.g. ':') between prefix and list
% L =  SHORTLISTSTR(..,'quotes',C) - attach quote-char C (e.g. '"') around every element
%           NOTE: using a quote-char sets 'mingroup' to Inf (i.e. no integer-grouping)
% L =  SHORTLISTSTR(..,'ellipsis',S) - use string S (e.g. '...') instead of 'and others'
% L =  SHORTLISTSTR(..,'delim',S) - use character(s) S instead of ', ' for joining.
% L =  SHORTLISTSTR(..,'and',S) - use character(s) S instead of ' and ' for joining last element.
% L =  SHORTLISTSTR(..,'to',S) - use character(s) S instead of ' to ' for integer ranges.
% L =  SHORTLISTSTR(..,'printlast',true) - print ELLIPSIS as the second-last element of the list
%           when skipping elements, for listings of the type 'a','k',...,'w'
%
% L =  SHORTLISTSTR(..,'newlines',true) - print elements in new lines. Actually a shortcut for:
%           'colon',':\n','ellipsis','\n...\n','joinchar','\n','and','\n','printlast',true
%
% Examples:
%     shortliststr({'j'},'row') -> 'row j'
%     shortliststr({'j','k'},'row') -> 'rows j and k'
%     shortliststr([7,11],{'vertex','vertices'}) -> 'vertices 7 and 11'
%     shortliststr({'reproduces','grows'},{'kipple'}) -> 'kipple reproduces and grows'
%     shortliststr([1,4,11,7],'',3) ->  '1, 4, 7 and others'
%     shortliststr({},'idea') -> 'no ideas'
%     shortliststr([],'number','empty','') -> ''
%     shortliststr({},'','empty','nothing at all') -> 'nothing at all'
%     shortliststr([-7, 1:10],'number') ->  'numbers -7 and 1 to 10'
%     shortliststr({'spam','spam','eggs'},'option') -> 'options spam and eggs'
% 
%     shortliststr({'spam','spam','spam'},'unique',false) -> 'spam, spam, and spam'
%     shortliststr(1:10,'number',5,'mingroup',Inf) -> 'numbers 1, 2, 3, 4, 5, and others'
%     shortliststr(1:3,'option','colon',':') -> 'options: 1, 2, and 3'
%     shortliststr({'foo','bar'},'quotes','''') -> ''foo' and 'bar''
%     shortliststr([1,4,11,7],'',2,'ellipsis','...') -> '1, 4, ...'
%     shortliststr(cellstr(['a':'z']'),'letter',4,'newlines',true) ->  'letters:
%                                                                       a
%                                                                       b
%                                                                       ...
%                                                                       z'
%
% See also: NTHINGS, STRJOIN, GETREPORT

    narginchk(1,Inf);
    if nargin < 2, prefix = ''; end
    if isempty(x) && isempty(prefix) && nargin < 3, thelist = ''; return; end
    
    % Scan for 'newlines' shortcut first, and set defaults accordingly
    [opt,remargs] = getflagoptions(varargin,{'-newlines','-printlast'});
    [opt,remargs] = getpairedoptions(remargs,opt);
    
    opt.unique = true;
    opt.mingroup = 3;
    opt.quotes = '';
    opt.maxn = 5;
    if isempty(prefix)
        opt.empty = '';
    else
        opt.empty = 'no';
    end
    
    if opt.newlines
        opt.colon = ':\n';
        opt.ellipsis = '\n...\n';
        opt.delim = '\n';
        opt.and = '\n';
        opt.to = '\n...\n';
    else
        opt.colon = '';
        opt.ellipsis = ' and others';
        opt.delim = ', ';
        opt.and = ' and ';
        opt.to = ' to ';
    end
    
    % Override defaults with any propperty-value pairs
    [opt,remargs] = getpairedoptions(remargs,opt);
    assert(numel(remargs) < 2,'Unrecognized option-value pair(s)');
    if ~isempty(remargs), opt.maxn = remargs{1}; end
    
    % sprintf() string inputs: let them take \n, \t, etc.
    for f = {'colon','quotes','ellipsis','delim','and','to'}, opt.(f{1}) = sprintf(opt.(f{1})); end
    
    opt.newlines = any(opt.delim == newline());
    if ~isempty(opt.quotes), opt.mingroup = Inf; end
    
    isint = @(x) isnumeric(x) && isscalar(x) && (isinf(x) || mod(x,1) == 0);
    assert(iscellstr(x) || isnumeric(x),'Expecting numeric vector or cell-array of strings'); %#ok<ISCLSTR>
    assert(isint(opt.maxn),'MAXN must be an integer');
    assert(isint(opt.mingroup),'''mingroup'' must be an integer');
    for f = {'colon','quotes','ellipsis','delim','and'}
        assert(ischar(opt.(f{1})),'''%s'' must be a string',f{1}); 
    end
    
    % Get plural of prefix, just in case
    if isempty(prefix), prefix = ''; plural = '';
    elseif iscellstr(prefix) && numel(prefix) <= 2
        if isscalar(prefix), prefix(2) = prefix; end    % {'kipple'}
        plural = prefix{2}; prefix = prefix{1};         % {'vertex','vertices'}
    elseif ischar(prefix)
        plural = prefix;
        if any(lower(plural(end)) == ['a':'r','t':'z']), plural(end+1)='s'; end
        % should work in most cases... excepts fors stuffs likes indexs, datas, etcs.
    else, error('THING must be a string, or 1-2 cell-array of strings'); 
    end
    
    if isempty(x)
        if isempty(prefix) || isempty(opt.empty), thelist = opt.empty;
        else, thelist = [opt.empty ' ' plural];
        end
        return;
    end
    
    % Attach colon, if required
    if ~isempty(prefix)
        prefix = [prefix,opt.colon];
        if ~any(prefix(end) == sprintf('\t\n\r ')), prefix(end+1)=' '; end 
    end
    if ~isempty(plural)
        plural = [plural,opt.colon];
        if ~any(plural(end) == sprintf('\t\n\r ')), plural(end+1)=' '; end
    end

    % Remove duplicates...
    if opt.unique, x = unique(x,'stable'); end
    n = numel(x);
    
    if isnumeric(x) 
        if opt.mingroup < Inf
        % ... and group integer-ranges
            x = findseries(x,opt.mingroup,opt.to);  
            if numel(x) < n, prefix = plural; end
        else, x = cellstr(num2str(x(:))); 
        end
    end
    if ~isempty(opt.quotes)
    % ... or add quotes
        if isnumeric(x), x = cellstr(num2str(x(:))); end
        x = cellfun(@(s) [opt.quotes,s,opt.quotes],x,'unif',0);
    end
    
    if opt.printlast && numel(x) > opt.maxn
        opt.and = opt.ellipsis;
        x(opt.maxn-1:end-1) = [];
    end
    
    switch numel(x)
       case 0, thelist = '';
       case 1
           thelist = [prefix,x{1}];
       otherwise
       if numel(x) <= opt.maxn
           thelist = strjoin(x(1:end-1)',opt.delim);
           thelist = [plural,thelist,opt.and,x{end}];
       else
           thelist = strjoin(x(1:opt.maxn)',opt.delim);
           thelist = [plural,thelist,opt.ellipsis];
       end
    end
end

function s = findseries(x,mingroupsize,TO)
% Return a (possibly reduced) cell-array of strings from a series of numbers, 
% For floating-point, non-integer numbers, return str2num for each element
% For integer indices, try to find sequences, e.g. [1:5,7] -> {'1 to 5', '7'}
    
    % for sets of integers that cannot be listed individually....
    if ~any(mod(x,1)) && numel(x) > 1
        x = sort(x(:));
        d = diff(x) == 1; d(end+1) = d(end);
        a = [d(1)==1;diff(d)==1];               % group starts
        b = [0;diff(d)==-1]; b(end) = d(end);   % group ends
        a = find(a); b = find(b);
        g = b - a + 1 >= mingroupsize;          % useful groups
        a = a(g); b = b(g);
    else
        a = []; b = [];
    end
    
    groups = cell(numel(a),1);
    ingroups = false(numel(x),1);
    for j = 1:numel(a)
        ingroups(a(j):b(j)) = true;
        groups{j} =[num2str(x(a(j))),TO,num2str(x(b(j)))];
    end
    
    singles = arrayfun(@(j) num2str(x(j)),find(~ingroups),'unif',0);
    
    [~,sortorder] = sort([x(a);x(~ingroups)]);
    s = [groups;singles];
    s = s(sortorder);
end
