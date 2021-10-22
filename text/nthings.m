function s = nthings(n,thing,varargin)
% L = NTHINGS(N,THING,..) - Return a string of the form '1 thing'/'N things', or additional
%   special forms (e.g. 'no things' / 'all things') depending on N.
%
%   THING: the name of the stuff represented on the list (e.g. 'row', 'element', 'input', ...)
%       the plural form (for irregular plurals) can be passed as the second element of a cell-
%       array, eg. {'radius','radii'}, and an uncountable THING can be specified as a single 
%       cell-array, e.g. {'stuff'}.
%       For char-arrays, 's' will be added to words not already ending in 's'. i.e. using
%       'row' is equivalent to {'row','rows'}.
%
% L = NTHINGS(..,'rep',{A,'a',B,'b',@(x) x > C,'c',..}) - define rules for replacement of the
%   default (NUM2STR) representation of N. Numeric values are tested for equality. Function
%   handles must return a boolean value when applied to N.*
%
% L =  NTHINGS(..,'zero',Z) - shorthand for 'rep',{..,0,Z}
% L = NTHINGS(..,'all',M) - if M > 1, shorthand for 'rep',{..,M,'all'}, ignored otherwise
% L = NTHINGS(..,'-noN') - shorthand for 'rep',{..,@(x) true,''} - print only 'thing(s)'*
%
%   (*) NOTE: rules are checked in the order provided, and shorthands ('zero','all','-noN') are
%   added to the end of the list. Only the first matching rule is applied.
%
% EXAMPLES:
%     nthings(randi(4)-1,'bored programmer','all',3,'zero','no')
%     nthings(randi(4)-1,'dog','rep',{0,'no',1,'a',@(x) x > 2,'many'})
%     nthings(randi(3),{'one','not one'},'-noN')
%     nthings(randi(3),{'one','not one'},'rep',{2,'two,'},'-noN')
%     nthings(1.5,'stuff')
%
% See also: SHORTLISTSTR

    narginchk(2,Inf);
    assert(isnumeric(n) && isscalar(n) && isreal(n) && isfinite(n),'Invalid number N');
    
    [opt,varargin] = getflagoptions(varargin,{{'-noN','-non'}});
    opt = completestruct(opt,getpairedoptions(varargin,{'rep','zero','all'},'restchk'));
    if isfield(opt,'rep')
        assert(iscell(opt.rep) && mod(numel(opt.rep),2) == 0,'Invalid set of replacement rules');
        opt.rep = reshape(opt.rep,2,[])';
    else
        opt.rep = cell(0,2);
    end
    if isfield(opt,'zero'), opt.rep(end+1,:) = {0,opt.zero}; end
    if isfield(opt,'all') && opt.all > 1, opt.rep(end+1,:) = {opt.all,'all'}; end
    if opt.noN, opt.rep(end+1,:) = {@(x) true,''}; end
    
    s = num2str(n);
    
    if ~isempty(opt.rep)
        rules = opt.rep(:,1);
        values = opt.rep(:,2);
        
        assert(iscellstr(values),'Invalid replacement strings'); %#ok<ISCLSTR>
        assert(all(cellfun(@(x) isnumeric(x) || isa(x,'function_handle'),opt.rep(:,1))),...
           'Invalid replacement rules');
       
        for j = 1:numel(rules)
           if isnumeric(rules{j}) 
               if rules{j} == n, s = values{j}; break; end
           elseif rules{j}(n)
               s = values{j}; break;
           end
        end
    end
    
    isplural = n ~= 1;
    if isempty(thing) % keep empty
    elseif iscellstr(thing) && numel(thing) <= 2 %#ok<ISCLSTR>
        [thing,plural] = deal(thing{:});
        if isplural, thing = plural; end
    elseif ischar(thing)
        if isplural && any(lower(thing(end)) == ['a':'r','t':'z']), thing(end+1)='s'; end
        % should work in most cases... excepts fors stuffs likes indexs, datas, etcs.
    else, error('THING must be a string, or 1-2 cell-array of strings'); 
    end
    
    s = strtrim([s ' ' thing]);
end
