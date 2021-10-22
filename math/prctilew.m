function y = prctilew(x,p,w)
% Returns the p percentile(s) of vector of values x with weights w

	n = numel(x);
	if nargin < 3, w = ones(n,1); end
	if numel(w) ~= n, error('x and w must match in size'); end
    if ~isvector(x), x = x(:); w = w(:); end

	[x,idx] = sort(x(:));
	w = cumsum(w(idx))/sum(w);
    [w,ic,ia] = unique(w);
    if numel(ic) < numel(ia), x = accumarray(ia,x); end
	y = interp1(w,x,p/100,'linear','extrap');
end
