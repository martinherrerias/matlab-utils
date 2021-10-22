function w = lambertwlog(logx)
% W = LAMBERTWLOG(LOGX) - Functional inverse of x = w·exp(w) for x > 0 (Upper/Real branch)
% Modified from: http://blogs.mathworks.com/cleve/2013/09/02/the-lambert-w-function/
% Uses Haley method up to MaxLog, then switches to direct iteration over log(x) = log(w) + w to
% avoid machine precision errors. 
% Tolerance is NEPS·eps(w) for all logx, so the best solution within numerical precision is granted.
%
% LOGX: log(x), a numeric array of any size
%    W: solution to exp(logx) = w.*exp(w) for every LOGX

    MaxIter = 1000;
    MaxLog = 680;       % Based on trial & error it should work fine up to 696
    NEPS = 3;           % maximum precision limited to NEPS·eps(w)
        
    finite = isfinite(logx) & isreal(logx);
    
    w = zeros(size(logx));
    big = (logx > MaxLog) & finite;
    small = ~big & finite;
    w(small) = lambertw_haley(logx(small));
    
    w(logx == -Inf) = 0;
    w(logx == Inf) = Inf;
    w(isnan(logx)) = NaN;
    
    if ~any(big), return; end

    for j = find(big)'
    % if the number is too big, do it without calculating th, but solving log(th) = W + log(W) iteratively
        w(j) = logx(j)-log(logx(j));   % use log(th)-log(log(th)) as seed
        err = Inf; 
        k=0;
        while (abs(err) > NEPS*eps(w(j)))
            err = w(j)+log(w(j))-logx(j);
            w(j)= w(j)-err;
            k=k+1;
            if (k >= MaxIter)
                w(j) = NaN;
            end
        end
    end
    
    function w = lambertw_haley(logx)
    % Haley's method (modified tolerance and starting guess)

        iter = 0;
        x = exp(logx);
        w = max(1,logx); % fancier starting guesses don't seem to improve speed
        v = Inf;
        f = Inf;

        while any(abs(w - v) > NEPS*eps(w) & abs(f) > NEPS*eps(x))
           v = w;
           e = exp(w);
           f = w.*e - x;  % Iterate to make this quantity zero
           w = w - f./((e.*(w+1) - (w+2).*f./(2*w+2)));
           iter = iter+1;
           assert(iter < 100, 'lambertwlog:iter','Iteration limit reached');
        end
    end
end

