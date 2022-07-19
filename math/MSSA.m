function [Xr,Sx,RC,L,opt] = MSSA(X,varargin)
% [XR,Sx,RC,L,OPT] = MSSA(X,..) - Gap-filling & filtering via Multi-Channel Singular Spectrum 
%   Analysis. Also called Empirical Orthogonal Function (EOF) Analysis, SSA is a form of Principal
%   Component Regression. It works by finding the eigenvectors of the autocovariance matrix of X, 
%   (i.e. the covariance of X and M-1 lagged copies of X), selecting a reduced set of them, and
%   expressing the signal as a sum of (lagged and weighted) orthogonal "Principal Components".
%
%   NOTE: currently uses a simple white-noise threshold to discriminate "significant" modes,
%       which is not very robust.
%
%   TODO: improve methods for determination of optimal KMAX and M, e.g. MC / MTM.
%
% INPUT: X is an [N,D] matrix, with rows representing observations of D covariates. NaN/Inf will 
%   be considered missing values.
%
%   ..,'maxlag',M - Use M-1 lagged copies of X, instead of the default M = 3. (!) NOTE that the
%       algorithm involves the eigenfactorization of an M·D square covariance matrix,
%       so things can quickly scale out of proportion.
%
%   ..,'kmax',K - Calculate only the first K <= M·D reconstructed components. This returns a
%       filtered (cleaner) XR, and improves performance (convergence of the last RC's when gap-
%       filling long time series can be a pain).
%
%   ..,'-periodic' - uses circular-shifts instead of ignoring missing values at the edges
%
%   ..,'-symmetric' - [experimental!] use lags 1-M:M-1 instead of 0:M-1 to introduce non-
%       causal correlation components.
%
%   ..,'lags',U - [even more experimental!] provide your own vector of lags, e.g. to try things 
%       like seasonal correlation components.
%
%   ..,'tol',T - (defaul 1e-3) NRMSE tolerance of reconstructed signal (for gap-filling)
%   ..,'maxiter',N - max. number of iterations (for each component)
%
%   ..,'exclude',F - Don't attempt to fill gaps flagged in F (boolean Nx1 vector). NOTE: this 
%       is NOT a means to specify invalid/missing values, these are allways set to ~isfinite(X).
%
%   ...,'-plot' - meant for debugging: plot iterative solution.
%   
% OUTPUT:
%   XR: reconstructed X, with any gaps filled, and potentially filtered (see KMAX, below).
%
%   Sx: estimated model uncertainty on XR, (for missing samples, Sx is prediction uncertainty).
%
%   CAUTION: this is NOT taken straight from [1-5] but adapted starting from [6] and using the
%   simplifying assumption that noise in the sample covariance propagates to the EOFs as
%   uncorrelated uncertainty in the eigenvalues. This produces very conservative estimates
%   (very high model uncertainty) when eigenvalues are close. The problem is explained in [6].
%
%   TODO: covariance shrinkage?, EOF rotation?
%
%   RC: [N,D,K] array of reconstructed components. Each plane RC(:,:,k) contains the k'th component
%       of the signal, in order of decreasing eigenvalue (i.e. decreasing explained variance).
%       Partial sums Y = sum(RC(:,:,1:k),3) can be used to get 'cleaner' versions of X.
%
%   L: K-vector of eigenvalues, in decreasing order. The variance explained by each RC is given
%       by L/sum(L).
%
%   OPT: Structure of used options, possibly completed with defaults.
%
% REFERENCES:
%  [1] M. Ghil et al., “Advanced Spectral Methods for Climatic Time Series,” Reviews of Geophysics,
%       vol. 40, no. 1, pp. 3-1-3–41, 2002.
%  [2] D. Kondrashov and M. Ghil, “Spatio-temporal filling of missing points in geophysical data 
%       sets,” Nonlinear Processes in Geophysics, vol. 13, no. 2, pp. 151–159, May 2006.
%  [3] H. Björnsson and S. A. Venegas, “A manual for EOF and SVD analyses of climate data.” 1997.
%  [4] M. R. Allen and L. A. Smith, “Monte Carlo SSA: Detecting irregular oscillations in the 
%       Presence of Colored Noise,” Journal of Climate, vol. 9, no. 12, pp. 3373–3404, Dec. 1996
%  [5] Y. S. Unal and M. Ghil, “Interannual and interdecadal oscillation patterns in sea level”
%  [6] G. R. North, T. L. Bell, R. F. Cahalan, and F. J. Moeng, “Sampling Errors in the Estimation 
%      of Empirical Orthogonal Functions,” Monthly Weather Review, vol. 110, no. 7, Jul. 1982
%  [7] Y. Zhang, “Quantiﬁcation of Prediction Uncertainty for Principal Components Regression and
%      Partial Least Squares Regression,” University College London.

% See also: FILLGAPS, SVD 

    if nargin == 0, test(); return; end

    N = size(X,1);
    D = size(X,2);
    
    [opt,varargin] = getflagoptions(varargin,{'-plot','-symmetric','-periodic','-robustcov'});
    
    opt.maxiter = 1000;
    opt.tol = 1e-2;
    opt.exclude = false(N,1);
    opt.threshold = 1; % norminv(0.975,1,1/sqrt(N-2));
    
    opt.lags = [];
    opt.maxlag = [];
    opt.Kmax = [];
    opt = getpairedoptions(varargin,opt,'restchk');
    
    validateattributes(opt.exclude,{'numeric','logical'},{'vector','binary','numel',N},...
        '','exclude');
         
    tofill = ~isfinite(X) & ~opt.exclude(:);
    missing = tofill | opt.exclude(:);
    X(missing) = NaN;
    
    anytofill = any(tofill,'all');
    % anymissing = any(missing,'all');
    nf = nnz(any(tofill,2));
    
    if isempty(opt.maxlag)
        
        r = lag1corr(X);
        tau = max(-1./log(abs(diag(r)))); % Unal and Ghil (1995)
        
        if anytofill
            [r,~,x] = find(diff([false(1,D);tofill;false(1,D)],1));
            maxgap = max(r(x < 0) - r(x > 0));
        else
            maxgap = 0;
        end
        
        opt.maxlag = max([ceil(N/10),2.5*maxgap,ceil(tau)]);
        opt.maxlag = min([opt.maxlag,N/2,floor(maxarraysize()/(4*N*D))]);
    end
    
    if ~isempty(opt.lags)
        validateattributes(opt.lags,{'numeric'},{'vector','integer','<',N,'>',-N},'','lags');
    else
        M = opt.maxlag;
        if opt.symmetric
            % opt.lags = [0:2:M-1,-(1:2:M-1)];
            opt.lags = ceil(-M/2):ceil(M/2);
        else
            opt.lags = 0:M-1;
        end
    end
    % [~,ic] = unique((opt.lags(:)),'sorted'); 
    % opt.lags = opt.lags(ic)';
    opt.lags = unique(opt.lags(:)','sorted');
    if ~any(opt.lags == 0), opt.lags = [0,opt.lags]; end % lag 0 at the start!
    
    M = numel(opt.lags);

    if isempty(opt.Kmax)
        opt.Kmax = M*D;
        hardsetK = false;
    else
        validateattributes(opt.Kmax,{'numeric'},{'scalar','integer','positive'},'','Kmax');
        opt.Kmax = min(M*D,opt.Kmax);
        hardsetK = true;
    end

    Xr = zeros(N,D);             % Reconstructed signal
    
    % these are kept for re-scaling
    XM = mean(X,1,'omitnan');
    XS = std(X,1,'omitnan');

    XS = max(eps(1),XS); % avoid division by zero on 'flat' signals
    X = (X - XM)./XS;

    % F = trajectorymatrix(X,opt.lags,opt.periodic);
    % % C = robustcov(F,'method','ogk','rows','pairwise');
    % C = cov(F,'partialrows');
    % % r = corr(F,X,'rows','pairwise');

    [C,V] = nanxcov(X,opt.lags,opt.periodic);
    assert(~any(isnan(C),'all'),'Undefined covariance terms, signal might be flat or too short');

    [Q,L] = eig(C);
    L = diag(L);                  % extract the diagonal
    [L,ind] = sort(L,'descend');  % sort eigenvalues
    Q = Q(:,ind);                 % ... and eigenvectors
    
    % [lo,hi] = MCSSA(X,Ed,opt.lags,N,opt.periodic);  % Not working, needs FIX!
    % significant = L >= hi; % | L < lo;

    if hardsetK
        significant = (1:M*D)' <= opt.Kmax;
    else
        significant = L > opt.threshold;
        % significant(K+1:end) = false; % (§)
    end

    % Eigenvector uncertainty, taken from North et al. 1982.
    sL = diag(Q'*V*Q);
    sL = sqrt(max(0,sL));
    significant = significant & L > 1 + 1.96*sL;
    % sL(~significant) = 0;

    % % Uncertainty of eigenvalues estimated as sqrt(2/DOF)·L following (Ghil et al. 2002)*
    % % (2/N)lambda factor can be traced back to North et al. 1982.
    % % (*) NOTE: Effective tau from Allen & Smith results in very low DOF in many cases. Since
    % % according to Ghil et al. N/M is a conservative estimate, this is now set as bound.
    % 
    % DOF = N/(1.5*tau);
    % DOF = max(DOF,N/M);
    % sL = sqrt(2/DOF).*L.*significant;

    Q(:,~significant) = 0;                
    S = Q*Q'; % signal projection matrix

    % DEBUG
    if opt.plot
        GUIfigure('MSSA'); % clf();
        subplot(2,1,1); cla(); hold on;
        L = max(L,0);
        plot(L);

        % errorbar(1:M*D,L,sqrt(2/DOF).*L,'.-');
        % errorbar(1:M*D,(lo+hi)/2,(hi-lo)/2,'.');
        errorbar(find(significant),L(significant),sL(significant),'.');
        plot([0,M*D+1],[1,1].*opt.threshold,'k--');
        set(gca,'yscale','log');
        xlim([0,nnz(L)])
        ylim([0.1,Inf])
    end

    % ... then set missing values to zero (mean) (§)
    X(  missing) = 0;
    
    min_err = Inf;
    Xr(:) = 0;
    for iter = 1:opt.maxiter 

        % Fill the multi-channel trajectory matrix
        F = trajectorymatrix(X,opt.lags,opt.periodic);

        last_Xr = Xr; 

        RC = F*S;
        [Xr,RC] = invtrajectory(RC,opt.lags,opt.periodic);    

        % update partially reconstructed signal (§)
        X(missing) = Xr(missing);

        err = rms(Xr(tofill) - last_Xr(tofill));
        min_err = min(err,min_err);

        % (§) Upon convergence of the RC(k) signal, move to RC(k+1)
        if ~anytofill || err < opt.tol/sqrt(nf) || err > min_err, break; end

        if iter == opt.maxiter
            warning('Failed to converge, NRMS = %0.2e',err)
            break;
        end
    end

    % DEBUG
    if opt.plot
        X(missing) = NaN;
        GUIfigure('MSSA');
        title(sprintf('K = %s',shortliststr(find(significant))));
        % title(sprintf('K = %d',K));
        subplot(2,1,2);
        cla(); hold on;
        plot(X,'.');
        plot(Xr,'-');
        X(missing) = Xr(missing);
    end

    if any(significant(opt.Kmax+1:end))
        warning(' should use Kmax > %d',find(significant,1,'last'));
    end
    % significant(opt.Kmax+1:end) = [];
    % Xr = sum(RC(:,:,significant),3);
    
    % decay timescale of the AR(1) series should be < 1/10 of series length (Ghil et al. 2002)
    r = lag1corr(X);
    tau = max(-1./log(abs(diag(r)))); % Allen & Smith, 1996
    if any(tau > N/10)
       warning(['Time series seems short. Detected AR(1) decay timescale is < %0.1f samples, '...
                'at least %d points are recommended'],tau,ceil(tau)*10); 
    end

    X(missing) = NaN;
    Se = mean((X - Xr).^2,1,'omitnan');
    
    % covariance of the smoothed signal Sr comes from assuming variance "modes" in eigenvalue
    % matrix Q proportional to L*Q, and propagating into the reconstructed components and Xr
    Sr = zeros(N,D);
    for j = 0:D-1
        ix = (1:M)+M*j;
        % Sr(:,j+1) = (2/DOF).^2*mean((RC(:,ix).*L(ix)').^2,2,'omitnan');
        Sr(:,j+1) = mean((RC(:,ix).*sL(ix)').^2,2,'omitnan');
    end
    Sx = Sr + Se;

    % For missing samples, estimated uncertainty is propagated iteratively
    if anytofill
        
        s = sensitivity(S,opt.lags);
    
        [r,c] = find(missing);
        for i = 1:opt.maxiter
            
            last_Sx = Sx;
            Sa = zeros(N,D);
            for k = 1:D
                ink = (c == k);
                if ~any(ink), continue; end   

                ij = s.i(:) + r(ink)';
                ok = ij > 0 & ij <= N;                
                ij = ij + (s.j(:)-1)*N;       
                vv = s.v{k}(:).*Sx(missing(:,k))';

                Sa(:) = Sa(:) + accumarray(ij(ok),vv(ok),[N*D,1],@sum);
            end
            % Sx = Sa + Sr + Se;
            Sx = min(1,Sa + Sr + Se); % avoid blow up

            err = abs(last_Sx(tofill) - Sx(tofill));
            if all(err < opt.tol)
                break; 
            end
            
            if i == opt.maxiter, warning('Failed to converge for SE'); end
        end
    end
    Sx = sqrt(Sx);
    
    % DEBUG/TEST
    if opt.plot
        GUIfigure('MSSA');
        subplot(2,1,2); hold on;
        for j = 1:D
            patch([1:N,N:-1:1]',[Xr(:,j)+1.96*Sx(:,j);flipud(Xr(:,j)-1.96*Sx(:,j))],j,...
                'facealpha',0.2,'edgecolor','none');
        end
    end

    Sx(opt.exclude,:) = NaN;
    Xr(opt.exclude,:) = NaN;
    
    % Rescale
    % XS = XS./std(X,1,'omitnan');
    % XM = XM - mean(X.*XS,1,'omitnan');
    
    Xr = Xr.*XS + XM;
    Sx = Sx.*XS;
        
    RC = permute(reshape(RC,N,M,D),[1,3,2]);
    RC = RC.*XS;
    RC(:,:,1) = RC(:,:,1) + XM; % assign offset to main component
end

function [C,V] = nanxcov(A,lags,periodic)
% Estimate the covariance matrix F'F, where F = trajectorymatrix(X,lags,periodic)

    [N,D] = size(A);
    M = numel(lags);

    L = lags - lags';
    [lags,~,ic] = unique(L);
    u = numel(lags);
    
    r = zeros(D,D,M,M);
    n = zeros(D,D,M,M);
    for i = 1:D
        xi = A(:,i);
        for k = 1:i
            for j = 1:u
                if abs(lags(j)) > N - 3, continue; end  % n(i,k,jj) = 0
                jj = (ic == j);
                if i == k && lags(j) == 0
                    r(i,k,jj) = 1;
                    n(i,k,jj) = nnz(~isnan(xi));
                else
                    xkj = circshift(A(:,k),lags(j));
                    if ~periodic
                        if lags(j) > 0
                            xkj(1:lags(j)) = NaN;
                        else
                            xkj(end+lags(j)+1:end) = NaN;
                        end
                    end
                    ok = ~(isnan(xkj) | isnan(xi));
                    nok = nnz(ok);
                    n(i,k,jj) = nok;
                    if nok < 3
                        % r(i,k,jj) = 0;
                    else
                        r(i,k,jj) = corr(xi(ok),xkj(ok));
                    end
                end
            end
            r(k,i,:,:) = permute(r(i,k,:,:),[1,2,4,3]);
            n(k,i,:,:) = permute(n(i,k,:,:),[1,2,4,3]);
        end
    end
    r = reshape(permute(r,[3,1,4,2]),M*D,M*D);
    n = reshape(permute(n,[3,1,4,2]),M*D,M*D);
    
    s = repelem(std(A,0,1,'omitnan'),1,M);
    S = s.*s';
    C = r.*S;
    
    V = hypot(C,S).*sqrt(2./max(0,n-3));
    
%     % Confidence interval using Fisher-transform
%     z = atanh(r);
%     dz = 1.96./sqrt(max(0,n-3));
%     dr = (tanh(z+dz) - tanh(z-dz))/1.96;
    
end

function F = trajectorymatrix(X,lags,periodic)
% Fill the multi-channel trajectory matrix: F = [B{1}(X), ... B{m}(X)]
% where B{k} is a backward shift by lags{k} - min(lags).
%
% e.g. for lags = 0:M, F(k,:) = [X(1,k),X(1,k-1),..,X(1,k-M),X(2,k),X(2,k-1),..,X(D,k-M)]
%
% If PERIODIC, elements that don't exist (i.e. X(j,k-l) for k-l < 1) are set to NaN, 
% and the matrix F will be extended by max(lags) - min(lags) to accomodate the lagged series.
% Otherwise B{k}(X) === circshift(X,lags(k)) and size(F,1) = size(X,1).

    [N,D] = size(X);
    M = numel(lags);
    
    [a,b] = bounds(lags);
    
    if issparse(X)
        [r,c,v] = find(X);
        r = r + lags(:)' - a;
        c = (c-1)*M + (1:M);
        v = repmat(v,1,M);
        if periodic
            r = mod(r-1,N)+1;
        else
            out = r < 1 | r > N;
            r(out) = [];
            c(out) = [];
            v(out) = [];
        end
        F = sparse(r,c,v);
    else
        if periodic
            F = zeros(N,M*D);
            for m = 1:M
                F(:,(0:D-1)*M + m) = circshift(X,lags(m)-a,1);
            end
        else
            F = zeros(N+b-a,M*D);
            X(N+1:N+b-a,:) = NaN;
            for m = 1:M
                F(:,(0:D-1)*M + m) = circshift(X,lags(m)-a,1);
                % F(1:lags(m),(0:D-1)*M + m) = NaN;
            end
        end
    end
end

function [X,F] = invtrajectory(F,lags,periodic)
% Recover X from reconstructed components

    [N,D] = size(F);
    M = numel(lags);
    D = D/M;
    
    [a,b] = bounds(lags);

    for m = 1:M
        F(:,(0:D-1)*M + m) = circshift(F(:,(0:D-1)*M + m),a-lags(m),1);
    end
    if ~periodic
        N = N-b+a;
        F(N+1:end,:) = [];
    end
    
    X = zeros(N,D);
    for j = 0:D-1
        X(:,j+1) = mean(F(:,(1:M)+M*j),2,'omitnan');
    end
end

function s = sensitivity(S,lags)
% Returns values for dX(l,j)/dX(i,k) in the form of cell array s.v such that:
%
%   s.v{k}(l - i + b-a + 1,j) = | dX(l,j)/dX(i,k) |²
%
% In other words, s.v[k} is a "kernel" for dX/dX(i,k), with row b-a+1 corresponding to lag 0.

    M = numel(lags);
    D = size(S,1)/M;
    
    [a,b] = bounds(lags);
    ne = (b-a)*2+1;
      
    s.i = repmat((a-b:b-a)',1,D);
    s.j = repmat(1:D,ne,1);
    s.v = cell(D,1);
    for k = 1:D
        Sk = zeros(ne,M*D);
        Sk(lags + (b-a) + 1,:) = S((k-1)*M+(1:M),:);
        for m = 1:M
            Sk(:,(0:D-1)*M + m) = circshift(Sk(:,(0:D-1)*M + m),-lags(m),1);
        end
        s.v{k} = zeros(ne,D);
        for j = 0:D-1
            s.v{k}(:,j+1) = mean(Sk(:,(1:M)+M*j),2);
        end
        s.v{k} = s.v{k}/s.v{k}((b-a)+1,k);
        s.v{k} = s.v{k}.^2;
    end
end

function [r,c0] = lag1corr(X)
% ubiased lag-1 autocorrelation according to Allen & Smith 1996

    r = corr(X,circshift(X,1),'rows','pairwise');
    c0 = cov(X,'partialrows');
    [N,~] = size(X);
    % r = diag(r);
    
    r0 = r;
    last_r = r;
    
    for j = 1:10
        a = 2/N^2;
        rk = r;
        mu2 = 1/N + a*(N-1)*rk;
        for k = 2:N-1
            rk = rk.*r;
            mu2 = mu2 + a*(N-k)*rk;
        end

        r = r0.*(1 - mu2) + mu2;
        if all(abs(r - last_r) < 1e-5), break; end
        last_r = r;
    end
    
    c0 = c0./(1-mu2);
end

% function [lo,hi] = MCSSA(X,Ed,lags,N,periodic)
% % According to Ghil et al. 2002 and Allen & Smith, 1996. the white-noise threshold used above
% % is not very robust, and one has to account for coloured-noise by using a surrogate AR(1)
% % process and evaluating its projected spectral response on the calculated EOFs.
% %
% % TODO: see what's wrong... it seems to flag arbitrary components
% 
%     steps = 500;
%     P = 90;
% 
%     [N,D] = size(X);
%     M = numel(lags);
%     
%     [r,c] = lag1corr(X);
% 
%     % decay timescale  of the AR(1) noise
%     tau = max(-1./log(abs(diag(r))));
%     if any(tau > N/10)
%        warning('Time series seems short'); 
%     end
%     
%     % Mdl = estimate(varm(D,1),X);
%     Mdl = varm('Constant',zeros(D,1),'AR',{zeros(D)},'Covariance',eye(D));
%     % Mdl = varm('Constant',zeros(D,1),'AR',{r},'Covariance',diag(diag(c)-diag(r).^2));
%     Ln = zeros(M*D,steps);
%     
%     for j = 1:steps
%         X = simulate(Mdl,N);
% 
%         F = zeros(N,D*M);
%         for m = 1:M
%             F(:,(0:D-1)*M + m) = circshift(X,-lags(m),1);
%             if ~periodic
%                 F(1:-lags(m),(0:D-1)*M + m) = NaN;
%                 F(end-lags(m)+1:end,(0:D-1)*M + m) = NaN;
%             end
%         end
%         C = cov(F,'omitrows');
%         L = Ed'*C*Ed;
%         Ln(:,j) = diag(L);
%     end
%     Ln(Ln <= 0) = NaN;
%     lim = prctile(Ln,[50-P/2,50 + P/2],2);
%     lo = lim(:,1);
%     hi = lim(:,2);
% end


function test()
% Modified from:
% A beginner's guide to SSA(Singular Spectrum Analysis)
% by David Claessen (CERES-ERTI) and Andreas Groth (LMD)
% CERES-ERTI, Ecole Normale Supérieure, Paris
% http://environnement.ens.fr/IMG/file/DavidPDF/SSA_beginners_guide_v9.pdf
    % 
    % rng(0.123456789);
    % N = 500;
    % DIMS = 1;
    % NOISE = 0.05;
    % GAPS = 10;
    % GAP_LEN = ceil(0.8*N/GAPS);

    N = 100*(randi(10)+2);
    DIMS = randi(3);
    NOISE = rand(1)*0.2;
    
    GAPS = randi(floor(N/8));
    GAP_LEN = ceil(rand(1)*N/GAPS);

    k = 10;
    t = (1:N)';
    
    for d = DIMS:-1:1
        freq = rand(1,k)*20 + 1;
        phase = rand(1,k)*360;
        wt = rand(1,k); wt = wt/sum(wt);

        X(:,d) = sum(sind( 360*freq.*t/N + phase ).*wt,2);
    end

    noise = NOISE*randn(size(X)); % Gaussian noise
    X = X+noise;
    
    d = size(X,2);
    gaps = arrayfun(@(a,n) a + (0:n-1),randi(N*d-GAP_LEN,GAPS,1),randi(GAP_LEN,GAPS,1),'unif',0);
    X([gaps{:}]) = NaN;

    MSSA(X,'-plot'); % '-symmetric','periodic','lags',U'
end
