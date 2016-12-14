function ksstat = estKS(varargin)
% ESTKS(X, xmin, xmax, alpha, X_dattype) computes the KS distance
% between a power-law fit and data X

X = varargin{1};
xmin = varargin{2};
xmax = varargin{3};
alpha = varargin{4};
X_dattype = varargin{5};
% ------------------------------------------------------------------------
min_nr_data_pts_allowed = 10;
trunc_X = X( X>=xmin & X<=xmax );
% ------------------------------------------------------------------------
if length(trunc_X) < min_nr_data_pts_allowed
    ksstat = inf;
else
    % Empirical CDF
    [tX_CDF, tX] = cdfcalc(trunc_X);
    % Compute KS value between power law in (xmin,xmax) and empirical CDF
    if strcmp(X_dattype, 'INTS')
        S = sum(tX.^(-abs(alpha)));
        theor_pl_CDF = cumsum(tX.^(-abs(alpha))./S);
    elseif strcmp(X_dattype, 'REAL')
        if abs(alpha)==1
            theor_pl_CDF = (log(tX)-log(xmin))/(log(xmax)-log(xmin));
        else
            theor_pl_CDF = 1/(xmax^(1-abs(alpha))-xmin^(1-abs(alpha)))*...
                (tX.^(1-abs(alpha))-xmin^(1-abs(alpha)));
        end
    end
    % Compute KS distance
    ksstat = max(max(abs(tX_CDF(2:end)-theor_pl_CDF)),...
        max(abs(tX_CDF(1:end-1)-theor_pl_CDF)));
end
% ------------------------------------------------------------------------
end