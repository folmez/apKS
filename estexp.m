function alpha = estexp(varargin)
% ESTEXP(X, xmin, xmax, X_dattype) estimates the power-law exponent on a 
% given interval.

X = varargin{1};
xmin = varargin{2};
xmax = varargin{3};
X_dattype = varargin{4};
% ------------------------------------------------------------------------
soa = (1:0.01:3.50)';    % Set of candidate power law exponents (alphas)
tX = X( X>=xmin & X<=xmax );    % Truncate X
% ------------------------------------------------------------------------
% A term needed in the likelihood function
if strcmp(X_dattype, 'REAL')
    U = (xmax.^(1-soa)-xmin.^(1-soa))./(1-soa);
elseif strcmp(X_dattype, 'INTS')
    U = zeros(length(soa), 1);
    for i = 1:length(soa)
        U(i) = sum((xmin:xmax).^(-soa(i)));
    end
end
U(isnan(U)|isinf(U)) = log(xmax/xmin);
% Compute the likelihood function of discretized alpha
L_of_alphas = (-1)*soa*sum(log(tX)) - length(tX)*log(U);
% Maximize likelihood vector
[~,L_max_index] = max(L_of_alphas);
% Determine alpha which maximizes the likelihood function
alpha = abs((-1)*soa(L_max_index));
% ------------------------------------------------------------------------
end