function sbX = gsbd(varargin)
%   GSBD(X, alpha, xmin, xmax, X_dattype) generates bootstrap samples of
%   data with a semiparametric approach as described in section 4.1 of
%   Clauset et. al. 'power-law distributions in empirical data'. The
%   approach can be summarized as follows: Let n be the length of the data
%   and n_pl be the number of data points in the fitted power-law region.
%   Generate a new data set one data point at a by
%       (i) picking an element from non-power-law subset of the original
%           data set with probability 1 - n_pl/n
%       (ii) picking a power-law distributed data point (according to the
%           power-law-fit) with probability n_pl/n

X = varargin{1};
alpha = varargin{2};
xmin = varargin{3}(1);
xmax = varargin{3}(2);
X_dattype = varargin{4};
% ------------------------------------------------------------------------
n = length(X);
X_non_pl = X(X<xmin | X>xmax);
n_non_pl = length(X_non_pl);
% ------------------------------------------------------------------------
% Choose uniform randomly n1 data points from the non-power-law portion of
% the given data using randsample with replacement
n1 = sum(rand(n,1)<(n_non_pl/n));
sbX1 = randsample(X_non_pl, n1, 'true');
% sbX1 = X_non_pl(ceil(n_non_pl*rand(n1,1)));
% ------------------------------------------------------------------------
n2 = n-n1;
% Choose n2 power-law distributed data points
switch X_dattype
    case 'REAL'
        if alpha-1 ~= 0
            sbX2 = (rand(n2,1)*(xmax^(1-alpha)-xmin^(1-alpha))+ ...
                xmin^(1-alpha)).^(-1/(alpha-1));
        else
            sbX2 = exp(rand(n2,1)*(log(xmax)-log(xmin))+log(xmin));
        end
    case 'INTS'
        temp_xmin = xmin-0.5;
        temp_xmax = xmax+0.5;
        if alpha-1 ~= 0
            sbX2 = (rand(n2,1)* ...
                (temp_xmax^(1-alpha)-temp_xmin^(1-alpha))+ ...
                temp_xmin^(1-alpha)).^(-1/(alpha-1));
        else
            sbX2 = exp(rand(n2,1)*(log(temp_xmax)-log(temp_xmin))+ ...
                log(temp_xmin));
        end
        sbX2 = round(sbX2);
end
% ------------------------------------------------------------------------
sbX = sort([sbX1; sbX2]);
% ------------------------------------------------------------------------
end