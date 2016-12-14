function LmX = ELSPD(varargin)
% ELSPD (X, m) finds data points from the given data set X that are closest
% to a number of logarithmically equally spaced points so that there is at
% least m points in every decade of the data range.

X = varargin{1};    % Given data set
m = varargin{2};    % At least m data points in every decade
% ------------------------------------------------------------------------
nr_pts = ceil(log10(max(X)/min(X)))*m;
if nr_pts >= length(X)
    LmX = X;
else
    LmX = zeros(nr_pts,1);
    % Identify logarithmicall equally spaced points between the minimum 
    % and the maximum of the data
    eq_log_spaced_pts = logspace(log10(min(X)), log10(max(X)), nr_pts)';
    % Find closest data points to the set logarithmicall equally spaced 
    % points without repeating the same data point
    for i=1:nr_pts
        [~, ind_closest_nbor] = min(abs(X-eq_log_spaced_pts(i)));
        LmX(i) = X(ind_closest_nbor);
        % Avoid repeating
        X(ind_closest_nbor) = [];
    end
end
LmX = sort(LmX);
% ------------------------------------------------------------------------
if mean(abs(diff(log10(LmX))))> (1/m)*1.1 && ~isinf(m) && 0
    % 1.1 means we allow it the mean difference between successive LmX
    % terms to be off by 10 percent
    error('Average distance between successive points of L_m(X) is greater than the expected value 1/m');
end
% ------------------------------------------------------------------------
end