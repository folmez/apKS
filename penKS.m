function [results, eeb] = penKS(varargin)
% PENKS(X, X_dattype) calculates the bounder power-law fit to the given 
% data set X. The penalty slope is 0 by default, i.e. the default is the
% KS method for bounded power law fitting.

%% Input arguments
X = varargin{1};
X_dattype = varargin{2};

min_nr_trial_pts_in_a_decade = 10;
data_title = 'Untitled data';
pen_slope = 0;
methods_needed = 5;
interval_length_threshold = 10;

i = 3;
while i<=length(varargin),
    switch varargin{i},
        case 'methods_needed',          methods_needed = varargin{i+1};
        case 'xmin_vector',             xmin_vector = varargin{i+1};
        case 'xmax_vector',             xmax_vector = varargin{i+1};
        case 'data_title',              data_title = varargin{i+1};
        case 'pen_slope',               pen_slope = varargin{i+1};
        case 'interval_length_threshold',
            interval_length_threshold = varargin{i+1};
        case 'min_nr_trial_pts_in_a_decade'
            min_nr_trial_pts_in_a_decade = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end

% Define default candidate intervals unless they are provided
if ~exist('xmin_vector','var') && ~exist('xmax_vector','var')
    LmX = elspd(X, min_nr_trial_pts_in_a_decade);
    xmin_vector = LmX;
    xmax_vector = LmX;
end

% Two-step-detection requires KS test as well
methods_needed = sort(unique([1, methods_needed]));

nr_xmins = length(xmin_vector);
nr_xmaxs = length(xmax_vector);

gm = zeros(nr_xmins,nr_xmaxs,5);

% results =   alpha, lower-bnd, upper-bnd, min , KS quality-of-fit
results = zeros(length(methods_needed),5);
eeb = zeros(5,5); % results for all methods

%% Model 
% Estimate exponents and compute quality of fit values
for j = 1:nr_xmaxs
    xmax = xmax_vector(j);
    for i = 1:nr_xmins
        xmin = xmin_vector(i);
        if xmin < xmax && ...
                sum(X>=xmin&X<=xmax)>2 && ...
                (xmax/xmin >= interval_length_threshold)
            % 1- Left end point must be smaller
            % 2- Make sure candidate interval contains more than 2 points
            % 3- Make sure candidate interval is long enough if short
            % candidate intervals are left out.
            
            % Compute power-law exponent
            alpha = estexp(X, xmin, xmax, X_dattype);
            gm(i,j,5) = alpha;
            
            % Compute K-S value
            gm(i,j,1) = estKS(X, xmin, xmax, alpha, X_dattype);
        else
            % Set everything to inf if [xmin, xmax] is not a valid
            % interval
            gm(i,j,5) = inf;
            gm(i,j,1) = inf;
        end
    end
end

% Add penalty to KS-metric
KS_penalty = log((xmin_vector.^(-1))*(xmax_vector'))*...
    (-1)*pen_slope;
gm(:,:,2) = gm(:,:,1) + KS_penalty;

% Determine bounds and corresponding exponents
counter = 1;
for i = setdiff(methods_needed, 5)
    vals = gm(:,:,i);
    % min index R and min index C, min quality-of-fit value
    min_qof = min(vals(:));
    [miR,miC] = ind2sub(size(vals),find(vals==min_qof,1));
    
    results(counter,1) = gm(miR,miC,5);
    results(counter,2) = xmin_vector(miR);
    results(counter,3) = xmax_vector(miC);
    results(counter,4) = gm(miR,miC,i);
    results(counter,5) = gm(miR,miC,1);
    counter = counter + 1;
end

% Two step detection based on KS
if ismember(5, methods_needed)
    use_all_points_in_iterative_steps = 0;
    % 1) Use KS-detected lower-bound and upper-bound
    % 2) Fix xmax as KS-detected-upper-bound and detect lower-bound
    % 3) Fix xmin detected as in step 2 and detect upper-bound
    % 4) (updated on 4-1) Loop between 2 and 3 until an equilibrium is
    % determined
    iter_limit = 100;
    prev_pl_int = zeros(iter_limit, 2);
    
    iter_step = 1;
    prev_pl_int(iter_step, :) = results(1,2:3);
    while 1
        temp_xmax = prev_pl_int(iter_step,2);
        
        if use_all_points_in_iterative_steps
            xmins = equally_log_spaced_point_detector(X(X<=temp_xmax), 100);
        else
            xmins = LmX(LmX<=temp_xmax);
        end
        
        r = penKS(X, X_dattype, 'xmin_vector', xmins, ...
            'xmax_vector', temp_xmax, 'methods_needed', (1:2), ...
            'data_title', data_title, ...
            'pen_slope', pen_slope,...
            'min_nr_trial_pts_in_a_decade', ...
            min_nr_trial_pts_in_a_decade, ...
            'interval_length_threshold', interval_length_threshold);
        
        temp_xmin = r(2,2);

        if use_all_points_in_iterative_steps
            xmaxs = equally_log_spaced_point_detector(X(X>=temp_xmin), 100);
        else
            xmaxs = LmX(LmX>=temp_xmin);
        end
        
        r = penKS(X, X_dattype, 'xmin_vector', temp_xmin, ...
            'xmax_vector', xmaxs, 'methods_needed', (1:2), ...
            'data_title', data_title, ...
            'pen_slope', pen_slope, ...
            'min_nr_trial_pts_in_a_decade', ...
            min_nr_trial_pts_in_a_decade, ...
            'interval_length_threshold', interval_length_threshold);
        
        new_pl_int = r(2,2:3);
        results(counter,:) = r(2,:);
        % Describe when the iteration stops
        if ~isempty(intersect(new_pl_int, prev_pl_int, 'rows'))
            [limit_pl_int, ~, limit_pl_index] = ...
                intersect(new_pl_int, prev_pl_int, 'rows');
            
            if limit_pl_index == iter_step
                % Single limit interval
                break;
            elseif limit_pl_index < iter_step
                % Oscillating limit interval
                % 1 - Choose the largest interval, if tied
                % 2 - Choose the interval w/ most data points, if tied
                % 3 - Choose the one to the right
                
                if limit_pl_int(2)/limit_pl_int(1) > ...
                        max(prev_pl_int(limit_pl_index+1:iter_step,2)./ ...
                        prev_pl_int(limit_pl_index+1:iter_step,1))
                    break;
                else
                    iter_step = iter_step+1;
                    prev_pl_int(iter_step,:) = new_pl_int;
                end
            end
        elseif iter_step > iter_limit
            error('TWO-STEP DETECTION DIDN''T CONVERGE IN 10 LOOPS!!!');
        else
            iter_step = iter_step+1;
            prev_pl_int(iter_step,:) = new_pl_int;
        end
    end
end

% Estimated exponents and bounds
eeb(methods_needed, :) = results;

end