function out = apKS(varargin)
% APKS(X) calculates the bounded power-law fit to the given data set X
% using the apKS method described in (CITE EPL)
%
% e.g.
% EPL1:	X = gsdf('EPL1', 1.5, [1 100],    1e3, 1); apKS(X);
% EPL2: X = gsdf('EPL2', 1.5, [1 100],    1e3, 1); apKS(X);
% EPL3: X = gsdf('EPL3', 1.5, [-1 1 100], 1e3, 1); apKS(X);
% IAPL: X = gsdf('IAPL', 1.5, [5 1 1000], 1e3, 1); apKS(X);
% EXP:  X = gsdf('EXP',  1,   1,          1e3, 1); apKS(X);
% NQPL: X = gsdf('NQPL', 5,   [1 10 100], 1e3, 1); apKS(X);
%
% If you find our method useful, please cite the following papers: 
% [1] EPL
% [2] RBM
%
%   1. If KS method for bounded power-law fitting is sufficient, try:
%       apKS(X, 'need_only_KS', 1);
%
%   2. If you want method to treat your data as continuous even though all
%   its elements are integers, try:
%       apKS(X, 'assume_data_is_real', 1);
%
%   3. If you want to change the p-value threshold from 0.10 to 0.05, try:
%       apKS(X, 'p_val_threshold', 0.05);
%   
%   4. If you want to include more than 10 data points from each decade of
%   the data range when searching for a bounded power-law in the data, try:
%       apKS(X, 'min_nr_trial_pts_in_a_decade', 100);
%
% Copyright (C) Fatih Olmez (folmez@gmail.com)
% apKS comes with ABSOLUTELY NO WARRANTY
%

%% Input arguments
X = varargin{1};

data_title = 'Untitled data';
plot_best_pl_fit = 1;
display_stuff = 1;
display_p_val_stuff = 1;
nr_reps = 25;
slope_diff_tol = 1e-2;
min_nr_trial_pts_in_a_decade = 10;
need_only_KS = 0;
assume_data_is_real = 0;
interval_length_threshold = 10;

need_p_val = 1;
p_val_threshold = 0.10;

slopes = [1e-5 1e0];

i = 2;
while i <= length(varargin),
    switch varargin{i},
        case 'assume_data_is_real',     assume_data_is_real = varargin{i+1};
        case 'data_title',              data_title = varargin{i+1};
        case 'need_only_KS',            need_only_KS = varargin{i+1};
        case 'display_stuff',           display_stuff = varargin{i+1};
        case 'need_p_val',              need_p_val = varargin{i+1};
        case 'display_p_val_stuff',     display_p_val_stuff = varargin{i+1};
        case 'p_val_threshold',         p_val_threshold = varargin{i+1};
        case 'interval_length_threshold',
            interval_length_threshold = varargin{i+1};
        case 'min_nr_trial_pts_in_a_decade'
            min_nr_trial_pts_in_a_decade = varargin{i+1};
        case 'nr_reps'
            nr_reps = varargin{i+1};
        case 'slope bounds'
            slopes = varargin{i+1};
        case 'slope_diff_tol',          slope_diff_tol = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end

%% Check inputs
if any(X<=0|isnan(X))
    error('Data set must contain only positive numbers.');
elseif length(slopes)~=2
    error('Exactly two slopes (one min, one max) must be inputted.');
end

% Select method (discrete or continuous) for fitting
if assume_data_is_real
    X_dattype = 'REAL';
else
    if isempty(setdiff(X,floor(X)))
        X_dattype = 'INTS';
        fprintf('Data set is discrete.\n\n');
    elseif isreal(X)
        X_dattype = 'REAL';
    else
        X_dattype = 'UNKN';
    end
    if strcmp(X_dattype,'INTS') && min(X) > 1000 && length(X)>100,
        fprintf('Data is discrete but will be treated as continous.\n\n');
        X_dattype = 'REAL';
    end
end

%% Model
slopes = [0 slopes];
% We would like to compare results against no slope
% slopes(2) = always the good slope
% slopes(3) = always the bad slope

% Preallocate vectors sufficiently long in order to find the right slope
alpha_hat_all = zeros(100, 1);
qof_values_all = zeros(100, 1);
xmin_hat_all = zeros(100, 1);
xmax_hat_all = zeros(100, 1);
p_val_hat_all = zeros(100,1);
slopes_all = zeros(100,1);
why_broken = [];

% 1) slope=0, 2) minimum slope, 3) maximum slope, 4) ...
slope_count = 0;
while (slopes(3)-slopes(2))/slopes(3) > slope_diff_tol && ...
        (slope_count==0 || ~need_only_KS)
    slope_count = slope_count+1;
    if slope_count<=3, current_slope = slopes(slope_count);
    elseif slope_count>3, current_slope = sqrt(prod(slopes(2:3)));
    end
    tSIM = tic;
    
    if display_stuff
        fprintf('Slope #%i = %1.5f\n', slope_count-1, current_slope);
    end
    
    % Compute power-law fit with the new slope
    [~,eeb] = penKS(X, X_dattype, ...
        'data_title', data_title, ...
        'pen_slope', current_slope, ...
        'min_nr_trial_pts_in_a_decade', min_nr_trial_pts_in_a_decade, ...
        'interval_length_threshold', interval_length_threshold);
    
    alpha_hat_temp = eeb(5, 1);
    xmin_hat_temp = eeb(5, 2);
    xmax_hat_temp = eeb(5, 3);
    qof_values_temp = eeb(5, end);
    
    if display_stuff
        fprintf('Power-law fit:\talpha=%1.2f on (%3.2f,%3.2f)\n', ...
            alpha_hat_temp, xmin_hat_temp, xmax_hat_temp);
    end
    
    % Estimate p-value
    p_val_hat_temp = 0;
    if need_p_val
        % Estimate p-value to the power-law fit only if it's a new fit
        [is_it_repeating, iorf] = ismember([xmin_hat_temp xmax_hat_temp], ...
            [xmin_hat_all xmax_hat_all], 'rows');
        if is_it_repeating
            p_val_hat_temp = p_val_hat_all(iorf(1));
            if display_stuff
                if p_val_hat_temp <= p_val_threshold
                    fprintf('Rejected b/c same');
                elseif p_val_hat_temp > p_val_threshold
                    fprintf('Accepted b/c same');
                end
                fprintf(' as the one for slope #%i=%1.5f\n', ...
                    iorf(1)-1, slopes_all(iorf(1)));
            end
        else
            p_val_hat_temp = estpval(X, X_dattype, ...
                alpha_hat_temp, ...
                [xmin_hat_temp, xmax_hat_temp], qof_values_temp, ...
                nr_reps, 'pen_slope', 0, ...
                'display_p_val_stuff', display_p_val_stuff, ...
                'min_nr_trial_pts_in_a_decade', ...
                min_nr_trial_pts_in_a_decade, ...
                'interval_length_threshold', interval_length_threshold);
            
            if display_stuff
                if p_val_hat_temp <= p_val_threshold
                    fprintf('Rejected with');
                elseif p_val_hat_temp > p_val_threshold
                    fprintf('Accepted with');
                end
                fprintf(' p-value %1.2f (computed from %i reps)\n', ...
                    p_val_hat_temp, nr_reps);
            end
        end
    end
    
    % Archive current one
    slopes_all(slope_count) = current_slope;
    xmin_hat_all(slope_count) = xmin_hat_temp;
    xmax_hat_all(slope_count) = xmax_hat_temp;
    alpha_hat_all(slope_count) = alpha_hat_temp;
    qof_values_all(slope_count) = qof_values_temp;
    p_val_hat_all(slope_count) = p_val_hat_temp;
    
    if slope_count==3
        if p_val_hat_all(2) <= p_val_threshold
            warning(['Minimum slope doesn''t produce a power-law fit\n' ...
                ' with good p-value (>p_val_threshold)!!!']);
            why_broken = 'minimum slope';
            break;
        elseif p_val_hat_all(3) > p_val_threshold
            warning(['Maximum slope doesn''t produce a power-law fit\n' ...
                ' with bad p-value (<p_val_threshold)!!!']);
            why_broken = 'maximum slope';
            break;
        end
    elseif slope_count>3
        % Which slope is replaced?
        if p_val_hat_temp > p_val_threshold
            slope_index_to_be_replaced = 2;
        elseif p_val_hat_temp <= p_val_threshold
            slope_index_to_be_replaced = 3;
        end
        
        % Replace the slopes vector
        slopes(slope_index_to_be_replaced) = current_slope;
    end
    
    if display_stuff
        fprintf('... finished in %3.2f minutes\n\n', toc(tSIM)/60);
    end
end

slopes_all(slope_count+1:end) = [];
xmin_hat_all(slope_count+1:end) = [];
xmax_hat_all(slope_count+1:end) = [];
alpha_hat_all(slope_count+1:end) = [];
qof_values_all(slope_count+1:end) = [];
p_val_hat_all(slope_count+1:end) = [];

[~, new_index_set] = sort(slopes_all);
slopes_all = slopes_all(new_index_set);
xmin_hat_all = xmin_hat_all(new_index_set);
xmax_hat_all = xmax_hat_all(new_index_set);
alpha_hat_all = alpha_hat_all(new_index_set);
qof_values_all = qof_values_all(new_index_set);
p_val_hat_all = p_val_hat_all(new_index_set);


% Display and plot results if a power-law interval has been found
if slope_count>3 || strcmp(why_broken, 'maximum slope')
    % Last good p-value index
    lgpvi = find(p_val_hat_all<=p_val_threshold,1) - 1;
    if isempty(lgpvi)
        lgpvi = 3;
    end
    % Percentage improvement of the power-law interval using penalty
    improvement = round(xmax_hat_all(lgpvi)/xmin_hat_all(lgpvi) * ...
        xmin_hat_all(1)/xmax_hat_all(1) * 100) - 100;

    % Display all results
    if display_stuff
        fprintf('\tKS_slope\talpha\tx_min\tx_max\tqof_val\tp-value\n');
        fprintf('\t%1.5f\t\t%3.2f\t%3.2f\t%1.1e\t%3.4f\t%3.2f\n', ...
            [slopes_all alpha_hat_all xmin_hat_all xmax_hat_all ...
            qof_values_all p_val_hat_all]');
        fprintf('\n\tKS_slope\talpha\tx_min\tx_max\tqof_val\tp-value\tImprovement(Percentage)\n');
        fprintf('\t%1.5f\t\t%3.2f\t%3.2f\t%1.1e\t%3.4f\t%3.2f\t%3d\n', ...
            [slopes_all([1,lgpvi]) alpha_hat_all([1,lgpvi]) ...
            xmin_hat_all([1,lgpvi]) xmax_hat_all([1,lgpvi]) ...
            qof_values_all([1,lgpvi]) p_val_hat_all([1,lgpvi]) ...
            [NaN; improvement]]');
    end
    
    if plot_best_pl_fit
        x0 = xmin_hat_all(lgpvi);
        x1 = xmax_hat_all(lgpvi);
        alp = alpha_hat_all(lgpvi);
        fig_title = [data_title ', slope of KS = ' ...
            num2str(slopes_all(lgpvi),'%1.2e')];
        papod(X, 'data_title', fig_title, 'power-law fit', [alp x0 x1]);
    end

    out.nr_slopes = slope_count;
    out.slopes = slopes_all;
    out.xmin_hat = xmin_hat_all;
    out.xmax_hat = xmax_hat_all;
    out.alpha_hat = alpha_hat_all;
    out.qof_values = qof_values_all;
    out.p_val_hat = p_val_hat_all;
    out.improvement = improvement;
    out.lgpvi = lgpvi; % Last good p-value index
    out.xmin = xmin_hat_all(lgpvi);
    out.xmax = xmax_hat_all(lgpvi);
    out.alpha = alpha_hat_all(lgpvi);
    out.p_val = p_val_hat_all(lgpvi);
elseif slope_count==3
    % Display all results
    if display_stuff
        fprintf('\tKS_slope\talpha\tx_min\tx_max\tqof_val\tp-value\n');
        fprintf('\t%1.5f\t\t%3.2f\t%3.2f\t%1.1e\t%3.4f\t%3.2f\n', ...
            [slopes_all alpha_hat_all xmin_hat_all xmax_hat_all ...
            qof_values_all p_val_hat_all]');
    end
    
    % Plot pdf of data
    if plot_best_pl_fit
        x0 = xmin_hat_all(1);
        x1 = xmax_hat_all(1);
        alp = alpha_hat_all(1);
        p_val = p_val_hat_all(1);
        papod(X, 'data_title', [data_title ' (non-PL)'], ...
            'power-law fit', [alp x0 x1 p_val]);
    end
    
    out.nr_slopes = slope_count;
    out.slopes = slopes_all;
    out.xmin_hat = xmin_hat_all;
    out.xmax_hat = xmax_hat_all;
    out.alpha_hat = alpha_hat_all;
    out.qof_values = qof_values_all;
    out.p_val_hat = p_val_hat_all;
    out.improvement = NaN;
    out.lgpvi = NaN; % Last good p-value index
    out.xmin = 0;
    out.xmax = 0;
    out.alpha = 0;
    out.p_val = 0;
elseif slope_count==1
    % Display all results
    if display_stuff
        fprintf('\tKS_slope\talpha\tx_min\tx_max\tqof_val\tp-value\n');
        fprintf('\t%1.5f\t\t%3.2f\t%3.2f\t%1.1e\t%3.4f\t%3.2f\n', ...
            [slopes_all alpha_hat_all xmin_hat_all xmax_hat_all ...
            qof_values_all p_val_hat_all]');
    end
    
    out.nr_slopes = slope_count;
    out.slopes = slopes_all;
    out.xmin_hat = xmin_hat_all;
    out.xmax_hat = xmax_hat_all;
    out.alpha_hat = alpha_hat_all;
    out.qof_values = qof_values_all;
    out.p_val_hat = p_val_hat_all;
    out.improvement = NaN;
    out.lgpvi = NaN; % Last good p-value index
    out.xmin = xmin_hat_all(1);
    out.xmax = xmax_hat_all(1);
    out.alpha = alpha_hat_all(1);
    out.p_val = p_val_hat_all(1);
    
    % Plot pdf of data
    if plot_best_pl_fit
        fig_title = [data_title ' (PL-fitted only by KS)'];
        papod(X, 'data_title', fig_title, 'power-law fit', out);
    end
end

end