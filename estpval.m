function p_values =  estpval(varargin)
% ESTPVAL(X, X_dattype, alpha, [xmin xmax], KS, nr_reps) estimates the 
% p-value corresponding to the power-law fit with power-law exponent alpha
% over the interval [xmin, xmax] with the KS distance value KS using 
% nr_reps semiparametric bootstrap samples.

X = varargin{1};
X_dattype = varargin{2};
alpha = varargin{3}; % must be positive
bounds = varargin{4};
qof_val = varargin{5};
nr_reps = varargin{6};

display_p_val_stuff = 1;
pen_slope = 0;
min_nr_trial_pts_in_a_decade = 10;
interval_length_threshold = 10;
% ------------------------------------------------------------------------
if alpha<0
    fprintf('\nalpha = %3.2f\n',alpha);
    warning('Exponent alpha must be positive! Its sign is changed!!!'); 
    alpha = abs(alpha);
end
% ------------------------------------------------------------------------
i = 7;
while i<=length(varargin),
    switch varargin{i},
        case 'display_p_val_stuff', display_p_val_stuff = varargin{i+1};
        case 'pen_slope',           pen_slope = varargin{i+1};
        case 'interval_length_threshold'
            interval_length_threshold = varargin{i+1};
        case 'min_nr_trial_pts_in_a_decade'
            min_nr_trial_pts_in_a_decade = varargin{i+1};
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end
% ------------------------------------------------------------------------
% Deafult parameters
which_qof_method = 1;   % Use KS method when validating
% ------------------------------------------------------------------------
% Estimate p-value for the fit if the KS value of the fit is worth it
reasonable_KS_upper_limit = 20/sqrt(length(X(X>bounds(1) & X<bounds(2))));
if qof_val < reasonable_KS_upper_limit
    tic;
    % KS-values of power-law-fits to semiparametric bootstrap samples
    sbqof = zeros(nr_reps,1);
    
    if display_p_val_stuff
        fprintf('[0]\t');
        fprintf('p_KS');
        fprintf('\tTime\t\tBounds(%6.2f, %6.2f)', bounds);
        fprintf('\tK-S(%6.4f)\t', qof_val);
        fprintf('alpha(%1.2f)\n', alpha);
    end
    
    for i = 1:nr_reps
        % Construct a semiparametric bootstrap sample from X
        sbX = gsbd(X, alpha, bounds, X_dattype);
        
        % Determine a power-law-fit and KS-distance
        [~,eeb] = penKS(sbX, X_dattype, ...
            'methods_needed', which_qof_method, ...
            'pen_slope', pen_slope, ...
            'min_nr_trial_pts_in_a_decade', ...
            min_nr_trial_pts_in_a_decade, ...
            'interval_length_threshold', interval_length_threshold);
        
        sbX_bounds_hat = eeb(which_qof_method,2:3);
        sbX_pl_fit_KS_val = eeb(which_qof_method, end);
        sbX_pl_fit_alpha = eeb(which_qof_method, 1);
        
        sbqof(i) = eeb(which_qof_method,end);
        if display_p_val_stuff
            fprintf('[%i]\t%6.4f', i, sum(sbqof(1:i)>=qof_val)./i);
            fprintf('\t[%4.2fm]', toc/60)
            fprintf('\t\t%6.2f, %6.2f', sbX_bounds_hat);
            fprintf('\t\t%6.4f', sbX_pl_fit_KS_val)
            fprintf('\t\t%1.2f\n', sbX_pl_fit_alpha);
        end        
    end
    p_values = sum(sbqof>=qof_val)./nr_reps;
else
    if display_p_val_stuff
        fprintf('KS-value (%3.4f) ', qof_val);
        fprintf('of the power-law fit (%3.2f,%3.2f)', bounds);
        fprintf('is too high.');
        fprintf(' (>%3.4f)\n', reasonable_KS_upper_limit);
        fprintf('Rejected without estimating p-value\n');
    end
    p_values = -Inf;   
end
% ------------------------------------------------------------------------
end