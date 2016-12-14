function [Xn, Xx, bin_lengths, nr_bins] = papod(varargin)
% PAPOD(X) plots the approximate pdf of the given data set X using
% logarithmic binning and plots a bounded power-law fit on top if provided.
% 
%    Example:
%       X = (1-rand(1e4,1)).^(-1/(2.5-1));
%       papod(X);
%
%   1. If you want to draw a bounded power-law fit with alpha=2.5, xmin=1
%   and xmax = 1000, try:
%       papod(X, 'power-law fit', [2.5, 1, 1000]);
%
%   2. If you want the power-law fit appear with dashed lines to emphasize
%   lack of statistical evidence for its validity, try
%       papos(X, 'power-law fit', [2.5, 1, 1000, 0]);
%
%   3. If you want to change the number of bins, try:
%       papod(X, 'nr_bins', 100);
%
%   4. If you want to change the figure title, try:
%       papod(X, 'data_title', 'hello world!')        
%
% Version 1.0	(2016 July)
% Copyright (C) Fatih Olmez (folmez@gmail.com)
% PAPOD comes with ABSOLUTELY NO WARRANTY
% 

%% Input arguments
X = varargin{1};

data_title = 'Untitled data';
nr_bins = min(round(length(X)*0.1), 5e1);
density_color = 'b';

i = 2;
while i<=length(varargin),
    switch varargin{i},
        case 'data_title',              data_title = varargin{i+1};
        case 'nr_bins',                 nr_bins = varargin{i+1};
        case 'power-law fit'
            alpha = varargin{i+1}(1);
            xmin = varargin{i+1}(2);
            xmax = varargin{i+1}(3);
            if length(varargin{i+1})==3
                von = 1;                % assume statistically valid
            elseif length(varargin{i+1})>=4
                von = varargin{i+1}(4); % statistically valid or not
            end
        otherwise,
            display(varargin{i});
            error('Unexpected inputs!!!');
    end
    i = i+2;
end

%% Check inputs
if any(X==0)
    error('There are zeroes in your data set. Remove them.');
end

%% Model
n = length(X);
bins = logspace(log10(min(X)), log10(max(X)+1e-10), nr_bins);

if isempty(setdiff(X,floor(X)))
    if n < 500
        Xx = unique(X);
        Xn = histc(X, Xx);
        warning(['Data points are used as bins as n<1000 and' ...
            ' data set consists of integeres']);
    elseif n > 500
        [Xn, Xx] = hist(X, bins);
        % bin_lengths = diff([Xx(1)-(Xx(2)-Xx(1)) Xx]); % This is wrong!!!
        bin_lengths = [Xx(2)-Xx(1), Xx(3:end)-Xx(1:end-2), ...
            Xx(end)-Xx(end-1)]*0.5;
        Xn = Xn./bin_lengths;
    end
else
    [Xn, Xx] = hist(X, bins);
    % bin_lengths = diff([Xx(1)-(Xx(2)-Xx(1)) Xx]); % This is wrong!!!
    bin_lengths = [Xx(2)-Xx(1), Xx(3:end)-Xx(1:end-2), ...
        Xx(end)-Xx(end-1)]*0.5;
        Xn = Xn./bin_lengths;
end

% Normalize
Xn = Xn/n;

if exist('xmin', 'var')
    i1 = find(Xx>=xmin,1);    % first index in the power-law region
    i2 = find(Xx>xmax,1)-1;   % last index in the power-law region
    
    C_hat = mean(Xn(i1:i2)'.*(Xx(i1:i2)'.^(alpha)));
    pl_exp_for_plot = (-1)*alpha;
    y0 = C_hat*xmin^(pl_exp_for_plot);
    y1 = C_hat*xmax^(pl_exp_for_plot);    
    data_legend_title = ['Approx PDF (' num2str(n) ' pts)'];
end

% Approximate PDF of data
nis = find(Xn);
figure
loglog(Xx(nis), Xn(nis), [density_color '.'] , 'MarkerSize', 20);
title(data_title, 'FontSize', 10);
hold on;

%% Add power-law fit to plot
if exist('xmin', 'var')
    if exist('p_val', 'var')
        if von==0
            plot([xmin xmax], [y0,y1], ['r' 'o--'], 'LineWidth', 2);
        elseif von==1
            plot([xmin xmax], [y0,y1], ['r' 'o-'], 'LineWidth', 2);
        end
    else
        plot([xmin xmax], [y0,y1], ['r' 'o-'], 'LineWidth', 2);
    end
    pl_vs_data_percentage = round(sum(X>=xmin & X<=xmax)/length(X)*100);
    pl_fit_title = ['PL(%' num2str(pl_vs_data_percentage) ...
        ') :[' num2str(xmin,'%4.2f') ',' ...
        num2str(xmax,'%4.1f') '], alpha=' num2str(alpha, '%1.2f')];
    h_legend = legend(data_legend_title, ...
        pl_fit_title, 'Location','Best');
    set(h_legend,'FontSize',15);
end
axis tight;

end