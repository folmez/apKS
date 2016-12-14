function [T, data_title, PDFx, PDFn] = gsdf(varargin)
% GSDF generates synthetic data from CDF.
%   e.g.
%       X = gsdf('EPL1', 1.5, [1 100],    1e3, 1);
%       X = gsdf('EPL2', 1.5, [1 100],    1e3, 1);
%       X = gsdf('EPL3', 1.5, [-1 1 100], 1e3, 1);
%       X = gsdf('IAPL', 1.5, [5 1 1000], 1e3, 1);
%       X = gsdf('EXP',  1,   1,          1e3, 1);
%       X = gsdf('NQPL', 5,   [1 10 100], 1e3, 1);

%% Input arguments
type = varargin{1};
n = varargin{4};
plot_log_log_pdf = varargin{5};

%% Model
tSIM = tic;
switch type
    case 'EPL1'
        %% Exact power-law in an interval
        alpha = varargin{2};        % power-law exponent
        xmin = varargin{3}(1);     % power-law lower bound
        M = varargin{3}(2);         % power-law upper bound
        
        data_title = ['EPL1(' num2str(n,'%1.0e')  ...
            'pts): PL(' num2str(alpha) ') in [' num2str(xmin) ',' ...
            num2str(M) ']'];
        
        C = (1-alpha)/(M^(1-alpha)-xmin^(1-alpha));
        
        U = rand(n,1); 
        T = zeros(n,1);
        
        T = ((1-alpha)*U/C + xmin^(1-alpha)).^(1/(1-alpha));
    case 'EPL2'
        %% Exact power-law in an interval with sharp transitions to non-PL
        alpha = varargin{2};    % power-law exponent
                               	% absolute lower-bound
        if length(varargin{3})==2,      t0 = 0;
        elseif length(varargin{3})==3,  t0 = varargin{3}(1);
        end
        xmin = varargin{3}(end-1);     % power-law lower bound
        M = varargin{3}(end);           % power-law upper bound

        beta = alpha*log(M/xmin)/(M-xmin);
        
        data_title = ['EPL2(' num2str(n,'%1.0e')  ...
            'pts): PL(' num2str(alpha) ') in [' num2str(xmin) ',' ...
            num2str(M) '], exp(' num2str(beta) ') otherwise'];
        
        % Compute A and C
        temp = [exp(-beta*M), -M^(-alpha); ...
            (exp(-beta*t0)+exp(-beta*M)-exp(-beta*xmin))/beta , ...
            (M^(1-alpha)-xmin^(1-alpha))/(1-alpha)]\[0;1];
        % 0 - continuity condition, 1 - probability cond
        A = temp(1);
        C = temp(2);
        
        U = rand(n,1); 
        T = zeros(n,1);
        
        CDF_at_xmin = A/(-beta)*(exp(-beta*xmin)-exp(-beta*t0));
        CDF_at_M = CDF_at_xmin + C/(1-alpha)*(M^(1-alpha)-xmin^(1-alpha));
        
        % Construct synthetic data back from uniformly distributed U in
        % interval (0,1) by inverse sampling
        index_set = U<CDF_at_xmin;
        T(index_set) = -1/beta*log( -beta/A*U(index_set)+exp(-beta*t0) );
        index_set = U<=CDF_at_M & U>=CDF_at_xmin;
        T(index_set) = ( (1-alpha)/C *(U(index_set)-CDF_at_xmin) + ...
            xmin^(1-alpha) ).^(1/(1-alpha));
        index_set = U>CDF_at_M;
        T(index_set) = -1/beta * log(-beta/A*(U(index_set)-CDF_at_M) + ...
            exp(-beta*M));        
    case 'EPL3'
        %% Exact power-law in an interval with smooth transitions to non-PL
        alpha = varargin{2};    % power-law exponent
        mu = varargin{3}(1);    % log-normal(mu,sigma)
        xmin = varargin{3}(2); % power-law lower bound
        M = varargin{3}(3);     % power-law upper bound

        sigma = sqrt((log(xmin)-mu)/(alpha-1)); % continuous slope at xmin
        beta = alpha/M; % needed for continuous slope at M
        
        if ~isreal(sigma)
            error(['Sigma is complex, change mu for lognormal. ' ...
                'Currently mu=' num2str(mu) '. ' ...
                'Try choosing mu on the other side of ' ...
                num2str(log(xmin))]);
        end
        
        data_title = ['EPL3(' num2str(n,'%1.0e')  ...
            'pts): 0 < log-n(' num2str(mu) ',' ...
            num2str(sigma) ') < ' num2str(xmin) ...
            ' < PL(' num2str(alpha) ') < ' num2str(M) ...
            ' < exp(' num2str(beta) ')'];
        
        % Compute A, C and L
        temp = [exp(-beta*M) , -M^(-alpha) ; ...
            exp(-beta*M)/beta , ...
            (M^(1-alpha)-xmin^(1-alpha))/(1-alpha) + ...
            xmin^(-alpha)/lognpdf(xmin,mu,sigma)* ...
            logncdf(xmin,mu,sigma)]\[0;1];
        A = temp(1);
        C = temp(2);
        L = C*xmin^(-alpha)/lognpdf(xmin,mu,sigma); % continuity at xmin
        
        % Compute A, C and L
        temp = [exp(-beta*M) , -M^(-alpha) , 0 ; ...
            0 , xmin^(-alpha) , -lognpdf(xmin,mu,sigma) ; ...
            exp(-beta*M)/beta , ...
            (M^(1-alpha)-xmin^(1-alpha))/(1-alpha) , ...
            logncdf(xmin,mu,sigma)]\[0;0;1];
        A = temp(1);
        C = temp(2);
        L = temp(3);
        
        U = rand(n,1); 
        T = zeros(n,1);

        CDF_at_xmin = L*logncdf(xmin,mu,sigma);
        CDF_at_M = L*logncdf(xmin,mu,sigma) + ...
            C*(M^(1-alpha)-xmin^(1-alpha))/(1-alpha);
        
        % Construct synthetic data back from uniformly distributed U in
        % interval (0,1) by inverse transform sampling
        T(U<=CDF_at_xmin) = logninv(U(U<=CDF_at_xmin)*1/L,mu,sigma);
        T(U<=CDF_at_M & U>CDF_at_xmin) = ((U(U<=CDF_at_M & U>CDF_at_xmin)...
            -CDF_at_xmin)*(1-alpha)/C + xmin^(1-alpha)).^(1/(1-alpha));
        T(U>CDF_at_M) = -1/beta * ...
            log(-beta/A*(U(U>CDF_at_M) - CDF_at_M - A*exp(-beta*M)/beta));
    case 'IAPL'
        %% IAPL
        alpha = varargin{2};            % power-law exponent
        if length(varargin{3})==2
            t0 = 0;
        elseif length(varargin{3})==3
            t0 = varargin{3}(1);        % absolute lower-bound
        end
        xmin = varargin{3}(end-1);     % power-law lower bound
        M = varargin{3}(end);           % power-law upper bound

        beta = alpha/(M+xmin);
        
        data_title = ['IAPL(' num2str(n,'%1.0e') ...
            'pts): t>' num2str(t0) ...
            ', t raised to (-' num2str(alpha) ')' ...
            ' in [' num2str(xmin) ',' num2str(M) ']'];
        
        P_T_temp = @(t) (t+xmin).^(-alpha).*exp(-beta.*t);
        C = 1/integral(P_T_temp,t0,inf);

        T = zeros(n,1);
        data_point_counter = 0;
        while data_point_counter<n
            % Rejection sampling of 100 points at every turn
            sample_size = 1000;
            u = rand(sample_size,1);
            % Construct data points according to asym. pl using inverse
            % sampling
            UU = rand(sample_size,1);
            T_prime = ((1.-UU).^(1/(1-alpha))-1).*xmin ...
                + ((1.-UU).^(1/(1-alpha))).*t0;
            % Reject if u > C_prime/C * P_T/P_T_prime
            T_prime(u>=exp(-beta*T_prime)) = []; 
            T(data_point_counter+(1:length(T_prime))) = T_prime;
            data_point_counter = data_point_counter + length(T_prime);
            if data_point_counter > n
                T(n+1:data_point_counter) = [];
            end
        end  
    case 'EXP'
        %% Exponentially distributed data
        lambda = varargin{2};
        t0 = varargin{3};
        
        data_title = ['EXP(' num2str(lambda) ...
            '), t > ' num2str(t0)];
        
        T = exprnd(1/lambda, n, 1) + t0;
        
        C = 1/exp(-lambda*t0);   
    case 'NQPL'
        %% NQPL
        k = varargin{2};
        if length(varargin{3})==2
            t0 = 0;
        elseif length(varargin{3})==3
            t0 = varargin{3}(1);        % absolute lower-bound
        end
        xmin = varargin{3}(end-1);     % power-law lower bound
        M = varargin{3}(end);           % power-law upper bound
        
        % Exponent is set to 1.5 by default!!!!
        alpha = 1.5;
        
        beta = alpha/(xmin + M);
        data_title = ['NQPL(' num2str(1.5) ...
            ') syn. data (k=' num2str(k) ...
            ') > ' num2str(t0)];
        
        P_T_temp = @(t) (xmin.^(-alpha)) .* ...
            exp(-beta.*t -(k*alpha).*((1+t/xmin).^(1/k)-1));
        C = 1/integral(P_T_temp, t0, inf);
        
        T = zeros(n,1);
        data_point_counter = 0;
        while data_point_counter<n
            % Rejection sampling of 1000 points at every turn
            sample_size = 1000;
            u = rand(sample_size,1);
            % Construct exponentially distributed random values greater
            % than t0. Simply add t0 to exponentially distributed random
            % numbers
            T_prime = exprnd(1/beta, sample_size, 1) + t0;
            % Reject if u > C_prime/C * P_T/P_T_prime
            T_prime(u>=xmin^(-alpha)*exp(-(k*alpha).* ...
                ((1+T_prime/xmin).^(1/k)-1))) = []; 
            T(data_point_counter+(1:length(T_prime))) = T_prime;
            data_point_counter = data_point_counter + length(T_prime);
            if data_point_counter > n
                T(n+1:data_point_counter) = [];
            end
        end  
    otherwise
        error('Unexpected data name');
end

fprintf('Synthetic data generation completed in %3.2f minutes\n', ...
    toc(tSIM)/60);
fprintf([data_title '\n']);

%% Theoretical pdfs
%   Function that calculates pdf values at given points for synthetic
%   distributions
    function PDFn = calcpdfvals(type, PDFx)
        PDFn = zeros(length(PDFx),1);
        switch type
            case 'EPL1'
                PDFn(PDFx<xmin) = 0;
                PDFn(PDFx>=xmin & PDFx<=M) = C* ...
                    (PDFx(PDFx>=xmin & PDFx<M)).^(-alpha);
                PDFn(PDFx>M) = 0;
            case 'EPL2'
                PDFn(PDFx<xmin) = A*exp(-beta*PDFx(PDFx<xmin));
                PDFn(PDFx>=xmin & PDFx<M) = C* ...
                    (PDFx(PDFx>=xmin & PDFx<M)).^(-alpha);
                PDFn(PDFx>=M) = A*exp(-beta*PDFx(PDFx>=M));
            case 'EPL3'
                PDFn(PDFx<xmin) = L*lognpdf(PDFx(PDFx<xmin),mu,sigma);
                PDFn(PDFx>=xmin & PDFx<M) = C* ...
                    (PDFx(PDFx>=xmin & PDFx<M)).^(-alpha);
                PDFn(PDFx>=M) = A*exp(-beta*PDFx(PDFx>=M));
            case 'IAPL'
                PDFn = C*(PDFx+xmin).^(-alpha).*exp(-beta.*PDFx);
            case 'NQPL'
                PDFn = C*xmin^(-alpha)*exp( -beta.*PDFx - ...
                    (k*alpha).*((1+PDFx/xmin).^(1/k)-1) );
            case 'EXP'
                PDFn = C*lambda*exp(-lambda*PDFx);
            otherwise
                error('Unexpected data name!');
        end
    end

%% Plot
if plot_log_log_pdf
    % Calculate and plot empirical pdf
    [~, Tx] = papod(T, 'data_title', data_title);
    % Calculate true PDF
    PDFx = Tx;
    PDFn = calcpdfvals(type, PDFx);
    % Plot true PDF
    hold on;
    loglog(PDFx, PDFn, 'b', 'LineWidth', 2);
    h_legend = legend('Approx PDF', 'Theoretical PDF', ...
        'Location','Best');
    set(h_legend, 'FontSize', 15);
    hold off;
end

end