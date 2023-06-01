classdef GeneralizedGaussianDistribution < prob.ToolboxFittableParametricDistribution
    % GENERALIZEDGAUSSIANDISTRIBUTION or Generalized Normal Distribution.
    %   An object of the GeneralizedGaussianDistribution class.
    %
    %   GeneralizedGaussianDistribution methods:
    %       cdf     - Cumulative distribution function
    %       mean    - Mean
    %       median  - Median
    %       pdf     - Probability distribution function
    %       random  - Random number generation
    %       std     - standard deviation
    %       var     - variance
    %
    %
    %   GeneralizedGaussianDistribution properties:
    %       DistributionName    - Name of the distribution
    %       beta                - value of the shape parameter
    %       sigma               - value of the scale parameter
    %       mu                  - value of the location parameter
    %       NumParameters       - Number of parameters
    %       ParameterNames      - Names of the parameters
    %       ParameterValues     - Vector of values of parameters
    %       ParameterCovariance - Covariance matrix of estimated parameters
    %       InputData           - structure containing data used to fit the distribution
    %
    % See also: fitdist, makedist

    properties(GetAccess='public', Constant=true)
        DistributionName = 'Generalized Gaussian';
        ParameterNames = {'mu', 'alpha', 'beta'};
        ParameterDescription = {'location', 'scale', 'shape'};
        NumParameters = 3;
    end

    properties(Dependent)
        alpha
        beta
        mu
    end

    properties(GetAccess='public',SetAccess='protected')
        ParameterValues
    end

    methods(Hidden)
        function pd = GeneralizedGaussianDistribution(mu, alpha, beta)
            % default normal distribution if no additional parameters.
            if (nargin == 0)
                mu = 0;
                alpha = sqrt(2);
                beta = 2;
            end
            pd.ParameterValues = [mu, alpha, beta];
            pd.ParameterIsFixed = [true, true, true];
            pd.ParameterCovariance = zeros(pd.NumParameters);
        end
    end

    methods
        function b = get.beta(this)
            % shape parameter.
            b = this.ParameterValues(3);
        end

        function a = get.alpha(this)
            % scale parameter.
            a = this.ParameterValues(2);
        end

        function mu = get.mu(this)
            % mean/mu.
            mu = this.ParameterValues(1);
        end

        function m = mean(this)
            % same as the shape parameter.
            m = this.mu;
        end

        function m = median(this)
            % same as shape parameter.
            m = this.mu;
        end

        function sigma = std(this)
            % standard deviation. Based on scale and shape terms.
            sigma = this.alpha * ...
             sqrt(gamma(3 * (1 / this.beta) / gamma(1 / this.beta)));
        end

        function v = var(this)
            % variance. Based on scale and shape terms.
            v = this.alpha^2 * (gamma(3 * (1 / this.beta)) / gamma(1 / this.beta));
        end

        function k = kurt(this)
            % Excess kurtosis, (3 is already subtracted.)
            k = ggdkurt(this.beta);
        end

        function e = entropy(this)
            % entropy of the distribution.
            binv = 1 / this.beta;
            e = binv - log(this.beta / (2*this.alpha*gamma(binv)));
        end
    end

    methods(Static)
        function y = pdffunc(x, mu, alpha, beta)
            y = ( beta / ( 2.0 * alpha * gamma(1 / beta) ) ) * exp( -( abs(x - mu) / alpha )^beta );
        end

        function y = cdffunc(x, mu, alpha, beta)
            % CDF of gamma function provides the incomplete gamma function,
            % as well as the necessary 1 / gamma(a) normalization.
            binv = 1 / beta;
            v = gamcdf((abs(x - mu) / alpha)^beta, binv, 1) * 0.5;
            y = 0.5 + sign(x - mu) * v;

        end

        function y = invfunc(x, mu, alpha, beta)
            % TODO: is the CDF function invertable?
            y = NaN;
        end

        function y = randfunc(mu, alpha, beta, varargin)
            % Utilize the pre-written GGDRND function.
            sigma = alpha * sqrt(gamma(3 * (1/beta)) / gamma(1 / beta));
            y = ggdrnd(mu, sigma ,beta, varargin{:});
        end


        function obj = fit(x, varargin)
            % FIT distribution based on input observation.

            % get the mean value and the std.
            mu = mean(x); m1 = mean(abs(x));
            sigma = std(x);
            
            % using the Globally convergent method.
            beta = gcm2(x, m1/sigma+3, 100);
            alpha = sigma * sqrt(gamma(1/beta) / gamma(3/beta));

            % construct the "observed" data.
            obj = GeneralizedGaussianDistribution(mu, alpha, beta);
            obj.InputData = struct('data', x, 'cens', [], 'freq', []);
        end
    end
end
