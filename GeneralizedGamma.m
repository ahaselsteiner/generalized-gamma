% This software was written for the publication
% "Predicting wave heights for marine design by prioritizing extreme 
% events in a global model" by A.F. Haselsteiner and K-D. Thoben, see
% https://arxiv.org/pdf/1911.12835.pdf .

classdef GeneralizedGamma < handle
% We use the parameterization and variables names that are also used in 
% doi.org/10.1061/9780872629332.038
   properties
      Lambda
      C
      M
   end
   
   methods
      function obj = GeneralizedGamma(lambda, c, m)
         if nargin > 2
            obj.Lambda = lambda;
            obj.C = c;
            obj.M = m;
         end
      end
      
      function parmHat = fitDist(this, sample)
          % Estimate the parameters of the distribution using maximum 
          % likelihood estimation.
          start = [3 1 1.5];
          lb = [0 0 0];
          ub = [Inf Inf Inf];
          parmHat = mle(sample, 'pdf', @(x, lambda, c, m) ...
              this.pdf(sample, lambda, c, m), ...
              'start', start, 'lower', lb, 'upper', ub);
          this.Lambda = parmHat(1);
          this.C = parmHat(2);
          this.M = parmHat(3);
      end
      
      function f = pdf(this, x, lambda, c, m)
          pdf = @(x, lambda, c, m) c ./ gamma(m) .* lambda.^(c .* m) ...
              .* x.^(c .* m -1) .* exp(-1 .* (lambda .* x).^c);
          
          if nargin < 3
              f = pdf(x, this.Lambda, this.C, this.M);
          else
              f = pdf(x, lambda, c, m);
          end
      end
      
      function F = cdf(this, x)
          % Cumulative distribution function.
          F = nan(size(x));
          for i = 1:length(x)
              if x(i) > 0
                  F(i) = gammainc((this.Lambda .* x(i)).^this.C, this.M);
                  % Matlab's gammainc() already divides by the gamma
                  % function, see https://de.mathworks.com/help/matlab/ref/gammainc.html
                  % We could also write 
                  % F(i) = 1 - igamma(this.M, (this.Lambda .* x(i)).^this.C) / gamma(this.M)
              else
                  F(i) = 0;
              end
          end
      end
      
      function x = icdf(this, p)
          % Inverse cumulative distribution function.
          %fun = @this.cdf;
          %x0 = 2;
          %x = fzero(fun,x0) - p;
          x = 1 ./ this.Lambda .* gammaincinv(p, this.M).^(1/this.C);
      end
      
      function x = drawSample(this, n)
          if n < 2
              n = 1;
          end
          p = rand(n, 1);
          x = this.icdf(p);
      end
      
      function val = negativeloglikelihood(this, x)
          % Negative log-likelihood value (as a metric of goodness of fit).
          val = sum(-log(pdf(x, this.Lambda, this.C, this.M)));
      end
      
      function mae = meanabsoluteerror(this, sample, pi)
          % Mean absolute error (as a measure of goodness of fit).
          n = length(sample);
          if nargin < 3
              i = [1:n]';
              pi = (i - 0.5) ./ n;
          end
          xi = sort(sample);
          xhati = this.icdf(pi); % The prediction.
          mae = sum(abs(xi - xhati)) / n;
      end
      
      function ax = qqplot(this, sample, qqFig, qqAx, lineColor)
          if nargin > 2
              set(0, 'currentfigure', qqFig);
              set(qqFig, 'currentaxes', qqAx);
          else
              qqFig = figure();
          end
          if nargin < 4
              lineColor = [0 0 0];
          end
          n = length(sample);
          i = [1:n]';
          pi = (i - 0.5) ./ n;
          xi = sort(sample);
          xhati = this.icdf(pi); % The prediction.
          hold on
          plot(xhati, xi, 'kx'); 
          plot(xhati, xhati, 'color', lineColor, 'linewidth', 1.5);
          xlabel('Theoretical quantiles');
          ylabel('Ordered values');
      end
   end
end
