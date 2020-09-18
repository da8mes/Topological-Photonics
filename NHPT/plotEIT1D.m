function [h] = plotEIT1D(xdat, params, varargin)
% function [h] = plotEIT1D(xdat, params)
% goes with function [pout] = fitEIT1D(xdat, ydat, pguess)

if length(varargin)>0
    plotstyle = varargin{1};
else
    plotstyle = '-k';
end

generalParams; % at least Kappa in here

fitfun = @(fdata) curvedef_EIT(Kappa, GammaP, params(4), params(5), params(6),...
    params(1), params(2), params(3), ...
    params(7), fdata, params(8));

h=plot(xdat, fitfun(xdat), plotstyle);

end