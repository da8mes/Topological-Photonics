function [pout] = fitEIT1D(xdat, ydat, pguess)
% function [pout] = fitEIT1D(xdat, ydat, pguess)
% pguess = [dc, de, dr, GammaR, G, Omega, amp, offset]

generalParams; % Kappa, GammaP


fitfun = @(params, fdata) curvedef_EIT(Kappa, GammaP, params(4), params(5), params(6),...
    params(1), params(2), params(3), ...
    params(7), fdata, params(8));

% check array shapes =P
if size(xdat,2)>size(xdat,1)
    xdat = xdat';
end
if size(ydat,2)>size(ydat,1)
    ydat = ydat';
end

opt=optimoptions('lsqcurvefit');
opt.Display = 'none';
pout = lsqcurvefit(fitfun, pguess, xdat, ydat, [],[],opt);

end

