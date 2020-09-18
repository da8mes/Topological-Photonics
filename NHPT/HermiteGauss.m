function [field_strength] = HermiteGauss(XY, waists, modenums)
% function [field_strength] = HermiteGauss(XY, waists, modenums)
% Normalized such that the maximum field_strength for HG00 would be 1
% For now, assuming longitudinal mode structure is irrelevant (and only
% including "transverse" positions)

XY_norm = XY./repmat(waists, [size(XY,1), 1]);

field_strength = hermiteh(modenums(1), sqrt(2)*XY_norm(:,1)) ...
    .*hermiteh(modenums(2), sqrt(2)*XY_norm(:,2)) ...
    .*exp(-sum(XY_norm.^2,2)) ...
    ./sqrt(prod(2.^modenums).*prod(factorial(modenums))) ; % normalization; since we will use this to set appropriate "g" factors, we want HG00(0,0)=1; 
        % moreover, higher modes should be normalized such that the sum of 
        % squared field amplitudes is the same as the 00 mode (in general,
        % the field max will no longer be 1)

end

