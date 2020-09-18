function [ trans ] = curvedef_EIT(Kappa, GammaP, GammaR, G, Omega, ...
    dc, de, dr, ...
    amp, freq_data, varargin)
% Regular, boring EIT =P
% allows for anything to be "scanned"! if any variable (dc, de, G, whatever)
% has the same length as freq_data, then the each entry will be treated
% accordingly

if length(varargin)>0
   offset = varargin{1};
else
   offset = 0; 
end

Ec = dc - freq_data + 1i*Kappa/2;
Ee = de - freq_data + 1i*GammaP/2;
Er = dr - freq_data + 1i*GammaR/2;
A = amp./(Ec - (abs(G).^2)./(Ee-(abs(Omega).^2./Er)));
trans = abs(A).^2*Kappa+offset;

end

