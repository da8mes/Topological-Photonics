function [field_strength] = LaguerreGauss(XY, waists, modenums)
% function [field_strength] = LaguerreGauss(XY, waists, modenums)
% Normalized such that the maximum field_strength for LG00 would be 1
% For now, assuming the LLL [only one modenum matters]
% For now, assuming longitudinal mode structure is irrelevant (and only
% including "transverse" positions)
% 
% if length(modenums)>1
%     disp('Warning! Only LLL LG modes currently implemented; modenums(1)');
% end
if length(modenums)<2
    modenums(2)=0;
end

XY_norm = XY./repmat(waists, [size(XY,1), 1])*2;
theta = angle(XY_norm(:,1)+1i*XY_norm(:,2));
rr = sqrt(XY_norm(:,1).^2+XY_norm(:,2).^2);
m=modenums(1);
p=abs(modenums(2));
if modenums(2) <0
    disp('Warning! radial index should be >=0');
end

% original, for just the LLL
% field_strength = exp(1i*m*theta).*rr.^m.*exp(-rr.^2/4)/sqrt(2^m*factorial(m));

% all LG modes
field_strength = rr.^abs(m).*exp(-rr.^2/4).*polyval(LaguerreGen(p, abs(m)), rr.^2/2).*exp(-1i*m*theta);

% brute force normalization
field_strength = field_strength/sqrt(sum(sum(abs(field_strength).^2)));
% all modes



%         exp(1i*m*theta).*rr.^m.*exp(-rr.^2/4)/sqrt(2*pi*2^m*factorial(m))
        
end



% each radial excitation counts as much as two angular excitations