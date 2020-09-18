function [Ac, Bcc, M2, V2] = NHPT_MM(C, dc, de, d2, g, Omega, U, varargin)
% Non-Hermitian Perturbation Theory Solver, Multimode

%% input parser
p = inputParser;
% addRequired(p,'d',@istable);
% check if provided path exists
% validationFcn = @(x) exist(x,'dir');
% addParameter(p,'folder_name',d.folder_name{:},validationFcn);
addParameter(p,'TwoPhoton', true ,@islogical); % by default, include the two photon part
addParameter(p,'ShowTiming', false ,@islogical); % by default, include the two photon part

parse(p, varargin{:});

FLAG_2p= p.Results.TwoPhoton; % solve the two-photon part?
FLAG_timing = p.Results.ShowTiming;

%% Preliminaries
Ncav = length(dc);
Nat = size(g, 1);

if sum([size(Omega,1), size(U,1), size(U,2)]==size(g,1)) ~= 3
    disp('Warning! Atom parameter set sizes do not match');
end
if sum([length(C), size(g, 2)]==length(dc)) ~= 2
    disp('Warning! Cavity parameter set sizes do not match');
end

%% Single photon solve
% initialize state amplitudes
Ac = zeros(Ncav, 1);   % cavity mode amplitudes, single excitation
Ae = zeros(Nat, 1); % p-state amplitudes, single excitation
Ar = zeros(Nat, 1); % rydberg amplitudes, single excitation

% solve single photon part (1st order NHPT)
M=zeros(Ncav, Ncav);
V=C';
% assemble the matrix... do it with for loops, since time should be
% negligible anyway
for ix = 1:Ncav % which entry of C?
    for jx = 1:Ncav
        if ix==jx % if this is the mode corresponding to C...
            M(ix, jx) = dc(ix)+sum(abs(g(:, ix)).^2./(abs(Omega).^2/d2-de));
        else
            M(ix, jx) = sum(conj(g(:, ix)).*g(:, jx)./(abs(Omega).^2/d2-de));
        end
    end
end
% solve for single photon manifold cavity amplitudes
Ac = linsolve(M, V);

% keyboard
for ix = 1:Nat
    Ae(ix) = sum(g(ix, :)*Ac)/(abs(Omega(ix))^2/d2-de);
end
Ar = -Omega.*Ae/d2;

%% Two photon steady state solve
if FLAG_2p
    tic
    % % initialize state amplitudes (not actually necessary, but conceptually
    % % useful)
    % Bcc = zeros(Ncav*(Ncav+1)/2, 1); % two cavity modes occupied;
    %     % indexed (2*Ncav-mode1)*(mode1-1)/2+mode2
    % Bce = zeros(Ncav*Nat, 1); % one cavity, one p-state;
    %     % indexed [Ncav*(Ncav+1)/2+]  Nat*(mode#-1)+atom#
    % Bcr = zeros(Ncav*Nat, 1); % one cavity, one rydberg;
    %     % indexed [Ncav*(Ncav+1)/2+Ncav*Nat+]  Nat*(mode#-1)+atom#
    
    % % index for the ix, jx double cavity mode:
    % (2*Ncav-ix)*(ix-1)/2+jx
    
    M2 = zeros(Ncav*(Ncav+1)/2+2*Ncav*Nat, Ncav*(Ncav+1)/2+2*Ncav*Nat);
    V2 = zeros(Ncav*(Ncav+1)/2+2*Ncav*Nat, 1);
    
    % build V2 first:
    % double cavity
    for ix = 1:Ncav % mode 1
        for jx = 1:Ncav % mode 2 % include EITHER mode becoming excited
            if ix==jx
                mult = sqrt(2);
                totix =(2*Ncav-ix)*(ix-1)/2+jx;
            elseif ix>jx
                totix=(2*Ncav-jx)*(jx-1)/2+ix;
                mult = 1;
            elseif jx>ix               
                totix=(2*Ncav-ix)*(ix-1)/2+jx;
                mult = 1;
            end
            
            V2(totix) = V2(totix)+ mult*conj(C(ix))*Ac(jx);
        end
    end
    % cavity plus P-state
    offset = Ncav*(Ncav+1)/2;
    for ix = 1:Ncav
        for jx = 1:Nat
            V2(offset+Nat*(ix-1)+jx)=conj(C(ix))*Ae(jx);
        end
    end
    % cavity plus Rydberg
    offset = offset+ Ncav*Nat;
    for ix = 1:Ncav
        for jx = 1:Nat
            V2(offset+Nat*(ix-1)+jx)=conj(C(ix))*Ar(jx);
        end
    end
    
    
    % % Build the main matrix M2:
    % double cavity rows
    for ix = 1:Ncav % mode 1
        for jx = ix:Ncav % mode 2
            row = (2*Ncav-ix)*(ix-1)/2+jx;
            if ix==jx % same mode
                %             mult = sqrt(2);
                M2(row, row)=2*dc(ix);
                M2(row, Ncav*(Ncav+1)/2+Nat*(ix-1)+(1:Nat))=sqrt(2)*conj(transpose(g(:, ix)));
            else % different modes
                %             mult = 1;
                M2(row, row)=dc(ix)+dc(jx);
                M2(row, Ncav*(Ncav+1)/2+Nat*(ix-1)+(1:Nat))=conj(transpose(g(:, jx)));
                M2(row, Ncav*(Ncav+1)/2+Nat*(jx-1)+(1:Nat))=conj(transpose(g(:, ix)));
            end
        end
    end
    
    % helper coefficients
    Xee=zeros(Nat, Nat);
    Yee=zeros(Nat, Nat);
    Zee=zeros(Nat, Nat);
    Xer=zeros(Nat, Nat);
    Yer=zeros(Nat, Nat);
    Zer=zeros(Nat, Nat);
    for j = 1:Nat
        for k = 1:Nat
            % % % % % % % % % %         SLOW AND HAS ERROR:
            % %         % for Bee^jk, ORIGINAL
            % %         denom = -2*(d2+de)*conj(Omega(j))*conj(Omega(k))*(-Omega(j)*(U(j,k)+2*d2)*conj(Omega(k))-2*Omega(j)*de*conj(Omega(k)))...
            % %             -(-Omega(j)*conj(Omega(j))*conj(Omega(k))+conj(Omega(k))*((-U(j, k)-2d2)*(d2+de)+Omega(k)*conj(Omega(k))))...
            % %             *(-2*de*(d2+de)*conj(Omega(k))+conj(Omega(k))*(-Omega(j)*conj(Omega(j))+Omega(k)*conj(Omega(k))));
            % %         Xee(j,k) = (U(j,k)*d2^2*conj(Omega(k))^2+2*d2^3*conj(Omega(k))^2+2*U(j,k)*d2*de*conj(Omega(k))^2 ...
            % %             +4*d2^2*de*conj(Omega(k))^2+U(j,k)*de^2*conj(Omega(k))^2+2*d2*de^2*conj(Omega(k))^2 ...
            % %             -Omega(j)*(d2+de)*conj(Omega(j))*(conj(Omega(k))^2)-Omega(k)*(d2+de)*conj(Omega(k))^3) ...
            % %             /denom;
            % %         Zee(j,k) = (-U(j,k)*(d2+de)*conj(Omega(j))*conj(Omega(k))^2-2*d2*(d2+de)*conj(Omega(j))*conj(Omega(k))^2 ...
            % %             +Omega(j)*conj(Omega(j))^2*conj(Omega(k))^2-Omega(k)*conj(Omega(j))*conj(Omega(k))^3)...
            % %             /denom;
            % %         Yee(j,k) = (-U(j,k)*(d2+de)*conj(Omega(k))^3-2*d2*(d2+de)*conj(Omega(k))^3 ...
            % %             -Omega(j)*conj(Omega(j))*conj(Omega(k))^3+Omega(k)*conj(Omega(k))^4)...
            % %             /denom;
            
            % for Bee^jk, using FullSimplify
            denom = conj(Omega(k))^2*(2*abs(Omega(j))^2*(d2+de)*(U(j,k)+2*d2+2*de) ...
                -(-(U(j,k)+2*d2)*(d2+de)-abs(Omega(j))^2+abs(Omega(k))^2)*(-2*de*(d2+de)-abs(Omega(j))^2+abs(Omega(k))^2));
            Xee(j,k)= (d2+de)*conj(Omega(k))^2*((U(j,k)+2*d2)*(d2+de)-abs(Omega(j))^2-abs(Omega(k))^2)  ...
                /denom;
            Yee(j,k)= conj(Omega(k))^3*(-(U(j,k)+2*d2)*(d2+de)-abs(Omega(j))^2+abs(Omega(k))^2) ...
                /denom;
            Zee(j,k)= conj(Omega(j))*conj(Omega(k))^2*(-(U(j,k)+2*d2)*(d2+de)+abs(Omega(j))^2-abs(Omega(k))^2)  ...
                /denom;
            
            % for Ber^ (used FullSimplify =P Should go back and redo above)
            % this set is about 3x faster than above
            denom = ((U(j,k)+2*d2)*(d2+de)-abs(Omega(j))^2)*(2*de*(d2+de)-abs(Omega(j))^2)...
                - abs(Omega(k))^2*((d2+de)*(U(j,k)+2*d2+2*de)+2*abs(Omega(j))^2)+abs(Omega(k))^4;
            Xer(j,k) =  Omega(k)*((U(j,k)+2*d2)*(d2+de)+abs(Omega(j))^2-abs(Omega(k))^2) ...
                /denom;
            Yer(j,k) =  (-(U(j,k)+2*d2)*(2*de*(d2+de)-abs(Omega(j))^2)+2*de*abs(Omega(k))^2) ...
                /denom;
            Zer(j,k) =  (-Omega(k)*conj(Omega(j))*(U(j,k)+2*d2+2*de))...
                /denom;
        end
    end
    Xee(isnan(Xee))=0; % any couplings that evaluated to NaN can be set to zero
    Yee(isnan(Yee))=0; % any couplings that evaluated to NaN can be set to zero
    Zee(isnan(Zee))=0; % any couplings that evaluated to NaN can be set to zero
    Xer(isnan(Xer))=0; % any couplings that evaluated to NaN can be set to zero
    Yer(isnan(Yer))=0; % any couplings that evaluated to NaN can be set to zero
    Zer(isnan(Zer))=0; % any couplings that evaluated to NaN can be set to zero
    
    
    % one photon one P-state atom rows
    offsetCE = Ncav*(Ncav+1)/2;
    offsetCR = Ncav*(Ncav+1)/2+Ncav*Nat;
    for kx = 1:Ncav
        for mx = 1:Nat
            row = offsetCE+Nat*(kx-1)+mx;
            
            % cc
            indices = [];
            for px = 1:Ncav
                indices = [indices, (2*Ncav-min([kx, px]))*(min([kx, px])-1)/2+max([kx, px])];
            end
            M2(row, indices) = g(mx, :);
            M2(row, indices(kx)) = sqrt(2)*M2(row, indices(kx));
            
            % energy of this state
            M2(row, row) = dc(kx)+de;
            
            % direct coupling frm rydberg state
            M2(row, offsetCR+Nat*(kx-1)+mx) = conj(Omega(mx));
            
            otherAt = [1:(mx-1), (mx+1):Nat];
            for px = 1:Ncav
                % ce contribution from ee
                %             [Ncav*(Ncav+1)/2+]  Nat*(mode#-1)+atom#
                M2(row, offsetCE+Nat*(px-1)+mx)=M2(row, offsetCE+Nat*(px-1)+mx)+...
                    sum(conj(g(otherAt, kx)).*g(otherAt, px).*Xee(otherAt, mx));
                %
                M2(row, offsetCE+Nat*(px-1)+otherAt)=M2(row, offsetCE+Nat*(px-1)+otherAt)...
                    +transpose(conj(g(otherAt, kx)).*Xee(otherAt, mx)*g(mx, px));
                
                % cr contribution from ee
                M2(row, offsetCR+Nat*(px-1)+mx)=M2(row, offsetCR+Nat*(px-1)+mx) ...
                    +sum(conj(g(otherAt, kx)).*g(otherAt, px).*Yee(otherAt, mx));
                
                M2(row, offsetCR+Nat*(px-1)+otherAt)=M2(row, offsetCR+Nat*(px-1)+otherAt) ...
                    +transpose(conj(g(otherAt, kx)).*Zee(otherAt, mx)*g(mx, px));
            end
        end
    end
    
    
    % one photon one Rydberg atom rows
    for kx = 1:Ncav
        for mx = 1:Nat
            row = offsetCR+Nat*(kx-1)+mx;
            
            % excitation from CE by blue
            M2(row, offsetCE+Nat*(kx-1)+mx) = Omega(mx);
            
            % energy of this state
            M2(row, row) = dc(kx)+d2;
            
            % couplings from the B_er term
            otherAt = [1:(mx-1), (mx+1):Nat];
            for px = 1:Ncav
                % ce contribution from er
                %             [Ncav*(Ncav+1)/2+]  Nat*(mode#-1)+atom#
                M2(row, offsetCE+Nat*(px-1)+mx)=M2(row, offsetCE+Nat*(px-1)+mx)+...
                    sum(conj(g(otherAt, kx)).*g(otherAt, px).*Xer(otherAt, mx));
                
                M2(row, offsetCE+Nat*(px-1)+otherAt)=M2(row, offsetCE+Nat*(px-1)+otherAt)+...
                    transpose(conj(g(otherAt, kx)).*Xer(otherAt, mx)*g(mx, px));
                
                % cr contribution from er
                M2(row, offsetCR+Nat*(px-1)+mx)=M2(row, offsetCR+Nat*(px-1)+mx) ...
                    +sum(conj(g(otherAt, kx)).*g(otherAt, px).*Yer(otherAt, mx));
                
                M2(row, offsetCR+Nat*(px-1)+otherAt)=M2(row, offsetCR+Nat*(px-1)+otherAt) ...
                    +transpose(conj(g(otherAt, kx)).*Zer(otherAt, mx)*g(mx, px));
            end
        end
    end
    
    
    AssemblyTime=toc;
    if FLAG_timing
        disp(['Took ', num2str(AssemblyTime), ' seconds to assemble the matrix']);
    end
    tic
    
    % solve the linear system to obtain the second-order NHPT result
    % X = linsolve(M2, V2);
    % if Nat*Ncav > 2000 % use GPU [03/13/18 Logan: my GPU is not good enough to help]
    %     X = gather(gpuArray(M2)\gpuArray(V2));
    % else % use CPU
    X = linsolve(M2, V2);
    % end
    
    
    Bcc = X(1:Ncav*(Ncav+1)/2);
    Bce = X(offsetCE+1:offsetCE+Ncav*Nat);
    Bcr = X(offsetCR+1:offsetCR+Ncav*Nat);
    
    SolveTime=toc;
    if FLAG_timing
        disp(['Took ', num2str(SolveTime), ' seconds to solve the linear system']);
    end
else % don't solve the two-photon part
    Bcc = NaN(Ncav*(Ncav+1)/2,1);
    Bce = NaN(Ncav*Nat,1);
    Bcr = NaN(Ncav*Nat,1);
end


% keyboard
end

