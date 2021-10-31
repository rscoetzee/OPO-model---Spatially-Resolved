function[k1p,k1s,k1i] = RHSDiff(h,FAp,FAs,FAi,FMixP,FMixS,FMixI,z,d,C1p,C1s,C1i,C2p,C2s,C2i)
                                

alpha_s = 0.00;       % Absorption
alpha_i = 0.00;  
alpha_p = 0.00; 

% rho_p = 0.00;      % Walkoff angle Rho.
% rho_s = 0.00;
% rho_i = 0.00;

% k1p = d.*h.*(-1i.*((2*(pi^2)/kp).*(sx.^2 + sy.^2) -alpha_p + (2.*pi.*sy.*tan(rho_p))).*(FAp) + FMixP);
% k1s = d.*h.*(-1i.*((2*(pi^2)/ks).*(sx.^2 + sy.^2) -alpha_s + (2.*pi.*sy.*tan(rho_s))).*(FAs) + FMixS);
% k1i = d.*h.*(-1i.*((2*(pi^2)/ki).*(sx.^2 + sy.^2) -alpha_i + (2.*pi.*sy.*tan(rho_i))).*(FAi) + FMixI);

k1p = d.*h.*((C1p -alpha_p + (C2p)).*(FAp) + FMixP);
k1s = d.*h.*((C1s -alpha_s + (C2s)).*(FAs) + FMixS);
k1i = d.*h.*((C1i -alpha_i + (C2i)).*(FAi) + FMixI);

%%
