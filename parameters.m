function[h,deff,L,Lcav,c,lam_p,lam_s,lam_i,deltaK,w0x,w0y,w0sx,w0sy,w0ix,w0iy, ...
        AreaP,AreaS,AreaI,AreaCM,tp,PRF,vp,vs,vi,hplank,hbar,kp,ks,ki,n_p,n_s,n_i,e0, ...
        z,tRT,R1p,R2p,R1s,R2s,R1i,R2i,AR1p,AR2p,AR1s,AR2s,AR1i,AR2i,t,Pphase,  ...
        Sphase,Iphase,gg1,RT,X,Y,GridStepX,GridStepY,rho_p,rho_s,rho_i] = parameters()
    
clc    

lam_p = 1064.21*10^-9;                      % pump wavelength (m)
lam_s = 2140.*10^-9;                      % signal wavelength (m)
theta = 0;
invlam_i = (1./lam_p) - (cosd(theta)./lam_s);     
lam_i = 1/(invlam_i)                    % idler wavelength (m)

deltaK = 0;                              % wave vector mismatch (1/m)
c = 2.99792458*10^8;                          % speed of light. m/s
L = 0.012;                               % length of crystal (m).
Lcav = 0.020;                            % length of cavity (m)   
deff =  12*10^-12;                  % m/V    %SNLO is a good source for these values.
e0 = 8.854187817620*10^-12;                       % F/m
w0x = 1200*10^-6;                         % pump beam radius (1/e^2)  (m)
w0y = 1200*10^-6;
w0cmx = 0.01200;
w0cmy = 0.01200;
w0sx = 1200*10^-6;
w0sy = 1200*10^-6;
w0ix = w0sx;
w0iy = w0sy;
AreaP = (pi*w0x*w0y)/2;
AreaS = (pi*w0sx*w0sy)/2;
AreaI = (pi*w0ix*w0iy)/2;
AreaCM = (pi*w0cmx*w0cmy)/2;
tp = 10*10^-9;                           % pulse duration (s)
PRF = 100;                               % Pulse Repition rate H                           
vp = 2*pi*c/lam_p;                       % angular frequency
vs = 2*pi*c/lam_s;
vi = 2*pi*c/lam_i;
hplank = 6.62606957*10^-34;              % J.s
hbar = hplank/(2*pi);
lamp_mu = 1.06421;
lams_mu = 2.140;
lami_mu = ((1./lamp_mu) - (1./lams_mu))^-1;

n_p = sqrt(4.59423 + (0.06206/((lamp_mu^2) - 0.04763)) ...     %Sellmeier KTP
    + (110.80672/((lamp_mu^2) - 86.12171)));           %Refractive indicies

n_s = sqrt(4.59423 + (0.06206/((lams_mu^2) - 0.04763)) ...
    + (110.80672/((lams_mu^2) - 86.12171)));

n_i = sqrt(4.59423 + (0.06206/((lami_mu^2) - 0.04763)) ...
    + (110.80672/((lami_mu^2) - 86.12171)));

invlam_i = (1/cosd(theta))*(1/n_i)*((n_p/lam_p) - (n_s*cosd(theta)./lam_s));     
lam_i = 1/(invlam_i)    

kp = 2*pi*n_p/(lam_p);    % wave numbers
ks = 2*pi*n_s/(lam_s);
ki = 2*pi*n_i/(lam_i);        
h = L*0.01;    % Step size   
z = 0:h:L;
tRT = 2*(Lcav*n_i)/c;       % Calculate round trip time for given cavity L.

R1s = 0.99;     % Mirror reflectivites.
R2s = 0.5;
R1i = 0.99;   
R2i = 0.001;
R1p = 0.001; 
R2p = 0.001; 
AR1p = 0.001;                       % Reflectivity on 1st and 2nd side of crystal (if no AR-coating).
AR1s = 0.001;                       % Set to zero if crystal is coated.
AR1i = 0.001;                       % pump ~ 8.6 %, signal idler ~ 8.2 % reflectivity
AR2p = 0.001;
AR2s = 0.001;
AR2i = 0.001;

RT = 200;
gg1 = 1;   %(number of slices per round trip)
dq = 1/gg1; % 1/(number of slices per round trip)
t = (-RT/2:dq:(RT/2-dq))*tRT;

Pphase = 2*pi*rand;
Sphase = 2*pi*rand;
Iphase = 2*pi*rand;

rho_p = 0.001;      % Walkoff angle Rho.
rho_s = 0.002;
rho_i = 0.002;

% GridStep = 188*10^-6
% 
% X = [-3*w0x:GridStep:+3*w0x];      % Specify x and y intervals (m).
% Y = [-3*w0y:GridStep:+3*w0y];
% length(X)
% 1;

% X = linspace(-2.8*w0x,+2.8*w0x,64);
% Y = linspace(-2.8*w0y,+2.8*w0y,64);

X = linspace(-2.8*w0x,+2.8*w0x,64);
Y = linspace(-2.8*w0y,+2.8*w0y,64);

GridStepX = X(2)-X(1);
GridStepY = Y(2)-Y(1);
numel(X);








