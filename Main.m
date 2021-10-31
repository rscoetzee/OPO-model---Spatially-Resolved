clc
close all force
close all
clear 
tic

 [h,deff,L,Lcav,c,lam_p,lam_s,lam_i,deltaK,w0x,w0y,w0sx,w0sy,w0ix,w0iy, ...
 AreaP,AreaS,AreaI,AreaCM,tp,PRF,vp,vs,vi,hplank,hbar,kp,ks,ki,n_p,n_s,n_i,e0, ...
 z,tRT,R1p,R2p,R1s,R2s,R1i,R2i,AR1p,AR2p,AR1s,AR2s,AR1i,AR2i,t,Pphase,  ...
 Sphase,Iphase,gg1,RT,X,Y,GridStepX,GridStepY,rho_p,rho_s,rho_i] = parameters();


Pumpvalues = 10;
MeasuredOut = 1;

% Pumpvalues = linspace(1,20,5);
% MeasuredOut = linspace(1,20,5);

SimulatedOut = zeros(1,length(Pumpvalues));
OutOut = zeros(1,length(Pumpvalues));
SNLO = zeros(1,length(Pumpvalues));
SNLO2 = SNLO;

feature('jit', 'on'); feature('accel', 'on');

wait1 = waitbar(0,'Total Progress'); %'Position',[0.25 0.4 0.25 0.2]);
wait2 = waitbar(0,'Calculating','Position',[580 300 300 50]);


[x,y] = meshgrid(X,Y);

dx   = GridStepX;
nffx = length(X);
fx   = 1/dx;
dsx  = fx/(nffx);
ssx  = -fx/2:dsx:(fx/2)-dsx;    % spatial frequency X

dy   = GridStepY;
nffy = length(Y);
fy   = 1/dy;
dsy  = fy/(nffy);
ssy  = -fy/2:dsy:(fy/2)-dsy;    % spatial frequency Y

[sx,sy] = meshgrid(ssx,ssy);

% Preallocation stuff
Track8b = zeros(1,length(t));
Track8 = Track8b; Track8c = Track8b;
CC = zeros(1,length(X));
CCP = CC; CCS = CC;
IdOut = zeros(1,length(Pumpvalues));
SigOut = zeros(1,length(Pumpvalues));
ResPump = zeros(1,length(Pumpvalues));
PumpP = zeros(1,length(t));
SigP = zeros(1,length(t));
IdP = zeros(1,length(t));
CenterPhase = zeros(1,length(t));
% Preallocation stuff

for ii = 1:length(Pumpvalues)
    
    
% Eseed = 9.3330*10^-20; %J    % page 425 in Yariv, gives expression, for a single mode. (photon energy*bandwidth)
  Eseed = 2.9321*10^-19;       % Multimode - see page 432 Yariv.
% Eseed = 0.001;               % So assume bandwidth ~ cm^-1 ~ 3*10^8 Hz ~ 0.3 Ghz.
                               % E = hv = 9.3330e -20 J
                               % Therefore input power ~ 2.799*10^-9 W, 2.799 nW.

Ppseed = sqrt(4*log(2)/pi)*(Eseed/tp);
Ep(ii) = Pumpvalues(ii)/1000;                 % pulse energy J
Pp(ii) = sqrt(4*log(2)/pi)*(Ep(ii)/tp);


initial_pumpenvelope = sqrt(Pp(ii)).*exp(-1.*(t.^2)./((tp/sqrt(2*log(2))).^2));
Seed_env             = sqrt(Ppseed).*exp(-1.*(t.^2)./(tp/sqrt(2*log(2))).^2);
% Seed_env             = 1*10^-9;   %Eseed*PRF;

trapz(t,initial_pumpenvelope.^2);         % integrate under envelope = pulse energy.
trapz(t,initial_pumpenvelope.^2).*PRF;    % integrate under envelope * PRF = pulse power.

AA = exp((1i.*2.*pi.*c.*t./lam_p)).*initial_pumpenvelope.^2;
% figure; plot(t,real(AA));



Signal_power_out = cell(length(t),1);
Idler_power_out  = cell(length(t),1);
Total_output     = cell(length(t),1);
OutputPulse      = cell(length(t),1);
DepPump          = cell(length(t),1);
Theta1           = cell(length(t),1);
Theta2           = cell(length(t),1);
pumpphase = Theta1; sigphase = Theta1; idphase = Theta1;
A0p  = zeros(length(X),length(Y));       
A0s  = A0p; A0i = A0p;          
A0pr = A0p; 
A0sr = A0pr; A0ir = A0pr;

A0prgrid = cell(1,length(t));  
A0pin = A0prgrid; A0sin = A0pin; A0iin = A0pin;

d = 1;
Noshots = 1;
Outvalue = zeros(1,Noshots);

E0 = zeros(1,length(t));
E0s = E0; E0i = E0;

Exyp = cell(1,length(t));
Exys = Exyp; Exyi = Exyp;  

Ixyp = cell(1,length(t));
Ixys = Ixyp; Ixyi = Ixyp;

   for mm = 1:length(t)

       Ixyp{mm} = (initial_pumpenvelope(mm).^2).*(1/(AreaP)).*exp(-(2).*((x.^2)+(y.^2))./(w0x*w0y));
       Ixys{mm} = (Seed_env(mm).^1).*(1/(AreaS)).*exp(-(2).*((x.^2)+(y.^2))./(w0sx*w0sy));
       Ixyi{mm} = (Seed_env(mm).^1).*(1/(AreaI)).*exp(-(2).*((x.^2)+(y.^2))./(w0ix*w0iy));
    
   end
      
   for mm = 1:length(t)  

        Exyp{mm} = sqrt(2.*(Ixyp{mm})./(e0.*n_p.*c)).*exp(-1i.*Pphase);
        Exyi{mm} = sqrt(2.*(Ixyi{mm})./(e0.*n_i.*c)).*exp(-1i.*Iphase);
        Exys{mm} = sqrt(2.*(Ixys{mm})./(e0.*n_s.*c)).*exp(-1i.*Sphase); 

   end
   

A0prgrid      = cell(1,gg1); 
A0prgrid(:,:) = {zeros(length(X),length(Y))};
A0srgrid      = A0prgrid; A0irgrid = A0prgrid;

Pumpout = Total_output;
Sigout  = Total_output;
Idout   = Total_output;
TotTot  = Total_output;

SigProfile = zeros(1,length(t));
DepProfile = zeros(1,length(t));
IdProfile  = zeros(1,length(t));

SignalField  = cell(1,length(t));
DepPumpField = cell(1,length(t));
IdlerField   = cell(1,length(t));

C1p = (-1i.*((2*(pi^2)/kp).*(sx.^2 + sy.^2)));
C1s = (-1i.*((2*(pi^2)/ks).*(sx.^2 + sy.^2)));
C1i = (-1i.*((2*(pi^2)/ki).*(sx.^2 + sy.^2)));

C2p = 2.*pi.*sy.*tan(rho_p);
C2s = 2.*pi.*sy.*tan(rho_s);
C2i = 2.*pi.*sy.*tan(rho_i);

for jj = 1:Noshots
    for kk = 0:RT-1 

        vv = 1;
        
            for nn = 1+kk*gg1:kk*gg1+gg1

                 A0pin{nn} = sqrt((1-AR1p)).*sqrt((1-R1p)).*(abs(Exyp{nn})).*exp(-1i.*angle(Exyp{nn})) + sqrt((1-AR1p)).*(abs(A0prgrid{vv})).*exp(+1i.*angle(A0prgrid{vv}));   
                 A0sin{nn} = sqrt((1-AR1s)).*sqrt((1-R1s)).*(abs(Exys{nn})).*exp(-1i.*angle(Exys{nn})) + sqrt((1-AR1s)).*(abs(A0srgrid{vv})).*exp(+1i.*angle(A0srgrid{vv})); 
                 A0iin{nn} = sqrt((1-AR1i)).*sqrt((1-R1i)).*(abs(Exyi{nn})).*exp(-1i.*angle(Exyi{nn})) + sqrt((1-AR1i)).*(abs(A0irgrid{vv})).*exp(+1i.*angle(A0irgrid{vv}));

                 vv = vv +1; 

            end        
        
                pp = 1;
                
            for dd = 1+kk*gg1:kk*gg1+gg1
                
                A0p =  A0pin{dd};    
                A0s =  A0sin{dd};   
                A0i =  A0iin{dd};   

                % Forward prop
                [As,Ai,Ap] = rk4(h,A0p,A0s,A0i,L,z,d,sx,sy,GridStepX,GridStepY,deltaK,c,deff,n_p,n_s,n_i,...
                                              vp,vs,vi,kp,ks,ki,C1p,C1s,C1i,C2p,C2s,C2i,X,Y);

                Theta2   {dd} = angle(Ap) - angle(As) - angle(Ai);
                pumpphase{dd} = angle(Ap);
                sigphase {dd} = angle(As);
                idphase  {dd} = angle(Ai);

                Signal_power_out{dd} = ((1-AR2s)).*(1-R2s).*((abs(As)).^2).*((GridStepX*GridStepY).*e0.*n_s.*c./2);
                Idler_power_out {dd} = ((1-AR2i)).*(1-R2i).*((abs(Ai)).^2).*((GridStepX*GridStepY).*e0.*n_i.*c./2);
                DepPump         {dd} = ((1-AR2p)).*(1-R2p).*((abs(Ap)).^2).*((GridStepX*GridStepY).*e0.*n_p.*c./2);
                
                Total_output{dd} = Signal_power_out{dd} + Idler_power_out{dd};
                
                SigProfile(dd) = sum(sum(Signal_power_out{dd}));
                DepProfile(dd) = sum(sum(DepPump{dd}));
                IdProfile (dd) = sum(sum(Idler_power_out{dd}));  
                
                SignalField {dd} = sqrt((1-AR2s)).*sqrt(1-R2s).*((abs(As).^1)).*exp(-1i.*angle(As)); % These used in M2 calcs.
                DepPumpField{dd} = sqrt((1-AR2p)).*sqrt(1-R2p).*((abs(Ap).^1)).*exp(-1i.*angle(Ap));
                IdlerField  {dd} = sqrt((1-AR2i)).*sqrt(1-R2i).*((abs(Ai).^1)).*exp(-1i.*angle(Ai));

                A0p = sqrt((1-AR2p)).*sqrt(R2p).*abs(Ap).*exp(-1i.*angle(Ap));        
                A0s = sqrt((1-AR2s)).*sqrt(R2s).*abs(As).*exp(-1i.*angle(As)); 
                A0i = sqrt((1-AR2i)).*sqrt(R2i).*abs(Ai).*exp(-1i.*angle(Ai)); 
                 
                d = d; 
                z = sort(z,'descend');

                % Backward Prop
                [As,Ai,Ap] = rk4(h,A0p,A0s,A0i,L,z,d,sx,sy,GridStepX,GridStepY,deltaK,c,deff,n_p,n_s,n_i,...
                                              vp,vs,vi,kp,ks,ki,C1p,C1s,C1i,C2p,C2s,C2i,X,Y);

                A0prgrid{pp} = sqrt((1-AR1p)).*sqrt(R1p).*abs(Ap).*exp(-1i.*angle(Ap));                             
                A0srgrid{pp} = sqrt((1-AR1s)).*sqrt(R1s).*abs(As).*exp(-1i.*angle(As));  
                A0irgrid{pp} = sqrt((1-AR1i)).*sqrt(R1i).*abs(Ai).*exp(-1i.*angle(Ai));  

                d = d;
                z = sort(z,'ascend');

                Theta1{dd} = angle(Ap) - angle(As) - angle(Ai);
                            
                pp = pp + 1;
                
            end  

        waitbar(kk/(RT-1),wait2);
        
    end
    
    OutputPulse = (Total_output); 

    waitbar(ii/length(Pumpvalues),wait1);

     for i = 1:length(X)
         for j = 1:length(Y)
             for bb = 1:length(t)  
         Track7  = Idler_power_out{bb};
         Track7b = DepPump{bb};
         Track7c = Signal_power_out{bb};
         Track8(bb)  = Track7(i,j);
         Track8b(bb) = Track7b(i,j);
         Track8c(bb) = Track7c(i,j);
             end
             CC(i,j)  = trapz(t,Track8);   % Idler
             CCP(i,j) = trapz(t,Track8b);  % DepPump
             CCS(i,j) = trapz(t,Track8c);  % Signal
         end
     end
   
   IdOut(ii)   = sum(sum(CC))*1000;
   ResPump(ii) = sum(sum(CCP))*1000;
   SigOut(ii)  = sum(sum(CCS))*1000;

   for bb = 1:length(t) 
       
     PHASE = Theta2{dd};
     PHASE1 = pumpphase{dd};
     PHASE2 = sigphase{dd};
     PHASE3 = idphase{dd};
     CenterPhase(bb) = PHASE(round(length(X)/2),round(length(Y)/2));
     PumpP(bb) = PHASE1(round(length(X)/2),round(length(Y)/2));
     SigP(bb) = PHASE2(round(length(X)/2),round(length(Y)/2));
     IdP(bb) = PHASE3(round(length(X)/2),round(length(Y)/2));
     
   end
end
end

   figure;
   plot(CenterPhase,'r'); hold on;
   plot(PumpP,'b'); 
   plot(SigP,'k');
   plot(IdP,'g');
   legend('thp-ths-thi','P','S','I');

   plot(Pumpvalues,MeasuredOut,'r+','Linewidth',2);
   hold on
   grid on
   plot(Pumpvalues,IdOut,'g+','Linewidth',2);
   hold on
   plot(Pumpvalues,SigOut,'k+','Linewidth',2);
   hold on
   plot(Pumpvalues,SigOut+IdOut,'b+','Linewidth',2);
   legend('Measured','Idler','Signal','Total');
   
   toc
   
   
figure; 
plot(Pumpvalues,IdOut,'g+','Linewidth',2);
hold on
% plot(Pumpvalues,MeasuredOutI,'go','Linewidth',2);
grid on
% hold on
% plot(Pumpvalues,SNLODiff,'ko','Linewidth',2);
% hold on
% plot(Pumpvalues,SNLOResPump,'bo','Linewidth',2);
hold on
plot(Pumpvalues,SigOut,'k+','Linewidth',2);
hold on
% plot(Pumpvalues,MeasuredOutS,'ko','Linewidth',2);
hold on

xlabel('Pump energy (mJ)');
ylabel('Total output (mJ)');  



 
     for i = 1:length(X)
         for j = 1:length(Y)
             for bb = 1:length(t)  
         Track7 = Total_output{bb};
         Track7b = DepPump{bb};
         Track8(bb) = Track7(i,j);
         Track8b(bb) = Track7b(i,j);
             end
%              plot(Track8b);
             CC(i,j) = trapz(t,Track8);
             CCP(i,j) = trapz(t,Track8b);
         end
     end
    
   OutOut = sum(sum(CC));
   ResPump = sum(sum(CCP));
     
toc
    
Output_envelope = OutputPulse;

figure;


plot(t.*10^9,abs(SigProfile)./max(abs(initial_pumpenvelope.^2)),'g','Linewidth',1.5)
hold on
plot(t.*10^9,abs(DepProfile)./max(abs(initial_pumpenvelope.^2)),'b-.','Linewidth',1.5) 
hold on
plot(t.*10^9,abs(initial_pumpenvelope).^2/max(abs(initial_pumpenvelope.^2)),'r-','Linewidth',1.5)
hold on
plot(t.*10^9,abs(IdProfile)./max(abs(initial_pumpenvelope.^2)),'k','Linewidth',1.5)
hold on
plot(t.*10^9,abs(SigProfile+IdProfile)./max(abs(initial_pumpenvelope.^2)),'m:','Linewidth',1.5)
xlabel('Time (ns)');
ylabel('Normalized Instantaneous Power (W)');
legend('Signal','DepPump','Pump','Idler')
grid on
axis([t(1).*10^9 t(end).*10^9 0 1.2.*10^9])

figure;
surf(X,Y,CCS);

figure;
surf(X,Y,CCP);
