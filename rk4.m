
function[As,Ai,Ap] = rk4(h,A0p,A0s,A0i,L,z,d,sx,sy,GridStepX,GridStepY,deltaK,c,deff,n_p,n_s,n_i,vp,vs,vi,kp,ks,ki,C1p,C1s,C1i,C2p,C2s,C2i,X,Y);

N = round(abs(L/h));
Ap = (A0p);  
As = (A0s);
Ai = (A0i);

for jj = 1:N+1 

        
        MixP = (1i.*deff.*vp./(c.*n_p)).*Ai.*As.*exp(-1i.*deltaK.*z(jj));
        MixS = (1i.*deff.*vs./(c.*n_s)).*Ap.*conj(Ai).*exp(1i.*deltaK.*z(jj));
        MixI = (1i.*deff.*vi./(c.*n_i)).*Ap.*conj(As).*exp(1i.*deltaK.*z(jj));

        FMixP = fftshift(fft2(ifftshift(MixP))); %.*GridStep^2;
        FMixS = fftshift(fft2(ifftshift(MixS))); %.*GridStep^2;
        FMixI = fftshift(fft2(ifftshift(MixI))); %.*GridStep^2;   
        
        FAp = fftshift(fft2(ifftshift(Ap))); %.*GridStep^2;         
        FAs = fftshift(fft2(ifftshift(As))); %.*GridStep^2;
        FAi = fftshift(fft2(ifftshift(Ai))); %.*GridStep^2;    
        
        [k1p,k1s,k1i] = RHSDiff(h,FAp,FAs,FAi,FMixP,FMixS,FMixI,z(jj),d,C1p,C1s,C1i,C2p,C2s,C2i); 
        
        MixP = (1i.*deff.*vp./(c.*n_p)).*Ai.*As.*exp(-1i.*deltaK.*(z(jj)+h/2));
        MixS = (1i.*deff.*vs./(c.*n_s)).*Ap.*conj(Ai).*exp(1i.*deltaK.*(z(jj)+h/2));
        MixI = (1i.*deff.*vi./(c.*n_i)).*Ap.*conj(As).*exp(1i.*deltaK.*(z(jj)+h/2));

        FMixP = fftshift(fft2(ifftshift(MixP))); %.*GridStep^2;
        FMixS = fftshift(fft2(ifftshift(MixS))); %.*GridStep^2;
        FMixI = fftshift(fft2(ifftshift(MixI))); %.*GridStep^2;  
              
        [k2p,k2s,k2i] = RHSDiff(h,FAp+k1p/2,FAs+k1s/2,FAi+k1i/2,FMixP,FMixS,FMixI,z(jj)+h/2,d,C1p,C1s,C1i,C2p,C2s,C2i);  
        
        [k3p,k3s,k3i] = RHSDiff(h,FAp+k2p/2,FAs+k2s/2,FAi+k2i/2,FMixP,FMixS,FMixI,z(jj)+h/2,d,C1p,C1s,C1i,C2p,C2s,C2i);

        MixP = (1i.*deff.*vp./(c.*n_p)).*Ai.*As.*exp(-1i.*deltaK.*(z(jj)+h));
        MixS = (1i.*deff.*vs./(c.*n_s)).*Ap.*conj(Ai).*exp(1i.*deltaK.*(z(jj)+h));
        MixI = (1i.*deff.*vi./(c.*n_i)).*Ap.*conj(As).*exp(1i.*deltaK.*(z(jj)+h));

        FMixP = fftshift(fft2(ifftshift(MixP))); %.*GridStep^2;
        FMixS = fftshift(fft2(ifftshift(MixS))); %.*GridStep^2;
        FMixI = fftshift(fft2(ifftshift(MixI))); %.*GridStep^2;  
        
        [k4p,k4s,k4i] = RHSDiff(h,FAp+k3p,FAs+k3s,FAi+k3i,FMixP,FMixS,FMixI,z(jj)+h,d,C1p,C1s,C1i,C2p,C2s,C2i);
        
        FAp = FAp + (k1p + 2*k2p +2*k3p + k4p)./6;
        FAs = FAs + (k1s + 2*k2s +2*k3s + k4s)./6;
        FAi = FAi + (k1i + 2*k2i +2*k3i + k4i)./6;
        
        Ap = fftshift(ifft2(ifftshift(FAp))); %.*(1/(GridStep^2));
        As = fftshift(ifft2(ifftshift(FAs))); %.*(1/(GridStep^2));
        Ai = fftshift(ifft2(ifftshift(FAi))); %.*(1/(GridStep^2));

end


%         [k1p,k1s,k1i] = RHSDiff(h,FAp,FAs,FAi,FMixP,FMixS,FMixI,z,sx,sy,d);
%         
%         [k2p,k2s,k2i] = RHSDiff(h,FAp+k1p/2,FAs+k1s/2,FAi+k1i/2,FMixP+k1p/2,FMixS+k1s/2,FMixI+k1i/2,z+h/2,sx,sy,d);
%         
%         [k3p,k3s,k3i] = RHSDiff(h,FAp+k2p/2,FAs+k2s/2,FAi+k2i/2,FMixP+k2p/2,FMixS+k2s/2,FMixI+k2i/2,z+h/2,sx,sy,d);
%         
%         [k4p,k4s,k4i] = RHSDiff(h,FAp+k3p,FAs+k3s,FAi+k3i,FMixP+k3p,FMixS+k3s,FMixI+k3i,z+h,sx,sy,d);

% Algorithm for Diff model:

% (1) Start with Ap,As,Ai -> Calculate Mixing term(s)-> MixP, MixS, MixI.
% (2) FT(MixP,MixS,MixI) -> FMixP, FMixS, FmixI; FT(Ap,As,Ai) -> FAp,FAs,FAi
% (3) Use FMixP, FMixS, FmixI & FAp,FAs,FAi in RK4 to calculate FAp, etc at
% next step.
% (4) ifft(FAp(z+h)) -> Ap(z+h), etc.



% Old Rk4 for just mixing terms, i.e. no diffraction terms etc.

% N = floor(abs(L/h));
% 
% As = zeros(length(A0p),N);      % 2D matrix, rows -> Pulse values, coloumns -> N position.
% Ai = zeros(length(A0p),N);
% Ap = zeros(length(A0p),N);
% 
% Ap(:,1) = transpose(A0p);  
% As(:,1) = transpose(A0s);
% Ai(:,1) = transpose(A0i);
% 
% for j = 1:N+1
%         
%     [k1p,k1s,k1i] = RHS(h,As(:,j),Ap(:,j),Ai(:,j),z(j),d);  
%     
%     [k2p,k2s,k2i] = RHS(h,As(:,j)+k1s/2,Ap(:,j)+k1p/2,Ai(:,j)+k1i/2,z(j)+h/2,d);
%     
%     [k3p,k3s,k3i] = RHS(h,As(:,j)+k2s/2,Ap(:,j)+k2p/2,Ai(:,j)+k2i/2,z(j)+h/2,d);
%     
%     [k4p,k4s,k4i] = RHS(h,As(:,j)+k3s,Ap(:,j)+k3p,Ai(:,j)+k3i,z(j)+h,d);
%  
%     As(:,j+1) = As(:,j) + (k1s + 2*k2s +2*k3s + k4s)/6;
%     Ai(:,j+1) = Ai(:,j) + (k1i + 2*k2i +2*k3i + k4i)/6;
%     Ap(:,j+1) = Ap(:,j) + (k1p + 2*k2p +2*k3p + k4p)/6;   
%     
% end
% 

%% Note about fftshift etc:

        %  y=exp(-t.^2);
        %  figure;plot(abs(fftshift(fft(y))))
        %  figure;plot(phase(fftshift(fft(y))))
        %  figure;plot(phase(fftshift(fft(ifftshift(y)))))
        
