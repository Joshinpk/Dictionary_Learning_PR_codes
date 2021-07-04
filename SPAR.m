%V. Katkovnik,  April, 2016, TUT
function [RMSE_phi_SPAR,RMSE_ampl_SPAR,RMSE_ampl_TWF_set,...
    RMSE_phi_TWF_set,RMSE_phi_GS1,RMSE_ampl_GS1] = SPAR(varphi,Bo,BM3Dparam,ALGparam,PHI_dlpr,RMSE_phi1)
% addpath .\DATA
%addpath .\PUMA
J=sqrt(-1);
%==================================================================================================================
% Parameter reading
%==================================================================================================================
L=ALGparam.L;  
Fienup=ALGparam.Fienup; 
Pois_Gaus_sel=ALGparam.Pois_Gaus_sel;
IterNum=ALGparam.IterNum;        
unwrapping=ALGparam.unwrapping; 
filtering=ALGparam.ISfiltering;  
support_part=ALGparam.support_part;    
KAPPA_set=ALGparam.KAPPA_set;
SNR_Gaus_set=ALGparam.SNR_Gaus_set;
sig_num=ALGparam.signum;
TWF=ALGparam.TWF;


Nstep=BM3Dparam.Nstep;  
N11=BM3Dparam.N11;  
N22=BM3Dparam.N22; 
threshType=BM3Dparam.threshType; 
threshold_ampl=BM3Dparam.threshold_ampl;
threshold_phi=BM3Dparam.threshold_phi;


x=Bo.*exp(1j*varphi);     % TRUE COMPLEX-VALUED OBJECT
[yN,xN]=size(varphi); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% OBSERVATION MODEL and OPERATORS %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Masks   = zeros(yN,xN,L);        %% Storage for L masks, each of dim yN x xN
m       = yN*xN*L;                %% Data length

rand('seed',10001), randn('seed',1001);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GENERATION of RAMDOM PHASE MODULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 1
    %% Four random values for Phase Mask
for ll = 1:L, Masks(:,:,ll) = randsrc2(yN, xN,[1i -1i 1 -1]); end  %% As in TWF
% Sample phases: each symbol in alphabet {1, -1, i , -i} has equal prob.
else
    
    %% Gaussian Phase Mask
for ll = 1:L, Masks(:,:,ll) = randn(yN, xN); end
Masks=exp(1j*Masks*pi);
end
 
if 0                          %% small magnitude integer random mask
maxim=5;
for ll = 1:L, Masks1(:,:,ll) = randi((2*maxim)+1,[yN xN])-(1+maxim); end

Masks=exp(1j*Masks1*2*pi/(2*maxim+1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Linear operators for As 
%% Forward propagation
A = @(I)  fft2(conj(Masks) .* reshape(repmat(I,[1 L]), size(I,1), size(I,2), L));  % Input is yN x xN image, output is yN x xN x L array

%% Backward propagation
At = @(Y) sum(Masks .* ifft2(Y), 3)/L;                        %% Input is yN x xN x L array, output is yN x xN image

At0 = @(Y) sum(Masks .* ifft2(Y), 3) * size(Y,1) * size(Y,2); %% Backward for initialization only

y = abs(A(x)).^2;                                             %% True intensity at the sensor plane
%%%%

s_noise_level=0;
if Pois_Gaus_sel==1
    noise_set=KAPPA_set;
else
    noise_set=SNR_Gaus_set;
end

for NL=noise_set      %% Loop for KAPPA values
s_noise_level=s_noise_level+1;
% gamma_1=1/KAPPA;         %% Parameter of the iterative algorithm
if Pois_Gaus_sel==1
    KAPPA=NL;
    gamma_1=1/KAPPA;         %% Parameter of the iterative algorithm
    gamma_2=gamma_1;
else
   KAPPA=NL;
   SNR=NL;
   ssigma=sqrt((sum(y(:).^2)/length(y(:)))/(10^(0.1*SNR)));
   gamma_1=ssigma*ssigma/10;
   gamma_2=1;
end


% SNR=0;
% sigma= sqrt((sum(y(:).^2)/length(y(:)))/(10^(0.1*SNR)));
% gamma_1=sigma*sigma/10;

%% Scaled Poissonian data modeling z
% sigma=0.9;
rand('seed',100001), randn('seed',1001);
if Pois_Gaus_sel==1
    z  = poissrnd(y*KAPPA);                %% Poissonain Intensity Observations
else
    Gnoise=ssigma*randn(size(y));
    z=y+Gnoise;
end

% SNR=40;
% sigma=sqrt((sum(y(:).^2)/length(y(:)))/(10^(0.1*SNR)));
% Gnoise=sigma*randn(size(y));
% z=y+Gnoise;


% z  = poissrnd(y*KAPPA);                %% Poissonain Intensity Observations

Number_Photons_per_Pixel(s_noise_level)=mean(z(:));                      %% Mean umber of Photons per Pixel
SNR(s_noise_level)=10*log10(KAPPA^2*sum(y(:).^2)/sum(y(:)*KAPPA-z(:)).^2); %% SNR for Poisson Data
%%%%%%%%%%%%%% SUPPORT SENSOR %%%%%%%%%%%
Masks_Support=zeros(size(varphi));
support_y=[max(-floor(yN/2*support_part)+yN/2,1):min(+floor(yN/2*support_part)+yN/2,yN)];
support_x=[max(-floor(xN/2*support_part)+xN/2,1):min(+floor(xN/2*support_part)+xN/2,yN)];

for ll = 1:L, Masks_Support(support_y,support_x,ll) = 1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% TWF algorithm %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if TWF==1
    
    %% TWF algorithm from 
    %% 'Solving Quadratic Systems of Equations via Truncated Wirtinger Flow,' Yuxin Chen and Emmanuel J. Cand?s, 2015
    %% codes available at http://web.stanford.edu/~yxchen/TWF/

%% INITIALIZATION of ALGORITHM
if exist('Params')                == 0,  Params.n1          = 100;  end
if isfield(Params, 'n2')          == 0,  Params.n2          = 100;  end             % signal dimension
if isfield(Params, 'L')           == 0,  Params.L           = 12;   end             % number of measurements
if isfield(Params, 'grad_type')   == 0,  Params.grad_type   = 'TWF_Poiss';  end     % 'TWF_Poiss': Poisson likelihood

if isfield(Params, 'alpha_lb')    == 0,  Params.alpha_lb    = 0.3;  end
if isfield(Params, 'alpha_ub')    == 0,  Params.alpha_ub    = 5;    end
if isfield(Params, 'alpha_h')     == 0,  Params.alpha_h     = 5;    end
if isfield(Params, 'alpha_y')     == 0,  Params.alpha_y     = 3;    end 
if isfield(Params, 'T')           == 0,  Params.T           = 100/2; end    	% number of iterations
if isfield(Params, 'mu')          == 0,  Params.mu          = 0.2;  end		% step size / learning parameter
if isfield(Params, 'npower_iter') == 0,  Params.npower_iter = 50;   end		% number of power iterations

Params.n1=yN;
Params.n2=xN; Params.L=L; Params.m=m;

tic       %% for TWF

[Relerrs,RMSE_phase_TWF, RMSE_ampl_TWF,x_hat,x_hat_0,bb1_TWF] = TWF_my(z.*Masks_Support,x, Params, A, At0,KAPPA); %% Truncated WF algorithm
RMSE_phi_TWF_set(s_noise_level)=RMSE_phase_TWF(end);
RMSE_ampl_TWF_set(s_noise_level)=RMSE_ampl_TWF(end);


time_TWF(s_noise_level)=toc;  %% for TWF


if 0
fprintf('Relative error after initialization: %f\n', Relerrs(1))
fprintf('Relative error after %d iterations: %f\n', T, Relerrs(T+1))
 
figure(1), semilogy(0:Params.T,Relerrs) 
xlabel('Iteration'), ylabel('Relative error (log10)'), ...
     title('Relative error vs. iteration count'), grid on
 
 figure(2), plot(0:Params.T,RMSE_phase_TWF) 
xlabel('Iteration'), ylabel('Relative RMSE'), ...
     title('Phase RMSE vs. iteration count'), grid on
end
%% Visualization for TWF

if 0
figure(1), plot(0:Params.T,Relerrs),... %% semilogy(0:Params.T,Relerrs) 
xlabel('Iteration'), ylabel('Relative Error'),grid on, ...
title('Relative error vs. iteration count')
figure(2), subplot(2,2,1), imshow(angle(x),[]),title('phase-x'),...
subplot(2,2,2), imshow(angle(x_hat),[]),title('phase-x-TWF'),
subplot(2,2,3), imshow(abs((x)),[]),title('ampl-x'),...
subplot(2,2,4), imshow(abs(x_hat),[]),title('ampl-x-TWF'),

end

end   %% end of TWF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% SPAR algorithm %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic                 %% for SPAR
%% INITIALIZATION FOR SPAR
V=varphi;                                         %% for size only
V=ones(size(V)).* exp(1j*pi*randn(size(V))*0.1);  %% RANDOM PHASE INITIALIZATION
V_init=V;
RMSE_ampl(s_noise_level,1)=norm(abs(x)-abs(V_init),'fro')/sqrt(yN*xN);  %% RMSE Initial ampl error

%% Phase unwrapping
if unwrapping==0
potential.quantized = 'no';
potential.threshold = pi;
verbose='no';
phi_Uo_init=puma_ho(angle(V_init),.5,'verbose',verbose); %% for demo unwrapping from raw data
bb_init=mean(varphi(:) - phi_Uo_init(:));
RMSE_phi(s_noise_level,1) = norm(varphi - phi_Uo_init-bb_init, 'fro')/sqrt(yN*xN); %% RMSE Initial rel. phase error

else
    V_cor_init=exp(-1i*angle(trace(x'*V_init))) * V_init;   %% PHASE CORRECTED INITIAL GUESS 
    bb1_init=-angle(trace(x'*V_init));
    phi_Uo_init=angle(V_init);
    Uo_corr_init=exp(-1i*angle(trace(x'*V_init))) * V_init;   %% CORRECTED INITIAL GUESS 
   %% RMSE_phi(s_noise_level,1) = norm(angle(x) - angle(Uo_corr_init), 'fro')/sqrt(yN*xN); % EMSE Initial rel. phase error
    RMSE_phi(s_noise_level,1) = norm(wrap(angle(x) - angle(Uo_corr_init)), 'fro')/sqrt(yN*xN); % EMSE Initial rel. phase error
end

%% RMSE of the phase and amplitude estimates for use in BM3D filtered estimates
sigma_phi(1)=function_stdEst2D(phi_Uo_init,2);
sigma_ampl(1)=function_stdEst2D(abs(V),2);

%% MAIN ITERATIONS of SPAR PHASE RETRIEVAL
gamma_10=gamma_1;
for ss=1:IterNum
% display(strcat('Sig No:',num2str(sig_num),', Kappa',num2str(KAPPA),', Main Iteration :',num2str(ss)))
ss0=ss;       %% ss0 is iteration when BM3D grouping is updated
%% STEP 1:     %% Forward propagation 
Us1=A(V);       

%% STEP 2:     %% Poissonian Filtering
if Fienup==1   %% No Filtering at Sensor Plane
    Us=ifftshift(sqrt((fftshift(z).*Masks_Support/KAPPA+(fftshift(abs(Us1)).^2).*(1-Masks_Support)))).*exp(j*angle(Us1));
else
%     [Us]=NoiseFiltering_3D(Us,z,Masks_Support,KAPPA,gamma_1,L);  %% Filtering at Sensor Plane
%     [Us]=GaussianNoiseFilterCardan(Us,z,ssigma, Masks_Support,gamma_1,L);
    if Pois_Gaus_sel==1
        [Us]=NoiseFiltering_3D(Us1,z,Masks_Support,KAPPA,gamma_1,L);  %% Filtering at Sensor Plane for Poisson Observation
    else
        [Us]=GaussianNoiseFilterCardan(Us1,z,ssigma, Masks_Support,gamma_1,L);
    end
end
% perc=sqrt(sum(abs(Us(:)-Us1(:)).^2))/sqrt(sum(abs(Us1(:)).^2))*100


%% STEP 3: Backward propagation

Uo=At(Us);  

%% STEP 4: Phase unwrapping ;
if unwrapping==1
potential.quantized = 'no';
potential.threshold = pi;
verbose='no';
phi_Uo=puma_ho(angle(Uo),.5,'verbose',verbose); %% for demo unwrapping from raw data

else
    
phi_Uo=double(angle(Uo));
end   
    
%% STEP 5: BM3D phase and amplitude filtering
%% BM3D filtering of phase

sigma_phi(ss+1)=function_stdEst2D(phi_Uo,2);  %% STD of phase
phi_u0_SPAR=BM3D_SPAR_UNWRAP_PHIp(phi_Uo,threshType, sigma_phi(ss+1)*threshold_phi, N11, N22, Nstep,filtering,ss,ss0);


%% BM3D filtering of amplitude           
sigma_ampl(ss+1)=function_stdEst2D(abs(Uo),2);  %% STD of amplitude

B_u0_SPAR=BM3D_SPAR_ABSp(abs(Uo),threshType, sigma_ampl(ss+1)*threshold_ampl, N11, N22, Nstep,filtering,ss,ss0);

 %% STEP 6: UPDATE of x

V=B_u0_SPAR.*exp(1j*phi_u0_SPAR);

V_cor=exp(-1i*(angle(trace(x'*V)))) * V;  %% Phase correction by an invariant shift according to x

if unwrapping==1
bb=mean(varphi(:) - phi_u0_SPAR(:));                                       %% absolute phase shift
% RMSE_phi(s_noise_level,ss+1) = norm(varphi - phi_u0_SPAR-bb, 'fro')/sqrt(yN*xN); %% RMSE absolute phase error
RMSE_phi(s_noise_level,ss+1) = norm(wrap(angle(x) - (angle(V_cor))), 'fro')/sqrt(yN*xN);
else
    
RMSE_phi(s_noise_level,ss+1) = norm(wrap(angle(x) - (angle(V_cor))), 'fro')/sqrt(yN*xN); %% RMSE phase error
end

RMSE_ampl(s_noise_level,ss+1)=norm(abs(x)-abs(V),'fro')/sqrt(yN*xN);                      %% RMSE amplitude error

% figure(3), subplot(2,2,1), imshow(angle(V_cor),[]), title(['IterNum=', num2str(ss+1),'; RMSE-phase=', num2str(RMSE_phi(ss+1),3)]),...
%     subplot(2,2,2), imshow(abs(V_cor),[]), title(['IterNum=', num2str(ss+1),'; RMSE-ampl=', num2str(RMSE_ampl(ss+1),3)])
%     subplot(2,2,3), imshow(angle(V_cor_init),[]), title(['IterNum=', num2str(0),'; RMSE-phase=', num2str(RMSE_phi(1),3)]),...
%     subplot(2,2,4), imshow(abs(V_cor_init),[]), title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(1),3)])
% 
if 0   %% ITERATION SHOW
if unwrapping==1
figure(3), subplot(2,2,1), imshow(phi_Uo+bb,[]), title(['IterNum=', num2str(ss+1),'; RMSE-abs.phase=', num2str(RMSE_phi(ss+1),3)]),...
    subplot(2,2,2), imshow(abs(V),[]), title(['IterNum=', num2str(ss+1),'; RMSE-ampl=', num2str(RMSE_ampl(ss+1),3)])
    subplot(2,2,3), imshow(phi_Uo_init+bb_init,[]), title(['IterNum=', num2str(0),'; RMSE-abs.phase=', num2str(RMSE_phi(1),3)]),...
    subplot(2,2,4), imshow(abs(V_init),[]), title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(1),3)])
    else
 figure(3), subplot(2,2,1), imshow(angle(V_cor),[]), title(['IterNum=', num2str(ss+1),'; RMSE-abs.phase=', num2str(RMSE_phi(ss+1),3)]),...
    subplot(2,2,2), imshow(abs(V_cor),[]), title(['IterNum=', num2str(ss+1),'; RMSE-ampl=', num2str(RMSE_ampl(ss+1),3)])
    subplot(2,2,3), imshow(angle(V_cor_init),[]), title(['IterNum=', num2str(0),'; RMSE-abs.phase=', num2str(RMSE_phi(1),3)]),...
    subplot(2,2,4), imshow(abs(V_cor_init),[]), title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(1),3)])
    end
end
%%[ss RMSE_phi(s_noise_level,ss) RMSE_ampl(s_noise_level,ss)]


end
SPAR_time(s_noise_level)=toc ;   %% for SPAR

if 0
if unwrapping==1
figure(3), subplot(2,2,1), imshow(phi_Uo+bb,[]), title(['IterNum=', num2str(ss+1),'; RMSE-abs.phase=', num2str(RMSE_phi(ss+1),3)]),...
    subplot(2,2,2), imshow(abs(V),[]), title(['IterNum=', num2str(ss+1),'; RMSE-ampl=', num2str(RMSE_ampl(ss+1),3)])
    subplot(2,2,3), imshow(phi_Uo_init+bb_init,[]), title(['IterNum=', num2str(0),'; RMSE-abs.phase=', num2str(RMSE_phi(1),3)]),...
    subplot(2,2,4), imshow(abs(V_init),[]), title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(1),3)])
    else
 figure(3), subplot(2,2,1), imshow(angle(V_cor),[]), title(['IterNum=', num2str(ss+1),'; RMSE-abs.phase=', num2str(RMSE_phi(ss+1),3)]),...
    subplot(2,2,2), imshow(abs(V_cor),[]), title(['IterNum=', num2str(ss+1),'; RMSE-ampl=', num2str(RMSE_ampl(ss+1),3)])
    subplot(2,2,3), imshow(angle(V_cor_init),[]), title(['IterNum=', num2str(0),'; RMSE-abs.phase=', num2str(RMSE_phi(1),3)]),...
    subplot(2,2,4), imshow(abs(V_cor_init),[]), title(['IterNum=', num2str(0),'; RMSE-ampl=', num2str(RMSE_ampl(1),3)])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% GS ALGORITHM %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic                         %% for GS
PHASE_RETRIEVAL_CODED_FOURIER_GS_ALGORITHM
time_GS(s_noise_level)=toc;        %% for GS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Images for results %%%%%%
if 0
if TWF==1
figure,
    colormap gray; fontsize = 12;
    subplot (2,2,2), surfl( puma_ho(angle(V_cor_GS),.5,'verbose',verbose)),shading interp, title(['\fontsize{20} GS-F, RMSE =', num2str(RMSE_phi_GS(s_noise_level,end),2)]),...
axis tight
subplot(2,2,1), surfl(puma_ho(angle(x_hat),.5,'verbose',verbose)),shading interp, title(['\fontsize{20} TWF, RMSE = ', num2str(RMSE_phi_TWF_set(end),2)]),...
axis tight    
subplot(2,2,3), surfl(puma_ho(angle(V_cor),.5,'verbose',verbose)),shading interp, title(['\fontsize{20} SPAR, RMSE = ', num2str(RMSE_phi(s_noise_level,end),2)]) 
axis tight
subplot(2,2,4), surfl(puma_ho(PHI_dlpr,.5,'verbose',verbose)), ,shading interp, title(['\fontsize{20} DLPR, RMSE = ', num2str(RMSE_phi1(s_noise_level,end),2)]) 
figure
    colormap gray; fontsize = 12;
    subplot (2,2,2), imagesc(angle(V_cor_GS)),shading interp, title(['GS-F, RMSE =', num2str(RMSE_phi_GS(s_noise_level,end),2)]),...
    subplot(2,2,1), imagesc(angle(x_hat)),shading interp, title(['TWF, RMSE = ', num2str(RMSE_phi_TWF_set(end),2)]),...
    subplot(2,2,3), imagesc(angle(V_cor)),shading interp, title(['SPAR, RMSE = ', num2str(RMSE_phi(s_noise_level,end),2)]) 
    subplot(2,2,4), imagesc(wrap(PHI_dlpr)), title(['DLPR, RMSE = ', num2str(RMSE_phi1(s_noise_level,end),2)]) ,shading interp
    
    
% figure, subplot (1,3,2), imshow(angle(V_cor_GS),[]), title(['GS, RMSE =', num2str(RMSE_phi_GS(s_noise_level,end),2)]),...
% subplot(1,3,1), imshow(angle(x_hat),[]), title(['TWF, RMSE = ', num2str(RMSE_phi_TWF_set(end),2)]),...
% subplot(1,3,3), imshow(wrap(phi_u0_SPAR),[]), title(['SPAR, RMSE = ', num2str(RMSE_phi(s_noise_level,end),2)]) 
% subplot(1,4,4), imshow(wrap(varphi),[]), title('TRUE') 
% pause(.2)

if 0
figure, subplot (1,3,2), imshow(abs(V_GS),[]), title(['GS AMPL, RMSE = ', num2str(RMSE_ampl_GS(s_noise_level,end),2)]),...
    subplot(1,3,1), imshow(abs(x_hat),[]), title(['TWF AMPL,  \chi =  ', num2str(KAPPA,2),', RMSE = ', num2str(RMSE_ampl_TWF_set(end),2)]),...
    subplot(1,3,3), imshow(abs(V),[]), title(['SPAR AMPL, RMSE = ', num2str(RMSE_ampl(s_noise_level,end),2)]) 
% pause(.2)
end
else
    figure, subplot (1,2,1), imshow(angle(V_cor_GS),[]), title(['GS PHASE, \chi =  ', num2str(KAPPA,2),', RMSE =', num2str(RMSE_phi_GS(s_noise_level,end),2)]),...
   %% subplot(1,3,1), imshow(angle(x_hat),[]), title(['TWF PHASE, \chi =  ', num2str(KAPPA,2),', RMSE = ', num2str(RMSE_phi_TWF_set(end),2)]),...
    subplot(1,2,2), imshow(angle(V_cor),[]), title(['SPAR PHASE, RMSE = ', num2str(RMSE_phi(s_noise_level,end),2)]) 
% pause(.2)
end
if 0
figure, subplot (1,2,1), imshow(abs(V_GS),[]), title(['GS AMPL,  \chi =  ', num2str(KAPPA,2),', RMSE = ', num2str(RMSE_ampl_GS(s_noise_level,end),2)]),...
 %%   subplot(1,3,1), imshow(abs(x_hat),[]), title(['TWF AMPL, \chi =  ', num2str(KAPPA,2),', RMSE = ', num2str(RMSE_ampl_TWF_set(end),2)]),...
        subplot(1,2,2), imshow(abs(V),[]), title(['SPAR AMPL, RMSE = ', num2str(RMSE_ampl(s_noise_level,end),2)]) 
% pause(.2)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
% PLOTS_Combined
end  %% KAPPA_set
RMSE_phi_SPAR=RMSE_phi(:,end);
RMSE_ampl_SPAR=RMSE_ampl(:,end);
RMSE_phi_GS1=RMSE_phi_GS(:,end);
RMSE_ampl_GS1=RMSE_ampl_GS(:,end);

% RMSE_BM3D=RMSE_phi';
% savefile = sprintf('RMSE_BM3D_K01_L12.mat');
% save(savefile,'RMSE_BM3D');
end
%% Time_total=toc