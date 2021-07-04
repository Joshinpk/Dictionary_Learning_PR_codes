% Phase retrieval framework from V. Katkovnik,  April, 2016, TUT
function [RMSE_phi1,RMSE_ampl1,PHI_dlpr] = DLPR(varphi,amplitude,CDSFparam,ALGparam)

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
ISfiltering=ALGparam.ISfiltering;  
support_part=ALGparam.support_part;    
KAPPA_set=ALGparam.KAPPA_set;
SNR_Gaus_set=ALGparam.SNR_Gaus_set;
sig_num=ALGparam.signum;
ext_dic=CDSFparam.ext_dic;
x=amplitude.*exp(1j*varphi);     % TRUE COMPLEX-VALUED OBJECT
[yN,xN]=size(varphi); 
patsize = CDSFparam.patsize;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% GENERATION of RAMDOM PHASE MODULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Masks   = zeros(yN,xN,L);        %% Storage  for L masks, each of dim yN x xN
m       = yN*xN*L;                %% Data length
rand('seed',10001), randn('seed',1001);
for ll = 1:L, Masks(:,:,ll) = randsrc2(yN, xN,[1i -1i 1 -1]); end  %% As in TWF
% Sample phases: each symbol in alphabet {1, -1, i , -i} has equal prob.

%%%%%%%%%%%%%%%%%%%%%%%%%
% GEneration of Ri
%%%%%%%%%%%%%%%%%%%%%%%%%

Sum_Ri_H_Ri = RiGen(CDSFparam);



%% Linear operators for As 
% Forward propagation
A = @(I)  fft2(conj(Masks) .* reshape(repmat(I,[1 L]), size(I,1), size(I,2), L));  % Input is yN x xN image, output is yN x xN x L array
y = abs(A(x)).^2; % True intensity at the sensor plane
s_noise_level=0;

if Pois_Gaus_sel==1
    noise_set=KAPPA_set;
else
    noise_set=SNR_Gaus_set;
end

for NL=noise_set      %% Loop for KAPPA values
s_noise_level=s_noise_level+1;
%% Backward Propagation
if Pois_Gaus_sel==1
    KAPPA=NL;
    gamma_1=1/(KAPPA);         %% Parameter of the iterative algorithm
    Beta=0.0001/(gamma_1);
else
   KAPPA=NL;
   SNR=NL
   ssigma=sqrt((sum(y(:).^2)/length(y(:)))/(10^(0.1*SNR)));
   gamma_1=ssigma*ssigma/10;
   Beta=0.0001/(gamma_1);
%    gamma_2=1;
end

At = @(Y1,Y2,scale) (sum(Masks .* ifft2(Y1), 3)+Beta*scale*gamma_1*Y2)./(L+Beta*scale*gamma_1*Sum_Ri_H_Ri);                        %% Input is yN x xN x L array, output is yN x xN image

%% Scaled Poissonian data modeling z
rand('seed',100001), randn('seed',1001);
if Pois_Gaus_sel==1
    z  = poissrnd(y*KAPPA);                %% Poissonain Intensity Observations
else
    Gnoise=ssigma*randn(size(y));
    z=y+Gnoise;
end


%%%%%%%%%%%%%% SUPPORT SENSOR %%%%%%%%%%%
Masks_Support=zeros(size(varphi));
support_y=[max(-floor(yN/2*support_part)+yN/2,1):min(+floor(yN/2*support_part)+yN/2,yN)];
support_x=[max(-floor(xN/2*support_part)+xN/2,1):min(+floor(xN/2*support_part)+xN/2,yN)];

for ll = 1:L, Masks_Support(support_y,support_x,ll) = 1; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% SPAR algorithm %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic                 %% for DLPR
%% INITIALIZATION FOR DLPR
x_t=varphi;                                         %% for size only
x_t=ones(size(x_t)).* exp(1j*pi*randn(size(x_t))*0.1);  %% RANDOM PHASE INITIALIZATION
% x_t=ones(size(x_t)).* exp(1j*x_t);  %% RANDOM PHASE INITIALIZATION

x_init=x_t;
RMSE_ampl(s_noise_level,1)=norm(abs(x)-abs(x_init),'fro')/sqrt(yN*xN);  %% RMSE Initial ampl error
% % x_t=x_t.*Sum_Ri_H_Ri;
%% Phase unwrapping
if unwrapping==0
    potential.quantized = 'no';
    potential.threshold = pi;
    verbose='no';
    phi_Uo_init=puma_ho(angle(x_init),.5,'verbose',verbose); %% for demo unwrapping from raw data
    bb_init=mean(varphi(:) - phi_Uo_init(:));
    RMSE_phi(s_noise_level,1) = norm(varphi - phi_Uo_init-bb_init, 'fro')/sqrt(yN*xN); %% RMSE Initial rel. phase error
else
    V_cor_init=exp(-1i*angle(trace(x'*x_init))) * x_init;   %% PHASE CORRECTED INITIAL GUESS 
    bb1_init=-angle(trace(x'*x_init));
    phi_Uo_init=angle(x_init);
    Uo_corr_init=exp(-1i*angle(trace(x'*x_init))) * x_init;   %% CORRECTED INITIAL GUESS 
   %% RMSE_phi(s_noise_level,1) = norm(angle(x) - angle(Uo_corr_init), 'fro')/sqrt(yN*xN); % EMSE Initial rel. phase error
    RMSE_phi(s_noise_level,1) = norm(wrap(angle(x) - angle(Uo_corr_init)), 'fro')/sqrt(yN*xN); % EMSE Initial rel. phase error
end



%% MAIN ITERATIONS of SPAR PHASE RETRIEVAL
% gamma_10=gamma_1;
tic
for ss=1:IterNum
gsc=10/ss;%(IterNum+1-ss)/IterNum;
display(strcat('Sig No:',num2str(sig_num),', Kappa',num2str(KAPPA),', Main Iteration :',num2str(ss)))
%%            
ss0=ss;       %% ss0 is iteration when BM3D grouping is updated
%% STEP 1:     %% Forward propagation 
Vs=A(x_t);       

%% STEP 2:     %% Poissonian Filtering

if Fienup==1   %% No Filtering at Sensor Plane

    Us=ifftshift(sqrt((fftshift(z).*Masks_Support/KAPPA+(fftshift(abs(Vs)).^2).*(1-Masks_Support)))).*exp(j*angle(Vs));
else
    if Pois_Gaus_sel==1
        [Us]=NoiseFiltering_3D(Vs,z,Masks_Support,KAPPA,gsc*gamma_1,L);  %% Filtering at Sensor Plane for Poisson Observation
%         [Us]=NoiseFiltering_3D(Vs,z,Masks_Support,KAPPA,gamma_1,L);  %% Filtering at Sensor Plane for Poisson Observation

    else
        [Us]=GaussianNoiseFilterCardan(Vs,z,ssigma, Masks_Support,gsc*gamma_1,L);
    end
end

%% STEP 3: Backward propagation
% perc=sqrt(sum(abs(Us(:)-Vs(:)).^2))/sqrt(sum(abs(Vs(:)).^2))*100
% perc=norm(Us(:)-Vs(:))/norm(Vs(:))*100

Sum_Ri_H_D_alpha=x_t.*Sum_Ri_H_Ri;
x_t_half=At(Us,Sum_Ri_H_D_alpha,gsc);  
 
    

    %% STEP 4: Sparse Dictionery Learning and filtering
    %Dic Learning
    addpath(strcat(pwd,'/SpaRSAL'));
    if ss~=1
        CDSFparam.Dinitial=D;
    end
    if ext_dic==1
        if ss==1
            strD=strcat('Prior_Dic_Num_',num2str(sig_num),'.mat');
            load (strD); 
        end
    else
        clear D
        D= CD_Dic_Learning(x_t_half,CDSFparam);
    end
    %Sparse Filtering
    CDSFparam.Omp_ctrl_flag=0;  %make it 0
    CDSFparam.main_itern=ss;
    CDSFparam.KAPPA=KAPPA;
    %OMP
    x_t= CD_Sparse_Filtering(x_t_half,D,CDSFparam);


phi_u0_SPAR=angle(x_t);
V_cor=exp(-1i*(angle(trace(x'*x_t)))) * x_t;  %% Phase correction by an invariant shift according to x 

if unwrapping==1
bb=mean(varphi(:) - phi_u0_SPAR(:));                                       %% absolute phase shift
% RMSE_phi(s_noise_level,ss+1) = norm(varphi - phi_u0_SPAR-bb, 'fro')/sqrt(yN*xN); %% RMSE absolute phase error
RMSE_phi(s_noise_level,ss+1) = norm(wrap(angle(x) - (angle(V_cor))), 'fro')/sqrt(yN*xN);
% jpk(s_noise_level,ss+1)= norm(wrap(angle(x) - (angle(x_t))), 'fro')/sqrt(yN*xN)
else    
RMSE_phi(s_noise_level,ss+1) = norm(wrap(angle(x) - (angle(V_cor))), 'fro')/sqrt(yN*xN);
%% RMSE phase error
end
RMSE_phi
RMSE_ampl(s_noise_level,ss+1)=norm(abs(x)-abs(x_t),'fro')/sqrt(yN*xN);                      %% RMSE amplitude errorRMSE_phi
end
DLPR_time(s_noise_level)=toc;    %% for DLPR


end%% KAPPA_set
RMSE_phi1=RMSE_phi(:,end);
RMSE_ampl1=RMSE_ampl(:,end);
PHI_dlpr=angle(V_cor);
end
