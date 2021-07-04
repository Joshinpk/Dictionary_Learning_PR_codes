clc;close all;clear all;
I=sqrt(-1);

addpath ./SpaRSAL
addpath ./PUMA

%================================
%Test Signal Generation 
%Group-1: sig-1,2  || Group-2: sig-3,4,5   ||  Group-3: sig-6,7  || Group-4: sig-8,9
% Please refer to the paper https://www.mdpi.com/1424-8220/18/11/4006 for
% the details of the test case groups
%================================
SIG_set=[1  ] %enter the test case numbers. Eg: SIG_set=[3]. Multiple test cases are allowed, Eg: SIG_set=[3 4 5]

%================================
%Observation model selection
%================================
ALGparam.Pois_Gaus_sel=0;    % 1 for Poissonian observation and 2 for Gaussian Observation

%================================
%enter the Noise level
%================================
%for Poissonian noise
ALGparam.KAPPA_set=[.00001  ]; %Eg: ALGparam.KAPPA_set=[.00001, .0001, .001]
%for Gaussian noise
ALGparam.SNR_Gaus_set=[1]; % Eg: ALGparam.SNR_Gaus_set=[1 3 7 10];

%====================================================================================================================
% Parameters
%====================================================================================================================
CDSFparam.ext_dic=0;
% CD_Sparse_Filtering
CDSFparam.K=256;%256
CDSFparam.T=30;
CDSFparam.patchnumDic=64;
CDSFparam.lamdaDic=0.11;
CDSFparam.err=1e-3;%
CDSFparam.IsSubMean=1;
CDSFparam.iterndisplay=0;
CDSFparam.patsize=10; 
CDSFparam.scale=1.5;
%BM3D Filtering
BM3Dparam.Nstep = 3;  
BM3Dparam.N11=8;  
BM3Dparam.N22=8; 
BM3Dparam.threshType='h'; 
BM3Dparam.threshold_ampl=1.4*4;
BM3Dparam.threshold_phi=1.4*4;
% Algorithm Parameters
ALGparam.IterNum=20;         % Iteration number 
ALGparam.unwrapping=1;       % 1 -means unwrapping on before BM3D phase filtering, 0- means unwrapping off
ALGparam.Fienup=0;           % 1 means that the Fienup rules is used instead of the optimal solution for data denoising
ALGparam.ISfiltering=1;      % 1 for BM3D filtering
ALGparam.TWF=1;              % Truncated Wirtinger Flow, TWF=1, acitivation, 0 -deactivation
ALGparam.support_part=1;     % Percentage of active pixels is support_part^2. support_part=.5 gives 25% 
ALGparam.L = 12;             % Number of Experiments 
ALGparam.SIG_set=SIG_set;
%====================================================================================================================
% Dictionery Learning Phase Retrieval
%====================================================================================================================
if 1
SIG_set=ALGparam.SIG_set;
for sig_indx =1: length(SIG_set)
    sig_num=SIG_set(sig_indx);
    ALGparam.signum=sig_num;
    [phase,amplitude] = sig_generator_DLPR(sig_num);
    CDSFparam.image_size=size(phase);
    [RMSE_phi1,RMSE_ampl1,PHI_dlpr] = DLPR(phase,amplitude,CDSFparam,ALGparam);
    RMSE_phi(sig_indx,:)=RMSE_phi1';
    RMSE_ampl(sig_indx,:)=RMSE_ampl1';
end
[RMSE_DLPR ~]=tablegenDL(RMSE_phi,RMSE_ampl,ALGparam)
end



