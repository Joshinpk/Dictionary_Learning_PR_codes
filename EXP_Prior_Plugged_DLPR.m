clc; close all; clear all;
I=sqrt(-1);
addpath(strcat(pwd,'/SpaRSAL'));
%================================
%Test Signal Generation 
% 1: Truncated Gaussian, 2: MRI Side View, 3: MRI Top View, 4: MRI Front View
%================================
testcase_no=1;
%================================
%enter the Noise level
%================================
%for Poissonian noise
ALGparam.KAPPA_set=[.00001]; %Eg: ALGparam.KAPPA_set=[.00001, .0001, .001]


%====================================================================================================================
% Parameters
%====================================================================================================================
% CD_Sparse_Filtering
CDSFparam.K=256
CDSFparam.T=30;
CDSFparam.patchnumDic=64;
CDSFparam.lamdaDic=0.11;%0.110
CDSFparam.err=1e-3;%1e-3
CDSFparam.patsize=7;
CDSFparam.IsSubMean=1;
CDSFparam.ext_dic=0;
dic_pre_learn_en=1;
CDSFparam.iterndisplay=1;
CDSFparam.scale=1.2;

%BM3D Filtering
BM3Dparam.Nstep = 3;  
BM3Dparam.N11=8;  
BM3Dparam.N22=8; 
BM3Dparam.threshType='h'; 
BM3Dparam.threshold_ampl=1.4*4;
BM3Dparam.threshold_phi=1.4*4;

% Algorithm Parameters
ALGparam.Pois_Gaus_sel=1;    % 1 for Poissonian observation and 2 for Gaussian Observation
ALGparam.IterNum=25;         % Iteration number 
ALGparam.unwrapping=1;       % 1 means unwrapping on before BM3D phase filtering, 0- means unwrapping off
ALGparam.Fienup=0;           % 1 means that the Fienup rules is used instead of the optimal solution for data denoising
ALGparam.ISfiltering=1;      % 1 for BM3D filtering
ALGparam.TWF=1;              % Truncated Wirtinger Flow, TWF=1, acitivation, 0 -deactivation
ALGparam.support_part=1;     % Percentage of active pixels is support_part^2. support_part=.5 gives 25% 
ALGparam.L = 12;             % Number of Experiments 
ALGparam.SNR_Gaus_set=[10];

%====================================================================================================================
% Prior (Dictionery) Learning
%====================================================================================================================
%loading
%load Cs_F_89_ts2_P3;load Cs_F_120_ts2_P1; load Cs_F_120_ts2_P4; load Cs_F_120_ts4_P2;
%load Cs_S_95_ts1_P1;load Cs_S_108_ts3_P2; load Cs_S_108_ts1_P3; load Cs_S_105_ts3_P4;
%load Cs_T_110_ts4_P1; load Cs_T_97_ts4_P2; load Cs_T_110_ts4_P3; load Cs_T_100_ts4_P4;


switch(testcase_no)
    case 1
        DictNo=4;
        varphiDict = [trunc_Gauss_gen('G') trunc_Gauss_gen('QTL') trunc_Gauss_gen('QTR') trunc_Gauss_gen('QBL') trunc_Gauss_gen('QBR')];
        varphi = trunc_Gauss_gen('CAKE'); CDSFparam.patsize=10;CDSFparam.scale=1.5;

    case 2
        %MRI Side view
        DictNo=1;
        varphiDict = [Cs_S_95_ts1_P1  Cs_S_108_ts1_P3 Cs_S_105_ts3_P4];
        varphi=Cs_S_108_ts3_P2; 
    case 3
        %MRI Top view
        DictNo=2;
        varphiDict = [Cs_T_110_ts4_P1 Cs_T_110_ts4_P3 Cs_T_100_ts4_P4];
        varphi=Cs_T_97_ts4_P2; 
    case 4
        %MRI Front view 
        DictNo=3;
        varphiDict = [Cs_F_120_ts2_P1  Cs_F_89_ts2_P3  Cs_F_120_ts2_P4 ];
        varphi=Cs_F_120_ts4_P2;
           
    otherwise
           error('Wrong test case number. Please enter numbers from 1 to 4')
    
end

ALGparam.signum=DictNo;

%====================================================================================================================
% prior plugged DLPR
%====================================================================================================================
amplitude=ones(size(varphi));

if dic_pre_learn_en==1
Vdic=1.*exp(I*varphiDict);
D= CD_Dic_Learning(Vdic,CDSFparam);
strD=strcat('Prior_Dic_Num_',num2str(DictNo),'.mat')
% strD='mri_front_P134'
savefile = sprintf(strD,CDSFparam.patsize,CDSFparam.patsize);
save(savefile,'D');
end
CDSFparam.image_size=size(varphi);

[RMSE_DLPR,~,PHI_dlpr] = DLPR(varphi,amplitude,CDSFparam,ALGparam);

    
RMSE_DLPR
