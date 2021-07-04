clear all
%% CD_Sparse_Filtering
CDSFparam.K=256;%256
CDSFparam.T=500;%500
CDSFparam.patchnumDic=64;
CDSFparam.lamdaDic=0.11;%0.11
CDSFparam.err=1e-3;%1e-3
CDSFparam.patsize=10;
CDSFparam.IsSubMean=0;
CDSFparam.Omp_ctrl_flag=0;
CDSFparam.main_itern=100;
CDSFparam.KAPPA=1;
patchsize=CDSFparam.patsize;
CDSFparam.iterndisplay=1;

sig_set=3;

%%

if 1
for sig_indx=1:length(sig_set)
sig_num=sig_set(sig_indx);
[varphi,amplitude] = sig_generator(sig_num);

[yN,xN]=size(varphi); 


  
x=amplitude.*exp(1j*varphi);                           %% TRUE COMPLEX-VALUED OBJECT
sigma=0.9;
CDSFparam.sigma2=sigma;
I=sqrt(-1);
noise = sigma*(randn(size(x))+I*randn(size(x)))/sqrt(2);
z=x+noise;

%% Dic Learning
Vdic=x;
% addpath(strcat(pwd,'\SpaRSAL'));
D= CD_Dic_Learning(Vdic,CDSFparam);

strD=strcat('DicSIG_',num2str(sig_num),'.mat')
savefile = sprintf(strD,patchsize,patchsize);
save(savefile,'D');
end
end
% load TGussDic.mat;
% %% Filter
% x_hat= CD_Sparse_Filtering(z,D,CDSFparam);
% varphi_hat=angle(x_hat);
% 
% %%
% % colormap gray
% % subplot(311)
% % imagesc(angle(x))
% % subplot(312)
% % imagesc(angle(z))
% % subplot(313)
% % imagesc(varphi_hat)
% patchsize=CDSFparam.patsize;
% [M N]=size(varphi);
% wraperr_norm = norm(wrap(varphi_hat - varphi),'fro')^2;
% PSNR1 = 10*log10(4*M*N*pi^2/wraperr_norm)

%%
if 0

I=sqrt(-1);
[varphi_SP,amplitude] = sig_generator(2);
[varphi_Quad,amplitude] = sig_generator(30);
[varphi_Gaus,amplitude] = sig_generator(31);
[varphi_Mou,amplitude] = sig_generator(32);
[varphi_Hsin,amplitude] = sig_generator(33);
[varphi_Vsin,amplitude] = sig_generator(34);

x1=amplitude.*exp(I*varphi_SP);
x2=amplitude.*exp(I*varphi_Quad);
x3=amplitude.*exp(I*varphi_Gaus);
x4=amplitude.*exp(I*varphi_Mou);
x5=amplitude.*exp(I*varphi_Hsin);
x6=amplitude.*exp(I*varphi_Vsin);

x=[x1 x2 x3 x4 x5 x6];

Vdic=x;
% addpath(strcat(pwd,'\SpaRSAL'));
D= CD_Dic_Learning(Vdic,CDSFparam);

strD=strcat('','.mat')
savefile = sprintf(strD,patchsize,patchsize);
save(savefile,'D');

surf(angle(x))

mask=zeros(size(varphi_SP));
mask11=mask;mask12=mask;mask13=mask;
mask21=mask;mask22=mask;mask32=mask;
mask31=mask;mask23=mask;mask33=mask;

mask11(1:33,1:33)=1; mask12(1:33,34:66)=1; mask13(1:33,67:100)=1;
mask21(34:66,1:33)=1;mask22(34:66,34:66)=1; mask23(34:66,67:100)=1;
mask31(67:100,1:33)=1;mask32(67:100,34:66)=1;mask33(67:100,67:100)=1;

varphi=mask11.*varphi_Vsin+mask12.*varphi_Mou+mask13.*varphi_Hsin+...
       mask21.*varphi_Quad+mask22.*varphi_Gaus+mask23.*varphi_SP+...
       mask31.*varphi_Hsin+mask32.*varphi_Vsin+mask33.*varphi_Mou;


% colormap gray; 
surf(wrap(varphi))

end





% figure
% subplot(321);imagesc(varphi2);
% subplot(322);imagesc(varphi30);
% subplot(323);imagesc(varphi31);
% subplot(324);imagesc(varphi32);
% subplot(325);imagesc(varphi33);
% subplot(326);imagesc(varphi34);
% 
% figure
% subplot(331);imagesc(mask11);
% subplot(334);imagesc(mask21);
% subplot(337);imagesc(mask31);
% subplot(332);imagesc(mask12);
% subplot(335);imagesc(mask22);
% subplot(338);imagesc(mask32);
% subplot(333);imagesc(mask13);
% subplot(336);imagesc(mask23);
% subplot(339);imagesc(mask33);















