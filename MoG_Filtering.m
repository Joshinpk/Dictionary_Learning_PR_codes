function [image_esti_c,SIGMA_MxMxK] = MoG_Filtering(Z,MoGparam,SIGMA_MxMxK,wt)

% sigma=0.1;

sigmar=function_stdEst2D_phase(real(Z(1:100,1:100)),1); 
sigmai=function_stdEst2D_phase(imag(Z(1:100,1:100)),1);
sigma = max(0.01,sqrt((sigmar^2+sigmai^2)));

% sigma=(std(real(Z(:))).^2+std(imag(Z(:))).^2)^0.5;


sigma_correction=0.1;

K=MoGparam.K; %10
patchsize=MoGparam.patchsize; %10
img_size=size(Z);
patch_vector_z = im2col(Z,[patchsize,patchsize]);

noise_dig_correction = (sigma_correction)^2;
for k=1:K
    [U,S] = eig(SIGMA_MxMxK(:,:,k));
    s_diag_new=max((diag(S)'-noise_dig_correction),0);
    S_new=diag(s_diag_new);
    SIGMA_MxMxK_denoise(:,:,k)=U*S_new*U';
    SIGMA_MxMxK_denoise(:,:,k)=(SIGMA_MxMxK_denoise(:,:,k)+SIGMA_MxMxK_denoise(:,:,k)')/2;
end

[patch_vector_mmse,~]=mmse_estimator(patch_vector_z,SIGMA_MxMxK_denoise,wt,sigma);
image_esti_c= depatch(patch_vector_mmse,ones(size(patch_vector_mmse)),patchsize,img_size);



end