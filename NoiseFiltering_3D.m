%%
function [Us]=NoiseFiltering_3D(V,I,Masks_Support,KAPPA,gamma_1,L)

Us=zeros(size(I));
for s=1:L
    temp1=V(:,:,s); temp2=I(:,:,s);
    temp1=fftshift(temp1); temp2=fftshift(temp2); %% for proper zeroing the central part of the sensor
    W=abs(temp1); gamma_11=gamma_1*KAPPA;

    W=(W/gamma_11+sqrt(W.^2/gamma_11^2+4*temp2*(1+1/gamma_11)/KAPPA))/2/(1+1/gamma_11);

    WW=abs(temp1).*(1-Masks_Support(:,:,s))+W.*Masks_Support(:,:,s);
    WW=ifftshift(WW); temp1=ifftshift(temp1);    %% return fft2 location of frequencies
    Us(:,:,s)=WW.*exp(1j*angle(temp1));
end