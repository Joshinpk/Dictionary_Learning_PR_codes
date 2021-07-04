function [Us]=GaussianNoiseFilterCardan(V,Z,sigma, Masks_Support,gamma_1,L)

Us=zeros(size(Z));
for s=1:L
  
    temp1=V(:,:,s); temp2=Z(:,:,s);
   %% temp1=fftshift(temp1); temp2=fftshift(temp2); %% for proper zeroing the central part of the sensor
  [Filtered]= GaussianNoiseFilter(temp1,temp2,gamma_1,sigma);
  
   %% [Filtered]=GaussianNoiseFilter(temp1,temp2,gamma_1,sigma);
  
    WW=abs(temp1).*(1-Masks_Support(:,:,s))+abs(Filtered).*Masks_Support(:,:,s);
   %% WW=ifftshift(WW); temp1=ifftshift(temp1);    %% return fft2 location of frequencies
    
    Us(:,:,s)=WW.*exp(1j*angle(temp1));

% temp1=V(:,:,s); temp2=Z(:,:,s);
% %    temp1=fftshift(temp1); temp2=fftshift(temp2); %% for proper zeroing the central part of the sensor
%      
%     [Us]=GaussianNoiseFilter(temp1,temp2,gamma_1,sigma);
%   
    %%WW=abs(temp1).*(1-Masks_Support(:,:,s))+abs(Us).*Masks_Support(:,:,s);
   %% WW=ifftshift(WW); 
%     temp1=ifftshift(temp1);    %% return fft2 location of frequencies
%     
%     Us(:,:,s)=WW.*exp(1j*angle(temp1));
   
end
    
