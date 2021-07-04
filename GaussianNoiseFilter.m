function [Us]=GaussianNoiseFilter(V,Z,gamma_1,sigma)
% https://proofwiki.org/wiki/Cardano%27s_Formula. Coefficents : a=1, b=0, c, d
%V : Signal after Forward propagation (Us in the main code)
%Z: Noisy observation

I=sqrt(-1);
C=sigma^2/2/gamma_1-Z;
D=-(sigma^2/2/gamma_1)*abs(V);

Q=C/3;
R=-D/2;
Disc=Q.^3+R.^2; % Discriminant

Spos=sign((R+Disc.^0.5)).*(abs(R+Disc.^0.5)).^(1/3);
Tpos=sign(R-Disc.^0.5).*(abs(R-Disc.^0.5)).^(1/3);

Sneg=((R+Disc.^0.5)).^(1/3);
Tneg=((R-Disc.^0.5)).^(1/3);

S=Spos.*(Disc>=0)+Sneg.*(Disc<0);
T=Tpos.*(Disc>=0)+Tneg.*(Disc<0);

% S=Sneg;
% T=Tneg;

root3D(:,:,1)=S+T;
root3D(:,:,2)=-(S+T)/2+I*(S-T)*sqrt(3)/2;
root3D(:,:,3)=-(S+T)/2-I*(S-T)*sqrt(3)/2;

real_mask=(imag(root3D)==0);
real_roots=root3D.*real_mask;
bsl=max(real_roots,[],3); %
%So max operation at this stage will automatically choose the positive root


%% based on the assumption that the discriminant Disc is positive. ...

%% NO,  we assume that there is a non-negative solution. The sign of the discriminant is not important.

% bsl=bsl.*(Z>=0);

Us=bsl.*exp(I*angle(V));

end



% %     temp1=fftshift(temp1); temp2=fftshift(temp2); %% for proper zeroing the central part of the sensor
% % WW=ifftshift(WW); temp1=ifftshift(temp1);    %% return fft2 location of frequencies

