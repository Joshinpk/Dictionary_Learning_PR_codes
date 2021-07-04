I=sqrt(-1);
a=1;
b=0;
c=255;
d=-1;

Q=(3*a*c-b*b)/(9*a*a);
R=(9*a*b*c-27*a^2*d-2*b^3)/(54*a^3);

Disc=Q^3+R^2

S=sign((R+Disc.^0.5))*(abs(R+Disc.^0.5)).^(1/3);
T=sign(R-Disc.^0.5)*(abs(R-Disc.^0.5)).^(1/3);

x1=S+T-b/3/a;
x2=-(S+T)/2+I*(S-T)*sqrt(3)/2-b/3/a;
x3=-(S+T)/2-I*(S-T)*sqrt(3)/2-b/3/a;

card=[x1; x2 ;x3]
root_real=card(card==real(card));
root_real_pos=root_real(root_real>=0)
% roots([a b c d])
%%
c=5;
d=-6;


Q=c/3;
R=-d/2;
Disc=Q^3+R^2
S=sign((R+Disc.^0.5))*(abs(R+Disc.^0.5)).^(1/3);
T=sign(R-Disc.^0.5)*(abs(R-Disc.^0.5)).^(1/3);

x1=S+T;
x2=-(S+T)/2+I*(S-T)*sqrt(3)/2;
x3=-(S+T)/2-I*(S-T)*sqrt(3)/2;

card=[x1; x2 ;x3]
root_real=card(card==real(card));
root_real_pos=max(root_real)
%%


card(:,:,1)=[1+3i 3;   -1-2i 4-5i];
card(:,:,2)=[1    2+3i;-1    2-5i];
card(:,:,3)=[1+3i 2;   +1       4];

real_mask=(imag(card)==0);
real_roots=card.*real_mask;
real_roots_positive=max(real_roots,[],3)





