function [phase,amplitude] = sig_generator(sig_num)
%% choose_image(#,image_param) # denotes the image number
% 0: Constant | 2: Gaussian          | 4:Shear Plane | 6: Camera Man | 8: Cheese 
% 1: Qudratic | 3:Truncated Gaussian | 5: Mountains  | 7: Lena       | 9: Test pat
%%


image_param.size=[100,100];
switch(sig_num)
    %Group 1: Invariant Amplitude 
    case 1    %const,Truncated Gauss
       phase= choose_image(3,image_param); 
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 2  %const,Shear Plane
       phase= choose_image(4,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);    
    case 3  %const,Sinusoid
       phase= choose_image(10,image_param);
%        phase=phase;
%        phase=phase-min(phase(:))+0.1;
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 4 
       phase=3*choose_image(6,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 5
       phase= choose_image(7,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 6
       phase= choose_image(9,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 7
       phase1= choose_image(7,image_param);
       image_param.size=size(phase1);
       phase2= choose_image(1,image_param);
       phase=phase1+phase2;
       amplitude=choose_image(0,image_param);
       %Group 2: Completely different amplitude
    case 8
       phase= choose_image(7,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(6,image_param); 
    case 9   %Mount, Shear Plane
       phase= choose_image(4,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(5,image_param);
       amplitude=(amplitude-min(amplitude(:)))/8+0.1;
    case 10
       phase= choose_image(7,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(5,image_param);
       amplitude=(amplitude-min(amplitude(:)))/4+0.1;
    case 11  %Quad, TG
       phase= choose_image(3,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(1,image_param);
       amplitude=amplitude/8;
    case 12
       phase= choose_image(9,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(3,image_param);
    case 13 %Gauss, Shear Plane
       phase= choose_image(4,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(2,image_param);
       amplitude=amplitude/4+1;
       %Group 3 A: Large amplitude very similar to phase
    case 14
       phase= choose_image(6,image_param);
       image_param.size=size(phase);
       type='HSim';k0=1;k1=5;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 15
       phase= choose_image(7,image_param);
       image_param.size=size(phase);
       type='HSim';k0=1;k1=5;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 16  % TG
       phase= choose_image(3,image_param);
       image_param.size=size(phase);
       type='HSim';k0=1;k1=5;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 17  % SP
       phase= choose_image(4,image_param);
       image_param.size=size(phase);
       type='HSim';k0=1;k1=5;
       amplitude= amplitude_generator(k0,k1,type,phase);
       %Group 3 B: Small amplitude very similar to phase
    case 18
       phase= choose_image(6,image_param);
       image_param.size=size(phase);
       type='HSim';k0=1;k1=1;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 19
       phase= choose_image(7,image_param);
       image_param.size=size(phase);
       type='HSim';k0=1;k1=1;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 20 % TG
       phase= choose_image(3,image_param);
       image_param.size=size(phase);
       type='HSim';k0=1;k1=1;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 21 % SP
       phase= choose_image(4,image_param);
       image_param.size=size(phase);
       type='HSim';k0=1;k1=1;
       amplitude= amplitude_generator(k0,k1,type,phase);
       %Group 4 A: Large amplitude low similarity to phase
    case 22
       phase= choose_image(6,image_param);
       image_param.size=size(phase);
       type='LSim';k0=1;k1=5;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 23
       phase= choose_image(7,image_param);
       image_param.size=size(phase);
       type='LSim';k0=1;k1=5;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 24 % TG
       phase= choose_image(3,image_param);
       image_param.size=size(phase);
       type='LSim';k0=1;k1=5;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 25 % SP
       phase= choose_image(4,image_param);
       image_param.size=size(phase);
       type='LSim';k0=1;k1=5;
       amplitude= amplitude_generator(k0,k1,type,phase);
%        Group 4 B: Small amplitude low similarity to phase
    case 26
       phase= choose_image(6,image_param);
       image_param.size=size(phase);
       type='LSim';k0=1;k1=1;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 27
       phase= choose_image(7,image_param);
       image_param.size=size(phase);
       type='LSim';k0=1;k1=1;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 28 % TG
       phase= choose_image(3,image_param);
       image_param.size=size(phase);
       type='LSim';k0=1;k1=1;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 29 % SP
       phase= choose_image(4,image_param);
       image_param.size=size(phase);
       type='LSim';k0=1;k1=1;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 30
       phase= choose_image(1,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 31
       phase= choose_image(2,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 32
       phase= choose_image(5,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 33
       phase= choose_image(10,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 34
       phase= choose_image(11,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 35
       phase= choose_image(14,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 36
       image_param.orientation='side';
       phase= choose_image(15,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 37
       image_param.data_sel='longs';
       phase= choose_image(16,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 38
       image_param.data_sel='isola';
       phase= choose_image(16,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
end
end