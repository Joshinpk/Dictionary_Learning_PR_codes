function [phase,amplitude] = sig_generator_DLPR(sig_num)
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
       % Gropu-2 different amplitude
    case 3   %Mount, Shear Plane
       phase= choose_image(4,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(5,image_param);
       amplitude=(amplitude-min(amplitude(:)));
       amplitude=amplitude/max(amplitude(:))+1;
    case 4  %Quad, TG
       phase= choose_image(3,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(1,image_param);
       amplitude=amplitude/max(amplitude(:))+1;
    case 5 %Gauss, Shear Plane
       phase= choose_image(4,image_param);
       image_param.size=size(phase);
       amplitude=choose_image(2,image_param);
       amplitude=amplitude/max(amplitude(:))+1;
%      amplitude=amplitude/4+1;
       
       %Group 3 : Small amplitude very similar to phase
    case 6 % TG
       phase= choose_image(3,image_param);
       image_param.size=size(phase);
       type='HSim';k0=1;k1=1;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 7 % SP
       phase= choose_image(4,image_param);
       image_param.size=size(phase);
       type='HSim';k0=1;k1=1;
       amplitude= amplitude_generator(k0,k1,type,phase);
       
%        Group 4 B: Small amplitude low similarity to phase
    case 8 % TG
       phase= choose_image(3,image_param);
       image_param.size=size(phase);
       type='LSim';k0=1;k1=1;
       amplitude= amplitude_generator(k0,k1,type,phase);
    case 9 % SP
       phase= choose_image(4,image_param);
       image_param.size=size(phase);
       type='LSim';k0=1;k1=1;
       amplitude= amplitude_generator(k0,k1,type,phase);

    case 10    %MRI side view
        MRI_orientation='side';
        param.orientation=MRI_orientation;
       [varphi, mask]=MRIdataRead(param);
       phase=varphi.*mask;     
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 11    %MRI top view
        MRI_orientation='top';
        param.orientation=MRI_orientation;
       [varphi, mask]=MRIdataRead(param);
       phase=varphi.*mask;     
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
    case 12    %MRI front view
       MRI_orientation='front';
       param.orientation=MRI_orientation;
       [varphi, mask]=MRIdataRead(param);
       phase=varphi.*mask;     
       image_param.size=size(phase);
       amplitude=choose_image(0,image_param);
       
    otherwise
           error('Wrong test case number. Please enter numbers from 1 to 12')
end
end