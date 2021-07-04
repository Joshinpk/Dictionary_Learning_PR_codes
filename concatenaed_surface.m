function [phase] = concatenaed_surface(sig_num)
%% choose_image(#,image_param) # denotes the image number
% 0: Constant | 2: Gaussian          | 4:Shear Plane | 6: Camera Man | 8: Cheese 
% 1: Qudratic | 3:Truncated Gaussian | 5: Mountains  | 7: Lena       | 9: Test pat
%%
% image_param.size=[100,100];
% phase1= choose_image(3,image_param); 
% phase2= choose_image(4,image_param); 
% phase3= choose_image(5,image_param); 
% phase4= choose_image(10,image_param); 
% phase5= choose_image(1,image_param);
% % am=6; frq=6; mode='Disc'; phase6= sine_wave_generator(am,frq,mode,grid_size);
% phase6= choose_image(13,image_param);
% phase=[phase1 phase3 phase2 phase4 phase5 phase6];

phase1= trunc_Gauss_gen('G');
phase2= trunc_Gauss_gen('QTL');
phase3= trunc_Gauss_gen('QTR');
phase4= trunc_Gauss_gen('QBL');
phase5= trunc_Gauss_gen('QBR');
phase=[phase1 phase2 phase4 phase5 phase3];

end