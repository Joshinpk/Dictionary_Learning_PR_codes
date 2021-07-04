function varphi = shear_plane_gen(mode)
xN=100;
yN=100;
vone=ones(xN,yN);

switch mode
    case 'pyr'
        x = linspace(-50,49,100);
        y = linspace(-50,49,100);
        [X,Y] = meshgrid(x,y);
        varphi = (50-abs(X)) + (50-abs(Y)); 
    case 'step'
        [X,Y] = meshgrid(0:xN-1,0:yN-1);
        mask=ones(xN,yN); mask(:,[1:10 21:30 41:50 61:70 81:90])=0;
        varphi=Y.*(1-mask)+X.*mask;
    case 'mpyr'
        x = linspace(-50,49,100);
        y = linspace(-50,49,100);
        [X,Y] = meshgrid(x,y);
        varphi = (50-abs(X)) + (50-abs(Y));
        mask=ones(xN,yN);mask([1:15 86:100],:)=0;mask(:,[1:15 86:100])=0;
        varphi=varphi.*mask;
end

end
% 
% x = linspace(-50,49,100);
% y = linspace(-50,49,100);
% [X,Y] = meshgrid(x,y);
% Z = (50-abs(X)) + (50-abs(Y));
% % Z(Z < 0) = NaN;
% imagesc(wrap(Z))
% 
% 
% xN=100;yN=100;
% [X,Y] = meshgrid(0:xN-1,0:yN-1);
% mask=ones(xN,yN); mask(:,[1:10 21:30 41:50 61:70 81:90])=0;
% imagesc(mask)
% surfl(-X.*mask)
% figure
% varphi=Y.*(1-mask)+X.*mask;
% 
% 
% 
% Y = Y-xN/5;
% varphi=Y.*mask.*(Y>0);