function [varphi] = choose_image(image_numbr,image_param)

if isfield(image_param, 'size')
    siz=image_param.size;
    xN=siz(1);
    yN=siz(2);
end


if isfield(image_param, 'xN')
    xN=image_param.xN;
end

switch(image_numbr)
   case 0
        varphi=ones(xN,yN);
   case 1
        [X,Y] = meshgrid(-xN/2:xN/2-1,-yN/2:yN/2-1);
        A1=2;A2=13
        varphi=(A1-((X/(xN/2)).^2+(Y/(yN/2)).^2))*A2+0.1;     % Quadratic amplitude
   case 2
        c=floor(xN/20);
        varphi=gaussele(xN,yN,14*pi,2*c,3*c);  %Gaussian abs phase
%       load DEMO_1                     
%       varphi=varphi;
   case 3
        c=floor(xN/20);
        varphi=gaussele(xN,yN,14*pi,2*c,3*c);
        varphi(1:xN/2,1:yN/2) = 0;
%       load DEMO_3                      %Truncated Gaussian abs phase
%       varphi=varphi;
   case 4
        [X,Y] = meshgrid(0:xN-1,0:yN-1);
        %make discontinuous
        mask=ones(xN,yN); mask(:,1:yN/2)=0;Y = Y-xN/5;
        varphi=Y.*mask.*(Y>0);
%       load DEMO_4_small                %Shear plane abs phase
%       varphi=varphi;

    case 14
        varphi=pvs_gen(xN,yN);
    case 15
        param.orientation=image_param.orientation;
        [varphi, mask] = MRIdataRead(param);
end

end