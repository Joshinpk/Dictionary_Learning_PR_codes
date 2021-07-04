function varphi = trunc_Gauss_gen(mode)
xN=100;
yN=100;
vone=ones(xN,yN);
uvone=triu(vone);rotuvone=rot90(uvone);Upper=rotuvone.*uvone;
Upper(:,1:yN/2)=0;
cake=Upper+rot90(Upper)+rot90(rot90(Upper))+rot90(rot90(rot90(Upper)));

c=floor(xN/20);
varphi=gaussele(xN,yN,14*pi,2*c,3*c);
varphistore=varphi;
switch mode
    case 'G'
        varphi = varphi;
    case 'TL'
        varphi(1:xN/2,1:yN/2) = 0;
    case 'BL'
       varphi(xN/2+1:xN,1:yN/2) = 0;
    case 'TR'
        varphi(1:xN/2,1+yN/2:yN) = 0;
    case 'BR'
        varphi(xN/2+1:xN,1+yN/2:yN) = 0;
    case 'QTL'
        varphi(1:xN/2,1:yN/2) = 0;
        varphi=varphistore-varphi;
    case 'QBL'
       varphi(xN/2+1:xN,1:yN/2) = 0;
        varphi=varphistore-varphi;
    case 'QTR'
        varphi(1:xN/2,1+yN/2:yN) = 0;
        varphi=varphistore-varphi;
    case 'QBR'
        varphi(xN/2+1:xN,1+yN/2:yN) = 0;
        varphi=varphistore-varphi;
    case 'CAKE'
        varphi=varphistore.*cake;
    case 'SYNC' 
        [X,Y] = meshgrid(-50:49,-50:49);
        varphi=50*sinc(sqrt((X/(2*pi)).^2+(Y/(2*pi)).^2));
        imagesc(wrap(varphi))
        
end

end