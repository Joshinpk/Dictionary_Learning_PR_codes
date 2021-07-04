function [D] = CD_Dic_Learning(Z,CDSFparam)
if isfield(CDSFparam, 'Dinitial')
    paramDL.Dinitial = CDSFparam.Dinitial;
end
% mask=CDSFparam.mask;
patsize = CDSFparam.patsize;
paramDL.K = CDSFparam.K;%256
paramDL.T = CDSFparam.T;%500
paramDL.patchnum = CDSFparam.patchnumDic;
paramDL.lamda = CDSFparam.lamdaDic;%0.11
paramDL.err = CDSFparam.err;%1e-3
paramDL.iterndisplay =CDSFparam.iterndisplay;
IsSubMean=CDSFparam.IsSubMean;
% mask=CDSFparam.mask;
X = im2col(Z,[patsize patsize],'sliding');
% maskC = im2col(mask,[patsize patsize],'sliding');
% indexm = find(sum(abs(maskC))~=0);
% clear maskC;
% Xmask = X(:,indexm);
Xmask=X;
if IsSubMean == 1
    meancolx = angle(sum(Xmask));
    Xmask = Xmask.*(ones(size(Xmask,1),1)*exp(-1i*meancolx));

end


D = DicLearningM(Xmask,paramDL);

end