function [image_esti_c] = CD_Sparse_Filtering(Z,D,CDSFparam)
main_itern=CDSFparam.main_itern;
patsize = CDSFparam.patsize;
IsSubMean = CDSFparam.IsSubMean;
KAPPA=CDSFparam.KAPPA;
% sigma=CDSFparam.sigma;
scale=CDSFparam.scale;
sigmar=function_stdEst2D(real(Z),1); 
sigmai=function_stdEst2D(imag(Z),1);
sigma =scale*sqrt(sigmar^2+sigmai^2);
sigma=max(sigma,0.05);
% 
X = im2col(Z,[patsize patsize],'sliding');
if IsSubMean == 1
    meancolx = angle(sum(X));
    X = X.*(ones(size(X,1),1)*exp(-1i*meancolx));
end

% fprintf('Evaluating cost function...\n');
paramOMP.err = chi2inv(0.96,size(X,1)*2)*sigma^2/2;
paramOMP.Omp_ctrl_flag=CDSFparam.Omp_ctrl_flag;
% tic
alpha=OMP_C(D,X,paramOMP);
% omp_time=toc

X = D*alpha;
if IsSubMean == 1
    X=X.*(ones(size(X,1),1)*exp(1i*meancolx));
end


blockx = patsize;
blocky = patsize;
final_numestimate = zeros(size(Z));
final_extentestimate = zeros(size(Z));
for indexi = 1:blocky
    for indexj = 1:blockx
        tempesti = reshape(X((indexi-1)*blockx+indexj,:),size(Z)-[blockx,blocky]+1);
        numestimate = zeros(size(Z));
        extentestimate = zeros(size(Z));
        extentestimate(1:size(tempesti,1),1:size(tempesti,2)) = tempesti;
        numestimate(1:size(tempesti,1),1:size(tempesti,2)) = 1;
        
        extentestimate = circshift(extentestimate, [indexj,indexi]-1);
        numestimate = circshift(numestimate, [indexj,indexi]-1);
        
        final_numestimate = final_numestimate+numestimate;
        final_extentestimate = final_extentestimate+extentestimate;
    end
end
image_esti_c = final_extentestimate./final_numestimate;

end %fn CD_Sparse_Filtering end