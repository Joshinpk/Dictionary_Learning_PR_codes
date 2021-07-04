function [Ri] = RiGen(CDSFparam)
patsize = CDSFparam.patsize;
image_size=CDSFparam.image_size;
Z=ones(image_size);
X = im2col(Z,[patsize patsize],'sliding');


blockx = patsize;
blocky = patsize;
final_numestimate = zeros(size(Z));
final_extentestimate = zeros(size(Z));
for indexi = 1:blocky
    for indexj = 1:blockx
        tempesti = reshape(X((indexi-1)*blockx+indexj,:),size(Z)-[blockx,blocky]+1);
        numestimate = zeros(size(Z));
        numestimate(1:size(tempesti,1),1:size(tempesti,2)) = 1;
        numestimate = circshift(numestimate, [indexj,indexi]-1);
        final_numestimate = final_numestimate+numestimate;
    end
end
Ri = final_numestimate;
end 