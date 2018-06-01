
addpath(genpath('./cvx'));
addpath(genpath('./ColorFeatures'));
thresh = 26; 
rgb_im = imread('image.jpg');
%% obtain the interesting points
corner_im2 = getsalientpoints(rgb_im);
%% elimate the points too closing to the boundary of images
corner_im = elimatepoint(corner_im2,thresh);
%% obtain the convex hull
[row,col] = size(corner_im);
[y,x] = ind2sub([row,col],find(corner_im == 1));%find positiion of the corner points
dt = DelaunayTri(x,y);
if(~size(dt,1))
    return;
end
[k,av] = convexHull(dt);%find the points to plot the convex hull
sulabel_im = ReadDAT(size(corner_im),'image.dat');% obtain the superpixel labels
BW = roipoly(corner_im,x(k),y(k));%obtain the pixels inside the convex hull
pixel = regionprops(BW,'all');
ind = pixel.PixelIdxList;
out_ind = setdiff(1:row*col,ind);%obtain the pixels outside the convex hull
%% obtain the prior saliency map
psal = sscpriormap(rgb_im, ind, sulabel_im);

%% compute the likelihood probability
[PrI_sal, PrI_bk,PrO_sal,PrO_bk] = likelihoodprob(rgb_im, ind,out_ind);
%% compute final saliency map using the Bayes formula
psal_I = psal(ind);
psal_O = psal(out_ind');
Pr_0=(PrI_sal.*psal_I)./(PrI_sal.*psal_I+PrI_bk.*(1 - psal_I));%so called saliency 窗内的saliency
Pr_B=(PrO_sal.*psal_O)./(PrO_sal.*psal_O+PrO_bk.*(1-psal_O));%so called saliency 窗外的saliency
saliencymap = zeros(row,col);
saliencymap(ind) = Pr_0;
saliencymap(out_ind) = Pr_B;
saliencymap = (saliencymap - min(saliencymap(:)))/(max(saliencymap(:)) - min(saliencymap(:)));
figure
imshow(saliencymap);
