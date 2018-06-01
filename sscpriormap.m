function sal = sscpriormap(input_im,ind,superlabel)
%% function sal = sscpriormap(input_im,ind,superlabel)
%%compute prior saliency map
%%Input:
%       input_im      : RGB color image
%       ind           : pixels inside the convex hull 
%%Output:
%       sal           : prior saliency map 


%  read superpixel label
[row,col] = size(superlabel);
super_prop = regionprops(superlabel, 'all');
super_num = numel(super_prop);
% compute the clustering result of the superpixel
aa = ssccluster(input_im,superlabel);
% find cluster label of each pixel 
im_label(row,col) = 0;
for m = 1: super_num
    im_label(super_prop(m).PixelIdxList) = aa(m);
end 
% define the saliency map by the cluster_ind of the ind inside the hull
cluster_prop = regionprops(im_label, 'all');
sal(row,col) = 0;
for m = 1:7
    cind1 = cluster_prop(1).PixelIdxList;
    cind2 = cluster_prop(2).PixelIdxList;
    cind3 = cluster_prop(3).PixelIdxList;
    cind4 = cluster_prop(4).PixelIdxList;
    cind5 = cluster_prop(5).PixelIdxList;
    cind6 = cluster_prop(6).PixelIdxList;
    cind7 = cluster_prop(7).PixelIdxList;
    sal(cind1) = numel(intersect(cind1,ind))/numel(cind1);
    sal(cind2) = numel(intersect(cind2,ind))/numel(cind2);
    sal(cind3) = numel(intersect(cind3,ind))/numel(cind3);
    sal(cind4) = numel(intersect(cind4,ind))/numel(cind4);
    sal(cind5) = numel(intersect(cind5,ind))/numel(cind5);
    sal(cind6) = numel(intersect(cind6,ind))/numel(cind6);
    sal(cind7) = numel(intersect(cind7,ind))/numel(cind7);
end
