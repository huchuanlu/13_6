% function clusterhullsuper
% function sal = superhull(input_im, innerharrispoint, superlabel)
% input_im:         r g b color image
% innerharrispoint: 2*m ,first row for x and second row for y
% superlabel      : superpixel label
%% analogy the input parameters
clear all
clc
addpath('k_means');
sigma_g=1.5;
sigma_a=5;
nPoints=30;
thresh = 26;
im_dir = dir('E:\ICIP2011result\binarymaskimage\*.jpg')
% for i = 4:length(im_dir)
 for i = 2:10
    imname = im_dir(i).name;
% imname = '0_2_2721.jpg';
super_im=(imread([ 'E:\ICIP2011result\super 200\' imname(1:end-4) '_slic.jpg'])); 

input_im = double(imread(['E:\ICIP2011result\binarymaskimage\' imname]));
Mboost = BoostMatrix(input_im);
boost_im= BoostImage(input_im,Mboost);
[EnIm]= ColorHarris(boost_im,sigma_g,sigma_a,0.04);
[x_max,y_max,corner_im2,num_max]=getmaxpoints(EnIm,nPoints);
corner_im2 = elimatepoint(corner_im2,thresh);
[row,col] = size(corner_im2);
corner_imind = find(corner_im2 == 1);
siz = numel(corner_imind);
[y,x] = ind2sub([row,col],find(corner_im2 == 1));%y表示行数，x表示列数
dt = DelaunayTri(x,y);
if(~size(dt,1))
    return;
end
[k,av] = convexHull(dt);%求得外圈的点所在的位置
% innerharrispoint= [setdiff(x,x(k));setdiff(y,y(k))];
 %% 得到圈内的坐标点
 BW = roipoly(corner_im2,x(k),y(k));
 pixel = regionprops(BW,'all');
  ind = pixel.PixelIdxList;
  out_ind = setdiff(1:row*col,ind);
  input_im = RGB2Lab(input_im);
  R = input_im(:,:,1);
 G = input_im(:,:,2);
 B = input_im(:,:,3); % 3 color channels
 out_mean = [mean(R(out_ind)),mean(G(out_ind)),mean(B(out_ind))];
 superlabel = ReadDAT(size(corner_im2),['E:\ICIP2011result\super 200\' imname(1:end - 4) '.dat']);
 %% 聚成两类，分别求其于背景的聚离从而在圈内得到前景点背景点
STATS = regionprops(superlabel, 'all');
sup_num = numel(STATS);
innersuper = [];
insup_mean = [];
for r = 1:sup_num %% check the superpixel not along the image sides
    indxy = STATS(r).PixelIdxList;
    if (numel(intersect(indxy,ind)) > 0.6 * numel(indxy))
        innersuper = [innersuper,r];
        insup_mean = [insup_mean;mean(R(indxy)),mean(G(indxy)), mean(B(indxy))];
    end   
end
[cluster,Centroid] = cvKmeans(insup_mean',2);
cluster = cluster';
Centroid = Centroid';
dis = sum((Centroid - repmat(out_mean,2,1)).*(Centroid - repmat(out_mean,2,1)),2);
if(dis(1) > dis(2))
    one_object = 1; %1类是前景
else
    one_object = 2; %否则2类是前景
end
objinnersup_ind = find(cluster == one_object);% 内部点中属于前景概率比较大的superpixel
innersuper = innersuper(objinnersup_ind);

super_im=(imread([ 'E:\ICIP2011result\super 200\' imname(1:end-4) '_slic.jpg']));
R = super_im(:,:,1);
G = super_im(:,:,2);
B = super_im(:,:,3);

for m = 1:sup_num
    if(sum(find(innersuper == m)))
    pixelind = STATS(m).PixelIdxList;
    R(pixelind) = 255;
    G(pixelind) = 255;
    B(pixelind) = 0; 
    end
end 
super_im(:,:,1)= R;
super_im(:,:,2) = G;
super_im(:,:,3) = B;

figure(1)
imshow(super_im); hold on

  
 plot(x(k),y(k),'LineWidth',4,'Color','y');
 hold off
  F = getframe;
   close all
  imwrite(F.cdata,['E:\ICIP2011result\clusterhullsuper\',imname])
  clear F x y dt R G B row col 
end

