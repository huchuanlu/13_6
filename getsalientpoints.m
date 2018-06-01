function corner_im2 = getsalientpoints(rgb_im)
%% function corner_im2 = getsalientpoints(rgb_im)
%%compute salient points
%%Input:
%       rgb_im        : RGB color image
%%Output:
%       corner_im2    : detected salient points 

sigma_g=1.5; % parameters for computing harris points
sigma_a=5;  % parameters for computing harris points
nPoints=30; % number of salient points 
input_im=double(rgb_im); 
% input_im = RGB2Lab(input_im1);
Mboost = BoostMatrix(input_im);
boost_im= BoostImage(input_im,Mboost);
[EnIm]= ColorHarris(boost_im,sigma_g,sigma_a,0.04,1);
[x_max,y_max,corner_im2,num_max]=getmaxpoints(EnIm,nPoints);
