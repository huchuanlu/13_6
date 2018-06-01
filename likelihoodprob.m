function [PrI_sal, PrI_bk,PrO_sal,PrO_bk] = likelihoodprob(input_im, ind,out_ind)
%% function [PrI_sal, PrI_bk,PrO_sal,PrO_bk] = likelihoodprob(rgb_im, ind,out_ind);
%%Input:
%       input_im      : RGB color image
%       ind           : pixels inside the convex hull
%%Output:
%       PrI_sal       : observation likelihood probability of pixel inside the convex hull be salient 
%       PrI_bk        : observation likelihood probability of pixel inside the convex hull belongs to background
%       PrO_sal       : observation likelihood probability of pixel outside the convex hull be salient
%       PrO_bk        : observation likelihood probability of pixel outside the convex hull belongs to background
%% compute iner and outer histogram
input_im = RGB2Lab(input_im);
imsize = numel((ind)) +  numel((out_ind));
mat_im = [reshape(input_im(:,:,1),1,imsize);reshape(input_im(:,:,1),1,imsize);reshape(input_im(:,:,1),1,imsize)];
maxValO = max(mat_im(:,ind),[],2);
minVal0 = min(mat_im(:,ind),[],2);
maxValB = max(mat_im(:,out_ind),[],2);
minValB = min(mat_im(:,out_ind),[],2);
numBin=[60,60,60]; % Number of bins in histogram 
smoothFactor=[5,6,6]; % Smoothing factor
smoothingKernel=cell(1,3);
PrI_sal = 1;%numel(ind)/(row*col)
PrI_bk = 1;
PrO_sal =1;
PrO_bk = 1;
for i = 1:3
    cur_im = input_im(:,:,i);
    dataMat = cur_im(ind);
    [innerHist,innerBin] = ComputeHistogram(dataMat,numBin(i),minVal0(i),maxValO(i));
    smoothingKernel{i}=getSmoothKernel_(smoothFactor(i));
    innerHist=filterDistribution_(smoothingKernel{i},innerHist',numBin(i));
    
    dataMat = cur_im(out_ind);
    [outerHist,outerBin] = ComputeHistogram(dataMat,numBin(i),minValB(i),maxValB(i));
    smoothingKernel{i}=getSmoothKernel_(smoothFactor(i));
    outerHist=filterDistribution_(smoothingKernel{i},outerHist',numBin(i));
    % compute bin values at the pixels inside the convex hull 
    PrO_H1 = innerHist(innerBin);
    PrO_H0 = outerHist(innerBin);
    PrI_sal=PrI_sal.*PrO_H1;
    PrI_bk=PrI_bk.*PrO_H0; 
    % compute bin values at the pixels outside the convex hull  
    PrB_H1 = innerHist(outerBin);
    PrB_H0 = outerHist(outerBin);
    PrO_sal=PrO_sal.*PrB_H1;
    PrO_bk=PrO_bk.*PrB_H0;    
end

function [intHist,binInd]=ComputeHistogram(dataMat,numBin,minVal,maxVal)
%% function [intHist,binInd]=ComputeHistogram(dataMat,numBin,minVal,maxVal)
%input:
%  dataMat:L or A orB with inner or outer index
%  numBin: number of bins
%  minVal: mimimal value of dataMat
%  maxVal: maximal value of dataMat
%output:
%  intHist: histogram of dataMat
%  binInd : regularized image

binInd=max( min(ceil(numBin*(double(dataMat-minVal)/(maxVal-minVal))),numBin),1);% 将像素值 转换到[1 - 60]之间，大小中原图等同
intHist=zeros(numBin,1);
for i = 1:length(dataMat)
    intHist(binInd(i))=intHist(binInd(i))+1;
end
%% Get smoothing filter
function [smKer]=getSmoothKernel_(sigma)

if sigma==0
    smKer=1;
    return;
end

dim=length(sigma); % row, column, third dimension
sz=max(ceil(sigma*2),1);
sigma=2*sigma.^2;

if dim==1
    d1=-sz(1):sz(1);
    
    smKer=exp(-((d1.^2)/sigma));
    
elseif dim==2
    [d2,d1]=meshgrid(-sz(2):sz(2),-sz(1):sz(1));
    
    smKer=exp(-((d1.^2)/sigma(1)+(d2.^2)/sigma(2)));
    
elseif dim==3
    [d2,d1,d3]=meshgrid(-sz(2):sz(2),-sz(1):sz(1),-sz(3):sz(3));
    
    smKer=exp(-((d1.^2)/sigma(1)+(d2.^2)/sigma(2)+(d3.^2)/sigma(3)));
    
else
    error('Not implemented');
end

smKer=smKer/sum(smKer(:));





%% Smooth distribution
function dist=filterDistribution_(filterKernel,dist,numBin)

if numel(filterKernel)==1
    dist=dist(:)/sum(dist(:));
    return;
end

numDim=length(numBin);

if numDim==1
    %smoothDist=conv(dist,filterKernel,'same');
    
    lenDist=length(dist);
    hlenKernel=(length(filterKernel)-1)/2;

    dist=[dist(1)*ones(1,hlenKernel),dist,dist(end)*ones(1,hlenKernel)];
    dist=conv(dist,filterKernel);
    lenSmoothDist=length(dist);
    offset=(lenSmoothDist-lenDist)/2;
    dist=dist((offset+1):(lenSmoothDist-offset));
    
elseif numDim==2
    dist=reshape(dist,numBin);
    
    dist=conv2(filterKernel,filterKernel,dist,'same');
    
else
    dist=reshape(dist,numBin);
    
    for i=1:numDim
        fker=ones(1,numDim);
        fker(i)=length(filterKernel);
        fker=zeros(fker);
        fker(:)=filterKernel(:);
        
        dist=convn(dist,fker,'same');
        
    end
    
end

dist=dist(:)/sum(dist(:));
    





 