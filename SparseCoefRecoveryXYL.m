%--------------------------------------------------------------------------
% This function takes the D x N matrix of N data points and the constrained
% matrix, and write every
% point as a sparse linear combination of other points.
% Xp: D x N matrix of N data points
% cst: 1 if using the affine constraint sum(c)=1, else 0
% Opt: type of optimization, {'L1Perfect','L1Noisy','Lasso','L1ED'}
% lambda: regularizartion parameter of LASSO, typically between 0.001 and 
% 0.1 or the noise level for 'L1Noise'
% W: the constrained matxi
% CMat: N x N matrix of coefficients, column i correspond to the sparse
% coefficients of data point in column i of Xp


function CMat = SparseCoefRecoveryXYL(Xp,cst,Opt,lambda,W,beta)

D = size(Xp,1);
N = size(Xp,2);
CMat = zeros(N);
param.lambda = lambda;
param.mode = 2;
param.lambda2 = 0;
for i = 1:N
      l = W(:,i);
    y = Xp(:,i);
    if i == 1
        Y = Xp(:,i+1:end);
        CMat_temp = CMat(2:end,:); % get rid of the ith row
        ll = (CMat_temp*l);
    elseif ( (i > 1) && (i < N) )
        Y = [Xp(:,1:i-1) Xp(:,i+1:N)]; 
        CMat_temp = [CMat(1:i-1,:); CMat(i+1:end,:)];
        ll = (CMat_temp*l);
         
    else
        Y = Xp(:,1:N-1);
        CMat_temp = CMat(1:end-1,:);
        ll = (CMat_temp*l);
    end
  
    % L1 initial optimization using CVX
        if ( strcmp(Opt , 'Lasso') )          
            cvx_begin;
            cvx_precision medium
            variable c(N-1,1);
            minimize( norm(c,1) + lambda * norm(Y * c  - y)) ;
            subject to
            sum(c) == 1;
            cvx_end;      
        end
    
    % place 0's in the diagonals of the coefficient matrix
    if i == 1   
        CMat(1,1) = 0;
        CMat(2:N,1) = c(1:N-1);       
    elseif ( (i > 1) && (i < N) )
        CMat(1:i-1,i) = c(1:i-1);
        CMat(i,i) = 0;
        CMat(i+1:N,i) = c(i:N-1);
    else
        CMat(1:N-1,N) = c(1:N-1);
        CMat(N,N) = 0;
    end
end
%% Solving LSSC
for iter = 1:10
for i = 1:N
      l = W(:,i);
    y = Xp(:,i);
    if i == 1
        Y = Xp(:,i+1:end);
        CMat_temp = CMat(2:end,:); % get rid of the ith row
        ll = (CMat_temp*l);
    elseif ( (i > 1) && (i < N) )
        Y = [Xp(:,1:i-1) Xp(:,i+1:N)]; 
        CMat_temp = [CMat(1:i-1,:); CMat(i+1:end,:)];
        ll = (CMat_temp*l);
    else
        Y = Xp(:,1:N-1);
        CMat_temp = CMat(1:end-1,:);
        ll = (CMat_temp*l);
    end
  
    % L1 optimization using CVX utilize the code of CVPR 10: Gao et al.
        if ( strcmp(Opt , 'Lasso') )
           c = zeros(N-1,1);
           c = FeatureSign_Laplacian_gao(Y,y,l,CMat,lambda,beta,i,CMat(:,i));
        end
    
    % place 0's in the diagonals of the coefficient matrix
    if i == 1   
        CMat(1,1) = 0;
        CMat(2:N,1) = c(1:N-1);       
    elseif ( (i > 1) && (i < N) )
        CMat(1:i-1,i) = c(1:i-1);
        CMat(i,i) = 0;
        CMat(i+1:N,i) = c(i:N-1);
    else
        CMat(1:N-1,N) = c(1:N-1);
        CMat(N,N) = 0;
    end
end
end


