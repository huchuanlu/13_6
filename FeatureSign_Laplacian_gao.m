%% V is the sparse codes of all the feature. L_i is the i^{th} column of
%% laplacian matrix;
%% X is coeffcient matrix for all the instances.
%% min_x_i ||y_i-Bx_i||^2 +\gamma||x_i||_1 + \beta(2*v_i^T(VL_i)^T-v_i^TL_{ii}v_i).
%% The code is written by Shenghua Gao, SCE, NTU, SG. (email: gaos0004@ntu.edu.sg)
%  Feature sign search
%  reference
%  Efficient sparse coding algorithms
%  Honglak Lee Alexis Battle Rajat Raina Andrew Y. Ng
%       Computer Science Department
%       Stanford University
%       Stanford, CA 94305

% S(:,baseNum) = FeatureSign_Laplacian_gao(B,X(:,baseNum),...
%               L(:,baseNum),S, pars.gamma, pars.beta, baseNum,init_x);

function [x]= FeatureSign_Laplacian_gao(B,y,Li,V,gamma,beta,baseNum,init_x1)

Lii = Li(baseNum);

tmpLi = Li;
tmpLi(baseNum) = 0;
tmpVLi1 = V*tmpLi;
nbases=size(B,2);
init_x = zeros(nbases,1);
tmpVLi = zeros(nbases,1);
if baseNum == 1   
        tmpVLi = tmpVLi1(2:nbases+1); 
        init_x = init_x1(2:nbases+1);
    elseif ( (baseNum > 1) && (baseNum < nbases) )
        tmpVLi(1:baseNum-1) = tmpVLi1(1:baseNum-1);
        tmpVLi(baseNum:nbases) = tmpVLi1(baseNum+1:nbases+1);
        init_x(1:baseNum-1) = init_x1(1:baseNum-1);
        init_x(baseNum:nbases) = init_x1(baseNum+1:nbases+1);
    else
        tmpVLi(1:nbases) = tmpVLi1(1:nbases);
         init_x = init_x1(1:nbases);
end
% VLi = V*Li;


betaLii = diag(Lii*beta*ones(1,nbases));

OptTol = 1e-5;

if (nargin<8)
    x=zeros(nbases, 1);
else
    x = init_x;
end;

theta=sign(x);          %sign flag
a = (x~=0);             %active set

optc = 0;

By=B'*y;
B_h=B(:,a);
x_h=x(a);
Bx_h=B_h*x_h;

all_d=2*(B'*Bx_h-By + beta*tmpVLi+beta*Lii*x);

[ma mi]=max(abs(all_d).*(~a));



while optc==0,

    optc=1;
	
    if all_d(mi)>gamma+1e-10,
        theta(mi)=-1;
        a(mi)=1;
        %b=B(:,mi);
        %x(mi)=(gamma-all_d(mi))/(b'*b*2+beta*Lii*2); %% is this part
        %necessary?
    elseif all_d(mi)<-gamma-1e-10,
        theta(mi)=1;
        a(mi)=1;
        %b=B(:,mi);
        %x(mi)=(-gamma-all_d(mi))/(b'*b*2+beta*Lii*2);%% is this part
        %necessary?
    else
        if sum(a)==0,
%             gamma=ma-2*1e-10;
%             optc=0;
            b=B(:,mi);
            x(mi)=(By(mi)-beta*tmpVLi(mi))/(b'*b+beta*Lii);
            break;
        end
    end

    opts=0;
    B_h=B(:,a);
    x_h=x(a);
    theta_h=theta(a);
%     tmpVLi_h = tmpVLi(a);
%     VLi_h = VLi(a);
    while opts==0,
        opts=1;
        
        BB=B_h'*B_h;
        
        M = (B_h'*y-gamma*theta_h/2-beta*tmpVLi(a));
        N = (BB+betaLii(a,a));

        x_new=N\M;


        o_new = L1_cost(y,B_h,x_new,gamma,beta,tmpVLi(a)+Lii*x_new,Lii);

        

        %cost based on changing sign
        s=find(sign(x_new)~=theta_h);
        x_min=x_new;
        o_min=o_new;
        for j=1:length(s),
            zd=s(j);
            x_s=x_h-x_h(zd)*(x_new-x_h)/(x_new(zd)-x_h(zd));
            x_s(zd)=0;  %make sure it's zero
            o_s=L1_cost(y,B_h,x_s,gamma,beta, tmpVLi(a)+Lii*x_s,Lii);
            if o_s<o_min,
                x_min=x_s;
                o_min=o_s;
            end
        end

        x(a)=x_min;

        a=(x~=0);
        theta=sign(x);

        B_h=B(:,a);
        x_h=x(a);
        theta_h=theta(a);
        Bx_h=B_h*x_h;
		
        active_d=2*(B_h'*(Bx_h-y))+2*beta*tmpVLi(a)+2*beta*Lii*x_h+gamma*theta_h;

        if ~isempty(find(abs(active_d)>OptTol))
            opts=0;
        end
    end

    all_d=2*(B'*Bx_h-By + beta*tmpVLi + beta*Lii*x);

    [ma mi]=max(abs(all_d).*(~a));
    if ma>gamma+OptTol,
        optc=0;
    end
end

return;

function cost=L1_cost(y,B,x,gamma,beta, VLi, Lii)

tmp = y-B*x;
cost = tmp'*tmp+gamma*norm(x,1)+2*beta*x'*VLi-beta*x'*Lii*x;

return
