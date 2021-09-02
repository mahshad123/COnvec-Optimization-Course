clc 
clear all
A=[1 -1 1 0; -1 -1 0 1];
b=[1 ; -1];
H=[8 2 0 0 ; 2 2 0 0; 0 0 0 0 ; 0 0 0 0 ];
p=[-26; -9; 0; 0 ];
x0=[ 1 1 0.1 0.1]';
mu0=[1 1 .5 .6]';
lambda0=[1 .1]';
% H=[4 0 0; 0 1 -1; 0 -1 1];
% p=[-8 -6 -6]';
% A=[1 1 1];
% b=3;
% x0=[1;2;2];
% mu0=[0.2 ;0.2 ;0.2 ];
% lambda0=-1;
x=x0;
mu=mu0;
lambda=lambda0;
k=0;
while mu'*x >= 0.000001
    rd=H*x+p+A'*lambda-mu;
    rp=b-A*x;
    gamma=(diag(mu)+diag(x)*H)^-1;
    Y0= (A*gamma*diag(x)*A')^-1;
    T= (x'*mu)/(size(x,1)+ 2*sqrt(size(x,1))+size(x,1));
    yd= gamma* [ diag(x)*(mu+rd)-T*(ones(size(x,1),1))];
    delLambda=-Y0*( A*yd+rp);
    delX=-gamma*diag(x)*A'*delLambda-yd;
    delmu=H*delX+A'*delLambda+rd;
    alphap= 100000;
    for i=1 : size(delX,1)
        if delX(i)<0
            alphap= min (alphap ,  -x(i)/delX(i));
        end
    end
    alphad=100000;
    for i=1 : size(delmu,1)
        if delmu(i)<0
            alphad= min (alphad ,  -mu(i)/delmu(i));
        end
    end

    alpha=(1-10^(-6))*min( alphad, alphap);
    x=x+alpha*delX;
    mu=mu+alpha*delmu;
    lambda=lambda+alpha*delLambda;
    k=k+1;
end

3-x(1)
-x(2)
