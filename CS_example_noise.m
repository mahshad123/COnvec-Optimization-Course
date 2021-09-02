clear all
close all
clc

%% Load a sample signal %%
load('HeaviSine_2048.mat');
x=HeaviSine';
n=length(x)


figure(1)
plot(x,'b')


%% Creat a Sparse signal
for k=1:length(x)
    flag=rand(1);
    if flag>0.3
        x(k)=0;
    end
end
figure(2)
plot(x,'r*')
hold on


%% Calculate the Sparsity and the sample size
K=length(find(abs(x)>0.0001))
c=2
m=ceil(c*K*log(n/K))


%% Create the sampling matrix Phi
Phi=sqrt(1/m)*randn(m,n);

SD=0.1;
omega=SD*randn(m,1);
y=Phi*x+omega;
figure(3)
plot(y)

% %% Solve the l1 minimization problem to recover x from y
% cvx_begin
% cvx_solver gurobi
% variable x_rec(2048)
% minimize (norm(x_rec,1))
% subject to
% norm(Phi*x_rec-y) <= 1*sqrt(m)*SD 
% cvx_end
% figure(2)
% plot(x_rec,'gO')


%% Solve the unconstrained l1-LS problem to recover x from y
mu=0.000001; %0.000001  0.00001   0.0001    0.001
cvx_begin
cvx_solver gurobi
variable x_rec_2(2048)
minimize (0.5*x_rec_2'*Phi'*Phi*x_rec_2-y'*Phi*x_rec_2+0.5*y'*y+mu*norm(x_rec_2,1))
cvx_end

figure(2)
plot(x_rec_2,'gO')





