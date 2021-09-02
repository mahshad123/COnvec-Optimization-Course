clear all
close all
clc


%% Load the original signal of length 2048
load('HeaviSine_2048.mat');
x=HeaviSine';
n=length(x)


figure(1)
subplot(2,2,1)
plot(x,'r')
% hold on

%% Generate the TV transformation
T=-eye(n);
for i=1:n-1
    T(i,i+1)=1;
end

theta=T*x;

figure(1)
subplot(2,2,3)
plot(theta,'b*')


%% Generate a noise and the noisy measurement
SD=0.5;
omega=SD*randn(n,1);
u=x+omega;
figure(1)
subplot(2,2,2)
plot(u,'r')


% %% Run the l1 minimization to recover the sparse theta
% cvx_begin
% cvx_solver gurobi
% variable x_rec(n)
% minimize (norm(T*x_rec,1))
% subject to
% norm(x_rec-u)<=0.91*sqrt(n)*SD  % The epsilon needs to be very carefully chosen
% cvx_end
% 
% %% Recover the x
% % x_rec=T'*theta_rec;
% figure(1)
% subplot(2,2,4)
% plot(x_rec,'b')


% Run the unconstrained composite minimization to recover the theta and x
mu=1; 
cvx_begin
cvx_solver gurobi
variable x_rec_2(n)
minimize (0.5*x_rec_2'*x_rec_2-u'*x_rec_2+0.5*u'*u+mu*norm(T*x_rec_2,1))
cvx_end

% Recover the x
% x_rec_2=T'*theta_rec_2;
figure(1)
subplot(2,2,4)
plot(x_rec_2,'b')

