clear all
close all
clc
GraphicsSmoothing = 'off'

%% Load the original signal of length 2048
load('HeaviSine_2048.mat');
x=HeaviSine';
n=length(x)


figure(1)
subplot(2,2,1)
plot(x,'r')
% hold on

%% Generate the orthonormal wavelet transformation of length 10
% T=gen_wave(n,10,2);

load('Wavelet_T_2048_N_10'); % This gives you the wavelet transformation matrix T

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
% variable theta_rec(n)
% minimize (norm(theta_rec,1))
% subject to
% norm(T'*theta_rec-u)<=sqrt(n)*SD 
% cvx_end
% 
% %% Recover the x
% x_rec=T'*theta_rec;
% figure(1)
% subplot(2,2,4)
% plot(x_rec,'b')


%% Run the unconstrained composite minimization to recover the theta and x
mu=1; 
cvx_begin
cvx_solver gurobi
variable theta_rec_2(n)
minimize (0.5*theta_rec_2'*T*T'*theta_rec_2-u'*T'*theta_rec_2+0.5*u'*u+mu*norm(theta_rec_2,1))
cvx_end

%% Recover the x
x_rec_2=T'*theta_rec_2;
figure(1)
subplot(2,2,4)
plot(x_rec_2,'b')


% %% Use soft shrinkage
% mu=sqrt(2*log(n))*SD;
% theta_star=s_shrinkage_b_mu(T*u,0.6*mu);
% x_star=T'*theta_star;
% figure(1)
% subplot(2,2,4)
% plot(x_star,'b')

