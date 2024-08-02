clearvars
clear,clc
dbstop if error

%% Parameters of the problem

% scale of fluctuation of D (=2*correlation length)
theta = 2;
% number of elements
nele = 200;
% number of KL-modes
mkl = 200;
% length of domain
L = 2;
% element length
l = L/nele;
% a) Solve Integral eigenvalue problem for the exponential kernel
[lambda,phi] = EigFcnKL(mkl,L,theta/2);

% get quadrature nodes
%xGL=[-sqrt(3/5) 0 sqrt(3/5)];wGL=[5 8 5]/9; % 3-point gauss, l.77: use febar_gauss3.m
xGL=0; wGL=1;                               % midpoint-rule, l.76: use febar.m

XQ = l/2*(repmat(xGL',1,nele)+repmat(1:2:2*nele,length(xGL),1));
xq = XQ(:);
% b) plot eigenfunctions and average variance error
figure(1); for i=1:5,plot(phi{i}(xq));hold on; end; title('Eigenfunctions');legend('1','2','3','4','5');

err_var_avg = 1-cumsum(lambda)/L;

figure(2); plot(err_var_avg); xlabel('# of terms'); title('average variance error');

phi_at_xq = cell2mat(cellfun(@(c) c(xq'),phi,'un',0));

% axial rigidity
% mean, coefficient of variation and auto-correlation coefficient
mu_D = 100;
delta_D = 0.2;
sigma_D = delta_D*mu_D;
rho_D = @(z1,z2) exp(-abs(z1-z2));

% distributed load
q = 10;

% RF realization
xi = randn(mkl,1); %
Dsam = mu_D(ones(length(xq),1));
for i=1:mkl
    Dsam = Dsam+sqrt(lambda(i))*phi{i}(xq)*sigma_D*xi(i);
end
    
qsam = q*ones(length(xq),1);

usam = febar(Dsam,qsam,L,nele)';

figure(3)
plot(0:l:L,usam)
title('Displacement of bar')