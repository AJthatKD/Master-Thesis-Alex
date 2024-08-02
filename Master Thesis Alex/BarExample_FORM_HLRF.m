clear 
close all 
clc

nrv = 5; % Number of random variables

rho = 0.7;
mu_Xi = 100;
sigma_Xi =10;
mu_X4 = 300;
sigma_X4 = 40;

dist(1) = ERADist('normal','MOM',[mu_Xi, sigma_Xi]);      %R1
dist(2) = ERADist('normal','MOM',[mu_Xi, sigma_Xi]);      %R2
dist(3) = ERADist('normal','MOM',[mu_Xi, sigma_Xi]);      %R3
dist(4) = ERADist('normal','MOM',[mu_Xi, sigma_Xi]);      %R4
dist(5) = ERADist('gumbel','MOM',[mu_X4, sigma_X4]);      %L

Rho = [1 rho rho rho 0;
       rho 1 rho rho 0;
       rho rho 1 rho 0;
       rho rho rho 1 0;
       0   0   0   0 1]; %correlation matrix

% if ~issymmetric(Rho) && ~isreal(Rho)
% error ('Rho is not symmetrical/Real'):
%end

%% applying transformation   
T_Nataf = ERANataf(dist,Rho);

%% limit state function
%limit state function in the original space
g = @(x) x(1) + x(2) + x(3) + x(4) - x(5); % I only have the LSF in physical space
Dg = [1 1 1 1 -1]; % I can also get the gradient vector in physical spaces

% LSF in the standard normal space 
%G = @(u) g(T_Nataf.U2X(u)); % [G(u) = g(T(u))] 

%% FORM using the HLRF

maxit = 50;     % Maximum number of iterations
tol = 1e-5;     % Tolerance, for sufficient convergence
u = zeros(nrv, maxit); % start from [0 0 0 0]
beta = zeros(1, maxit); % initialize vector for storing beta evaluated at the different iterations

%HLRF method
k = 1;
while true
    % 1. evaluate LSF at point u_k
    [x_k, J] = T_Nataf.U2X(u(:,k), 'Jac');

    % But i cannot evaluate G(u), because i dont have the expression 
    % ... 
    H_uk = g(x_k);
    
    % 2. evaluate LSF gradient at point u_k and direction cosines
    DH_uk = Dg/J;
    norm_DH_uk = norm(DH_uk);
    alpha = -DH_uk / norm_DH_uk;
    
    % 3. calculate d
    %d = (H_uk/norm_DH_uk + alpha*u(:,k))*alpha' - u(:,k);
 
    % 4. calculate u_(k+1)
    %u(:,k+1) = u(:,k) + d;
    u(:,k+1) = (H_uk/norm_DH_uk + alpha*u(:,k))*alpha';

    % 5. calculate beta
    beta(k) = norm(u(:,k));
       
    if (norm(u(:,k+1) - u(:,k)) <= tol) || (k == maxit)
       break;
    else
       k = k+1;
    end
end

u_star = u(:,k);

Pf = normcdf(-beta(k));