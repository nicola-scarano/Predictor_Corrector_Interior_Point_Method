%% Test of the IMP method on given A,b and c

clear 
clc
format long

% IMP paramenters loading
load test_iplp.mat  

% 
% eta = 0.95;
% eps = 1e-7;
% kmax = 

% Setting of n and a and generation of A, b and c 
%n = [10^4, 10^6]
%a = [2,20,200]
n = 1e6;
a = 200;
[A,b,c] = Abc_setup(n,a);
%% IMP infeasible 

% set of x0, lambda0, s0 to random value
% we set the random generator seed to repeatability 
rng(5); 
x0 = rand(n,1);
rng(10);
lambda0 = rand(3,1);
rng(20);
s0 = rand(n,1);

disp('**** IMP infeasible: OPTIONS *****')
disp(['eta: ', eta])
disp(['eps: ', eps])

disp('**** IMP : START *****')
tic
[xk, lambdak, sk, fk, muk, k_inf, mu_inf, x_inf, lambda_inf, s_inf] = ...
ipm_lp(A, b, c, eta, eps, kmax, x0, lambda0, s0);
t = toc;
disp('**** IMP : FINISHED *****')
disp('**** IMP : RESULTS *****')
disp('************************************')
% disp(['xk: ', vec2str(xk), ';'])
% disp(['sk: ', vec2str(sk), ';'])
% disp(['lambdak: ', vec2str(lambdak), ';'])
disp(['muk: ', num2str(muk), ';'])
disp(['fk: ', num2str(fk), ';'])
disp(['N. of Iterations: ', num2str(k_inf),'/',num2str(kmax), ';'])
disp(['time of execution: ', num2str(t), ';'])

disp('************************************')

%% IMP feasible

% set of x0, lambda0, s0 to feasible points using the function
% lp_pdfeasible()
[x0, lambda0, s0] = lp_pdfeasible(A, b, c);

% disp('*** new x0, s0 and lambda0')
% disp(['x0: ', vec2str(x0), ';'])
% disp(['s0: ', vec2str(s0), ';'])
% disp(['lambda0: ', num2str(lambda0), ';'])

disp('**** IMP feasible: OPTIONS *****')
disp(['eta: ', eta])
disp(['eps: ', eps])

disp('**** IMP : START *****')
tic
[xk, lambdak, sk, fk, muk, k_f, mu_f, x_f, lambda_f, s_f] = ...
ipm_lp(A, b, c, eta, eps, kmax, x0, lambda0, s0);
t = toc;
disp('**** IMP : FINISHED *****')
disp('**** IMP : RESULTS *****')
disp('************************************')
% disp(['xk: ', vec2str(xk), ';'])
% disp(['sk: ', vec2str(sk), ';'])
% disp(['lambdak: ', vec2str(lambdak), ';'])
disp(['muk: ', num2str(muk), ';'])
disp(['fk: ', num2str(fk), ';'])
disp(['N. of Iterations: ', num2str(k_f),'/',num2str(kmax), ';'])
disp(['time of execution: ', num2str(t), ';'])

disp('************************************')

%% Computation of the series of value of dual function and primal function

pri_inf = zeros(1,k_inf);
dual_inf = zeros(1,k_inf);
pri_f = zeros(1,k_inf);
dual_f = zeros(1,k_inf);

for i = 1:k_inf
    pri_inf(i) = c'*x_inf(:,i);
    dual_inf(i) = b'*lambda_inf(:,i);
end

for i = 1:k_f
    pri_f(i) = c'*x_f(:,i);
    dual_f(i) = b'*lambda_f(:,i);
end

M = zeros(k_inf,4);
M(:,1) = pri_inf';
M(:,2) = dual_inf';
M(:,3) = pri_f';
M(:,4) = dual_f';

writematrix(M,'M_tab.csv','Delimiter',',')
type 'M_tab.csv'
