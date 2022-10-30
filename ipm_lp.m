function [xk, lambdak, sk, fk, muk, k, mu_series, x_seq, lambda_seq, s_seq] = ...
    ipm_lp(A, b, c, eta, eps, kmax, x0, lambda0, s0)
%
% function [xk, lambdak, sk, fk, muk, k] = ...
%     ipm_lp(A, b, c, eta, eps, kmax, x0, lambda0, s0)
%
% Function that performs the Interior Point Method for a LP 
% 
% INPUTS:
% c = vector of the objective function of the LP;
% A = matrix of the equality constraints in the LP;
% b = vector of the equality constraints in the LP;
% eta = parameter for the correction of the sequence.
% kmax = maximum number of iterations;
% eps = parameter for the stopping criterion;
% x0 = starting element of the xk sequence;
% lambda0 = starting element of the lambdak sequence;
% s0 = starting element of the sk sequence;
%
% OUTPUTS:
% xk = the last x computed by the function;
% lambdak = the last lambda computed by the function;
% sk = the last s computed by the function;
% fxk = objective function computed in xk;
% muk = the mu computed w.r.t. yk and lambdak;
% k = last iteration.
% mu_series = vector with the value of mu at each iteration
% x_seq = matrix with in the columns the value of x at each iteration
% lambda_seq = matrix with in the columns the value of lambda at each iteration
% s_seq = matrix with in the columns the value of s at each iteration 
%


% INITIALIZATION - PREPARING THE ITERATIONS
n = length(x0);
mu0 = (x0' * s0) / n; 

mu_series = zeros(1, kmax);
x_seq = zeros(n, kmax);
lambda_seq = zeros(3, kmax);
s_seq = zeros(n, kmax);


xk = x0;
lambdak = lambda0;
sk = s0;

muk = mu0;
k = 0;

% START OF The WHILE CYCLE
while k < kmax && muk > eps * mu0
    
    % INITIALIZE -F(xk, lambdak, sk)
    r1k = -A' * lambdak - sk + c;
    r2k = b - A * xk;
    r3k = - xk .* sk;
        
    % PREDICTION: COMPUTE (dxk_aff, dlambdak_aff, dsk_aff), i.e. solve the
    % Newton step represented by the lin. syst. 
    % JF(xk, lambdak, sk) * (dx, dlambda, ds)' = -F(xk, lambdak, sk)
    AD = A .* (xk ./ sk)';
    ADAt = AD * A';
    dlambdak_aff = ADAt\(AD * r1k - A * (r3k ./ sk) + r2k);
    dsk_aff = r1k - A' * dlambdak_aff;
    dxk_aff = (r3k - dsk_aff .* xk) ./ sk;
    
    % COMPUTE aP_aff, aD_aff, muk_aff, sigmak
    iidxk_aff = (dxk_aff < 0);
    iidsk_aff = (dsk_aff < 0);
    
    aP_aff = min([1, ...
        min(-xk(iidxk_aff)./dxk_aff(iidxk_aff))]);
    aD_aff = min([1, ...
        min(-sk(iidsk_aff)./dsk_aff(iidsk_aff))]);
    
    muk_aff = ((xk + aP_aff * dxk_aff)' * (sk + aD_aff * dsk_aff)) / n;
    sigmak = (muk_aff/muk)^3;
    
    % CORRECTION: COMPUTE (dxk, dlambdak, dsk)
    % REMEMBER: modify r3k
    r3k = r3k - dxk_aff .* dsk_aff + sigmak * muk;
    % Newton step represented by the lin. syst. 
    % JF(xk, lambdak, sk) * (dx, dlambda, ds)' = -F(xk, lambdak, sk)+Corr.
    dlambdak = ADAt\(AD * r1k - A * (r3k ./ sk) + r2k);
    dsk = r1k - A' * dlambdak;
    dxk = (r3k - dsk .* xk) ./ sk;
    
    % COMPUTE aP and aD
    iidxk = (dxk < 0);
    iidsk = (dsk < 0);
    
    aP_hat = min([1, ...
        min(-xk(iidxk)./dxk(iidxk))]);
    aD_hat = min([1, ...
        min(-sk(iidsk)./dsk(iidsk))]);
    
    aP = min([1, eta * aP_hat]);
    aD = min([1, eta * aD_hat]);
    
    % UPDATE
    xk = xk + aP * dxk;
    lambdak = lambdak + aD * dlambdak;
    sk = sk + aD * dsk;
    
    muk = (xk' * sk) / n;
    k = k + 1;
    
    mu_series(k) = muk;
    x_seq(:,k) = xk;
    lambda_seq(:,k) = lambdak;
    s_seq(:,k) = sk;
end

fk = c' * xk;
mu_series = mu_series(1:k);
x_seq = x_seq(:,1:k);
lambda_seq = lambda_seq(:,1:k);
s_seq = s_seq(:,1:k);


end