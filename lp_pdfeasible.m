function [x0, lambda0, s0] = lp_pdfeasible(A, b, c)

% "TILDE"-STEP) (very high) does not allow us to solve the linear system 
y_tilde = (A * A')\(b);
x_tilde = A' * y_tilde;

lambda_tilde = (A * A')\(A * c); % (A)
s_tilde = c - A' * lambda_tilde;

dx_tilde = max([0, ...
    -1.5 * min(x_tilde)]);

ds_tilde = max([0, ...
    -1.5 * min(s_tilde)]);

% "HAT"-STEP
x_hat = x_tilde + dx_tilde;
s_hat = s_tilde + ds_tilde;

dx_hat = 0.5 * (x_hat' * s_hat) / sum(s_hat);
ds_hat = 0.5 * (x_hat' * s_hat) / sum(x_hat);

% FINAL STEP
x0 = x_hat + dx_hat;
s0 = s_hat + ds_hat;
lambda0 = lambda_tilde;

end