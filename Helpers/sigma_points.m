function sig_x = sigma_points(alpha, kappa, x_mean, x_var)
n = size(x_var,1);

%sigma points
c = (alpha^2)*(n + kappa);
deltaX = sqrt(c)*chol(x_var);
sig_x = zeros(n,2*n+1);
sig_x(:,1) = x_mean;
sig_x(:,2:n+1) = x_mean + deltaX.';
sig_x(:,n+2:end) = x_mean - deltaX.';
end