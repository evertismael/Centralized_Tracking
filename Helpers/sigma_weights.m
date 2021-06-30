function [sig_wm, sig_wc] = sigma_weights(n,alpha, beta, kappa)
% weights:
lambda = (alpha^2)*(n + kappa) - n;
sig_wm = zeros(1,2*n+1);
sig_wm(:,1) = lambda/(n+lambda);
sig_wm(:,2:end) = 1/(2*(n+lambda));

sig_wc = zeros(1,2*n+1);
sig_wc(:,1) = lambda/(n+lambda) + 1 - alpha^2 + beta;
sig_wc(:,2:end) = 1/(2*(n+lambda));
end