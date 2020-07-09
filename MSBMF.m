function [X, Y, iter] = MSBMF(M, D, R, lambda1, lambda2, lambda3, k, tol1, tol2, maxiter)
% MSBMF: Drug repositioning based on multi-similarity bilinear matrix factorization.
% Usage: [X, Y, iter] = MSBMF(M, D, R, lambda1, lambda2, lambda3, k, tol1, tol2, maxiter)
%
% Inputs:
%        M                  - the target matrix with only known entries and the unobserved entries are 0.
%        D                  - disease similarity matrix
%        R                  - drug similarity matrix
%        lambda1,2,3        - parameters needed to give.
%        k                  - the latent dimension of matrix factorization.
%        tol1, tol2         - tolerance of termination conditions.
%        maxiter            - maximum number of iterations.
%
% Outputs:
%        X, Y               - two latent low-rank matrices of the completed matrix.
%        iter               - the number of iterations.
%
% Written by: Mengyun Yang
% Email: mengyunyang@csu.edu.cn
% Created: December 16, 2019

rand('state', 2019); % fix random seed
omega = double(M ~= 0);
omega_ = ones(size(omega)) - omega;
U = rand(size(M, 1), k);
V = rand(size(M, 2), k);
P = rand(size(D, 2), k);
Q = rand(size(R, 2), k);
X = U;
Y = V;
Z = M;
W1 = zeros(size(U));
W2 = zeros(size(V));
XY = M;

rho = 1.05;
mu = 1e-4;
max_mu = 1e20;

stop1 = 1;
stop2 = 1;

for i = 1: maxiter
    U = (Z * V + lambda2 * D * P - W1 + mu * X) * inv(V' * V + lambda2 * P' * P + (lambda1 + mu) * eye(k));
    
    V = (Z' * U + lambda2 * R * Q - W2 + mu * Y) * inv(U' * U + lambda2 * Q' * Q + (lambda1 + mu) * eye(k));
    
    P = (lambda2 * D' * U) * inv(lambda2 * U' * U + lambda3 * eye(k));
    
    Q = (lambda2 * R' * V) * inv(lambda2 * V' * V + lambda3 * eye(k));
    
    X = U + (1 / mu) * W1;
    X(X < 0) = 0;
    
    Y = V + (1 / mu) * W2;
    Y(Y < 0) = 0;
    
    Z = M .* omega + (U * V') .* omega_;
    
    W1 = W1 + mu * (U - X);
    
    W2 = W2 + mu * (V - Y);
    
    stop1_0 = stop1;
    stop1 = norm(X * Y' - XY, 'fro') / norm(XY, 'fro');
    stop2 = abs(stop1 - stop1_0)/ max(1, abs(stop1_0));
    XY = X * Y';
    if stop1 < tol1 && stop2 < tol2
        iter = i;
        break
    else
        iter = i;
        mu = min(mu * rho, max_mu);
    end
    
end

end