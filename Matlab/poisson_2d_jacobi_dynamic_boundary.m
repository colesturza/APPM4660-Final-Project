clear
close all

% -u_xx - u_yy = sin(x) + 4*sin(2*y)
% u = sin(x) + sin(2*y)

u = @(x, y) sin(x) + sin(2.*y);  % actual solution
f = @(x, y) sin(x) + 4.*sin(2.*y); % rhs of PDE

N = 50;

a = 0;
b = 2*pi;
c = 0;
d = 2*pi;

[x, y] = meshgrid(linspace(a, b, N), linspace(c, d, N));
h = x(1,1) - x(1, 2);  % h is the same in both the y and x direction

alpha = u(x(:, 1), y(:, 1));
beta = u(x(:, end), y(:, end));
gamma =  u(x(1, :), y(1, :));
delta =  u(x(end, :), y(end, :));

B = (h^2).*f(x(2:end-1, 2:end-1), y(2:end-1, 2:end-1));
B(:, 1) = B(:, 1) + alpha(2:end-1);
B(:, end) = B(:, end) + beta(2:end-1);
B(1, :) = B(1, :) + gamma(2:end-1);
B(end, :) = B(end, :) + delta(2:end-1);

n = size(B, 1);  % size of the K2D U = B system

B = reshape(B, [n^2, 1]);

e = ones(n, 1);
K = diag(-e(1:end-1), -1) + diag(2*e) + diag(-e(1:end-1), 1);
K2D = kron(eye(n), K) + kron(K, eye(n));

uapp = jacobi(K2D, B, 750);

uapp = reshape(uapp, [n, n]);

uapp = [alpha(2:end-1) uapp beta(2:end-1)];
uapp = [gamma; uapp; delta];

uex = u(x, y);

figure(1)
mesh(x, y, uapp)
figure(2)
mesh(x, y, abs(uapp-uex))

function x = jacobi(A, b, iters)
    
    N = length(b);

    % decompose A into A = D + L + U
    D = diag(A);  % Diagonal components of A
    L = tril(A, -1);  % Lower triangular components of A
    U = triu(A, 1);  % Upper triangular components of A
    
    LU = L + U;
    
    x = zeros(N, 1);
    
    for i = 1:iters
        
        x = (b - LU * x) ./ D;
        
    end

end