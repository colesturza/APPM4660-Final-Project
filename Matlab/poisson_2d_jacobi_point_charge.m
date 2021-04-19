clear
close all

N = 50;

% -u_xx - u_yy = 0

u = @(x, y) sin(x) + sin(2.*y);  % actual solution
f = @(x, y) 0 * ones(N-2, N-2); % rhs of PDE

a = 0;
b = 1;
c = 0;
d = 1;

[x, y] = meshgrid(linspace(a, b, N), linspace(c, d, N));
h = x(1,1) - x(1, 2);  % h is the same in both the y and x direction


x_src = [1.25 1.75];

points = [x(:), y(:)];
x_trg = points.';

% evaluate the solution at the target
uex = make_refsoln(x_trg, x_src);

uex = reshape(uex, [N, N]);

alpha = uex(:, 1);
beta = uex(:, end);
gamma =  uex(1, :);
delta =  uex(end, :);

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

% uapp = K2D \ B;

uapp = jacobi(K2D, B, 2500);

uapp = reshape(uapp, [n, n]);

uapp = [alpha(2:end-1) uapp beta(2:end-1)];
uapp = [gamma; uapp; delta];

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

function uex = make_refsoln(x_trg,x_src)

    N = size(x_trg,2);

    dd1 = (x_src(1).*ones(1, N)-x_trg(1, :)).^2;
    dd2 = (x_src(2).*ones(1, N)-x_trg(2, :)).^2;

    uex = -1/(2*pi)*log(sqrt(dd1+dd2)).';

end