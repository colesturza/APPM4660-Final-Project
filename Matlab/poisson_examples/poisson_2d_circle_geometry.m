clear
close all

% -u_xx - u_yy = sin(x) + 4*sin(2*y)
% u = sin(x) + sin(2*y)

u = @(x, y) exact_solution(x, y);  % actual solution
f = @(x, y) density(x, y); % rhs of PDE

N = 100;

a = 0;
b = 5;
c = 0;
d = 5;

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

% slow way (works)
% K = diag(2*ones(N-2,1))-diag(ones(N-3,1),-1)-diag(ones(N-3,1),1);
% K2D = kron(eye(N-2), K) + kron(K, eye(N-2));
% rhs = reshape(rhs, [(N-2)^2, 1]);
% sol = K2D\rhs;
% sol = reshape(sol, [N-2, N-2]);

n = N-2;

% generate the eigenvalues
lambda = zeros(n,1);
for i=1:n
   lambda(i) = 2*(1-cos(pi*i/(n+1)));
end
L = lambda + lambda';

% U = Q * inv(Λ) * Q * B
B_prime = fast_sine_transform(B);
B_prime = fast_sine_transform(B_prime');  % need to perform fft on rows and cols
U_prime = B_prime ./ L;
U_prime = fast_sine_transform(U_prime);
U = fast_sine_transform(U_prime');

% add back in the boundary conditions
uapp = [alpha(2:end-1) U beta(2:end-1)];
uapp = [gamma; uapp; delta];

uex = u(x, y);

figure(1)
mesh(x, y, uapp)
figure(2)
mesh(x, y, abs(uapp-uex))
figure(3)
mesh(x, y, uex)

% fast sine transform for a matrix V
function Y = fast_sine_transform(V)
    [m, n] = size(V);
    % extra term makes eigenvectors orthonormal so we only need one
    % transform for both ifft and fft (Q)
    const = -sqrt(2/(n+1));
    % need to shift the vales of the matrix
    V_ext = [zeros(1,m); V; zeros(n+1,m)];
    V_ext = imag(fft(V_ext));
    Y = const.*V_ext(2:n+1, :);
end

function Z = density(X, Y)
    Z = sin(X) + 4.*sin(2.*Y);
    Z((X - 3).^2 + (Y - 3).^2 > 1) = 0;
end

function Z = exact_solution(X, Y)
    Z = sin(X) + sin(2.*Y);
    Z((X - 3).^2 + (Y - 3).^2 > 1) = 0;
end