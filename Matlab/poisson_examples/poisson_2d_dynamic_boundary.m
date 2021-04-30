clear
close all

% -u_xx - u_yy = sin(x) + 4*sin(2*y)
% u = sin(x) + sin(2*y)

u = @(x, y) exp(x) + (1/2).*sin(6*pi.*y);  % actual solution
f = @(x, y) -exp(x) + 18*pi^2.*sin(6*pi.*y); % rhs of PDE

N = 102;

a = 0;
b = 1;
c = 0;
d = 1;

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

n = N-2;

% generate the eigenvalues
lambda = 2*(1-cos(pi*(1:n)/(n+1)));
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

% fast sine transform for a matrix V
function Y = fast_sine_transform(V)
    [m, n] = size(V);
    % extra term makes eigenvectors orthonormal so we only need one
    % transform for both ifft and fft (Q)
    const = sqrt(2/(n+1));
    % need to shift the vales of the matrix
    V_ext = [zeros(1,n); V; zeros(m+1,n)];
    V_ext = imag(fft(V_ext));
    Y = const.*V_ext(2:m+1, :);
end