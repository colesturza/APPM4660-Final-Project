clear
close all

N = 102;

% -u_xx - u_yy = 2sin(x)sin(y)
% u = sin(x)sin(y)

u = @(x, y) sin(pi.*x).*sin(pi.*y);  % actual solution
f = @(x, y) 2 * pi^2 .* sin(pi.*x).*sin(pi.*y); % rhs of PDE

a = 1;
b = 1;

[x, y] = meshgrid(linspace(0, a, N), linspace(0, b, N));
h = x(1,1) - x(1, 2);  % h is the same in both the y and x direction

% evaluate the solution at the target
uex = u(x, y);

left_boundary = uex(:, 1);
right_boundary = uex(:, end);
bottom_boundary = uex(1, :);
top_boundary = uex(end, :);

B = (h^2).*f(x(2:end-1, 2:end-1), y(2:end-1, 2:end-1));
B(:, 1) = B(:, 1) + left_boundary(2:end-1);
B(:, end) = B(:, end) + right_boundary(2:end-1);
B(1, :) = B(1, :) + bottom_boundary(2:end-1);
B(end, :) = B(end, :) + top_boundary(2:end-1);

n = N-2;

% generate the eigenvalues
lambda = 2*(1-cos(pi*(1:n)/(n+1)));
L = lambda + lambda';

% U = Q * inv(Λ) * Q^{-1} * B
B_prime = fast_sine_transform(B);
B_prime = fast_sine_transform(B_prime');
U_prime = B_prime ./ L;
U_prime = fast_sine_transform(U_prime);
U = fast_sine_transform(U_prime');

% add back in the boundary conditions
uapp = [left_boundary(2:end-1) U right_boundary(2:end-1)];
uapp = [bottom_boundary; uapp; top_boundary];

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