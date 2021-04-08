clear
close all

% -u_xx - u_yy = 2sin(x)sin(y)
% u = sin(x)sin(y)

u = @(x, y) sin(x).*sin(y);  % actual solution
f = @(x, y) 2 * sin(x).*sin(y); % rhs of ODE

a = 0;
b = 2*pi;
c = 0;
d = 2*pi;

alpha = 0;
beta = 0;
gamma = 0;
delta = 0;

N = 100;
[x, y] = meshgrid(linspace(a, b, N), linspace(a, b, N));
h = x(1,1) - x(1, 2);  % h is the same in both the y and x direction

% K = diag(2*ones(N-2,1))-diag(ones(N-3,1),-1)-diag(ones(N-3,1),1);
% K2D = kron(eye(N), K) + kron(K, eye(N));

rhs = h^2 * f(x(2:end-1, 2:end-1), y(2:end-1, 2:end-1));

% generate the eigenvalues
[k, l] = size(rhs);
L = 2*(1-cos((1:k)*pi/(k+1))) + 2*(1-cos((1:l)*pi/(l+1))); 
LL = diag(1./L);

% U = Q * inv(Λ) * Q * B
rhs1 = fast_sine_transform(rhs);
rhs1 = fast_sine_transform(rhs1');  % need to perform fft on rows and cols
rhs2 = LL * rhs1;
rhs2 = fast_sine_transform(rhs2);
sol = fast_sine_transform(rhs2');

% add back in the boundary conditions
uapp = [alpha * ones(1, l); sol; beta * ones(1, l)];
uapp = [gamma * ones(k+2, 1) uapp delta * ones(k+2, 1)];

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
    const = -sqrt(2/(n+1));
    % need to shift the vales of the matrix
    V1 = [zeros(1,m); V; zeros(n+1,m)];
    V2 = const * imag(fft(V1));
    Y = V2(2:n+1, :);
end