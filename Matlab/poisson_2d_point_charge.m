clear
close all

N = 100;

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

figure(1)
mesh(x, y, uex)

alpha = uex(:, 1);
beta = uex(:, end);
gamma =  uex(1, :);
delta =  uex(end, :);

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

figure(2)
mesh(x, y, uapp)
figure(3)
mesh(x, y, abs(uapp-uex))

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

function uex = make_refsoln(x_trg,x_src)

N = size(x_trg,2);

dd1 = (x_src(1).*ones(1, N)-x_trg(1, :)).^2;
dd2 = (x_src(2).*ones(1, N)-x_trg(2, :)).^2;

uex = -1/(2*pi)*log(sqrt(dd1+dd2)).';

end