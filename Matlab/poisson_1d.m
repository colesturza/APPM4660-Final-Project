clear
close all

% -u_xx = -sin(x)
% u = -sin(x)

u = @(x) -sin(x);  % actual solution
f = @(x) -sin(x); % rhs of ODE

a = 0; b = 2*pi;
alpha = 0; beta = -sin(2*pi);

N = 100;
x = linspace(a, b, N);
h = x(2) - x(1);

% K = diag(2*ones(N-2,1))-diag(ones(N-3,1),-1)-diag(ones(N-3,1),1);

rhs = h^2 * f(x(2:end-1));
rhs(1) = rhs(1) + alpha;
rhs(end) = rhs(end) + beta;

% sol = K\b';  % slow method

% [S, Lam] = eig(K);  % remove later
% sol = S * inv(Lam) * inv(S) * b';  % inefficient way fast solver

n = N-2;

% generate the eigenvalues
lambda = zeros(n,1);
for i=1:n
   lambda(i) = 2*(1-cos(pi*i/(n+1)));
end
L = lambda;

b_prime = fast_sine_transform(rhs');
u_prime = b_prime ./ L;
sol = fast_sine_transform(u_prime);

uapp =[alpha; sol; beta];

uex = u(x);

figure(1)
plot(x,uapp,x,uex)
figure(2)
plot(x,abs(uapp-uex'),'o-')
err = abs(uapp-uex');

% fast sine transform for a 1D vector v
function y = fast_sine_transform(v)
    n = length(v);
    % extra term makes eigenvectors orthonormal so we only need one
    % transform for both ifft and fft (Q)
    const = -sqrt(2/(n+1));
    % need to shift the vales of the vector
    v1 = [0; v; zeros(n+1, 1)];
    v2 = const * imag(fft(v1));
    y = v2(2:n+1);
end