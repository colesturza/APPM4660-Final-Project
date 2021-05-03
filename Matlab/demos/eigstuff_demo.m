N = 5;
e = ones(N, 1);


K = diag(-e(1:N-1), -1) + diag(2*e) + diag(-e(1:N-1), 1);

K2D = kron(eye(N), K) + kron(K, eye(N));

k = 1; l = 3;

j = (1:5)'; i = (1:5)';

S = sin(k*i*pi/(N+1)) * sin(l*j*pi/(N+1))';

S = reshape(S, [N^2, 1]);

lam = 2 - 2*cos(k*pi/(N+1)) + 2 - 2*cos(l*pi/(N+1));

lam_S_1 = K2D * S;
lam_S_2 = lam * S;