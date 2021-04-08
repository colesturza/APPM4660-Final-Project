N = 100;
e = ones(N, 1);


K = diag(-e(1:N-1), -1) + diag(2*e) + diag(-e(1:N-1), 1);

K2D = kron(eye(N), K) + kron(K, eye(N));

a = 2*pi;
b = 0;

h = (a-b)/(N -1);

[x, y] = meshgrid(b:h:a, b:h:a);


f = sin(x).*sin(2*y);




f = reshape(f, [N^2, 1]);




%%Method using just inverse not good change this to fft later.

U = K2D\f;

U = reshape(U, [N, N]);

Ucir = U;
Utri = U;

%%Circular Geometry. 
Ucir((x - pi).^2 + (y - pi).^2 > pi) = NaN;



%%Triangle Geometry
Utri((y) - x > 0) = NaN;


mesh(x, y, Utri)