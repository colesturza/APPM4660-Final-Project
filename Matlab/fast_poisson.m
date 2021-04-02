%                         d^2 u(x,y)   d^2 u(x,y)
%    2D-Laplacian(u)  =   ---------- + ----------  = f(x,y)   
%                           d x^2          d y^2
%    
%    u(x,y) = 0 if (x,y) is on the boundary of Omega

xmin = 0;
xmax = 1;
ymin = 0;
ymax = 1;

n = 45;
h = 1/(n+1);

[x,y] = meshgrid(xmin:h:xmax, ymin:h:ymax);

F = func(x, y);

B = F;

F = Poisson_FFT(B);

surf(x,y,F)

% RHS of 2D-Laplacian(u)
function f = func(x, y)
    f = sin(x).*cos(x.*y) - sqrt(y);
end

% Fast Sine Transform on input matrix V 
% of dimension n-by-m
function Y = fast_sine_transform(V)
    [n,m]=size(V); 
    V1=[zeros(1,m);V;zeros(n+1,m)];
    V2=imag(fft(V1)); 
    % In Matlab vectors and matrices are indexed 
    % starting at 1, not 0
    Y=V2(2:n+1,:);
end

% Solve the discrete Poisson equation 
% on an n-by-n grid with right hand side b
function X = Poisson_FFT(b)
    [n,m]=size(b);
    % Form eigenvalues of matrix T(nxn)
    L=2*(1-cos((1:n)*pi/(n+1))); 
    % Form reciprocal sums of eigenvalues
    % Include scale factor 2/(n+1)
    LL=(2/(n+1))*ones(n,n)./(L'*ones(1,n)+ones(n,1)*L);
    % Solve, using the Fast Sine Transform 
    X = fast_sine_transform(b');
    X = fast_sine_transform(X');
    X = LL.*X;
    X = fast_sine_transform(X');
    X = fast_sine_transform(X');
end