N = 250;

a = 2;
b = 2;

[x, y] = meshgrid(linspace(0, a, N), linspace(0, b, N));

% point charge location
x_src = [0.95 1.1];

points = [x(:), y(:)];
x_trg = points.';

% evaluate the solution at the target
uex = make_refsoln(x_trg, x_src);

uex = reshape(uex, [N, N]);

figure(1)
mesh(x, y, uex)

% evaluate solution for a point charge
function uex = make_refsoln(x_trg,x_src)
    N = size(x_trg,2);

    dd1 = (x_src(1).*ones(1, N)-x_trg(1, :)).^2;
    dd2 = (x_src(2).*ones(1, N)-x_trg(2, :)).^2;

    uex = -1/(2*pi)*log(sqrt(dd1+dd2)).';
end