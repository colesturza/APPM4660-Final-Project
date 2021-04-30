


%plotFundamentalSolution();

plotDoubleLayer();



function plotFundamentalSolution()
    x = -5:0.01:5;
    y = -5:0.01:5;
    [X, Y] = meshgrid(x, y);
    points = [X(:), Y(:)];
    x_src = [-1, -1; 1, 1; -1, 1; 1, -1; 0, -1; -1, 0; 1, 0; 0, 1;];
    Z = fSolTar(points, x_src);
    Z = reshape(Z, length(x), length(x));
    %contourf(Z, 'EdgeColor', 'interp' )
    image(points(:, 1), points(:, 2), Z,'CDataMapping','scaled')
    axis square
    title("Fundamental Solution of Laplace Equation")
end

function plotDoubleLayer()
    a = 5;
    x = -a:a*1e-3:a;
    y = -a:a*1e-3:a;
    [X, Y] = meshgrid(x, y);
    points = [X(:), Y(:)];
    Z = 10000*([1, 0]*points.')./(2*pi*vecnorm(points.').^2);
    Z = reshape(Z, length(x), length(x));
    F = zeros(length(x), length(x));
    for i = round(3*length(x)/4):round(1.25*length(x))
        Temp = [Z(mod(i, length(x))+1:end, :); Z((1:(mod(i, length(x)))), :)];
        F = F + Temp;
    end
    
    
    
    image(points(:, 1), points(:, 2), F,'CDataMapping','scaled')
    axis square;
    title("Fundamental Solution of Laplace Equation")
end


function output = fSolTar(x_trg, x_src)

    output = 0;
    for i = 1:length(x_src)
        output = output + sign(rand(1,1) - .5)*log(vecnorm((x_trg - x_src(i, :)).'));
    end
end
