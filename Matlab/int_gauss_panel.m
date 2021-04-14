function test_int_fixed2

% test integral equation discretization for Laplace 
% interior problem

%%% Set the geometry.
% flag_geom    = 'clover';
flag_geom = 'square';


% tt = linspace(0,2*pi);
% tt = tt(1:end-1);
N = 100;
%I want to use the gauss lengendre nodes, but scaled to go form 0 to 2pi


% list = createInterval([-.5, -.45, -.4, -.2, 0, .2, .4, .45, .5], 30);
% disp(length(list))
% [tt, w] = makeGaussianNodes(list, 10);

tt = 0:.01:4;
%%Otherwise this will end up counting the last point twice.
tt = tt(1:end-1);

%%Gaussian Panels

%[tt, w] = lgwt(20, 0, 4);
keyboard;
[C,ww] = make_geom(tt,flag_geom);
ww = ww*.01;

%ww = w.'.*ww;

keyboard;


%  validates arclength of circle
%abs(sum(ww)-2*pi)

x_src = [1.4;.5];



[x, y] = meshgrid(0.01:.01:1-.01, 0.01:.01:1-.01);

points = [x(:), y(:)];
keyboard;

x_trg = points.';


 A = LOCAL_construct_A_diag(C,ww);
keyboard;

g = make_bdry(C,x_src);

% solve for boundary charge density

sigma = A\g;

% evaluate the solution at the target
uex = make_refsoln(x_trg,x_src);

uapp = LOCAL_evalpot(x_trg,C,sigma,ww);

clf;
figure(1)
hold on;
scatter(C(1, :), C(4, :), 4, 'filled');
axis equal
disp(norm(uapp - uex))
mesh(x, y, reshape(uapp, length(x), length(x)))

keyboard
return
end

function vv = LOCAL_evalpot(x_trg,C,sigma,ww)

n  = size(C,2);
m  = size(x_trg,2);

X_g1  = x_trg(1,:)'*ones(1,n);
X_g2  = x_trg(2,:)'*ones(1,n);
Y_g1  = ones(m,1)*C(1,:);
Y_g2  = ones(m,1)*C(4,:);
Y_dg1 = ones(m,1)*C(2,:);
Y_dg2 = ones(m,1)*C(5,:);
ww    = ones(m,1)*ww;

EVAL = LOCAL_eval_kernel(X_g1,Y_g1,X_g2,Y_g2,Y_dg1,Y_dg2);
EVAL = EVAL.*ww;      
vv   = EVAL*sigma;

return
end

function list = createInterval(startList, numRefinement)
    N = length(startList);
    returnList = zeros(1, N+2);
    
    if numRefinement <= 0
       list = startList;
       return;
    end
    returnList(2) = (startList(1) + startList(2))/2;
    returnList(N+1) = (startList(N-1) + startList(N))/2;
    returnList(3:N) = startList(2:N-1);
    returnList(1) = startList(1);
    returnList(N+2) = startList(N);
    list = createInterval(returnList, numRefinement - 1);
end





function A = LOCAL_eval_kernel(X_g1,Y_g1,X_g2,Y_g2,Y_dg1,Y_dg2)

  nn1  = ( Y_dg2./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
  nn2  = (-Y_dg1./sqrt(Y_dg1.*Y_dg1 + Y_dg2.*Y_dg2));
  ddsq = (Y_g1 - X_g1).^2 + (Y_g2 - X_g2).^2;
  A    = -1/(2*pi)*(nn1.*(Y_g1 - X_g1) + nn2.*(Y_g2 - X_g2))./ddsq;

 
return
end

function uex = make_refsoln(x_trg,x_src)

N = size(x_trg,2);


dd1 = (x_src(1)-x_trg(1)).^2;
dd2 = (x_src(2)-x_trg(2)).^2;
uex = -1/(2*pi)*log(sqrt(dd1+dd2)).';

return

end

%Takes in a list, and spits out a bunch on nodes of order n:
function [tNodes, wNodes] = makeGaussianNodes(list, numNodes)
   Nnodes = length(list) - 1;
   nodesAtList = zeros(Nnodes);
   if Nnodes == length(numNodes)
       nodesAtList = numNodes;
   else
      nodesAtList(:) = numNodes(1,1); 
   end
    [tNodes, wNodes] = lgwt(nodesAtList(1), list(1), list(2));
   for i = 2:Nnodes
        [a, b] = lgwt(nodesAtList(i), list(i), list(i + 1));
        tNodes = [tNodes; a];
        wNodes = [wNodes; b];
   end
end

function g = make_bdry(C,x_src)

N = size(C,2);


dd1 = (x_src(1)-C(1,:)).^2;
dd2 = (x_src(2)-C(4,:)).^2;
g = -1/(2*pi)*log(sqrt(dd1+dd2)).';



return

end

function [C,ww] = make_geom(t,flag_geom)

% N = number of nodes
N = length(t);
C = zeros(6,N);

if(strcmp(flag_geom,'circle'))
%added second derivatives to calcluate the middle diagonal
C(1,:) = cos(t);
C(2,:) = -sin(t);
C(3,:) = -cos(t);
C(4,:) = sin(t);
C(5,:) = cos(t);
C(6,:) = -sin(t);

elseif(strcmp(flag_geom,'clover'))

C(1,:) = (cos(3*t) +2).*cos(t);
C(2,:) = -(3*sin(3*t).*cos(t)+(cos(3*t)+2).*sin(t));
C(4,:) = (cos(3*t)+2).*sin(t);
C(5,:) = (-3*sin(3*t).*sin(t)+(cos(3*t)+2).*cos(t));


elseif(strcmp(flag_geom,'unit-square'))
X = 2;
Y = 2;
x0 = 0;
y0 = 0;
    
    
C(1,:) = cos(t).*sec(t - pi/2*floor((4*t + pi)/(2*pi)))*X/2 + x0;
C(4,:) = sin(t).*sec(t - pi/2*floor((4*t + pi)/(2*pi)))*Y/2 + y0;

elseif(strcmp(flag_geom,'square'))
    
C(1, :) = [t(0 <= t & t < 1), ones(1, length(t(1 <= t & t < 2))), 3 - t(2 <= t & t < 3), zeros(1, length(t(3 <= t & t <= 4)))]; 
C(2, :) = [ones(1, length(t(0 <= t & t < 1))), zeros(1, length(t(1 <= t & t < 2))), -ones(1, length(t(2 <= t & t < 3))), zeros(1, length(t(3 <= t & t <= 4)))];
C(3, :) = zeros(1, length(C(1, :)));
C(4, :) = [zeros(1, length(t(0 <= t & t < 1))), t(1 <= t & t < 2) - 1, ones(1, length(t(2 <= t & t < 3))), 4 - t(3 <= t & t <= 4)];
C(5, :) = [zeros(1, length(t(0 <= t & t < 1))), ones(1, length(t(1 <= t & t < 2))), zeros(1, length(t(2 <= t & t < 3))), -ones(1, length(t(3 <= t & t <= 4)))];
C(6, :) = zeros(1, length(C(1, :)));
elseif(strcmp(flag_geom,'star'))
% starfish
tt = t;
r             = 0.5;  % Parameter controlling the discretization.
k             = 5;    % Parameter controlling the discretization. 

C(1,:)        =   1.5*cos(tt) + (r/2)*            cos((k+1)*tt) + (r/2)*            cos((k-1)*tt);
C(2,:)        = - 1.5*sin(tt) - (r/2)*(k+1)*      sin((k+1)*tt) - (r/2)*(k-1)*      sin((k-1)*tt);
C(3,:)        = - 1.5*cos(tt) - (r/2)*(k+1)*(k+1)*cos((k+1)*tt) - (r/2)*(k-1)*(k-1)*cos((k-1)*tt);
C(4,:)        =       sin(tt) + (r/2)*            sin((k+1)*tt) - (r/2)*            sin((k-1)*tt);
C(5,:)        =       cos(tt) + (r/2)*(k+1)*      cos((k+1)*tt) - (r/2)*(k-1)*      cos((k-1)*tt);
C(6,:)        = -     sin(tt) - (r/2)*(k+1)*(k+1)*sin((k+1)*tt) + (r/2)*(k-1)*(k-1)*sin((k-1)*tt);

elseif(strcmp(flag_geom,'pacman'))

  C(1,:)        = -     sign(-pi+2*pi*t).*sin(1.5*(-pi+2*pi*t)); % x-coordinate x(t)
  C(2,:)        = -  1.5*sign(-pi+2*pi*t).*cos(1.5*(-pi+2*pi*t))*2*pi; % x'(t)
  C(3,:)        =  2.25*sign(-pi+2*pi*t).*sin(1.5*(-pi+2*pi*t))*4*pi^2; % x''(t)
  C(4,:)        =     tan(3*pi/4)*sin(-pi+2*pi*t); % y(t)
  C(5,:)        =     tan(3*pi/4)*cos(-pi+2*pi*t)*2*pi; % y'(t)
  C(6,:)        = -   tan(3*pi/4)*sin(-pi+2*pi*t)*4*pi^2; %y''(t)

elseif(strcmp(flag_geom,'tear'))

  C(1,:)        =      2*sign(-pi+2*pi*t).*sin((-pi+2*pi*t)/2);
  C(2,:)        =      2*pi*sign(-pi+2*pi*t).*cos((-pi+2*pi*t)/2);
  C(3,:)        =      (2*pi)^2*-(1/2)*sign(-pi+2*pi*t).*sin((-pi+2*pi*t)/2);
  C(4,:)        =     -tan(pi/4)*sin(-pi+2*pi*t);
  C(5,:)        =     -2*pi*tan(pi/4)*cos(-pi+2*pi*t);
  C(6,:)        =     (2*pi)^2*tan(pi/4)*sin(-pi+2*pi*t);



else
  fprintf(1,'This option for the geometry is not implemented.\n');
  keyboard

end

ww = sqrt(C(2,:).*C(2,:) + C(5,:).*C(5,:));
% arclength = sqrt(C(2,:).*C(2,:) + C(5,:).*C(5,:));


return

end


function B = LOCAL_construct_A_diag(C,ww)

N = size(C,2);
dd1   = C(1,:)' * ones(1,N) - ones(N,1) * C(1,:);
dd2   = C(4,:)' * ones(1,N) - ones(N,1) * C(4,:);
nn1   = ones(N,1)*(C(5,:)./sqrt(C(2,:).*C(2,:) + C(5,:).*C(5,:)));
nn2   = ones(N,1)*(-C(2,:)./sqrt(C(2,:).*C(2,:) + C(5,:).*C(5,:)));
ddsq  = (dd1.*dd1 + dd2.*dd2) + eye(N);
B     = (1/(2*pi))*(nn1.*dd1 + nn2.*dd2)./ddsq;
d_vec = (-C(6,:).*C(2,:) + C(3,:).*C(5,:))./...
        ((4*pi)*(C(2,:).^2 + C(5,:).^2).^1.5); %rmk 12.2
B     = B + diag(d_vec);
B     = B.*(ones(N,1)*ww) - 0.5*eye(N);
return

end