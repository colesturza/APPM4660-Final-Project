
pseudoSpectral
%finiteDif
    %Interval 
%     a = 0;
%     b = 2*pi;
%     N = 100;
%     K's in the fft
%     k = (2*pi/(b-a)*[(-N/2):(N/2-1)]).';
%     k = fftshift(k);
%     dx = (b-a)/N;
%     x = (a + dx*(0:N-1)).';
%     
%     u0 = sin(x);
%     utt = real(ifft(1i*k.*fft(alpha(x).*ifft(1i*k.*fft(u0)))));
% 
%     plot(x, utt)
%     
%     NTime = 5000;
%     Time = 2*pi;
%     dt = Time/NTime;
%     
% 
% 
% 
%     uSol2 = zeros(N, NTime);
%     uSol2(:,1) = u0;
% 
%     uSol2(:, 2) = u0 + utt*dt^2 - 2*dt*phi(x); 
% 
%     for index = 3:NTime 
%         utt = real(ifft(1i*k.*fft(alpha(x).*ifft(1i*k.*fft(uSol2(:, index - 1))))));
%         uSol2(:, index) = 2*uSol2(:, index - 1) + utt*(dt^2) - uSol2(:, index - 2); 
%     end
% 
%     X = zeros(N, NTime);
%     T = zeros(N, NTime);
%     for i = 1:NTime
% 
%         X(:, i) = x;
%     end
%     for i = 1:N
%        T(i, :) = [0:dt:Time-dt]; 
%     end
% 
%   xlim([0 2*pi])
%   ylim([-1 1])
% 
%   mesh(X, T, uSol2);
%     xlabel('x')
%   ylabel('time (t)')
%   zlabel('u(x, t)')
 
  
  




function output = gaussianBump(x)
    output = exp(-(x-pi).^2/2);
end
function output = special(x)
    output = 1./(1 + (x-pi).^2);
end
function output = alpha(x)
    a = x(1:50);
    b = x(51:100);
    output = [a - a + 1; b - b + 1];
    
end
function output = phi(x)
    output = 0;
end
function output = triangle(x)
    a = x(1:50);
    b = x(51:100);
    output = [a/pi; (2*pi - b)/pi]; 
end



function finiteDif
    a = 0;
    b = 2*pi;
    N = 100;
    x = (b-a)/N*[1:N].';
    dx = (b-a)/N;
    A = sparse(N, N);
    for i = 1:N-1
        A(i, i) = -(alpha((i-1/2)*(b-a)/N) + alpha((i+1/2)*(b-a)/N));
        A(i + 1, i) = alpha((i + 1/2)*(b-a)/N);
        A(i, i + 1) = alpha((i + 1/2)*(b-a)/N);
    end
    A(N, N) = A(i, i);
    A(N, 1) = alpha((b-a)/(2*N));
    A(1, N) = A(N, 1);
    
    u0 = sin(x);
    
    
    utt = A*u0/(dx^2);
    
    NTime = 5000;
    Time = 2*pi;
    dt = Time/NTime;
    
    uSol = zeros(N, NTime);
    uSol(:,1) = u0;
    uSol(:, 2) = u0 + utt*dt^2 - 2*dt*phi(x); 
    
     for index = 3:NTime 
        utt = A*uSol(:, index - 1)/(dx^2);
        uSol(:, index) = 2*uSol(:, index - 1) + utt*dt^2 - uSol(:, index - 2); 
     end
     
    X = zeros(N, NTime);
    T = zeros(N, NTime);
    for i = 1:NTime

        X(:, i) = x;
    end
    for i = 1:N
       T(i, :) = [0:dt:Time-dt]; 
    end

  
  contour(X, T, uSol);
  xlim([0 2*pi])
  ylim([0 Time])
  xlabel('x')
  ylabel('time (t)')
  zlabel('u(x, t)')
    
    
%     h = figure;
%     axis tight manual
%     ax = gca;
%     ax.NextPlot = 'replaceChildren';
%     
%     M(NTime) = struct('cdata',[],'colormap',[]);
%     
%     for j = 1:NTime
%         xlim([0 2*pi])
%         ylim([-1 1])
%         plot(x, uSol(:, j));
%         drawnow
%         M(j) = getframe;
%     end
%     
%     h.Visible = 'on';
    
    
%         for j = 1:NTime
% 
%             % draw stuff
%             
%             
%             plot(x, uSol(:, j));
%             xlim([0 2*pi])
%             ylim([-1 1])
%             
%             frame = getframe(gcf);
%             im{j} = frame2im(frame);
%             [A,map] = rgb2ind(im{j},256);
%             if j == 1
%                 imwrite(A,map,'animation.gif','gif','LoopCount',Inf,'DelayTime',0);
%             else
%                 imwrite(A,map,'animation.gif','gif','WriteMode','append','DelayTime',0);
%             end
%         end
%     
end






function pseudoSpectral
    
    a = 0;
    b = 2*pi;
    N = 100;
    k = (2*pi/(b-a)*[(-N/2):(N/2-1)]).';
    k = fftshift(k);
    dx = (b-a)/N;
    x = (a + dx*(0:N-1)).';
    
    u0 = gaussianBump(x);
    utt = real(ifft(1i*k.*fft(alpha(x).*ifft(1i*k.*fft(u0)))));

    %plot(x, utt)
    
    NTime = 5000;
    Time = 2*pi;
    dt = Time/NTime;
    



    uSol = zeros(N, NTime);
    uSol(:,1) = u0;

    uSol(:, 2) = u0 + utt*dt^2 - 2*dt*phi(x); 

    for index = 3:NTime 
        utt = real(ifft(1i*k.*fft(alpha(x).*ifft(1i*k.*fft(uSol(:, index - 1))))));
        uSol(:, index) = 2*uSol(:, index - 1) + utt*(dt^2) - uSol(:, index - 2); 
    end

    X = zeros(N, NTime);
    T = zeros(N, NTime);
    for i = 1:NTime

        X(:, i) = x;
    end
    for i = 1:N
       T(i, :) = [0:dt:Time-dt]; 
    end

  xlim([0 2*pi])
  ylim([-1 1])
  xlabel('x')
  ylabel('time (t)')
  zlabel('u(x, t)')
  mesh(X, T, uSol);
 
  
  
  
  
  
% 
%     h = figure;
%     axis tight manual
%     ax = gca;
%     ax.NextPlot = 'replaceChildren';
%     
%     M(NTime) = struct('cdata',[],'colormap',[]);
%     
%     for j = 1:NTime
%         xlim([0 2*pi])
%         ylim([-1 1])
%         plot(x, uSol(:, j));
%         drawnow
%         M(j) = getframe;
%     end
%     
%     h.Visible = 'on';
    
%         for j = 1:NTime
% 
%             % draw stuff
%             
%             
%             plot(x, uSol(:, j));
%             xlim([-2*pi 2*pi])
%             ylim([-1 1])
%             
%             frame = getframe(gcf);
%             im{j} = frame2im(frame);
%             [A,map] = rgb2ind(im{j},256);
%             if j == 1
%                 imwrite(A,map,'animation.gif','gif','LoopCount',Inf,'DelayTime',0);
%             else
%                 imwrite(A,map,'animation.gif','gif','WriteMode','append','DelayTime',0);
%             end
%         end


end





