% slow way (works)
N = 5;
K = diag(2*ones(N,1))-diag(ones(N-1,1),-1)-diag(ones(N-1,1),1);
figure(1)
hold on
spy(K==2,'g')
spy(K==-1,'b')
hold off
K2D = kron(eye(N), K) + kron(K, eye(N));
figure(2)
hold on
spy(K2D==4,'r')
spy(K2D==-1,'b')
hold off