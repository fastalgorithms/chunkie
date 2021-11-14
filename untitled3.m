nout = 400;
tol  = 10^(-12);
tic; atmp = rand_fft_transf(A12,nout); toc;
tic; [SK,RD,T] = id(atmp,tol); toc;

nskel = numel(SK);
Ttot = zeros(nskel,size(A12,2));
Ttot(1:nskel,SK) = eye(nskel);
Ttot(1:nskel,RD) = T;
T12 = Ttot;

skel12 = SK;
R12 = A12(:,skel12);

%%%%%%%%%%%%%%%%%%%%%%%%
%% now the 2,1 block

tic; atmp = rand_fft_transf(A21,nout); toc;
tic; [SK,RD,T] = id(atmp,tol); toc;

nskel = numel(SK);
Ttot = zeros(nskel,size(A21,2));
Ttot(1:nskel,SK) = eye(nskel);
Ttot(1:nskel,RD) = T;
T21 = Ttot;

skel21 = SK;
R21 = A21(:,SK);

%norm(A21(:,SK)*T21 - A21)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% creating F1 in factored form

tic; V1 = transpose(transpose(A11)\transpose(T21)); toc;
tic; U1 = transpose(transpose(A22)\transpose(T12)); toc;
tic; U = U1*R21; toc;
tic; U = R12*U; toc;
U = U;
V = V1;

tic;
S1 = eye(size(V,1))-V*U;
S1 = S1\V;
V = S1;
toc;

F11 = eye(size(U,1))+U*V;

tic; U = A11\U; toc;

%%%%%%%%%%%%  1st block done:
%%%%%%%%%%%%  ZZ11 = B11 + U*V;
tic;
W1 = -(A11\R12);
C2 = V*R12;
W2 = -U*C2;
toc;

W  = W1 + W2;

%%%%%%%%%%%%  2nd block done:
%%%%%%%%%%%%  ZZ12 = W*U1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% creating F2 in factored form

tic; V1 = transpose(transpose(A11)\transpose(T21)); toc;
tic; U1 = transpose(transpose(A22)\transpose(T12)); toc;

tic; P = V1*R12; toc;
tic; P = R21*P; toc;
Q = U1;

tic;
K1 = eye(size(Q,1))-Q*P;

K1 = K1\Q;
Q = K1;
toc;

F22 = eye(size(P,1))+P*Q;


tic; P = A22\P; toc;

%%%%%%%%%%%%  1st block done:
%%%%%%  ZZ22 = B22 + P*Q;


tic;
X1 = -(A22\R21);
C2 = Q*R21;
X2 = -P*C2;
toc;

X  = X1 + X2;

%%%%%%%%%%%%  2nd block done:
%%%%  ZZ21 = X*V1;



%F11 = eye(size(U1,1))-U1*V1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


