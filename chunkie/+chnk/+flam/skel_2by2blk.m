function [skel_struct] = skel_2by2blk(FS1,FS2,A12,A21,nmax,tol)

skel_struct = [];

tic; atmp = chnk.flam.rand_fft_transf(A12,nmax); toc;
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

atmp = chnk.flam.rand_fft_transf(A21,nmax);
[SK,RD,T] = id(atmp,tol);

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

%tic; V1 = transpose(transpose(B11)*transpose(T21)); toc;
V1 = transpose(rskelf_sv(FS1,transpose(T21),'t'));
%tic; U1 = transpose(transpose(B22)*transpose(T12)); toc;
U1 = transpose(rskelf_sv(FS2,transpose(T12),'t'));
U = U1*R21;
U = R12*U;
V = V1;


S1 = eye(size(V,1))-V*U;
S1 = S1\V;
V = S1;


%F11 = eye(size(U,1))+U*V;

%tic; U = B11*U; toc;
U = rskelf_sv(FS1,U);
%%%%%%%%%%%%  1st block done:
%%%%%%%%%%%%  ZZ11 = B11 + U*V;
skel_struct.U11 = U;
skel_struct.V11 = V;



W1 = -rskelf_sv(FS1,R12);
C2 = V*R12;
W2 = -U*C2;


W  = W1 + W2;

skel_struct.U12 = W;
skel_struct.V12 = U1;
%%%%%%%%%%%%  2nd block done:
%%%%%%%%%%%%  ZZ12 = W*U1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% creating F2 in factored form

%tic; V1 = transpose(transpose(B11)*transpose(T21)); toc;
%tic; U1 = transpose(transpose(A22)\transpose(T12)); toc;

P = V1*R12; 
P = R21*P; 
Q = U1;


K1 = eye(size(Q,1))-Q*P;

K1 = K1\Q;
Q = K1;


%F22 = eye(size(P,1))+P*Q;


%tic; P = B22*P; toc;
P = rskelf_sv(FS2,P); 

skel_struct.U22 = P;
skel_struct.V22 = Q;
%%%%%%%%%%%%  1st block done:
%%%%%%  ZZ22 = B22 + P*Q;



X1 = -rskelf_sv(FS2,R21);
C2 = Q*R21;
X2 = -P*C2;


X  = X1 + X2;

skel_struct.U21 = X;
skel_struct.V21 = V1;

end

