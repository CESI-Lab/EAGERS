function [u]=MPCsolve(u,Ap,Bp,Cp,W,Np,Xf)
%% Augment matrices
%Definition of augmented model for MPC computation
m=size(Cp);
[~,ni]=size(Bp);
A = [Ap zeros(m(2),m(1)); Cp*Ap eye(m(1));];
B = [Bp;Cp*Bp;];

%Weight matrices and Laguerre's model description parameters
Q=W;                        %Weight matrix on states
R=1*eye(ni,ni);         %Weight matrix on input
a=0.9*ones(1,ni);         %Parameter for Laguerre computation
N=6*ones(1,ni);           %Parameter for Laguerre computation

%% Computate optimization matrices
% [Omega,Psi]=dmpc(A,B,a,N,Np,Q,R);
%A_e;B_e define the extended state-space model when integrator is used
% a contains the Laguerre pole locations for each input
%N the number of terms for each input; sum(N) is the dimension of eta
%Np prediction horizon
%Q weight on the state variables
%R weight on the input variables assumed to be diagonal.
% The cost function is J= eta ^T E eta +2 eta ^T H x(k_i)
[n,na]=size(B);
c=1;
Rdiag = zeros(sum(N),1);
for i=1:na;
    Rdiag(c:(c+N(i)-1),1)=R(i,i);
    c=c+N(i);
end
R_para = diag(Rdiag);

S_in=zeros(n,sum(N));
s=1;
for j=1:na;
    [~,L0]=lagd(a(j),N(j));
    S_in(:,s:s+N(j)-1)=B(:,j)*L0';
    s=s+N(j);
end
phi=S_in;
Omega=(phi)'*Q*(phi);
Psi=phi'*Q*A;
for i=2:Np;%calculate the finite sum S for each input (phi)
    c=1;
    for k=1:na;
        [Al,~]=lagd(a(k),N(k));%specify the state matrix Al
        phi(:,c:(c+N(k)-1))=A*phi(:,c:(c+N(k)-1))+S_in(:,c:(c+N(k)-1))*(Al^(i-1))';% Laguerre function
        c = c+N(k);
    end
    Omega=Omega+phi'*Q*phi;
    Psi=Psi+phi'*Q*(A^i);
end
Omega=Omega+R_para;

%% Compute Lzerot
Lzerot=zeros(ni,sum(N));
c=0;
for i=1:ni;
    [~,L0]=lagd(a(i),N(i));
    Lzerot(i,c+1:c+N(i))=L0';
    c=c+N(i);
end
%% Compute u
eta=-(Omega\Psi)*Xf;
deltau=Lzerot*eta;
u=u+deltau;