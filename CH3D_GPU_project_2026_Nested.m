% This code applies BHM method to BCP model in 3D - Nested Version
% Memory saving Code : Written on 4-1-2026 with no energy computation
function [Stability]=CH3D_GPU_project_2026_Nested(dt,M1,tfinal,N)
alpha=20;
a=0;b=M*pi;
%number of grid points N
%N=256/4;
%uniform mesh thickness
h=(b-a)/N;
%xgrid formtation (a b] and eventually (a b]^2
x=(a:h:b-h);
%[X,Y,Z] = meshgrid(x,x,x);
%wave number generation (same as in 1d)
k=[[0:N/2] [-N/2+1:-1]]./(M/2);
%k=single(k);      %------------------------------Single
%k=gpuArray(k);    %------------------------------GPU vs CPU   
%now in 2-d, here k2 means k^2  and k=(k1,k2) 
 [k1x k1y k1z]=meshgrid(k.^1,k.^1,k.^1); [kx ky kz]=meshgrid(k.^2,k.^2,k.^2);
 k2=kx+ky+kz; k4=k2.^2;
 k2=gpuArray(k2);k4=gpuArray(k4);
 clear kx ky kz k x
%Initial Condition----------------------
rng(1527,'twister'); 
U=0.01*rand(N,N,N)+.25;
%U=single(U);              %---------------single
U=gpuArray(U);             %---------------GPU vs CPU   3
%parameters
epsilon=.05*2; eps2=epsilon^2;
%The LHS, left hand side of problem----------
lhs=1+dt*M1*k4*eps2/epsilon+dt*alpha;   % CH lhs
%Left hand side-------------------------------
hat_U = fftn(U); 
it = 0; j=0; nn=0;   t = 0.0;
%BCP
TheUbar=((1/(2*pi))^3)*h*h*h*sum(sum(sum(U)));
%constant mobility 
M=1;
while (t <  tfinal-dt/2 )
U1 = U;  %Just Eyres scheme - no initial guess        
for i=1:1
      hat_U1 = (hat_U + dt.*fftn((eps2*(M1-M)*ifftn(k4.*fftn(U1))+2*ifftn(-1*k2.*fftn(U1.^3-U1)))/epsilon+alpha*(TheUbar-U1*0)))./lhs;
      U1 = real(ifftn(hat_U1));
end
U=U1;    %update 
hat_U=hat_U1;
it = it+1;
t = t+dt;
end  %main loop
save('U_numerical','U')
end
