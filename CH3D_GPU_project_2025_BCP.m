%This code applies BHM method to the BCP equation in 3D
function [Stability]=CH3D_GPU_project_2026_BCP(dt,M1,tfinal,N)

M=2;
a=0;
b=M*pi;
%number of grid points N
h=(b-a)/N;

%xgrid formtation (a b] and eventually (a b]^2
x=(a:h:b-h);
x=gpuArray(x);       %-----------------------%GPU vs CPU      1
%x=single(x);

[X,Y,Z] = meshgrid(x,x,x);

%wave number generation (same as in 1D)
k=[[0:N/2] [-N/2+1:-1]]./(M/2);
k=gpuArray(k);       %-----------------------%GPU vs CPU     2
%k=single(k);

%now in 2-d, here k2 means k^2  and k=(k1,k2) 
 [k1x k1y k1z]=meshgrid(k.^1,k.^1,k.^1);
 [kx ky kz]=meshgrid(k.^2,k.^2,k.^2);
 k2=kx+ky+kz;k4=k2.^2;

%Initial Condition----------------------
U=0.01*rand(N,N,N)+.5*0;
U=gpuArray(U);        %------------------------%GPU vs CPU    3
%U=single(U);

%parameters
epsilon=.1; 
eps2=epsilon^2;

%The LHS, left hand side of problem----------
lhs=1+dt*M1*k4*eps2/epsilon+dt*15;   % CH LHS
hat_U = fftn(U); 
it = 0; t = 0.0;

M=1;
%BCP avg value
TheUbar=((1/(2*pi))^3)*h*h*h*sum(sum(sum(U)));
while (t <  tfinal-dt/2)

U1 = U;         
for i=1:1
      RHS=eps2*(M1-M)*ifftn(k4.*fftn(U1))+2*ifftn(-1*k2.*fftn(U1.^3-U1));
      RHS=RHS/epsilon+15*(TheUbar-U1*0);  %BCP ?
      hat_rhs =hat_U + dt.*fftn(RHS);
      hat_U1 = hat_rhs./lhs;
      U1 = real(ifftn(hat_U1));
end
U=U1;    %update 
hat_U=hat_U1;
it = it+1;
t = t+dt;  
end  %main loop
save('U_numerical','U');
end
