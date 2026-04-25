%This code applies CH equation in 3D
%4-8-2026   with no energy computation

function [Stability]=CH3D_GPU_project_2025_BCP0(dt,M1,iter,tfinal,N)
%tic; 
clear Energy time volume
%clc;
Stability=1;

alpha=20;


M=2;
a=0;
b=M*pi;
% for proper spacing   [0 2pi]
%number of grid points N
%N=256/4;
%uniform mesh thickness
h=(b-a)/N;
%(Periodic bdy conditions)
n=N;

%xgrid formtation (a b] and eventually (a b]^2
x=(a:h:b-h);
%x=single(x);   %------------------------------Single
%x=gpuArray(x);       %-----------------------%GPU vs CPU      1


%x=gpuArray([a:h:b-h]);
%x=single(x);

%[X,Y,Z] = meshgrid(x,x,x);

%wave number generation (same as in 1d)
k=[[0:N/2] [-N/2+1:-1]]./(M/2);
%k=single(k);      %------------------------------Single
%k=gpuArray(k); %-------------------------------GPU vs CPU     2


%k=gpuArray([[0:N/2] [-N/2+1:-1]]./((M+4*0)/2));
%k=single(k);

%now in 2-d, here k2 means k^2  and k=(k1,k2) 
 [k1x k1y k1z]=meshgrid(k.^1,k.^1,k.^1);
 
 [kx ky kz]=meshgrid(k.^2,k.^2,k.^2);
 k2=kx+ky+kz;
 k4=k2.^2;

 k2=gpuArray(k2);
 k4=gpuArray(k4);
 clear kx ky kz k x

%Initial Condition----------------------
 %U=0.5*rand(N,N,N)+.4*0;
 %U=0.01*rand(N,N,N)+.35;
 %save('U_init','U');          %to use the same I.C just load init
 %load('U_init128','U');
 %load('U_init','U');
rng(1527,'twister'); 
%U=0.01*rand(N,N,N,'gpuArray')+.0;
%U=0.05*rand(N,N,N)+.1;  %single vs double
U=0.01*rand(N,N,N)+.25;
%U=single(U);              %single
U=gpuArray(U);          %-------------------------GPU vs CPU   3

%understanding mass
InitialSum = gather(sum(U(:)))*h^3/(2*pi)^3;

%U=single(U);

%2drops
% % U=2*ones(N,N,N)+(-1.0)*((b/9)^2<(X-b/2.8).^2+(Y-b/2).^2+(Z-b/2).^2)+...
% %     (-1.0)*((b/9)^2<(X-b/1.7).^2+(Y-b/2).^2+(Z-b/2).^2);
% % U(U==0)=-1;
% % U=-U*.99;
% % U=gpuArray(U);  

% % figure(1)
% % isosurface(X,Y,Z,U)
% % isosurface(X,Y,Z,U,-1)
% % isosurface(X,Y,Z,U,1)
% % ax = gca; 
% % ax.FontSize = 14;
% % camlight; lighting phong
% % axis([a b a b a b])
% % %whos U

%parameters
epsilon=.05*2; 
eps2=epsilon^2;

%The LHS, left hand side of problem----------
lhs=1+dt*M1*k4*eps2/epsilon+dt*alpha;   % CH lhs
%Left hand side-------------------------------
hat_U = fftn(U); 
it = 0; j=0; nn=0;   t = 0.0;
%Energy Computations
% % Ue1=real(ifftn(-1i*k1x.*fftn(U)));
% % Ue2=real(ifftn(-1i*k1y.*fftn(U)));
% % Ue3=real(ifftn(-1i*k1z.*fftn(U)));
% % energy=eps2/2*(Ue1.^2+Ue2.^2+Ue3.^2)+(1/4)*U.^4-(1/2)*U.^2;


% % energy=eps2/2*( real(ifftn(-1i*k1x.*fftn(U))).^2+...
% %     real(ifftn(-1i*k1y.*fftn(U))).^2+...
% %     real(ifftn(-1i*k1z.*fftn(U))).^2)+...
% %     (1/4)*U.^4-(1/2)*U.^2;

% BCP nonlocal energy term: 
%BCP
TheUbar=((1/(2*pi))^3)*h*h*h*sum(sum(sum(U)));
%(eps*alpha/2) * integral |grad psi|^2 % Poisson solve: -Delta*psi = u - ubar => psi_hat = (u-ubar)_hat / k2 

%k2_safe = k2 + (k2==0); % avoid division by zero at k=0 

%k2_safe=gpuArray(k2_safe);  %-----------------------------GPUarray4

%psi_hat = fftn(U - TheUbar) ./ (k2 + (k2==0));


% grad psi components -- same style as Ue1, Ue2, Ue3 
% % psi_x = real(ifftn(-1i*k1x.*psi_hat)); 
% % psi_y = real(ifftn(-1i*k1y.*psi_hat));
% % psi_z = real(ifftn(-1i*k1z.*psi_hat)); 

% Add nonlocal term to energy 

% energy = eps2/2*(Ue1.^2+Ue2.^2+Ue3.^2) + (1/4)*U.^4 - (1/2)*U.^2 ... 
%     + (epsilon*15/2)*(psi_x.^2+psi_y.^2+psi_z.^2); 

% energy = (epsilon/2)*(Ue1.^2+Ue2.^2+Ue3.^2) + ...
%     (2/epsilon)*((1/4)*U.^4-(1/2)*U.^2) + ...
%     (alpha/2)*(psi_x.^2+psi_y.^2+psi_z.^2);

% energy = (epsilon/2)*((real(ifftn(-1i*k1x.*fftn(U)))).^2+...
%     (real(ifftn(-1i*k1y.*fftn(U)))).^2+...
%     (real(ifftn(-1i*k1z.*fftn(U)))).^2) +...
%     (2/epsilon)*((1/4)*U.^4-(1/2)*U.^2) + ...
%     (alpha/2)*(real(ifftn(-1i*k1x.*psi_hat)).^2 + ... 
%     real(ifftn(-1i*k1y.*psi_hat)).^2 + ...
%     real(ifftn(-1i*k1z.*psi_hat)).^2); 



%clear psi_hat % free ~1GB after each step! ✓

% Energy(it+1) = h*h*h*sum(sum(sum(energy)));




%%deltaU = ifftn(-1*k2.*fftn(U));  % ΔU in physical space
%%Menergy = energy + (eps2*(M1-1)/2)*(deltaU.^2);  removed for 512
%-Energy(1)=h*h*h*sum(sum(sum(energy)));   %4/2
%%MEnergy(1)=h*h*h*sum(sum(sum(Menergy)));      removed for 512
%time(1)=0;

%volume(1)=h*h*h*sum(sum(sum(U)));
%Energy=gpuArray(Energy);time=gpuArray(time);
%constant mobility 
M=1;

while (t <  tfinal-dt*1*0-.0000001 )

U1 = U;  %Just Eyres scheme - no initial guess        
for i=1:1
      
      %RHS=(eps2*(M1-M)*ifftn(k4.*fftn(U1))+2*ifftn(-1*k2.*fftn(U1.^3-U1)))/epsilon+alpha*(TheUbar-U1*0);

      %isgpuarray(RHS)
      %%RHS=RHS/epsilon+15*(TheUbar-U1*0);  %BCP ?
      %hat_rhs =hat_U + dt.*fftn(RHS);
      %hat_U1 = hat_rhs./lhs;

      %hat_U1 = (hat_U + dt.*fftn(RHS))./lhs;

      hat_U1 = (hat_U + dt.*fftn((eps2*(M1-M)*ifftn(k4.*fftn(U1))+2*ifftn(-1*k2.*fftn(U1.^3-U1)))/epsilon+alpha*(TheUbar-U1*0)))./lhs;
      U1 = real(ifftn(hat_U1));
end
%isgpuarray(U1)
U=U1;    %update 
hat_U=hat_U1;
it = it+1;
t = t+dt;

%Energy Computations
% Ue1=real(ifftn(-1i*k1x.*fftn(U)));
% Ue2=real(ifftn(-1i*k1y.*fftn(U)));
% Ue3=real(ifftn(-1i*k1z.*fftn(U)));

%time(it+1)=t;
% energy=eps2/2*( real(ifftn(-1i*k1x.*fftn(U))).^2+...
%     real(ifftn(-1i*k1y.*fftn(U))).^2+...
%     real(ifftn(-1i*k1z.*fftn(U))).^2)+...
%     (1/4)*U.^4-(1/2)*U.^2;


%psi_hat = fftn(U - TheUbar) ./ k2_safe; 
% grad psi components -- same style as Ue1, Ue2, Ue3 
% psi_x = real(ifftn(-1i*k1x.*psi_hat)); 
% psi_y = real(ifftn(-1i*k1y.*psi_hat));
% psi_z = real(ifftn(-1i*k1z.*psi_hat)); 

% Add nonlocal term to energy 

% energy = eps2/2*(Ue1.^2+Ue2.^2+Ue3.^2) + (1/4)*U.^4 - (1/2)*U.^2 ... 
%     + (epsilon*15/2)*(psi_x.^2+psi_y.^2+psi_z.^2); 

% energy = (epsilon/2)*(Ue1.^2+Ue2.^2+Ue3.^2) + ...
%     (2/epsilon)*((1/4)*U.^4-(1/2)*U.^2) + ...
%     (alpha/2)*(psi_x.^2+psi_y.^2+psi_z.^2);

% % psi_hat = fftn(U - TheUbar) ./ (k2 + (k2==0)); 
% % energy = (epsilon/2)*((real(ifftn(-1i*k1x.*fftn(U)))).^2+...
% %     (real(ifftn(-1i*k1y.*fftn(U)))).^2+...
% %     (real(ifftn(-1i*k1z.*fftn(U)))).^2) +...
% %     (2/epsilon)*((1/4)*U.^4-(1/2)*U.^2) + ...
% %     (alpha/2)*(real(ifftn(-1i*k1x.*psi_hat)).^2 + ... 
% %     real(ifftn(-1i*k1y.*psi_hat)).^2 + ...
% %     real(ifftn(-1i*k1z.*psi_hat)).^2); 


    %clear psi_hat % free ~1GB after each step! ✓

% % 
% % Energy(it+1) = h*h*h*sum(sum(sum(energy)));



%%deltaU = ifftn(-1*k2.*fftn(U));  % ΔU in physical space
%Menergy = energy + (eps2*(M1-1)/2)*(deltaU.^2);  512
%Energy(it+1)=h*h*h*sum(sum(sum(energy))); 
%MEnergy(it+1)=h*h*h*sum(sum(sum(Menergy)));        % 512

%volume(it+1)=h*h*h*sum(sum(sum(U)));

% figure(3);                  %if monitoring the solution is needed.
% mesh(X,Y,U)
% title(['Uexact time = ' num2str(t), ' dt = ' num2str(dt)] ); 

% % % figure(4); 
% % % mesh(X,Y, real(U) )%, shading interp, %axis('off'), %axis('equal');
% % % %pcolor(X,Y, U ), shading interp, axis('off'), axis('equal');
% % % title(['Eyre CH time = ' num2str(t), ' dt = ' num2str(dt)] ); 

%if Energy(it+1)>Energy(it)
%    Stability=0;
%    break                   % This break will stop the code if energy increases - unstable
%end
   
end  %main loop
save('U_numerical','U');

%%MEnergy=MEnergy/(.04*MEnergy(1))-55; %to normalize the modified energy

% % figure(1)                 %------------- Figure
% % %plot(time,Energy,'.-')
% % %semilogx(time,Energy,'b.-',time,MEnergy,'r.-')
% % semilogx(time,Energy,'b-','LineWidth',2)
% % %legend('E(t)','Modified E(t)')
% % legend('Energy(t)')
% % ax = gca; 
% % ax.FontSize = 14;
% % %axis([time(1) time(end) min(Energy) max(Energy)])
% % xlabel('time')
% % ylabel('Free Energy')

% % figure(111)
% % plot(time,Energy,'b.-',time,MEnergy,'r.-')
% % legend('E(t)','Modified E(t)')
% % %semilogx(time,Energy,'.-')
% % ax = gca; 
% % ax.FontSize = 14;
% % %axis([time(1) time(end) min(Energy) max(Energy)])
% % xlabel('time')
% % ylabel('Free Energy')

% a=toc;
% minutes=a/60;
% hours=a/60^2;
% minutes_hours=[minutes hours]

%save('U_init128','U'); % to avoid transient and spikes

U = gather(U);
FinalSum = sum(U(:))*h^3/(2*pi)^3;

MassError = abs(InitialSum - FinalSum);

fprintf('\n==========================================\n');
fprintf('     REPORTE DE VICTORIA (N=%d)           \n', N);
fprintf('==========================================\n');
fprintf('Masa Final:   %.16e\n', FinalSum);
fprintf('Error de Masa: %.16e\n', MassError);
fprintf('==========================================\n');


end
