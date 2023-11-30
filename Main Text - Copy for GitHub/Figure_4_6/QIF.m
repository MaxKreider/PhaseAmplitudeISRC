%% setup

%mr clean
clc
clf

%parameters
tau_m = 10;
Delta = 0.3;
J = 21;
Theta = 4;
tau_d = 5;

%for the ODEs
format long
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

%time
N = 2048*8+1;
Tmax = 27.58055;
dt = Tmax/(N-1);
t = 0:dt:Tmax;

%Jacobian
Jac = @(t,u) [2*u(1)/tau_m, -2*pi^2*tau_m*u(2), -J; ...
              2*u(2)/tau_m, 2*u(1)/tau_m, 0; ...
              0, 1/tau_d, -1/tau_d];


%% dynamics

%function
F = @(t,u) [1/tau_m*(u(1)^2 - (pi*tau_m*u(2))^2 - J*tau_m*u(3) + Theta); ...
            1/tau_m*(Delta/(pi*tau_m)+2*u(2)*u(1)); ...
            1/tau_d*(-u(3)+u(2))];

%solve it
x0 = [2.28564   0.0688145   0.0239971];
[Tc,U] = ode113(F,t,x0,opts);

%find the period
T = Tmax;


%% define phase and K_00(theta)

%K_00(theta)
K_00 = [U(1:end,1)'; U(1:end,2)'; U(1:end,3)'];


%% variational equation

%system
f_var = @(t,u) [1/tau_m*(u(1)^2 - (pi*tau_m*u(2))^2 - J*tau_m*u(3) + Theta); ...
                1/tau_m*(Delta/(pi*tau_m)+2*u(2)*u(1)); ...
                1/tau_d*(-u(3)+u(2)); ...
                (2*u(1)/tau_m)*u(4) + (-2*pi^2*tau_m*u(2))*u(7) + (-J)*u(10); ...
                (2*u(1)/tau_m)*u(5) + (-2*pi^2*tau_m*u(2))*u(8) + (-J)*u(11); ...
                (2*u(1)/tau_m)*u(6) + (-2*pi^2*tau_m*u(2))*u(9) + (-J)*u(12); ...
                (2*u(2)/tau_m)*u(4) + (2*u(1)/tau_m)*u(7) + (0)*u(10); ...
                (2*u(2)/tau_m)*u(5) + (2*u(1)/tau_m)*u(8) + (0)*u(11); ...
                (2*u(2)/tau_m)*u(6) + (2*u(1)/tau_m)*u(9) + (0)*u(12); ...
                (0)*u(4) + (1/tau_d)*u(7) + (-1/tau_d)*u(10); ...
                (0)*u(5) + (1/tau_d)*u(8) + (-1/tau_d)*u(11); ...
                (0)*u(6) + (1/tau_d)*u(9) + (-1/tau_d)*u(12)];

%solve it
[t,u] = ode113(f_var,t,[x0 1 0 0 0 1 0 0 0 1],opts);

%construct M
M = zeros(3,3);
M(1,1) = u(end,4);
M(1,2) = u(end,5);
M(1,3) = u(end,6);
M(2,1) = u(end,7);
M(2,2) = u(end,8);
M(2,3) = u(end,9);
M(3,1) = u(end,10);
M(3,2) = u(end,11);
M(3,3) = u(end,12);

%diagonalize it
[V,lambda] = eig(M);

%obtain eigenvalues mu and then exponents lambda
mu_1 = lambda(2,2);
mu_2 = lambda(3,3);
lambda_1 = 1/T*log(mu_1);
lambda_2 = 1/T*log(mu_2);

%obtain corresponding eigenvectors
v_1 = -V(:,2);
v_2 = V(:,3);


%% Floquet Matrices and Stuff (C and R and Q)

%diagonalize M (actually we already did it xD)
C = [V(:,1), -V(:,2), V(:,3)];

%define R
R = 1/T*C*[0 0 0; 0 log(mu_1) 0; 0 0 log(mu_2)]/C;

%compute Q(t) for each t
PHI=zeros(3,3);
Q = zeros(3,3,length(t));

for i=1:length(t)
    
    %fundamental matrix
    PHI(1,1) = u(i,4);
    PHI(1,2) = u(i,5);
    PHI(1,3) = u(i,6);
    PHI(2,1) = u(i,7);
    PHI(2,2) = u(i,8);
    PHI(2,3) = u(i,9);
    PHI(3,1) = u(i,10);
    PHI(3,2) = u(i,11);
    PHI(3,3) = u(i,12);
    
    %inner matrix
    D = [1 0 0; 0 exp(-t(i)*lambda_1) 0; 0 0 exp(-t(i)*lambda_2)];
    
    %compute
    Q(:,:,i) = (PHI*C*D)/C;

end


%% compute K_10 and K_01

%initialize
PHI = zeros(3,3);
K_10 = zeros(3,length(t));
K_01 = zeros(3,length(t));

%define numerical nicety coefficients
b1 = .2;
b2 = 1;

%loop
for i=1:length(t)
    
    %fundamental matrix
    PHI(1,1) = u(i,4);
    PHI(1,2) = u(i,5);
    PHI(1,3) = u(i,6);
    PHI(2,1) = u(i,7);
    PHI(2,2) = u(i,8);
    PHI(2,3) = u(i,9);
    PHI(3,1) = u(i,10);
    PHI(3,2) = u(i,11);
    PHI(3,3) = u(i,12);
    
    %compute K_10
    K_10(:,i) = PHI*exp(-lambda_2*t(i))*v_2*b1;
    
    %compute K_01
    K_01(:,i) = PHI*exp(-lambda_1*t(i))*v_1*b2;
    
end


%% automatic differentiation ... m = 2

%define order
L = 10;

%organize our (known) data
coeff_x = zeros(L,L,length(t));
coeff_y = zeros(L,L,length(t));
coeff_z = zeros(L,L,length(t));

%initialize with known values
coeff_x(1,1,:) = K_00(1,:);
coeff_y(1,1,:) = K_00(2,:);
coeff_z(1,1,:) = K_00(3,:);

coeff_x(1,2,:) = K_01(1,:);
coeff_y(1,2,:) = K_01(2,:);
coeff_z(1,2,:) = K_01(3,:);

coeff_x(2,1,:) = K_10(1,:);
coeff_y(2,1,:) = K_10(2,:);
coeff_z(2,1,:) = K_10(3,:);

%organize unknown data structure
comp_coeff_x = zeros(L,L,length(t));
comp_coeff_y = zeros(L,L,length(t));
comp_coeff_z = zeros(L,L,length(t));

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   z component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it
    temp = zeros(size(comp_coeff_y(:,:,1)));
    for m=2:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp(alpha+1,m-alpha+1) = temp(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    
    %other stuff
    comp_coeff_y(:,:,i) = 2/tau_m*temp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %do it (x^2)
    temp1 = zeros(size(comp_coeff_x(:,:,1)));
    for m=2:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp1(alpha+1,m-alpha+1) = temp1(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_x(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
            
    %do it (y^2)
    temp2 = zeros(size(comp_coeff_x(:,:,1)));
    for m=2:L-1
        for alpha=0:m

            for k=0:alpha
                for ell=0:m-alpha
                    temp2(alpha+1,m-alpha+1) = temp2(alpha+1,m-alpha+1) + coeff_y(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    temp2 = temp2*(-pi^2*tau_m^2);
        
    %all together
    comp_coeff_x(:,:,i) = (temp1 + temp2)/tau_m;
        
    %reset
    temp=[];
    temp1=[];
    temp2=[];
    
end

%beta's (for m=2)
beta_02 = [squeeze(comp_coeff_x(1,3,:))'; squeeze(comp_coeff_y(1,3,:))'; squeeze(comp_coeff_z(1,3,:))'];
beta_11 = [squeeze(comp_coeff_x(2,2,:))'; squeeze(comp_coeff_y(2,2,:))'; squeeze(comp_coeff_z(2,2,:))'];
beta_20 = [squeeze(comp_coeff_x(3,1,:))'; squeeze(comp_coeff_y(3,1,:))'; squeeze(comp_coeff_z(3,1,:))'];


%% Solve the Homological Equations with Floquet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_02 = zeros(3,length(t));
A_11 = zeros(3,length(t));
A_20 = zeros(3,length(t));

%solve
for i=1:length(t)
    A_02(:,i) = C\((Q(:,:,i))\beta_02(:,i));
    A_11(:,i) = C\((Q(:,:,i))\beta_11(:,i));
    A_20(:,i) = C\((Q(:,:,i))\beta_20(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_02_f = [fftshift(fft(A_02(1,1:end))); fftshift(fft(A_02(2,1:end))); fftshift(fft(A_02(3,1:end)))];
A_11_f = [fftshift(fft(A_11(1,1:end))); fftshift(fft(A_11(2,1:end))); fftshift(fft(A_11(3,1:end)))];
A_20_f = [fftshift(fft(A_20(1,1:end))); fftshift(fft(A_20(2,1:end))); fftshift(fft(A_20(3,1:end)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_02 = zeros(3,length(t));
u_11 = zeros(3,length(t));
u_20 = zeros(3,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_02(1,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(2-0)*lambda_1-0).*A_02_f(1,:);
u_02(2,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(2-0)*lambda_1-lambda_1).*A_02_f(2,:);
u_02(3,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(2-0)*lambda_1-lambda_2).*A_02_f(3,:);

u_11(1,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(2-1)*lambda_1-0).*A_11_f(1,:);
u_11(2,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(2-1)*lambda_1-lambda_1).*A_11_f(2,:);
u_11(3,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(2-1)*lambda_1-lambda_2).*A_11_f(3,:);

u_20(1,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(2-2)*lambda_1-0).*A_20_f(1,:);
u_20(2,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(2-2)*lambda_1-lambda_1).*A_20_f(2,:);
u_20(3,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(2-2)*lambda_1-lambda_2).*A_20_f(3,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_02 = ([ifft(ifftshift(u_02(1,:))); ifft(ifftshift(u_02(2,:))); ifft(ifftshift(u_02(3,:)))]);
u_11 = ([ifft(ifftshift(u_11(1,:))); ifft(ifftshift(u_11(2,:))); ifft(ifftshift(u_11(3,:)))]);
u_20 = ([ifft(ifftshift(u_20(1,:))); ifft(ifftshift(u_20(2,:))); ifft(ifftshift(u_20(3,:)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_02 = zeros(3,length(t));
K_11 = zeros(3,length(t));
K_20 = zeros(3,length(t));

%compute
for i=1:length(t)
    K_02(:,i) = Q(:,:,i)*C*u_02(:,i);
    K_11(:,i) = Q(:,:,i)*C*u_11(:,i);
    K_20(:,i) = Q(:,:,i)*C*u_20(:,i);
end


%% repeat auto-diff for m = 3 :)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - update the known coefficient books
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%K_02
coeff_x(1,3,:) = K_02(1,:);
coeff_y(1,3,:) = K_02(2,:);
coeff_z(1,3,:) = K_02(3,:);

%K_11
coeff_x(2,2,:) = K_11(1,:);
coeff_y(2,2,:) = K_11(2,:);
coeff_z(2,2,:) = K_11(3,:);

%K_20
coeff_x(3,1,:) = K_20(1,:);
coeff_y(3,1,:) = K_20(2,:);
coeff_z(3,1,:) = K_20(3,:);

%organize unknown data structure
comp_coeff_x = zeros(L,L,length(t));
comp_coeff_y = zeros(L,L,length(t));
comp_coeff_z = zeros(L,L,length(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           step 2 - autodifferentiate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   z component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it
    temp = zeros(size(comp_coeff_y(:,:,1)));
    for m=3:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp(alpha+1,m-alpha+1) = temp(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    
    %other stuff
    comp_coeff_y(:,:,i) = 2/tau_m*temp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %do it (x^2)
    temp1 = zeros(size(comp_coeff_x(:,:,1)));
    for m=3:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp1(alpha+1,m-alpha+1) = temp1(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_x(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
            
    %do it (y^2)
    temp2 = zeros(size(comp_coeff_x(:,:,1)));
    for m=3:L-1
        for alpha=0:m

            for k=0:alpha
                for ell=0:m-alpha
                    temp2(alpha+1,m-alpha+1) = temp2(alpha+1,m-alpha+1) + coeff_y(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    temp2 = temp2*(-pi^2*tau_m^2);
        
    %all together
    comp_coeff_x(:,:,i) = (temp1 + temp2)/tau_m;
        
    %reset
    temp=[];
    temp1=[];
    temp2=[];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               step 3 - save values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%beta's (for m=2)
beta_03 = [squeeze(comp_coeff_x(1,4,:))'; squeeze(comp_coeff_y(1,4,:))'; squeeze(comp_coeff_z(1,4,:))'];
beta_12 = [squeeze(comp_coeff_x(2,3,:))'; squeeze(comp_coeff_y(2,3,:))'; squeeze(comp_coeff_z(2,3,:))'];
beta_21 = [squeeze(comp_coeff_x(3,2,:))'; squeeze(comp_coeff_y(3,2,:))'; squeeze(comp_coeff_z(3,2,:))'];
beta_30 = [squeeze(comp_coeff_x(4,1,:))'; squeeze(comp_coeff_y(4,1,:))'; squeeze(comp_coeff_z(4,1,:))'];


%% Solve the homological equations m = 3

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_03 = zeros(3,length(t));
A_12 = zeros(3,length(t));
A_21 = zeros(3,length(t));
A_30 = zeros(3,length(t));

%solve
for i=1:length(t)
    A_03(:,i) = C\((Q(:,:,i))\beta_03(:,i));
    A_12(:,i) = C\((Q(:,:,i))\beta_12(:,i));
    A_21(:,i) = C\((Q(:,:,i))\beta_21(:,i));
    A_30(:,i) = C\((Q(:,:,i))\beta_30(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_03_f = [fftshift(fft(A_03(1,1:end))); fftshift(fft(A_03(2,1:end))); fftshift(fft(A_03(3,1:end)))];
A_12_f = [fftshift(fft(A_12(1,1:end))); fftshift(fft(A_12(2,1:end))); fftshift(fft(A_12(3,1:end)))];
A_21_f = [fftshift(fft(A_21(1,1:end))); fftshift(fft(A_21(2,1:end))); fftshift(fft(A_21(3,1:end)))];
A_30_f = [fftshift(fft(A_30(1,1:end))); fftshift(fft(A_30(2,1:end))); fftshift(fft(A_30(3,1:end)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_03 = zeros(3,length(t));
u_12 = zeros(3,length(t));
u_21 = zeros(3,length(t));
u_30 = zeros(3,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_03(1,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(3-0)*lambda_1-0).*A_03_f(1,:);
u_03(2,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(3-0)*lambda_1-lambda_1).*A_03_f(2,:);
u_03(3,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(3-0)*lambda_1-lambda_2).*A_03_f(3,:);

u_12(1,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(3-1)*lambda_1-0).*A_12_f(1,:);
u_12(2,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(3-1)*lambda_1-lambda_1).*A_12_f(2,:);
u_12(3,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(3-1)*lambda_1-lambda_2).*A_12_f(3,:);

u_21(1,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(3-2)*lambda_1-0).*A_21_f(1,:);
u_21(2,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(3-2)*lambda_1-lambda_1).*A_21_f(2,:);
u_21(3,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(3-2)*lambda_1-lambda_2).*A_21_f(3,:);

u_30(1,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(3-3)*lambda_1-0).*A_30_f(1,:);
u_30(2,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(3-3)*lambda_1-lambda_1).*A_30_f(2,:);
u_30(3,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(3-3)*lambda_1-lambda_2).*A_30_f(3,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_03 = ([ifft(ifftshift(u_03(1,:))); ifft(ifftshift(u_03(2,:))); ifft(ifftshift(u_03(3,:)))]);
u_12 = ([ifft(ifftshift(u_12(1,:))); ifft(ifftshift(u_12(2,:))); ifft(ifftshift(u_12(3,:)))]);
u_21 = ([ifft(ifftshift(u_21(1,:))); ifft(ifftshift(u_21(2,:))); ifft(ifftshift(u_21(3,:)))]);
u_30 = ([ifft(ifftshift(u_30(1,:))); ifft(ifftshift(u_30(2,:))); ifft(ifftshift(u_30(3,:)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_03 = zeros(3,length(t));
K_12 = zeros(3,length(t));
K_21 = zeros(3,length(t));
K_30 = zeros(3,length(t));

%compute
for i=1:length(t)
    K_03(:,i) = Q(:,:,i)*C*u_03(:,i);
    K_12(:,i) = Q(:,:,i)*C*u_12(:,i);
    K_21(:,i) = Q(:,:,i)*C*u_21(:,i);
    K_30(:,i) = Q(:,:,i)*C*u_30(:,i);
end


%% repeat auto-diff for m = 4 :)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - update the known coefficient books
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%K_03
coeff_x(1,4,:) = K_03(1,:);
coeff_y(1,4,:) = K_03(2,:);
coeff_z(1,4,:) = K_03(3,:);

%K_12
coeff_x(2,3,:) = K_12(1,:);
coeff_y(2,3,:) = K_12(2,:);
coeff_z(2,3,:) = K_12(3,:);

%K_21
coeff_x(3,2,:) = K_21(1,:);
coeff_y(3,2,:) = K_21(2,:);
coeff_z(3,2,:) = K_21(3,:);

%K_30
coeff_x(4,1,:) = K_30(1,:);
coeff_y(4,1,:) = K_30(2,:);
coeff_z(4,1,:) = K_30(3,:);

%organize unknown data structure
comp_coeff_x = zeros(L,L,length(t));
comp_coeff_y = zeros(L,L,length(t));
comp_coeff_z = zeros(L,L,length(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           step 2 - autodifferentiate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   z component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it
    temp = zeros(size(comp_coeff_y(:,:,1)));
    for m=4:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp(alpha+1,m-alpha+1) = temp(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    
    %other stuff
    comp_coeff_y(:,:,i) = 2/tau_m*temp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %do it (x^2)
    temp1 = zeros(size(comp_coeff_x(:,:,1)));
    for m=4:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp1(alpha+1,m-alpha+1) = temp1(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_x(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
            
    %do it (y^2)
    temp2 = zeros(size(comp_coeff_x(:,:,1)));
    for m=4:L-1
        for alpha=0:m

            for k=0:alpha
                for ell=0:m-alpha
                    temp2(alpha+1,m-alpha+1) = temp2(alpha+1,m-alpha+1) + coeff_y(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    temp2 = temp2*(-pi^2*tau_m^2);
        
    %all together
    comp_coeff_x(:,:,i) = (temp1 + temp2)/tau_m;
        
    %reset
    temp=[];
    temp1=[];
    temp2=[];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               step 3 - save values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%beta's (for m=4)
beta_04 = [squeeze(comp_coeff_x(1,5,:))'; squeeze(comp_coeff_y(1,5,:))'; squeeze(comp_coeff_z(1,5,:))'];
beta_13 = [squeeze(comp_coeff_x(2,4,:))'; squeeze(comp_coeff_y(2,4,:))'; squeeze(comp_coeff_z(2,4,:))'];
beta_22 = [squeeze(comp_coeff_x(3,3,:))'; squeeze(comp_coeff_y(3,3,:))'; squeeze(comp_coeff_z(3,3,:))'];
beta_31 = [squeeze(comp_coeff_x(4,2,:))'; squeeze(comp_coeff_y(4,2,:))'; squeeze(comp_coeff_z(4,2,:))'];
beta_40 = [squeeze(comp_coeff_x(5,1,:))'; squeeze(comp_coeff_y(5,1,:))'; squeeze(comp_coeff_z(5,1,:))'];


%% Solve the homological equations m = 4

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_04 = zeros(3,length(t));
A_13 = zeros(3,length(t));
A_22 = zeros(3,length(t));
A_31 = zeros(3,length(t));
A_40 = zeros(3,length(t));

%solve
for i=1:length(t)
    A_04(:,i) = C\((Q(:,:,i))\beta_04(:,i));
    A_13(:,i) = C\((Q(:,:,i))\beta_13(:,i));
    A_22(:,i) = C\((Q(:,:,i))\beta_22(:,i));
    A_31(:,i) = C\((Q(:,:,i))\beta_31(:,i));
    A_40(:,i) = C\((Q(:,:,i))\beta_40(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_04_f = [fftshift(fft(A_04(1,1:end))); fftshift(fft(A_04(2,1:end))); fftshift(fft(A_04(3,1:end)))];
A_13_f = [fftshift(fft(A_13(1,1:end))); fftshift(fft(A_13(2,1:end))); fftshift(fft(A_13(3,1:end)))];
A_22_f = [fftshift(fft(A_22(1,1:end))); fftshift(fft(A_22(2,1:end))); fftshift(fft(A_22(3,1:end)))];
A_31_f = [fftshift(fft(A_31(1,1:end))); fftshift(fft(A_31(2,1:end))); fftshift(fft(A_31(3,1:end)))];
A_40_f = [fftshift(fft(A_40(1,1:end))); fftshift(fft(A_40(2,1:end))); fftshift(fft(A_40(3,1:end)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_04 = zeros(3,length(t));
u_13 = zeros(3,length(t));
u_22 = zeros(3,length(t));
u_31 = zeros(3,length(t));
u_40 = zeros(3,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_04(1,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(4-0)*lambda_1-0).*A_04_f(1,:);
u_04(2,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(4-0)*lambda_1-lambda_1).*A_04_f(2,:);
u_04(3,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(4-0)*lambda_1-lambda_2).*A_04_f(3,:);

u_13(1,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(4-1)*lambda_1-0).*A_13_f(1,:);
u_13(2,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(4-1)*lambda_1-lambda_1).*A_13_f(2,:);
u_13(3,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(4-1)*lambda_1-lambda_2).*A_13_f(3,:);

u_22(1,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(4-2)*lambda_1-0).*A_22_f(1,:);
u_22(2,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(4-2)*lambda_1-lambda_1).*A_22_f(2,:);
u_22(3,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(4-2)*lambda_1-lambda_2).*A_22_f(3,:);

u_31(1,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(4-3)*lambda_1-0).*A_31_f(1,:);
u_31(2,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(4-3)*lambda_1-lambda_1).*A_31_f(2,:);
u_31(3,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(4-3)*lambda_1-lambda_2).*A_31_f(3,:);

u_40(1,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(4-4)*lambda_1-0).*A_40_f(1,:);
u_40(2,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(4-4)*lambda_1-lambda_1).*A_40_f(2,:);
u_40(3,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(4-4)*lambda_1-lambda_2).*A_40_f(3,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_04 = ([ifft(ifftshift(u_04(1,:))); ifft(ifftshift(u_04(2,:))); ifft(ifftshift(u_04(3,:)))]);
u_13 = ([ifft(ifftshift(u_13(1,:))); ifft(ifftshift(u_13(2,:))); ifft(ifftshift(u_13(3,:)))]);
u_22 = ([ifft(ifftshift(u_22(1,:))); ifft(ifftshift(u_22(2,:))); ifft(ifftshift(u_22(3,:)))]);
u_31 = ([ifft(ifftshift(u_31(1,:))); ifft(ifftshift(u_31(2,:))); ifft(ifftshift(u_31(3,:)))]);
u_40 = ([ifft(ifftshift(u_40(1,:))); ifft(ifftshift(u_40(2,:))); ifft(ifftshift(u_40(3,:)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_04 = zeros(3,length(t));
K_13 = zeros(3,length(t));
K_22 = zeros(3,length(t));
K_31 = zeros(3,length(t));
K_40 = zeros(3,length(t));

%compute
for i=1:length(t)
    K_04(:,i) = Q(:,:,i)*C*u_04(:,i);
    K_13(:,i) = Q(:,:,i)*C*u_13(:,i);
    K_22(:,i) = Q(:,:,i)*C*u_22(:,i);
    K_31(:,i) = Q(:,:,i)*C*u_31(:,i);
    K_40(:,i) = Q(:,:,i)*C*u_40(:,i);
end


%% repeat auto-diff for m = 5 :)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - update the known coefficient books
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%K_04
coeff_x(1,5,:) = K_04(1,:);
coeff_y(1,5,:) = K_04(2,:);
coeff_z(1,5,:) = K_04(3,:);

%K_13
coeff_x(2,4,:) = K_13(1,:);
coeff_y(2,4,:) = K_13(2,:);
coeff_z(2,4,:) = K_13(3,:);

%K_22
coeff_x(3,3,:) = K_22(1,:);
coeff_y(3,3,:) = K_22(2,:);
coeff_z(3,3,:) = K_22(3,:);

%K_31
coeff_x(4,2,:) = K_31(1,:);
coeff_y(4,2,:) = K_31(2,:);
coeff_z(4,2,:) = K_31(3,:);

%K_40
coeff_x(5,1,:) = K_40(1,:);
coeff_y(5,1,:) = K_40(2,:);
coeff_z(5,1,:) = K_40(3,:);

%organize unknown data structure
comp_coeff_x = zeros(L,L,length(t));
comp_coeff_y = zeros(L,L,length(t));
comp_coeff_z = zeros(L,L,length(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           step 2 - autodifferentiate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   z component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it
    temp = zeros(size(comp_coeff_y(:,:,1)));
    for m=5:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp(alpha+1,m-alpha+1) = temp(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    
    %other stuff
    comp_coeff_y(:,:,i) = 2/tau_m*temp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %do it (x^2)
    temp1 = zeros(size(comp_coeff_x(:,:,1)));
    for m=5:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp1(alpha+1,m-alpha+1) = temp1(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_x(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
            
    %do it (y^2)
    temp2 = zeros(size(comp_coeff_x(:,:,1)));
    for m=5:L-1
        for alpha=0:m

            for k=0:alpha
                for ell=0:m-alpha
                    temp2(alpha+1,m-alpha+1) = temp2(alpha+1,m-alpha+1) + coeff_y(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    temp2 = temp2*(-pi^2*tau_m^2);
        
    %all together
    comp_coeff_x(:,:,i) = (temp1 + temp2)/tau_m;
        
    %reset
    temp=[];
    temp1=[];
    temp2=[];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               step 3 - save values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%beta's (for m=4)
beta_05 = [squeeze(comp_coeff_x(1,6,:))'; squeeze(comp_coeff_y(1,6,:))'; squeeze(comp_coeff_z(1,6,:))'];
beta_14 = [squeeze(comp_coeff_x(2,5,:))'; squeeze(comp_coeff_y(2,5,:))'; squeeze(comp_coeff_z(2,5,:))'];
beta_23 = [squeeze(comp_coeff_x(3,4,:))'; squeeze(comp_coeff_y(3,4,:))'; squeeze(comp_coeff_z(3,4,:))'];
beta_32 = [squeeze(comp_coeff_x(4,3,:))'; squeeze(comp_coeff_y(4,3,:))'; squeeze(comp_coeff_z(4,3,:))'];
beta_41 = [squeeze(comp_coeff_x(5,2,:))'; squeeze(comp_coeff_y(5,2,:))'; squeeze(comp_coeff_z(5,2,:))'];
beta_50 = [squeeze(comp_coeff_x(6,1,:))'; squeeze(comp_coeff_y(6,1,:))'; squeeze(comp_coeff_z(6,1,:))'];


%% Solve the homological equations m = 5

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_05 = zeros(3,length(t));
A_14 = zeros(3,length(t));
A_23 = zeros(3,length(t));
A_32 = zeros(3,length(t));
A_41 = zeros(3,length(t));
A_50 = zeros(3,length(t));

%solve
for i=1:length(t)
    A_05(:,i) = C\((Q(:,:,i))\beta_05(:,i));
    A_14(:,i) = C\((Q(:,:,i))\beta_14(:,i));
    A_23(:,i) = C\((Q(:,:,i))\beta_23(:,i));
    A_32(:,i) = C\((Q(:,:,i))\beta_32(:,i));
    A_41(:,i) = C\((Q(:,:,i))\beta_41(:,i));
    A_50(:,i) = C\((Q(:,:,i))\beta_50(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_05_f = [fftshift(fft(A_05(1,1:end))); fftshift(fft(A_05(2,1:end))); fftshift(fft(A_05(3,1:end)))];
A_14_f = [fftshift(fft(A_14(1,1:end))); fftshift(fft(A_14(2,1:end))); fftshift(fft(A_14(3,1:end)))];
A_23_f = [fftshift(fft(A_23(1,1:end))); fftshift(fft(A_23(2,1:end))); fftshift(fft(A_23(3,1:end)))];
A_32_f = [fftshift(fft(A_32(1,1:end))); fftshift(fft(A_32(2,1:end))); fftshift(fft(A_32(3,1:end)))];
A_41_f = [fftshift(fft(A_41(1,1:end))); fftshift(fft(A_41(2,1:end))); fftshift(fft(A_41(3,1:end)))];
A_50_f = [fftshift(fft(A_50(1,1:end))); fftshift(fft(A_50(2,1:end))); fftshift(fft(A_50(3,1:end)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_05 = zeros(3,length(t));
u_14 = zeros(3,length(t));
u_23 = zeros(3,length(t));
u_32 = zeros(3,length(t));
u_41 = zeros(3,length(t));
u_50 = zeros(3,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_05(1,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(5-0)*lambda_1-0).*A_05_f(1,:);
u_05(2,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(5-0)*lambda_1-lambda_1).*A_05_f(2,:);
u_05(3,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(5-0)*lambda_1-lambda_2).*A_05_f(3,:);

u_14(1,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(5-1)*lambda_1-0).*A_14_f(1,:);
u_14(2,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(5-1)*lambda_1-lambda_1).*A_14_f(2,:);
u_14(3,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(5-1)*lambda_1-lambda_2).*A_14_f(3,:);

u_23(1,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(5-2)*lambda_1-0).*A_23_f(1,:);
u_23(2,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(5-2)*lambda_1-lambda_1).*A_23_f(2,:);
u_23(3,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(5-2)*lambda_1-lambda_2).*A_23_f(3,:);

u_32(1,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(5-3)*lambda_1-0).*A_32_f(1,:);
u_32(2,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(5-3)*lambda_1-lambda_1).*A_32_f(2,:);
u_32(3,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(5-3)*lambda_1-lambda_2).*A_32_f(3,:);

u_41(1,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(5-4)*lambda_1-0).*A_41_f(1,:);
u_41(2,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(5-4)*lambda_1-lambda_1).*A_41_f(2,:);
u_41(3,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(5-4)*lambda_1-lambda_2).*A_41_f(3,:);

u_50(1,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_2+(5-5)*lambda_1-0).*A_50_f(1,:);
u_50(2,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_2+(5-5)*lambda_1-lambda_1).*A_50_f(2,:);
u_50(3,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_2+(5-5)*lambda_1-lambda_2).*A_50_f(3,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_05 = ([ifft(ifftshift(u_05(1,:))); ifft(ifftshift(u_05(2,:))); ifft(ifftshift(u_05(3,:)))]);
u_14 = ([ifft(ifftshift(u_14(1,:))); ifft(ifftshift(u_14(2,:))); ifft(ifftshift(u_14(3,:)))]);
u_23 = ([ifft(ifftshift(u_23(1,:))); ifft(ifftshift(u_23(2,:))); ifft(ifftshift(u_23(3,:)))]);
u_32 = ([ifft(ifftshift(u_32(1,:))); ifft(ifftshift(u_32(2,:))); ifft(ifftshift(u_32(3,:)))]);
u_41 = ([ifft(ifftshift(u_41(1,:))); ifft(ifftshift(u_41(2,:))); ifft(ifftshift(u_41(3,:)))]);
u_50 = ([ifft(ifftshift(u_50(1,:))); ifft(ifftshift(u_50(2,:))); ifft(ifftshift(u_50(3,:)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_05 = zeros(3,length(t));
K_14 = zeros(3,length(t));
K_23 = zeros(3,length(t));
K_32 = zeros(3,length(t));
K_41 = zeros(3,length(t));
K_50 = zeros(3,length(t));

%compute
for i=1:length(t)
    K_05(:,i) = Q(:,:,i)*C*u_05(:,i);
    K_14(:,i) = Q(:,:,i)*C*u_14(:,i);
    K_23(:,i) = Q(:,:,i)*C*u_23(:,i);
    K_32(:,i) = Q(:,:,i)*C*u_32(:,i);
    K_41(:,i) = Q(:,:,i)*C*u_41(:,i);
    K_50(:,i) = Q(:,:,i)*C*u_50(:,i);
end


%% repeat auto-diff for m = 6 :)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - update the known coefficient books
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%K_05
coeff_x(1,6,:) = K_05(1,:);
coeff_y(1,6,:) = K_05(2,:);
coeff_z(1,6,:) = K_05(3,:);

%K_14
coeff_x(2,5,:) = K_14(1,:);
coeff_y(2,5,:) = K_14(2,:);
coeff_z(2,5,:) = K_14(3,:);

%K_23
coeff_x(3,4,:) = K_23(1,:);
coeff_y(3,4,:) = K_23(2,:);
coeff_z(3,4,:) = K_23(3,:);

%K_32
coeff_x(4,3,:) = K_32(1,:);
coeff_y(4,3,:) = K_32(2,:);
coeff_z(4,3,:) = K_32(3,:);

%K_41
coeff_x(5,2,:) = K_41(1,:);
coeff_y(5,2,:) = K_41(2,:);
coeff_z(5,2,:) = K_41(3,:);

%K_50
coeff_x(6,1,:) = K_50(1,:);
coeff_y(6,1,:) = K_50(2,:);
coeff_z(6,1,:) = K_50(3,:);

%organize unknown data structure
comp_coeff_x = zeros(L,L,length(t));
comp_coeff_y = zeros(L,L,length(t));
comp_coeff_z = zeros(L,L,length(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           step 2 - autodifferentiate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   z component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it
    temp = zeros(size(comp_coeff_y(:,:,1)));
    for m=6:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp(alpha+1,m-alpha+1) = temp(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    
    %other stuff
    comp_coeff_y(:,:,i) = 2/tau_m*temp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %do it (x^2)
    temp1 = zeros(size(comp_coeff_x(:,:,1)));
    for m=6:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp1(alpha+1,m-alpha+1) = temp1(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_x(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
            
    %do it (y^2)
    temp2 = zeros(size(comp_coeff_x(:,:,1)));
    for m=6:L-1
        for alpha=0:m

            for k=0:alpha
                for ell=0:m-alpha
                    temp2(alpha+1,m-alpha+1) = temp2(alpha+1,m-alpha+1) + coeff_y(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    temp2 = temp2*(-pi^2*tau_m^2);
        
    %all together
    comp_coeff_x(:,:,i) = (temp1 + temp2)/tau_m;
        
    %reset
    temp=[];
    temp1=[];
    temp2=[];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               step 3 - save values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%beta's (for m=4)
beta_06 = [squeeze(comp_coeff_x(1,7,:))'; squeeze(comp_coeff_y(1,7,:))'; squeeze(comp_coeff_z(1,7,:))'];
beta_15 = [squeeze(comp_coeff_x(2,6,:))'; squeeze(comp_coeff_y(2,6,:))'; squeeze(comp_coeff_z(2,6,:))'];
beta_24 = [squeeze(comp_coeff_x(3,5,:))'; squeeze(comp_coeff_y(3,5,:))'; squeeze(comp_coeff_z(3,5,:))'];
beta_33 = [squeeze(comp_coeff_x(4,4,:))'; squeeze(comp_coeff_y(4,4,:))'; squeeze(comp_coeff_z(4,4,:))'];
beta_42 = [squeeze(comp_coeff_x(5,3,:))'; squeeze(comp_coeff_y(5,3,:))'; squeeze(comp_coeff_z(5,3,:))'];
beta_51 = [squeeze(comp_coeff_x(6,2,:))'; squeeze(comp_coeff_y(6,2,:))'; squeeze(comp_coeff_z(6,2,:))'];
beta_60 = [squeeze(comp_coeff_x(7,1,:))'; squeeze(comp_coeff_y(7,1,:))'; squeeze(comp_coeff_z(7,1,:))'];


%% Solve the homological equations m = 6

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_06 = zeros(3,length(t));
A_15 = zeros(3,length(t));
A_24 = zeros(3,length(t));
A_33 = zeros(3,length(t));
A_42 = zeros(3,length(t));
A_51 = zeros(3,length(t));
A_60 = zeros(3,length(t));

%solve
for i=1:length(t)
    A_06(:,i) = C\((Q(:,:,i))\beta_06(:,i));
    A_15(:,i) = C\((Q(:,:,i))\beta_15(:,i));
    A_24(:,i) = C\((Q(:,:,i))\beta_24(:,i));
    A_33(:,i) = C\((Q(:,:,i))\beta_33(:,i));
    A_42(:,i) = C\((Q(:,:,i))\beta_42(:,i));
    A_51(:,i) = C\((Q(:,:,i))\beta_51(:,i));
    A_60(:,i) = C\((Q(:,:,i))\beta_60(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_06_f = [fftshift(fft(A_06(1,1:end))); fftshift(fft(A_06(2,1:end))); fftshift(fft(A_06(3,1:end)))];
A_15_f = [fftshift(fft(A_15(1,1:end))); fftshift(fft(A_15(2,1:end))); fftshift(fft(A_15(3,1:end)))];
A_24_f = [fftshift(fft(A_24(1,1:end))); fftshift(fft(A_24(2,1:end))); fftshift(fft(A_24(3,1:end)))];
A_33_f = [fftshift(fft(A_33(1,1:end))); fftshift(fft(A_33(2,1:end))); fftshift(fft(A_33(3,1:end)))];
A_42_f = [fftshift(fft(A_42(1,1:end))); fftshift(fft(A_42(2,1:end))); fftshift(fft(A_42(3,1:end)))];
A_51_f = [fftshift(fft(A_51(1,1:end))); fftshift(fft(A_51(2,1:end))); fftshift(fft(A_51(3,1:end)))];
A_60_f = [fftshift(fft(A_60(1,1:end))); fftshift(fft(A_60(2,1:end))); fftshift(fft(A_60(3,1:end)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_06 = zeros(3,length(t));
u_15 = zeros(3,length(t));
u_24 = zeros(3,length(t));
u_33 = zeros(3,length(t));
u_42 = zeros(3,length(t));
u_51 = zeros(3,length(t));
u_60 = zeros(3,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_06(1,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(6-0)*lambda_1-0).*A_06_f(1,:);
u_06(2,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(6-0)*lambda_1-lambda_1).*A_06_f(2,:);
u_06(3,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(6-0)*lambda_1-lambda_2).*A_06_f(3,:);

u_15(1,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(6-1)*lambda_1-0).*A_15_f(1,:);
u_15(2,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(6-1)*lambda_1-lambda_1).*A_15_f(2,:);
u_15(3,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(6-1)*lambda_1-lambda_2).*A_15_f(3,:);

u_24(1,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(6-2)*lambda_1-0).*A_24_f(1,:);
u_24(2,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(6-2)*lambda_1-lambda_1).*A_24_f(2,:);
u_24(3,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(6-2)*lambda_1-lambda_2).*A_24_f(3,:);

u_33(1,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(6-3)*lambda_1-0).*A_33_f(1,:);
u_33(2,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(6-3)*lambda_1-lambda_1).*A_33_f(2,:);
u_33(3,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(6-3)*lambda_1-lambda_2).*A_33_f(3,:);

u_42(1,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(6-4)*lambda_1-0).*A_42_f(1,:);
u_42(2,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(6-4)*lambda_1-lambda_1).*A_42_f(2,:);
u_42(3,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(6-4)*lambda_1-lambda_2).*A_42_f(3,:);

u_51(1,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_2+(6-5)*lambda_1-0).*A_51_f(1,:);
u_51(2,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_2+(6-5)*lambda_1-lambda_1).*A_51_f(2,:);
u_51(3,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_2+(6-5)*lambda_1-lambda_2).*A_51_f(3,:);

u_60(1,:) = 1./(2*pi*sqrt(-1).*k/T+6*lambda_2+(6-6)*lambda_1-0).*A_60_f(1,:);
u_60(2,:) = 1./(2*pi*sqrt(-1).*k/T+6*lambda_2+(6-6)*lambda_1-lambda_1).*A_60_f(2,:);
u_60(3,:) = 1./(2*pi*sqrt(-1).*k/T+6*lambda_2+(6-6)*lambda_1-lambda_2).*A_60_f(3,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_06 = ([ifft(ifftshift(u_06(1,:))); ifft(ifftshift(u_06(2,:))); ifft(ifftshift(u_06(3,:)))]);
u_15 = ([ifft(ifftshift(u_15(1,:))); ifft(ifftshift(u_15(2,:))); ifft(ifftshift(u_15(3,:)))]);
u_24 = ([ifft(ifftshift(u_24(1,:))); ifft(ifftshift(u_24(2,:))); ifft(ifftshift(u_24(3,:)))]);
u_33 = ([ifft(ifftshift(u_33(1,:))); ifft(ifftshift(u_33(2,:))); ifft(ifftshift(u_33(3,:)))]);
u_42 = ([ifft(ifftshift(u_42(1,:))); ifft(ifftshift(u_42(2,:))); ifft(ifftshift(u_42(3,:)))]);
u_51 = ([ifft(ifftshift(u_51(1,:))); ifft(ifftshift(u_51(2,:))); ifft(ifftshift(u_51(3,:)))]);
u_60 = ([ifft(ifftshift(u_60(1,:))); ifft(ifftshift(u_60(2,:))); ifft(ifftshift(u_60(3,:)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_06 = zeros(3,length(t));
K_15 = zeros(3,length(t));
K_24 = zeros(3,length(t));
K_33 = zeros(3,length(t));
K_42 = zeros(3,length(t));
K_51 = zeros(3,length(t));
K_60 = zeros(3,length(t));

%compute
for i=1:length(t)
    K_06(:,i) = Q(:,:,i)*C*u_06(:,i);
    K_15(:,i) = Q(:,:,i)*C*u_15(:,i);
    K_24(:,i) = Q(:,:,i)*C*u_24(:,i);
    K_33(:,i) = Q(:,:,i)*C*u_33(:,i);
    K_42(:,i) = Q(:,:,i)*C*u_42(:,i);
    K_51(:,i) = Q(:,:,i)*C*u_51(:,i);
    K_60(:,i) = Q(:,:,i)*C*u_60(:,i);
end


%% repeat auto-diff for m = 7 :)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - update the known coefficient books
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%K_06
coeff_x(1,7,:) = K_06(1,:);
coeff_y(1,7,:) = K_06(2,:);
coeff_z(1,7,:) = K_06(3,:);

%K_15
coeff_x(2,6,:) = K_15(1,:);
coeff_y(2,6,:) = K_15(2,:);
coeff_z(2,6,:) = K_15(3,:);

%K_24
coeff_x(3,5,:) = K_24(1,:);
coeff_y(3,5,:) = K_24(2,:);
coeff_z(3,5,:) = K_24(3,:);

%K_33
coeff_x(4,4,:) = K_33(1,:);
coeff_y(4,4,:) = K_33(2,:);
coeff_z(4,4,:) = K_33(3,:);

%K_42
coeff_x(5,3,:) = K_42(1,:);
coeff_y(5,3,:) = K_42(2,:);
coeff_z(5,3,:) = K_42(3,:);

%K_51
coeff_x(6,2,:) = K_51(1,:);
coeff_y(6,2,:) = K_51(2,:);
coeff_z(6,2,:) = K_51(3,:);

%K_60
coeff_x(7,1,:) = K_60(1,:);
coeff_y(7,1,:) = K_60(2,:);
coeff_z(7,1,:) = K_60(3,:);

%organize unknown data structure
comp_coeff_x = zeros(L,L,length(t));
comp_coeff_y = zeros(L,L,length(t));
comp_coeff_z = zeros(L,L,length(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           step 2 - autodifferentiate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   z component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it
    temp = zeros(size(comp_coeff_y(:,:,1)));
    for m=7:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp(alpha+1,m-alpha+1) = temp(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    
    %other stuff
    comp_coeff_y(:,:,i) = 2/tau_m*temp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %do it (x^2)
    temp1 = zeros(size(comp_coeff_x(:,:,1)));
    for m=7:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp1(alpha+1,m-alpha+1) = temp1(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_x(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
            
    %do it (y^2)
    temp2 = zeros(size(comp_coeff_x(:,:,1)));
    for m=7:L-1
        for alpha=0:m

            for k=0:alpha
                for ell=0:m-alpha
                    temp2(alpha+1,m-alpha+1) = temp2(alpha+1,m-alpha+1) + coeff_y(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    temp2 = temp2*(-pi^2*tau_m^2);
        
    %all together
    comp_coeff_x(:,:,i) = (temp1 + temp2)/tau_m;
        
    %reset
    temp=[];
    temp1=[];
    temp2=[];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               step 3 - save values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%beta's (for m=4)
beta_07 = [squeeze(comp_coeff_x(1,8,:))'; squeeze(comp_coeff_y(1,8,:))'; squeeze(comp_coeff_z(1,8,:))'];
beta_16 = [squeeze(comp_coeff_x(2,7,:))'; squeeze(comp_coeff_y(2,7,:))'; squeeze(comp_coeff_z(2,7,:))'];
beta_25 = [squeeze(comp_coeff_x(3,6,:))'; squeeze(comp_coeff_y(3,6,:))'; squeeze(comp_coeff_z(3,6,:))'];
beta_34 = [squeeze(comp_coeff_x(4,5,:))'; squeeze(comp_coeff_y(4,5,:))'; squeeze(comp_coeff_z(4,5,:))'];
beta_43 = [squeeze(comp_coeff_x(5,4,:))'; squeeze(comp_coeff_y(5,4,:))'; squeeze(comp_coeff_z(5,4,:))'];
beta_52 = [squeeze(comp_coeff_x(6,3,:))'; squeeze(comp_coeff_y(6,3,:))'; squeeze(comp_coeff_z(6,3,:))'];
beta_61 = [squeeze(comp_coeff_x(7,2,:))'; squeeze(comp_coeff_y(7,2,:))'; squeeze(comp_coeff_z(7,2,:))'];
beta_70 = [squeeze(comp_coeff_x(8,1,:))'; squeeze(comp_coeff_y(8,1,:))'; squeeze(comp_coeff_z(8,1,:))'];


%% Solve the homological equations m = 7

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_07 = zeros(3,length(t));
A_16 = zeros(3,length(t));
A_25 = zeros(3,length(t));
A_34 = zeros(3,length(t));
A_43 = zeros(3,length(t));
A_52 = zeros(3,length(t));
A_61 = zeros(3,length(t));
A_70 = zeros(3,length(t));

%solve
for i=1:length(t)
    A_07(:,i) = C\((Q(:,:,i))\beta_07(:,i));
    A_16(:,i) = C\((Q(:,:,i))\beta_16(:,i));
    A_25(:,i) = C\((Q(:,:,i))\beta_25(:,i));
    A_34(:,i) = C\((Q(:,:,i))\beta_34(:,i));
    A_43(:,i) = C\((Q(:,:,i))\beta_43(:,i));
    A_52(:,i) = C\((Q(:,:,i))\beta_52(:,i));
    A_61(:,i) = C\((Q(:,:,i))\beta_61(:,i));
    A_70(:,i) = C\((Q(:,:,i))\beta_70(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_07_f = [fftshift(fft(A_07(1,1:end))); fftshift(fft(A_07(2,1:end))); fftshift(fft(A_07(3,1:end)))];
A_16_f = [fftshift(fft(A_16(1,1:end))); fftshift(fft(A_16(2,1:end))); fftshift(fft(A_16(3,1:end)))];
A_25_f = [fftshift(fft(A_25(1,1:end))); fftshift(fft(A_25(2,1:end))); fftshift(fft(A_25(3,1:end)))];
A_34_f = [fftshift(fft(A_34(1,1:end))); fftshift(fft(A_34(2,1:end))); fftshift(fft(A_34(3,1:end)))];
A_43_f = [fftshift(fft(A_43(1,1:end))); fftshift(fft(A_43(2,1:end))); fftshift(fft(A_43(3,1:end)))];
A_52_f = [fftshift(fft(A_52(1,1:end))); fftshift(fft(A_52(2,1:end))); fftshift(fft(A_52(3,1:end)))];
A_61_f = [fftshift(fft(A_61(1,1:end))); fftshift(fft(A_61(2,1:end))); fftshift(fft(A_61(3,1:end)))];
A_70_f = [fftshift(fft(A_70(1,1:end))); fftshift(fft(A_70(2,1:end))); fftshift(fft(A_70(3,1:end)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_07 = zeros(3,length(t));
u_16 = zeros(3,length(t));
u_25 = zeros(3,length(t));
u_34 = zeros(3,length(t));
u_43 = zeros(3,length(t));
u_52 = zeros(3,length(t));
u_61 = zeros(3,length(t));
u_70 = zeros(3,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_07(1,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(7-0)*lambda_1-0).*A_07_f(1,:);
u_07(2,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(7-0)*lambda_1-lambda_1).*A_07_f(2,:);
u_07(3,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(7-0)*lambda_1-lambda_2).*A_07_f(3,:);

u_16(1,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(7-1)*lambda_1-0).*A_16_f(1,:);
u_16(2,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(7-1)*lambda_1-lambda_1).*A_16_f(2,:);
u_16(3,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(7-1)*lambda_1-lambda_2).*A_16_f(3,:);

u_25(1,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(7-2)*lambda_1-0).*A_25_f(1,:);
u_25(2,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(7-2)*lambda_1-lambda_1).*A_25_f(2,:);
u_25(3,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(7-2)*lambda_1-lambda_2).*A_25_f(3,:);

u_34(1,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(7-3)*lambda_1-0).*A_34_f(1,:);
u_34(2,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(7-3)*lambda_1-lambda_1).*A_34_f(2,:);
u_34(3,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(7-3)*lambda_1-lambda_2).*A_34_f(3,:);

u_43(1,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(7-4)*lambda_1-0).*A_43_f(1,:);
u_43(2,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(7-4)*lambda_1-lambda_1).*A_43_f(2,:);
u_43(3,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(7-4)*lambda_1-lambda_2).*A_43_f(3,:);

u_52(1,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_2+(7-5)*lambda_1-0).*A_52_f(1,:);
u_52(2,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_2+(7-5)*lambda_1-lambda_1).*A_52_f(2,:);
u_52(3,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_2+(7-5)*lambda_1-lambda_2).*A_52_f(3,:);

u_61(1,:) = 1./(2*pi*sqrt(-1).*k/T+6*lambda_2+(7-6)*lambda_1-0).*A_61_f(1,:);
u_61(2,:) = 1./(2*pi*sqrt(-1).*k/T+6*lambda_2+(7-6)*lambda_1-lambda_1).*A_61_f(2,:);
u_61(3,:) = 1./(2*pi*sqrt(-1).*k/T+6*lambda_2+(7-6)*lambda_1-lambda_2).*A_61_f(3,:);

u_70(1,:) = 1./(2*pi*sqrt(-1).*k/T+7*lambda_2+(7-7)*lambda_1-0).*A_70_f(1,:);
u_70(2,:) = 1./(2*pi*sqrt(-1).*k/T+7*lambda_2+(7-7)*lambda_1-lambda_1).*A_70_f(2,:);
u_70(3,:) = 1./(2*pi*sqrt(-1).*k/T+7*lambda_2+(7-7)*lambda_1-lambda_2).*A_70_f(3,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_07 = ([ifft(ifftshift(u_07(1,:))); ifft(ifftshift(u_07(2,:))); ifft(ifftshift(u_07(3,:)))]);
u_16 = ([ifft(ifftshift(u_16(1,:))); ifft(ifftshift(u_16(2,:))); ifft(ifftshift(u_16(3,:)))]);
u_25 = ([ifft(ifftshift(u_25(1,:))); ifft(ifftshift(u_25(2,:))); ifft(ifftshift(u_25(3,:)))]);
u_34 = ([ifft(ifftshift(u_34(1,:))); ifft(ifftshift(u_34(2,:))); ifft(ifftshift(u_34(3,:)))]);
u_43 = ([ifft(ifftshift(u_43(1,:))); ifft(ifftshift(u_43(2,:))); ifft(ifftshift(u_43(3,:)))]);
u_52 = ([ifft(ifftshift(u_52(1,:))); ifft(ifftshift(u_52(2,:))); ifft(ifftshift(u_52(3,:)))]);
u_61 = ([ifft(ifftshift(u_61(1,:))); ifft(ifftshift(u_61(2,:))); ifft(ifftshift(u_61(3,:)))]);
u_70 = ([ifft(ifftshift(u_70(1,:))); ifft(ifftshift(u_70(2,:))); ifft(ifftshift(u_70(3,:)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_07 = zeros(3,length(t));
K_16 = zeros(3,length(t));
K_25 = zeros(3,length(t));
K_34 = zeros(3,length(t));
K_43 = zeros(3,length(t));
K_52 = zeros(3,length(t));
K_61 = zeros(3,length(t));
K_70 = zeros(3,length(t));

%compute
for i=1:length(t)
    K_07(:,i) = Q(:,:,i)*C*u_07(:,i);
    K_16(:,i) = Q(:,:,i)*C*u_16(:,i);
    K_25(:,i) = Q(:,:,i)*C*u_25(:,i);
    K_34(:,i) = Q(:,:,i)*C*u_34(:,i);
    K_43(:,i) = Q(:,:,i)*C*u_43(:,i);
    K_52(:,i) = Q(:,:,i)*C*u_52(:,i);
    K_61(:,i) = Q(:,:,i)*C*u_61(:,i);
    K_70(:,i) = Q(:,:,i)*C*u_70(:,i);
end


%% repeat auto-diff for m = 8 :)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - update the known coefficient books
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%K_07
coeff_x(1,8,:) = K_07(1,:);
coeff_y(1,8,:) = K_07(2,:);
coeff_z(1,8,:) = K_07(3,:);

%K_16
coeff_x(2,7,:) = K_16(1,:);
coeff_y(2,7,:) = K_16(2,:);
coeff_z(2,7,:) = K_16(3,:);

%K_25
coeff_x(3,6,:) = K_25(1,:);
coeff_y(3,6,:) = K_25(2,:);
coeff_z(3,6,:) = K_25(3,:);

%K_34
coeff_x(4,5,:) = K_34(1,:);
coeff_y(4,5,:) = K_34(2,:);
coeff_z(4,5,:) = K_34(3,:);

%K_43
coeff_x(5,4,:) = K_43(1,:);
coeff_y(5,4,:) = K_43(2,:);
coeff_z(5,4,:) = K_43(3,:);

%K_52
coeff_x(6,3,:) = K_52(1,:);
coeff_y(6,3,:) = K_52(2,:);
coeff_z(6,3,:) = K_52(3,:);

%K_61
coeff_x(7,2,:) = K_61(1,:);
coeff_y(7,2,:) = K_61(2,:);
coeff_z(7,2,:) = K_61(3,:);

%K_70
coeff_x(8,1,:) = K_70(1,:);
coeff_y(8,1,:) = K_70(2,:);
coeff_z(8,1,:) = K_70(3,:);

%organize unknown data structure
comp_coeff_x = zeros(L,L,length(t));
comp_coeff_y = zeros(L,L,length(t));
comp_coeff_z = zeros(L,L,length(t));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           step 2 - autodifferentiate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   z component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it
    temp = zeros(size(comp_coeff_y(:,:,1)));
    for m=8:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp(alpha+1,m-alpha+1) = temp(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    
    %other stuff
    comp_coeff_y(:,:,i) = 2/tau_m*temp;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
        
    %do it (x^2)
    temp1 = zeros(size(comp_coeff_x(:,:,1)));
    for m=8:L-1
        for alpha=0:m
            
            for k=0:alpha
                for ell=0:m-alpha
                    temp1(alpha+1,m-alpha+1) = temp1(alpha+1,m-alpha+1) + coeff_x(k+1,ell+1,i)*coeff_x(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
            
    %do it (y^2)
    temp2 = zeros(size(comp_coeff_x(:,:,1)));
    for m=8:L-1
        for alpha=0:m

            for k=0:alpha
                for ell=0:m-alpha
                    temp2(alpha+1,m-alpha+1) = temp2(alpha+1,m-alpha+1) + coeff_y(k+1,ell+1,i)*coeff_y(alpha-k+1,m-alpha-ell+1,i);
                end
            end
            
        end
    end
    temp2 = temp2*(-pi^2*tau_m^2);
        
    %all together
    comp_coeff_x(:,:,i) = (temp1 + temp2)/tau_m;
        
    %reset
    temp=[];
    temp1=[];
    temp2=[];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               step 3 - save values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%beta's (for m=4)
beta_08 = [squeeze(comp_coeff_x(1,9,:))'; squeeze(comp_coeff_y(1,9,:))'; squeeze(comp_coeff_z(1,9,:))'];
beta_17 = [squeeze(comp_coeff_x(2,8,:))'; squeeze(comp_coeff_y(2,8,:))'; squeeze(comp_coeff_z(2,8,:))'];
beta_26 = [squeeze(comp_coeff_x(3,7,:))'; squeeze(comp_coeff_y(3,7,:))'; squeeze(comp_coeff_z(3,7,:))'];
beta_35 = [squeeze(comp_coeff_x(4,6,:))'; squeeze(comp_coeff_y(4,6,:))'; squeeze(comp_coeff_z(4,6,:))'];
beta_44 = [squeeze(comp_coeff_x(5,5,:))'; squeeze(comp_coeff_y(5,5,:))'; squeeze(comp_coeff_z(5,5,:))'];
beta_53 = [squeeze(comp_coeff_x(6,4,:))'; squeeze(comp_coeff_y(6,4,:))'; squeeze(comp_coeff_z(6,4,:))'];
beta_62 = [squeeze(comp_coeff_x(7,3,:))'; squeeze(comp_coeff_y(7,3,:))'; squeeze(comp_coeff_z(7,3,:))'];
beta_71 = [squeeze(comp_coeff_x(8,2,:))'; squeeze(comp_coeff_y(8,2,:))'; squeeze(comp_coeff_z(8,2,:))'];
beta_80 = [squeeze(comp_coeff_x(9,1,:))'; squeeze(comp_coeff_y(9,1,:))'; squeeze(comp_coeff_z(9,1,:))'];


%% Solve the homological equations m = 8

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_08 = zeros(3,length(t));
A_17 = zeros(3,length(t));
A_26 = zeros(3,length(t));
A_35 = zeros(3,length(t));
A_44 = zeros(3,length(t));
A_53 = zeros(3,length(t));
A_62 = zeros(3,length(t));
A_71 = zeros(3,length(t));
A_80 = zeros(3,length(t));

%solve
for i=1:length(t)
    A_08(:,i) = C\((Q(:,:,i))\beta_08(:,i));
    A_17(:,i) = C\((Q(:,:,i))\beta_17(:,i));
    A_26(:,i) = C\((Q(:,:,i))\beta_26(:,i));
    A_35(:,i) = C\((Q(:,:,i))\beta_35(:,i));
    A_44(:,i) = C\((Q(:,:,i))\beta_44(:,i));
    A_53(:,i) = C\((Q(:,:,i))\beta_53(:,i));
    A_62(:,i) = C\((Q(:,:,i))\beta_62(:,i));
    A_71(:,i) = C\((Q(:,:,i))\beta_71(:,i));
    A_80(:,i) = C\((Q(:,:,i))\beta_80(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_08_f = [fftshift(fft(A_08(1,1:end))); fftshift(fft(A_08(2,1:end))); fftshift(fft(A_08(3,1:end)))];
A_17_f = [fftshift(fft(A_17(1,1:end))); fftshift(fft(A_17(2,1:end))); fftshift(fft(A_17(3,1:end)))];
A_26_f = [fftshift(fft(A_26(1,1:end))); fftshift(fft(A_26(2,1:end))); fftshift(fft(A_26(3,1:end)))];
A_35_f = [fftshift(fft(A_35(1,1:end))); fftshift(fft(A_35(2,1:end))); fftshift(fft(A_35(3,1:end)))];
A_44_f = [fftshift(fft(A_44(1,1:end))); fftshift(fft(A_44(2,1:end))); fftshift(fft(A_44(3,1:end)))];
A_53_f = [fftshift(fft(A_53(1,1:end))); fftshift(fft(A_53(2,1:end))); fftshift(fft(A_53(3,1:end)))];
A_62_f = [fftshift(fft(A_62(1,1:end))); fftshift(fft(A_62(2,1:end))); fftshift(fft(A_62(3,1:end)))];
A_71_f = [fftshift(fft(A_71(1,1:end))); fftshift(fft(A_71(2,1:end))); fftshift(fft(A_71(3,1:end)))];
A_80_f = [fftshift(fft(A_80(1,1:end))); fftshift(fft(A_80(2,1:end))); fftshift(fft(A_80(3,1:end)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_08 = zeros(3,length(t));
u_17 = zeros(3,length(t));
u_26 = zeros(3,length(t));
u_35 = zeros(3,length(t));
u_44 = zeros(3,length(t));
u_53 = zeros(3,length(t));
u_62 = zeros(3,length(t));
u_71 = zeros(3,length(t));
u_80 = zeros(3,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_08(1,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(8-0)*lambda_1-0).*A_08_f(1,:);
u_08(2,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(8-0)*lambda_1-lambda_1).*A_08_f(2,:);
u_08(3,:) = 1./(2*pi*sqrt(-1).*k/T+0*lambda_2+(8-0)*lambda_1-lambda_2).*A_08_f(3,:);

u_17(1,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(8-1)*lambda_1-0).*A_17_f(1,:);
u_17(2,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(8-1)*lambda_1-lambda_1).*A_17_f(2,:);
u_17(3,:) = 1./(2*pi*sqrt(-1).*k/T+1*lambda_2+(8-1)*lambda_1-lambda_2).*A_17_f(3,:);

u_26(1,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(8-2)*lambda_1-0).*A_26_f(1,:);
u_26(2,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(8-2)*lambda_1-lambda_1).*A_26_f(2,:);
u_26(3,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_2+(8-2)*lambda_1-lambda_2).*A_26_f(3,:);

u_35(1,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(8-3)*lambda_1-0).*A_35_f(1,:);
u_35(2,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(8-3)*lambda_1-lambda_1).*A_35_f(2,:);
u_35(3,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_2+(8-3)*lambda_1-lambda_2).*A_35_f(3,:);

u_44(1,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(8-4)*lambda_1-0).*A_44_f(1,:);
u_44(2,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(8-4)*lambda_1-lambda_1).*A_44_f(2,:);
u_44(3,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_2+(8-4)*lambda_1-lambda_2).*A_44_f(3,:);

u_53(1,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_2+(8-5)*lambda_1-0).*A_53_f(1,:);
u_53(2,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_2+(8-5)*lambda_1-lambda_1).*A_53_f(2,:);
u_53(3,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_2+(8-5)*lambda_1-lambda_2).*A_53_f(3,:);

u_62(1,:) = 1./(2*pi*sqrt(-1).*k/T+6*lambda_2+(8-6)*lambda_1-0).*A_62_f(1,:);
u_62(2,:) = 1./(2*pi*sqrt(-1).*k/T+6*lambda_2+(8-6)*lambda_1-lambda_1).*A_62_f(2,:);
u_62(3,:) = 1./(2*pi*sqrt(-1).*k/T+6*lambda_2+(8-6)*lambda_1-lambda_2).*A_62_f(3,:);

u_71(1,:) = 1./(2*pi*sqrt(-1).*k/T+7*lambda_2+(8-7)*lambda_1-0).*A_71_f(1,:);
u_71(2,:) = 1./(2*pi*sqrt(-1).*k/T+7*lambda_2+(8-7)*lambda_1-lambda_1).*A_71_f(2,:);
u_71(3,:) = 1./(2*pi*sqrt(-1).*k/T+7*lambda_2+(8-7)*lambda_1-lambda_2).*A_71_f(3,:);

u_80(1,:) = 1./(2*pi*sqrt(-1).*k/T+8*lambda_2+(8-8)*lambda_1-0).*A_80_f(1,:);
u_80(2,:) = 1./(2*pi*sqrt(-1).*k/T+8*lambda_2+(8-8)*lambda_1-lambda_1).*A_80_f(2,:);
u_80(3,:) = 1./(2*pi*sqrt(-1).*k/T+8*lambda_2+(8-8)*lambda_1-lambda_2).*A_80_f(3,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_08 = ([ifft(ifftshift(u_08(1,:))); ifft(ifftshift(u_08(2,:))); ifft(ifftshift(u_08(3,:)))]);
u_17 = ([ifft(ifftshift(u_17(1,:))); ifft(ifftshift(u_17(2,:))); ifft(ifftshift(u_17(3,:)))]);
u_26 = ([ifft(ifftshift(u_26(1,:))); ifft(ifftshift(u_26(2,:))); ifft(ifftshift(u_26(3,:)))]);
u_35 = ([ifft(ifftshift(u_35(1,:))); ifft(ifftshift(u_35(2,:))); ifft(ifftshift(u_35(3,:)))]);
u_44 = ([ifft(ifftshift(u_44(1,:))); ifft(ifftshift(u_44(2,:))); ifft(ifftshift(u_44(3,:)))]);
u_53 = ([ifft(ifftshift(u_53(1,:))); ifft(ifftshift(u_53(2,:))); ifft(ifftshift(u_53(3,:)))]);
u_62 = ([ifft(ifftshift(u_62(1,:))); ifft(ifftshift(u_62(2,:))); ifft(ifftshift(u_62(3,:)))]);
u_71 = ([ifft(ifftshift(u_71(1,:))); ifft(ifftshift(u_71(2,:))); ifft(ifftshift(u_71(3,:)))]);
u_80 = ([ifft(ifftshift(u_80(1,:))); ifft(ifftshift(u_80(2,:))); ifft(ifftshift(u_80(3,:)))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_08 = zeros(3,length(t));
K_17 = zeros(3,length(t));
K_26 = zeros(3,length(t));
K_35 = zeros(3,length(t));
K_44 = zeros(3,length(t));
K_53 = zeros(3,length(t));
K_62 = zeros(3,length(t));
K_71 = zeros(3,length(t));
K_80 = zeros(3,length(t));

%compute
for i=1:length(t)
    K_08(:,i) = Q(:,:,i)*C*u_08(:,i);
    K_17(:,i) = Q(:,:,i)*C*u_17(:,i);
    K_26(:,i) = Q(:,:,i)*C*u_26(:,i);
    K_35(:,i) = Q(:,:,i)*C*u_35(:,i);
    K_44(:,i) = Q(:,:,i)*C*u_44(:,i);
    K_53(:,i) = Q(:,:,i)*C*u_53(:,i);
    K_62(:,i) = Q(:,:,i)*C*u_62(:,i);
    K_71(:,i) = Q(:,:,i)*C*u_71(:,i);
    K_80(:,i) = Q(:,:,i)*C*u_80(:,i);
end


%% save data

%m=0
save('K00.mat','K_00')

%m=1
save('K10.mat','K_10')
save('K01.mat','K_01')

%m=2
save('K20.mat','K_20')
save('K11.mat','K_11')
save('K02.mat','K_02')

%m=3
save('K30.mat','K_30')
save('K21.mat','K_21')
save('K12.mat','K_12')
save('K03.mat','K_03')

%m=4
save('K40.mat','K_40')
save('K31.mat','K_31')
save('K22.mat','K_22')
save('K13.mat','K_13')
save('K04.mat','K_04')

%m=5
save('K50.mat','K_50')
save('K41.mat','K_41')
save('K32.mat','K_32')
save('K23.mat','K_23')
save('K14.mat','K_14')
save('K05.mat','K_05')

%m=6
save('K60.mat','K_60')
save('K51.mat','K_51')
save('K42.mat','K_42')
save('K33.mat','K_33')
save('K24.mat','K_24')
save('K15.mat','K_15')
save('K06.mat','K_06')

%m=7
save('K70.mat','K_70')
save('K61.mat','K_61')
save('K52.mat','K_52')
save('K43.mat','K_43')
save('K34.mat','K_34')
save('K25.mat','K_25')
save('K16.mat','K_16')
save('K07.mat','K_07')

%m=8
save('K80.mat','K_80')
save('K71.mat','K_71')
save('K62.mat','K_62')
save('K53.mat','K_53')
save('K44.mat','K_44')
save('K35.mat','K_35')
save('K26.mat','K_26')
save('K17.mat','K_17')
save('K08.mat','K_08')


%% compute response curves 

%take gradient of K00
dKx = gradient(K_00(1,:),Tmax/(N-1));
dKy = gradient(K_00(2,:),Tmax/(N-1));
dKz = gradient(K_00(3,:),Tmax/(N-1));

%initialize solution vector
Response_Curves = zeros(3,3,length(t));

%solve
for j=1:length(t)
    
    %do it
    Response_Curves(:,:,j) = [[dKx(j);dKy(j);dKz(j)], K_01(:,j), K_10(:,j)]\eye(3,3);
end

%iPRC
Z0 = squeeze(Response_Curves(1,:,:))';

%iIRCs
I0_1 = squeeze(Response_Curves(2,:,:))';
I0_2 = squeeze(Response_Curves(3,:,:))';


%% save data

%save it
save('Z0.mat','Z0')
save('I0_1.mat','I0_1')
save('I0_2.mat','I0_2')
