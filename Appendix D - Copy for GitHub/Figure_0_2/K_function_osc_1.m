%% setup

%mr clean
clc
clf

%ODE options
format long
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

%time
N = 2048*4+1;
Tmax = 6.710088675031566;
t = 0:Tmax/(N-1):Tmax;

fac = 0.993064415014497;

%define order of Taylor expansion
L = 16;

%Jacobian
Jac = @(t,u) [fac-3*u(1)^2, -fac; fac, 0];


%% dynamics

%function
F = @(t,u) fac*[u(1)-u(1)^3-u(2); u(1)];

%solve it
x0 = [0.823222 -1.02267];
[~,U] = ode113(F,t,x0,opts);

%find the period
T = Tmax;


%% define phase and K_0(theta)

%K_0(theta)
K_0 = [U(1:end,1)'; U(1:end,2)'];


%% variational equation

%system
f_var = @(t,u) fac*[u(1)-u(1)^3-u(2); ...
    u(1); ...
    (1-3*u(1)^2)*u(3)-u(5); ...
    (1-3*u(1)^2)*u(4)-u(6); ...
    u(3); ...
    u(4)];

%solve it
[t,u] = ode113(f_var,t,[x0 1 0 0 1],opts);

%construct M
M=zeros(2,2);
M(1,1) = u(end,3);
M(1,2) = u(end,4);
M(2,1) = u(end,5);
M(2,2) = u(end,6);

%diagonalize it
[V,lambda] = eig(M);

%obtain eigenvalues mu and then exponents lambda
mu_1 = lambda(1,1);
lambda_1 = 1/T*log(mu_1);

%obtain corresponding eigenvectors
v_1 = V(:,1);


%% Floquet Matrices and Stuff (C and R and Q)

%diagonalize M (actually we already did it xD)
C = [V(:,2), V(:,1)];

%define R
R = 1/T*C*[0 0; 0 log(mu_1)]/C;

%compute Q(t) for each t
PHI=zeros(2,2);
Q = zeros(2,2,length(t));

for i=1:length(t)
    
    %fundamental matrix
    PHI(1,1) = u(i,3);
    PHI(1,2) = u(i,4);
    PHI(2,1) = u(i,5);
    PHI(2,2) = u(i,6);
    
    %inner matrix
    D = [1 0; 0 exp(-t(i)*lambda_1)];
    
    %compute
    Q(:,:,i) = (PHI*C*D)/C;
    
end


%% compute K_1

%initialize
PHI = zeros(2,2);
K_1 = zeros(2,length(t));

%define numerical nicety coefficients
b = -.33;

%loop
for i=1:length(t)
    
    %fundamental matrix
    PHI(1,1) = u(i,3);
    PHI(1,2) = u(i,4);
    PHI(2,1) = u(i,5);
    PHI(2,2) = u(i,6);
    
    %compute K_10
    K_1(:,i) = PHI*exp(-lambda_1*t(i))*v_1*b;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_1 = zeros(2,length(t));

%numerical derivatives
dK1 = [gradient(K_1(1,:),Tmax/(N)); gradient(K_1(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_01
    LHS_1 = dK1(:,i) + lambda_1*K_1(:,i);
    RHS_1 = Jac(t(i),K_0(:,i))*K_1(:,i);
    error_1(:,i) = (LHS_1-RHS_1);
    
end

%compute L1 norm
E_1_x = 1/N*sum(abs(error_1(1,:)));
E_1_y = 1/N*sum(abs(error_1(2,:)));


%% automatic differentiation ... m = 2

%organize our (known) data
coeff_x = zeros(L,length(t));
coeff_y = zeros(L,length(t));

%initialize with known values
coeff_x(1,:) = K_0(1,:);
coeff_y(1,:) = K_0(2,:);

coeff_x(2,:) = K_1(1,:);
coeff_y(2,:) = K_1(2,:);

%organize unknown data structure
comp_coeff_x = zeros(L,length(t));
comp_coeff_y = zeros(L,length(t));

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=2
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(3,i) = -temp_x_cubed;
    
end

%beta's (for m=2)
beta_2 = fac*[comp_coeff_x(3,:); comp_coeff_y(3,:)];


%% Solve the Homological Equations with Floquet (m = 2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_2 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_2(:,i) = C\((Q(:,:,i))\beta_2(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_2_f = [fftshift(fft(A_2(1,:))); fftshift(fft(A_2(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_2 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_2(1,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_1-0).*A_2_f(1,:);
u_2(2,:) = 1./(2*pi*sqrt(-1).*k/T+2*lambda_1-lambda_1).*A_2_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_2 = [ifft(ifftshift(u_2(1,:))); ifft(ifftshift(u_2(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_2 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_2(:,i) = Q(:,:,i)*C*u_2(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_2 = zeros(2,length(t));

%numerical derivatives
dK2 = [gradient(K_2(1,:),Tmax/(N)); gradient(K_2(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_2 = dK2(:,i) + 2*lambda_1*K_2(:,i);
    RHS_2 = Jac(t(i),K_0(:,i))*K_2(:,i) + beta_2(:,i);
    error_2(:,i) = (LHS_2-RHS_2);
    
end

%compute L1 norm
E_2_x = 1/N*sum(abs(error_2(1,:)));
E_2_y = 1/N*sum(abs(error_2(2,:)));


%% automatic differentiation ... m = 3

%initialize with known values
coeff_x(3,:) = K_2(1,:);
coeff_y(3,:) = K_2(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=3
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(4,i) = -temp_x_cubed;
    
end

%beta's (for m=2)
beta_3 = fac*[comp_coeff_x(4,:); comp_coeff_y(4,:)];


%% Solve the Homological Equations with Floquet (m = 3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_3 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_3(:,i) = C\((Q(:,:,i))\beta_3(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_3_f = [fftshift(fft(A_3(1,:))); fftshift(fft(A_3(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_3 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_3(1,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_1-0).*A_3_f(1,:);
u_3(2,:) = 1./(2*pi*sqrt(-1).*k/T+3*lambda_1-lambda_1).*A_3_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_3 = [ifft(ifftshift(u_3(1,:))); ifft(ifftshift(u_3(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_3 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_3(:,i) = Q(:,:,i)*C*u_3(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_3 = zeros(2,length(t));

%numerical derivatives
dK3 = [gradient(K_3(1,:),Tmax/(N)); gradient(K_3(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_3 = dK3(:,i) + 3*lambda_1*K_3(:,i);
    RHS_3 = Jac(t(i),K_0(:,i))*K_3(:,i) + beta_3(:,i);
    error_3(:,i) = (LHS_3-RHS_3);
    
end

%compute L1 norm
E_3_x = 1/N*sum(abs(error_3(1,:)));
E_3_y = 1/N*sum(abs(error_3(2,:)));


%% automatic differentiation ... m = 4

%initialize with known values
coeff_x(4,:) = K_3(1,:);
coeff_y(4,:) = K_3(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=4
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(5,i) = -temp_x_cubed;
    
end

%beta's (for m=2)
beta_4 = fac*[comp_coeff_x(5,:); comp_coeff_y(5,:)];


%% Solve the Homological Equations with Floquet (m = 4)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_4 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_4(:,i) = C\((Q(:,:,i))\beta_4(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_4_f = [fftshift(fft(A_4(1,:))); fftshift(fft(A_4(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_4 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_4(1,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_1-0).*A_4_f(1,:);
u_4(2,:) = 1./(2*pi*sqrt(-1).*k/T+4*lambda_1-lambda_1).*A_4_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_4 = [ifft(ifftshift(u_4(1,:))); ifft(ifftshift(u_4(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_4 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_4(:,i) = Q(:,:,i)*C*u_4(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_4 = zeros(2,length(t));

%numerical derivatives
dK4 = [gradient(K_4(1,:),Tmax/(N)); gradient(K_4(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_4 = dK4(:,i) + 4*lambda_1*K_4(:,i);
    RHS_4 = Jac(t(i),K_0(:,i))*K_4(:,i) + beta_4(:,i);
    error_4(:,i) = (LHS_4-RHS_4);
    
end

%compute L1 norm
E_4_x = 1/N*sum(abs(error_4(1,:)));
E_4_y = 1/N*sum(abs(error_4(2,:)));


%% automatic differentiation ... m = 5

%initialize with known values
coeff_x(5,:) = K_4(1,:);
coeff_y(5,:) = K_4(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=5
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(6,i) = -temp_x_cubed;
    
end

%beta's (for m=2)
beta_5 = fac*[comp_coeff_x(6,:); comp_coeff_y(6,:)];


%% Solve the Homological Equations with Floquet (m = 5)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_5 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_5(:,i) = C\((Q(:,:,i))\beta_5(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_5_f = [fftshift(fft(A_5(1,:))); fftshift(fft(A_5(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_5 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_5(1,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_1-0).*A_5_f(1,:);
u_5(2,:) = 1./(2*pi*sqrt(-1).*k/T+5*lambda_1-lambda_1).*A_5_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_5 = [ifft(ifftshift(u_5(1,:))); ifft(ifftshift(u_5(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_5 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_5(:,i) = Q(:,:,i)*C*u_5(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_5 = zeros(2,length(t));

%numerical derivatives
dK5 = [gradient(K_5(1,:),Tmax/(N)); gradient(K_5(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_5 = dK5(:,i) + 5*lambda_1*K_5(:,i);
    RHS_5 = Jac(t(i),K_0(:,i))*K_5(:,i) + beta_5(:,i);
    error_5(:,i) = (LHS_5-RHS_5);
    
end

%compute L1 norm
E_5_x = 1/N*sum(abs(error_5(1,:)));
E_5_y = 1/N*sum(abs(error_5(2,:)));


%% automatic differentiation ... m = 6

%initialize with known values
coeff_x(6,:) = K_5(1,:);
coeff_y(6,:) = K_5(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=6
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(7,i) = -temp_x_cubed;
    
end

%beta's
beta_6 = fac*[comp_coeff_x(7,:); comp_coeff_y(7,:)];


%% Solve the Homological Equations with Floquet (m = 6)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_6 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_6(:,i) = C\((Q(:,:,i))\beta_6(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_6_f = [fftshift(fft(A_6(1,:))); fftshift(fft(A_6(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_6 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_6(1,:) = 1./(2*pi*sqrt(-1).*k/T+6*lambda_1-0).*A_6_f(1,:);
u_6(2,:) = 1./(2*pi*sqrt(-1).*k/T+6*lambda_1-lambda_1).*A_6_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_6 = [ifft(ifftshift(u_6(1,:))); ifft(ifftshift(u_6(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_6 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_6(:,i) = Q(:,:,i)*C*u_6(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_6 = zeros(2,length(t));

%numerical derivatives
dK6 = [gradient(K_6(1,:),Tmax/(N)); gradient(K_6(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_6 = dK6(:,i) + 6*lambda_1*K_6(:,i);
    RHS_6 = Jac(t(i),K_0(:,i))*K_6(:,i) + beta_6(:,i);
    error_6(:,i) = (LHS_6-RHS_6);
    
end

%compute L1 norm
E_6_x = 1/N*sum(abs(error_6(1,:)));
E_6_y = 1/N*sum(abs(error_6(2,:)));


%% automatic differentiation ... m = 7

%initialize with known values
coeff_x(7,:) = K_6(1,:);
coeff_y(7,:) = K_6(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=7
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(8,i) = -temp_x_cubed;
    
end

%beta's
beta_7 = fac*[comp_coeff_x(8,:); comp_coeff_y(8,:)];


%% Solve the Homological Equations with Floquet (m = 7)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_7 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_7(:,i) = C\((Q(:,:,i))\beta_7(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_7_f = [fftshift(fft(A_7(1,:))); fftshift(fft(A_7(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_7 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_7(1,:) = 1./(2*pi*sqrt(-1).*k/T+7*lambda_1-0).*A_7_f(1,:);
u_7(2,:) = 1./(2*pi*sqrt(-1).*k/T+7*lambda_1-lambda_1).*A_7_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_7 = [ifft(ifftshift(u_7(1,:))); ifft(ifftshift(u_7(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_7 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_7(:,i) = Q(:,:,i)*C*u_7(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_7 = zeros(2,length(t));

%numerical derivatives
dK7 = [gradient(K_7(1,:),Tmax/(N)); gradient(K_7(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_7 = dK7(:,i) + 7*lambda_1*K_7(:,i);
    RHS_7 = Jac(t(i),K_0(:,i))*K_7(:,i) + beta_7(:,i);
    error_7(:,i) = (LHS_7-RHS_7);
    
end

%compute L1 norm
E_7_x = 1/N*sum(abs(error_7(1,:)));
E_7_y = 1/N*sum(abs(error_7(2,:)));


%% automatic differentiation ... m = 8

%initialize with known values
coeff_x(8,:) = K_7(1,:);
coeff_y(8,:) = K_7(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=8
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(9,i) = -temp_x_cubed;
    
end

%beta's
beta_8 = fac*[comp_coeff_x(9,:); comp_coeff_y(9,:)];


%% Solve the Homological Equations with Floquet (m = 8)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_8 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_8(:,i) = C\((Q(:,:,i))\beta_8(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_8_f = [fftshift(fft(A_8(1,:))); fftshift(fft(A_8(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_8 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_8(1,:) = 1./(2*pi*sqrt(-1).*k/T+8*lambda_1-0).*A_8_f(1,:);
u_8(2,:) = 1./(2*pi*sqrt(-1).*k/T+8*lambda_1-lambda_1).*A_8_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_8 = [ifft(ifftshift(u_8(1,:))); ifft(ifftshift(u_8(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_8 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_8(:,i) = Q(:,:,i)*C*u_8(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_8 = zeros(2,length(t));

%numerical derivatives
dK8 = [gradient(K_8(1,:),Tmax/(N)); gradient(K_8(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_8 = dK8(:,i) + 8*lambda_1*K_8(:,i);
    RHS_8 = Jac(t(i),K_0(:,i))*K_8(:,i) + beta_8(:,i);
    error_8(:,i) = (LHS_8-RHS_8);
    
end

%compute L1 norm
E_8_x = 1/N*sum(abs(error_8(1,:)));
E_8_y = 1/N*sum(abs(error_8(2,:)));


%% automatic differentiation ... m = 9

%initialize with known values
coeff_x(9,:) = K_8(1,:);
coeff_y(9,:) = K_8(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=9
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(10,i) = -temp_x_cubed;
    
end

%beta's
beta_9 = fac*[comp_coeff_x(10,:); comp_coeff_y(10,:)];


%% Solve the Homological Equations with Floquet (m = 9)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_9 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_9(:,i) = C\((Q(:,:,i))\beta_9(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_9_f = [fftshift(fft(A_9(1,:))); fftshift(fft(A_9(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_9 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_9(1,:) = 1./(2*pi*sqrt(-1).*k/T+9*lambda_1-0).*A_9_f(1,:);
u_9(2,:) = 1./(2*pi*sqrt(-1).*k/T+9*lambda_1-lambda_1).*A_9_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_9 = [ifft(ifftshift(u_9(1,:))); ifft(ifftshift(u_9(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_9 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_9(:,i) = Q(:,:,i)*C*u_9(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_9 = zeros(2,length(t));

%numerical derivatives
dK9 = [gradient(K_9(1,:),Tmax/(N)); gradient(K_9(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_9 = dK9(:,i) + 9*lambda_1*K_9(:,i);
    RHS_9 = Jac(t(i),K_0(:,i))*K_9(:,i) + beta_9(:,i);
    error_9(:,i) = (LHS_9-RHS_9);
    
end

%compute L1 norm
E_9_x = 1/N*sum(abs(error_9(1,:)));
E_9_y = 1/N*sum(abs(error_9(2,:)));


%% automatic differentiation ... m = 10

%initialize with known values
coeff_x(10,:) = K_9(1,:);
coeff_y(10,:) = K_9(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=10
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(11,i) = -temp_x_cubed;
    
end

%beta's
beta_10 = fac*[comp_coeff_x(11,:); comp_coeff_y(11,:)];


%% Solve the Homological Equations with Floquet (m = 10)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_10 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_10(:,i) = C\((Q(:,:,i))\beta_10(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_10_f = [fftshift(fft(A_10(1,:))); fftshift(fft(A_10(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_10 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_10(1,:) = 1./(2*pi*sqrt(-1).*k/T+10*lambda_1-0).*A_10_f(1,:);
u_10(2,:) = 1./(2*pi*sqrt(-1).*k/T+10*lambda_1-lambda_1).*A_10_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_10 = [ifft(ifftshift(u_10(1,:))); ifft(ifftshift(u_10(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_10 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_10(:,i) = Q(:,:,i)*C*u_10(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_10 = zeros(2,length(t));

%numerical derivatives
dK10 = [gradient(K_10(1,:),Tmax/(N)); gradient(K_10(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_10 = dK10(:,i) + 10*lambda_1*K_10(:,i);
    RHS_10 = Jac(t(i),K_0(:,i))*K_10(:,i) + beta_10(:,i);
    error_10(:,i) = (LHS_10-RHS_10);
    
end

%compute L1 norm
E_10_x = 1/N*sum(abs(error_10(1,:)));
E_10_y = 1/N*sum(abs(error_10(2,:)));


%% automatic differentiation ... m = 11

%initialize with known values
coeff_x(11,:) = K_10(1,:);
coeff_y(11,:) = K_10(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=11
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(12,i) = -temp_x_cubed;
    
end

%beta's
beta_11 = fac*[comp_coeff_x(12,:); comp_coeff_y(12,:)];


%% Solve the Homological Equations with Floquet (m = 11)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_11 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_11(:,i) = C\((Q(:,:,i))\beta_11(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_11_f = [fftshift(fft(A_11(1,:))); fftshift(fft(A_11(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_11 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_11(1,:) = 1./(2*pi*sqrt(-1).*k/T+11*lambda_1-0).*A_11_f(1,:);
u_11(2,:) = 1./(2*pi*sqrt(-1).*k/T+11*lambda_1-lambda_1).*A_11_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_11 = [ifft(ifftshift(u_11(1,:))); ifft(ifftshift(u_11(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_11 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_11(:,i) = Q(:,:,i)*C*u_11(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_11 = zeros(2,length(t));

%numerical derivatives
dK11 = [gradient(K_11(1,:),Tmax/(N)); gradient(K_11(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_11 = dK11(:,i) + 11*lambda_1*K_11(:,i);
    RHS_11 = Jac(t(i),K_0(:,i))*K_11(:,i) + beta_11(:,i);
    error_11(:,i) = (LHS_11-RHS_11);
    
end

%compute L1 norm
E_11_x = 1/N*sum(abs(error_11(1,:)));
E_11_y = 1/N*sum(abs(error_11(2,:)));


%% automatic differentiation ... m = 12

%initialize with known values
coeff_x(12,:) = K_11(1,:);
coeff_y(12,:) = K_11(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=12
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(13,i) = -temp_x_cubed;
    
end

%beta's
beta_12 = fac*[comp_coeff_x(13,:); comp_coeff_y(13,:)];


%% Solve the Homological Equations with Floquet (m = 12)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_12 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_12(:,i) = C\((Q(:,:,i))\beta_12(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_12_f = [fftshift(fft(A_12(1,:))); fftshift(fft(A_12(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_12 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_12(1,:) = 1./(2*pi*sqrt(-1).*k/T+12*lambda_1-0).*A_12_f(1,:);
u_12(2,:) = 1./(2*pi*sqrt(-1).*k/T+12*lambda_1-lambda_1).*A_12_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_12 = [ifft(ifftshift(u_12(1,:))); ifft(ifftshift(u_12(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_12 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_12(:,i) = Q(:,:,i)*C*u_12(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_12 = zeros(2,length(t));

%numerical derivatives
dK12 = [gradient(K_12(1,:),Tmax/(N)); gradient(K_12(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_12 = dK12(:,i) + 12*lambda_1*K_12(:,i);
    RHS_12 = Jac(t(i),K_0(:,i))*K_12(:,i) + beta_12(:,i);
    error_12(:,i) = (LHS_12-RHS_12);
    
end

%compute L1 norm
E_12_x = 1/N*sum(abs(error_12(1,:)));
E_12_y = 1/N*sum(abs(error_12(2,:)));


%% automatic differentiation ... m = 13

%initialize with known values
coeff_x(13,:) = K_12(1,:);
coeff_y(13,:) = K_12(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=13
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(14,i) = -temp_x_cubed;
    
end

%beta's
beta_13 = fac*[comp_coeff_x(14,:); comp_coeff_y(14,:)];


%% Solve the Homological Equations with Floquet (m = 13)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_13 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_13(:,i) = C\((Q(:,:,i))\beta_13(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_13_f = [fftshift(fft(A_13(1,:))); fftshift(fft(A_13(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_13 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_13(1,:) = 1./(2*pi*sqrt(-1).*k/T+13*lambda_1-0).*A_13_f(1,:);
u_13(2,:) = 1./(2*pi*sqrt(-1).*k/T+13*lambda_1-lambda_1).*A_13_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_13 = [ifft(ifftshift(u_13(1,:))); ifft(ifftshift(u_13(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_13 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_13(:,i) = Q(:,:,i)*C*u_13(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_13 = zeros(2,length(t));

%numerical derivatives
dK13 = [gradient(K_13(1,:),Tmax/(N)); gradient(K_13(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_13 = dK13(:,i) + 13*lambda_1*K_13(:,i);
    RHS_13 = Jac(t(i),K_0(:,i))*K_13(:,i) + beta_13(:,i);
    error_13(:,i) = (LHS_13-RHS_13);
    
end

%compute L1 norm
E_13_x = 1/N*sum(abs(error_13(1,:)));
E_13_y = 1/N*sum(abs(error_13(2,:)));


%% automatic differentiation ... m = 14

%initialize with known values
coeff_x(14,:) = K_13(1,:);
coeff_y(14,:) = K_13(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=14
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(15,i) = -temp_x_cubed;
    
end

%beta's
beta_14 = fac*[comp_coeff_x(15,:); comp_coeff_y(15,:)];


%% Solve the Homological Equations with Floquet (m = 14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_14 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_14(:,i) = C\((Q(:,:,i))\beta_14(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_14_f = [fftshift(fft(A_14(1,:))); fftshift(fft(A_14(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_14 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_14(1,:) = 1./(2*pi*sqrt(-1).*k/T+14*lambda_1-0).*A_14_f(1,:);
u_14(2,:) = 1./(2*pi*sqrt(-1).*k/T+14*lambda_1-lambda_1).*A_14_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_14 = [ifft(ifftshift(u_14(1,:))); ifft(ifftshift(u_14(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_14 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_14(:,i) = Q(:,:,i)*C*u_14(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_14 = zeros(2,length(t));

%numerical derivatives
dK14 = [gradient(K_14(1,:),Tmax/(N)); gradient(K_14(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_14 = dK14(:,i) + 14*lambda_1*K_14(:,i);
    RHS_14 = Jac(t(i),K_0(:,i))*K_14(:,i) + beta_14(:,i);
    error_14(:,i) = (LHS_14-RHS_14);
    
end

%compute L1 norm
E_14_x = 1/N*sum(abs(error_14(1,:)));
E_14_y = 1/N*sum(abs(error_14(2,:)));


%% automatic differentiation ... m = 15

%initialize with known values
coeff_x(15,:) = K_14(1,:);
coeff_y(15,:) = K_14(2,:);

%for each Theta value
for i=1:length(t)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   y component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %nothing to do :)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %   x component setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %do it (compute coefficients for x^2)
    temp_x_square = zeros(1,L);
    
    for k=1:L-1
        for j=0:k
            temp_x_square(k+1) = temp_x_square(k+1) + coeff_x(j+1,i)*coeff_x(k-j+1,i);
        end
    end
    
    %do it again (compute coefficients for x^2*x)
    temp_x_cubed = 0;
    
    for k=15
        for j=0:k
            temp_x_cubed = temp_x_cubed + coeff_x(j+1,i)*temp_x_square(k-j+1);
        end
    end
    
    %store it
    comp_coeff_x(16,i) = -temp_x_cubed;
    
end

%beta's
beta_15 = fac*[comp_coeff_x(16,:); comp_coeff_y(16,:)];


%% Solve the Homological Equations with Floquet (m = 15)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   step 1 - solve in real space for the function A_(a,m-a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize
A_15 = zeros(2,length(t));

%solve
for i=1:length(t)
    A_15(:,i) = C\((Q(:,:,i))\beta_15(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 2 - solve for the function A_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%solve
A_15_f = [fftshift(fft(A_15(1,:))); fftshift(fft(A_15(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 3 - solve for the function u_(a,m-a) in Fourier space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
u_15 = zeros(2,length(t));

%define k vector
k=0:N-1;
k=k-(N-1)/2;

%solve
u_15(1,:) = 1./(2*pi*sqrt(-1).*k/T+15*lambda_1-0).*A_15_f(1,:);
u_15(2,:) = 1./(2*pi*sqrt(-1).*k/T+15*lambda_1-lambda_1).*A_15_f(2,:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  step 4 - solve for the function u_(a,m-a) in real space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%invert it
u_15 = [ifft(ifftshift(u_15(1,:))); ifft(ifftshift(u_15(2,:)))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                   step 5 - solve for K
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
K_15 = zeros(2,length(t));

%compute
for i=1:length(t)
    K_15(:,i) = Q(:,:,i)*C*u_15(:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%        error
%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize
error_15 = zeros(2,length(t));

%numerical derivatives
dK15 = [gradient(K_15(1,:),Tmax/(N)); gradient(K_15(2,:),Tmax/(N))];

%compute difference E
for i=1:length(t)
    
    %K_02
    LHS_15 = dK15(:,i) + 15*lambda_1*K_15(:,i);
    RHS_15 = Jac(t(i),K_0(:,i))*K_15(:,i) + beta_15(:,i);
    error_15(:,i) = (LHS_15-RHS_15);
    
end

%compute L1 norm
E_15_x = 1/N*sum(abs(error_15(1,:)));
E_15_y = 1/N*sum(abs(error_15(2,:)));


%% plot
%
% %k0
% figure(1)
% hold on
% plot(t/Tmax,K_0(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_0(2,:),'m','LineWidth',4)
% title('K_{0}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k1
% figure(2)
% hold on
% plot(t/Tmax,K_1(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_1(2,:),'m','LineWidth',4)
% title('K_{1}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k2
% figure(3)
% hold on
% plot(t/Tmax,K_2(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_2(2,:),'m','LineWidth',4)
% title('K_{2}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k3
% figure(4)
% hold on
% plot(t/Tmax,K_3(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_3(2,:),'m','LineWidth',4)
% title('K_{3}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k4
% figure(5)
% hold on
% plot(t/Tmax,K_4(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_4(2,:),'m','LineWidth',4)
% title('K_{4}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k5
% figure(6)
% hold on
% plot(t/Tmax,K_5(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_5(2,:),'m','LineWidth',4)
% title('K_{5}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k6
% figure(7)
% hold on
% plot(t/Tmax,K_6(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_6(2,:),'m','LineWidth',4)
% title('K_{6}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k7
% figure(8)
% hold on
% plot(t/Tmax,K_7(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_7(2,:),'m','LineWidth',4)
% title('K_{7}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k8
% figure(9)
% hold on
% plot(t/Tmax,K_8(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_8(2,:),'m','LineWidth',4)
% title('K_{8}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k9
% figure(10)
% hold on
% plot(t/Tmax,K_9(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_9(2,:),'m','LineWidth',4)
% title('K_{9}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k10
% figure(11)
% hold on
% plot(t/Tmax,K_10(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_10(2,:),'m','LineWidth',4)
% title('K_{10}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k11
% figure(12)
% hold on
% plot(t/Tmax,K_11(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_11(2,:),'m','LineWidth',4)
% title('K_{11}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k12
% figure(13)
% hold on
% plot(t/Tmax,K_12(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_12(2,:),'m','LineWidth',4)
% title('K_{12}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k13
% figure(14)
% hold on
% plot(t/Tmax,K_13(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_13(2,:),'m','LineWidth',4)
% title('K_{13}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k14
% figure(15)
% hold on
% plot(t/Tmax,K_14(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_14(2,:),'m','LineWidth',4)
% title('K_{14}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')
%
%
% %k15
% figure(16)
% hold on
% plot(t/Tmax,K_15(1,:),'k','LineWidth',4)
% plot(t/Tmax,K_15(2,:),'m','LineWidth',4)
% title('K_{15}')
% box on
% axis square
% set(gca,'fontsize',15)
% ylabel('K')
% xlabel('Phase \theta')
% legend('x - comp','y - comp')

%all together to check with paper
figure(1)
hold on
plot(t/Tmax,circshift(K_1(1,:),0),'r','LineWidth',2)
plot(t/Tmax,circshift(K_5(1,:),0),'g','LineWidth',2)
plot(t/Tmax,circshift(K_10(1,:),0),'b','LineWidth',2)
plot(t/Tmax,circshift(K_15(1,:),0),'m','LineWidth',2)
title('Kn')
box on
axis square
set(gca,'fontsize',15)
ylabel('K')
xlabel('Phase \theta')
legend('K1','K5','K10','K15')


%% save data

save('K0_vdp_1.mat','K_0')
save('K1_vdp_1.mat','K_1')
save('K2_vdp_1.mat','K_2')
save('K3_vdp_1.mat','K_3')
save('K4_vdp_1.mat','K_4')
save('K5_vdp_1.mat','K_5')
save('K6_vdp_1.mat','K_6')
save('K7_vdp_1.mat','K_7')
save('K8_vdp_1.mat','K_8')
save('K9_vdp_1.mat','K_9')
save('K10_vdp_1.mat','K_10')
save('K11_vdp_1.mat','K_11')
save('K12_vdp_1.mat','K_12')
save('K13_vdp_1.mat','K_13')
save('K14_vdp_1.mat','K_14')
save('K15_vdp_1.mat','K_15')



