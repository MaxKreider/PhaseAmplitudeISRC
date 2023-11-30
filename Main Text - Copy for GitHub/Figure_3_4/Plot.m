%% setup

%mr clean
clc
clf

%for the ODEs
format long
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

%parameters
eps = 0.3;

%time
Tmax = 6.67;
N = 2048+1;
dt = Tmax/(N-1);
t = 0:dt:Tmax;

%initial conditions
x0 = [-0.496089165791201 1.171912243464496];


%% load data

Z0 = load('Z0K_vdp_1.mat');
Z0 = cell2mat(struct2cell(Z0));


%% unperturbed dynamics

%unperturbed vector field
F_u = @(t,u) [u(1)-u(1)^3-u(2); u(1)];

%solve
[~,U]=ode113(F_u,t,x0,opts);

%ease of notation
v_u = U(:,1);
n_u = U(:,2);

% %compute unperturbed period T_p
% [~,loc]=findpeaks(v_u,'MinPeakHeight',0);
% timerz = T(loc);
% T_u = mean(diff(timerz))


%% perturbed dynamics

%cheat
Tp = 7.792548033248546;

%perturbed vector field
F_p = @(t,u) Tp/Tmax*[u(1)-u(1)^3-u(2)+eps*u(1); u(1)+eps*u(2)];

%solve
x01 = [-0.759816621966552   1.942163942081037];
[~,U]=ode113(F_p,0:dt:1*Tmax,x01,opts);

%ease of notation
v_p = U(:,1);
n_p = U(:,2);

% %compute perturbed period T_p
% [~,loc]=findpeaks(v_p,'MinPeakHeight',0);
% timerz = T(loc);
% T_p = mean(diff(timerz(round(length(loc)/2):end)))


%% P section

%define variational equation to find floquet multipliers and vectors
F_v = @(t,u) [u(1)-u(1)^3-u(2); ...
    u(1); ...
    (1-3*u(1)^2)*u(3)-u(5); ...
    (1-3*u(1)^2)*u(4)-u(6); ...
    u(3); ...
    u(4)];

%solve variation equation
[T,U] = ode113(F_v,t,[x0 1 0 0 1],opts);

%Monodromy matrix
M=zeros(2,2);
M(1,1) = U(end,3);
M(1,2) = U(end,4);
M(2,1) = U(end,5);
M(2,2) = U(end,6);

%find floquet coordinates (eigenvectors)
[v,lambda] = eig(M);

%matlab formatting lol
if abs(lambda(4))>.5
    v=flip(v,2);
end

%define new coordinate system, q2 is poincare section
q1 = v(:,1);
q2=null(q1');

%choose several poincare sections for testing
theta = .75;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
q2 = R*q2;
n=null(q2');


%% nu1

%form integrand
integrand = zeros(1,length(t));
for i=1:length(t)
    integrand(i) = dot(Z0(:,i),[1*v_u(i);1*n_u(i)]);
end

%integrate
T1 = cumtrapz(t,integrand);
T1 = -T1(end);

%find nu
nu1 = T1/Tmax;


%% IC

%set IC
IC = [1 0 0 1];

%RHS
FM = [];
for j=1:length(U)
    
    %create df/deps
    dfdeps = [v_u(j);n_u(j)];
    
    %create fundamental matrix solution
    PHI = [U(j,3), U(j,4); U(j,5), U(j,6)];
    
    %create vector slot for FM
    temp = PHI\(dfdeps+nu1*F_u(0,[U(j,1), U(j,2)]));
    
    %FM
    FM = [FM temp];
    
end

%integrate
FM = [cumtrapz(dt,FM(1,:));cumtrapz(dt,FM(2,:))];

%pls god pls
term3=[];
for i=1:length(U)
    
    %define PHI
    PHI = [U(j,3), U(j,4); U(j,5), U(j,6)];
    
    %sol
    term3(:,i) = PHI*FM(:,i);
    
end

%end point
FM = [term3(1,end); term3(2,end)];

%setup RHS
test=M-eye(2);
test=[test; n(1) n(2)];
FM = -[FM; 0];

%solve
p1_lul = test\FM;


%% iSRC

%iSRC equation
ISRC = @(t,u) [u(1)-u(1)^3-u(2); ...
               u(1);
               (1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1);
               (1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2)];

%skip a few cycles
tran=1;

%initial cond
IC=p1_lul';

%solve
for k=1:tran
    
    %solve
    [T_src,U_src] = ode113(ISRC,0:dt:Tmax,[x0 IC],opts);
    
    %periodicity shooting thingy
    IC = U_src(end,3:4);
    
end


%% 1st order LC

%perturbed LC
gamma_eps_v1 = v_u + U_src(:,3)*eps;
gamma_eps_n1 = n_u + U_src(:,4)*eps;


%% nu2

%iSRC derivative
gamma_1_p1 = gradient(U_src(:,3),dt);
gamma_1_p2 = gradient(U_src(:,4),dt);

%H
H = zeros(2,length(t));
for i=1:length(t)
    H(1,i) = 1/2*kron(U_src(i,3:4),U_src(i,3:4))*[-6*v_u(i);0;0;0];
    H(2,i) = 1/2*kron(U_src(i,3:4),U_src(i,3:4))*[0;0;0;0];
end

%D
D = zeros(2,length(t));
for i=1:length(t)
    D(1,i) = U_src(i,3);
    D(2,i) = U_src(i,4);
end

%form lol
integrand = zeros(1,length(t));
lol = zeros(2,length(t));
for i=1:length(t)
    lol(:,i) = H(:,i)+D(:,i)+nu1*[gamma_1_p1(i);gamma_1_p2(i)];
    integrand(i) = dot(Z0(:,i),lol(:,i));
end

%integrate
T2 = cumtrapz(t,integrand);
T2 = -T2(end);

%find nu
nu2 = T2/Tmax;


%% IC 

%set IC
IC = [1 0 0 1];

%RHS
FM = [];
for j=1:length(U)
        
    %create fundamental matrix solution
    PHI = [U(j,3), U(j,4); U(j,5), U(j,6)];
    
    %create vector slot for FM
    temp = PHI\(lol(:,j)+nu2*F_u(0,[U(j,1), U(j,2)]));
    
    %FM
    FM = [FM temp];
    
end

%integrate
FM = [cumtrapz(dt,FM(1,:));cumtrapz(dt,FM(2,:))];

%pls god pls
term3=[];
for i=1:length(U)
    
    %define PHI
    PHI = [U(j,3), U(j,4); U(j,5), U(j,6)];
    
    %sol
    term3(:,i) = PHI*FM(:,i);
    
end

%end point
FM = [term3(1,end); term3(2,end)];

%setup RHS
test=M-eye(2);
test=[test; n(1) n(2)];
FM = -[FM; 0];

%solve
p2_lul = test\FM;


%% iSRC

%iSRC equation
ISRC = @(t,u) [u(1)-u(1)^3-u(2); ...
               u(1);
               
               (1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1);
               (1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2);
               
               (1-3*u(1)^2)*u(5)+(-1)*u(6) + nu1*((1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1)) + nu2*(u(1)-u(1)^3-u(2)) + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[-6*u(1);0;0;0] + u(3);
               (1)*u(5)+(0)*u(6)           + nu1*((1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2))                       + nu2*(u(1))             + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[0;0;0;0]       + u(4)];

%skip a few cycles
tran=1;

%initial cond
IC=[p1_lul' p2_lul'];

%solve
for k=1:tran
    
    %solve
    [T_src2,U_src2] = ode113(ISRC,0:dt:Tmax,[x0 IC],opts);
    
    %periodicity shooting thingy
    IC = U_src2(end,3:6);
    
end


%% 2st order LC

%perturbed LC
gamma_eps_v2 = v_u + U_src2(:,3)*eps + U_src2(:,5)*eps^2;
gamma_eps_n2 = n_u + U_src2(:,4)*eps + U_src2(:,6)*eps^2;


%% nu3

%iSRC derivative
gamma_2_p1 = gradient(U_src2(:,5),dt);
gamma_2_p2 = gradient(U_src2(:,6),dt);

%H
H = zeros(2,length(t));
for i=1:length(t)
    H(1,i) = 1/2*(kron(U_src(i,3:4),U_src2(i,5:6))*[-6*v_u(i);0;0;0]+kron(U_src2(i,5:6),U_src(i,3:4))*[-6*v_u(i);0;0;0]) + ...
             1/6*kron(kron(U_src(i,3:4),U_src(i,3:4)),U_src(i,3:4))*[-6, 0, 0, 0, 0, 0, 0, 0]';
    H(2,i) = 1/2*(kron(U_src(i,3:4),U_src2(i,5:6))*[0;0;0;0]+kron(U_src2(i,5:6),U_src(i,3:4))*[0;0;0;0]) + ...
             1/6*kron(kron(U_src(i,3:4),U_src(i,3:4)),U_src(i,3:4))*[0, 0, 0, 0, 0, 0, 0, 0]';     
end

%D
D = zeros(2,length(t));
for i=1:length(t)
    D(1,i) = U_src2(i,5);
    D(2,i) = U_src2(i,6);
end

%form lol
integrand = zeros(1,length(t));
lol = zeros(2,length(t));
for i=1:length(t)
    lol(:,i) = H(:,i)+D(:,i)+nu1*[gamma_2_p1(i);gamma_2_p2(i)]+nu2*[gamma_1_p1(i);gamma_1_p2(i)];
    integrand(i) = dot(Z0(:,i),lol(:,i));
end

%integrate
T3 = cumtrapz(t,integrand);
T3 = -T3(end);

%find nu
nu3 = T3/Tmax;


%% IC 

%set IC
IC = [1 0 0 1];

%RHS
FM = [];
for j=1:length(U)
        
    %create fundamental matrix solution
    PHI = [U(j,3), U(j,4); U(j,5), U(j,6)];
    
    %create vector slot for FM
    temp = PHI\(lol(:,j)+nu3*F_u(0,[U(j,1), U(j,2)]));
    
    %FM
    FM = [FM temp];
    
end

%integrate
FM = [cumtrapz(dt,FM(1,:));cumtrapz(dt,FM(2,:))];

%pls god pls
term3=[];
for i=1:length(U)
    
    %define PHI
    PHI = [U(j,3), U(j,4); U(j,5), U(j,6)];
    
    %sol
    term3(:,i) = PHI*FM(:,i);
    
end

%end point
FM = [term3(1,end); term3(2,end)];

%setup RHS
test=M-eye(2);
test=[test; n(1) n(2)];
FM = -[FM; 0];

%solve
p3_lul = test\FM;


%% iSRC

%iSRC equation
ISRC = @(t,u) [u(1)-u(1)^3-u(2); ...
               u(1);
               
               (1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1);
               (1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2);
               
               (1-3*u(1)^2)*u(5)+(-1)*u(6) + nu1*((1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1)) + nu2*(u(1)-u(1)^3-u(2)) + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[-6*u(1);0;0;0] + u(3);
               (1)*u(5)+(0)*u(6)           + nu1*((1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2))                       + nu2*(u(1))             + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[0;0;0;0]       + u(4);
               
               (1-3*u(1)^2)*u(7)+(-1)*u(8) + ...
               nu1*((1-3*u(1)^2)*u(5)+(-1)*u(6) + nu1*((1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1)) + nu2*(u(1)-u(1)^3-u(2)) + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[-6*u(1);0;0;0] + u(3)) + ...
               nu2*((1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1)) + ...
               nu3*(u(1)-u(1)^3-u(2)) + ...
               1/2*(kron([u(3);u(4)],[u(5);u(6)])'*[-6*u(1);0;0;0]+kron([u(5);u(6)],[u(3);u(4)])'*[-6*u(1);0;0;0]) + ...
               1/6*kron(kron([u(3);u(4)],[u(3);u(4)]),[u(3);u(4)])'*[-6, 0, 0, 0, 0, 0, 0, 0]' + ...
               u(5);
               (1)*u(7)+(0)*u(8) + ...
               nu1*((1)*u(5)+(0)*u(6)           + nu1*((1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2))                       + nu2*(u(1))             + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[0;0;0;0]       + u(4)) + ...
               nu2*((1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2)) + ...
               nu3*(u(1)) + ...
               1/2*(kron([u(3);u(4)],[u(5);u(6)])'*[0;0;0;0]+kron([u(5);u(6)],[u(3);u(4)])'*[0;0;0;0]) + ...
               1/6*kron(kron([u(3);u(4)],[u(3);u(4)]),[u(3);u(4)])'*[0, 0, 0, 0, 0, 0, 0, 0]' + ...
               u(6)];
           

%skip a few cycles
tran=1;

%initial cond
IC=[p1_lul' p2_lul' p3_lul'];

%solve
for k=1:tran
    
    %solve
    [T_src3,U_src3] = ode113(ISRC,0:dt:Tmax,[x0 IC],opts);
    
    %periodicity shooting thingy
    IC = U_src3(end,3:8);
    
end


%% 3st order LC

%perturbed LC
gamma_eps_v3 = v_u + U_src3(:,3)*eps + U_src3(:,5)*eps^2 + U_src3(:,7)*eps^3;
gamma_eps_n3 = n_u + U_src3(:,4)*eps + U_src3(:,6)*eps^2 + U_src3(:,8)*eps^3; 


%% nu4

%iSRC derivative
gamma_3_p1 = gradient(U_src3(:,7),dt);
gamma_3_p2 = gradient(U_src3(:,8),dt);

%H
H = zeros(2,length(t));
for i=1:length(t)
    H(1,i) = 1/2*(kron(U_src(i,3:4),U_src3(i,7:8))*[-6*v_u(i);0;0;0]+...
                  kron(U_src3(i,7:8),U_src(i,3:4))*[-6*v_u(i);0;0;0]+...
                  kron(U_src2(i,5:6),U_src2(i,5:6))*[-6*v_u(i);0;0;0])+...
             1/6*(kron(kron(U_src(i,3:4),U_src(i,3:4)),U_src2(i,5:6))*[-6, 0, 0, 0, 0, 0, 0, 0]'+...
                 kron(kron(U_src(i,3:4),U_src2(i,5:6)),U_src(i,3:4))*[-6, 0, 0, 0, 0, 0, 0, 0]'+...
                 kron(kron(U_src2(i,5:6),U_src(i,3:4)),U_src(i,3:4))*[-6, 0, 0, 0, 0, 0, 0, 0]');
    H(2,i) = 1/2*(kron(U_src(i,3:4),U_src3(i,7:8))*[0;0;0;0]+...
                 kron(U_src3(i,7:8),U_src(i,3:4))*[0;0;0;0]+ ...
                 kron(U_src2(i,5:6),U_src2(i,5:6))*[0;0;0;0])+ ...
             1/6*(kron(kron(U_src(i,3:4),U_src(i,3:4)),U_src2(i,5:6))*[0, 0, 0, 0, 0, 0, 0, 0]'+...
                 kron(kron(U_src(i,3:4),U_src2(i,5:6)),U_src(i,3:4))*[0, 0, 0, 0, 0, 0, 0, 0]'+...
                 kron(kron(U_src2(i,5:6),U_src(i,3:4)),U_src(i,3:4))*[0, 0, 0, 0, 0, 0, 0, 0]');     
end

%D
D = zeros(2,length(t));
for i=1:length(t)
    D(1,i) = U_src3(i,7);
    D(2,i) = U_src3(i,8);
end

%form lol
integrand = zeros(1,length(t));
lol = zeros(2,length(t));
for i=1:length(t)
    lol(:,i) = H(:,i)+D(:,i)+nu1*[gamma_3_p1(i);gamma_3_p2(i)]+nu2*[gamma_2_p1(i);gamma_2_p2(i)]+nu3*[gamma_1_p1(i);gamma_1_p2(i)];
    integrand(i) = dot(Z0(:,i),lol(:,i));
end

%integrate
T4 = cumtrapz(t,integrand);
T4 = -T4(end);

%find nu
nu4 = T4/Tmax;


%% IC 

%set IC
IC = [1 0 0 1];

%RHS
FM = [];
for j=1:length(U)
        
    %create fundamental matrix solution
    PHI = [U(j,3), U(j,4); U(j,5), U(j,6)];
    
    %create vector slot for FM
    temp = PHI\(lol(:,j)+nu4*F_u(0,[U(j,1), U(j,2)]));
    
    %FM
    FM = [FM temp];
    
end

%integrate
FM = [cumtrapz(dt,FM(1,:));cumtrapz(dt,FM(2,:))];

%pls god pls
term3=[];
for i=1:length(U)
    
    %define PHI
    PHI = [U(j,3), U(j,4); U(j,5), U(j,6)];
    
    %sol
    term3(:,i) = PHI*FM(:,i);
    
end

%end point
FM = [term3(1,end); term3(2,end)];

%setup RHS
test=M-eye(2);
test=[test; n(1) n(2)];
FM = -[FM; 0];

%solve
p4_lul = test\FM;


%% iSRC

%iSRC equation
ISRC = @(t,u) [u(1)-u(1)^3-u(2); ...
               u(1);
               
               (1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1);
               (1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2);
               
               (1-3*u(1)^2)*u(5)+(-1)*u(6) + nu1*((1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1)) + nu2*(u(1)-u(1)^3-u(2)) + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[-6*u(1);0;0;0] + u(3);
               (1)*u(5)+(0)*u(6)           + nu1*((1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2))                       + nu2*(u(1))             + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[0;0;0;0]       + u(4);
               
               (1-3*u(1)^2)*u(7)+(-1)*u(8) + ...
               nu1*((1-3*u(1)^2)*u(5)+(-1)*u(6) + nu1*((1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1)) + nu2*(u(1)-u(1)^3-u(2)) + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[-6*u(1);0;0;0] + u(3)) + ...
               nu2*((1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1)) + ...
               nu3*(u(1)-u(1)^3-u(2)) + ...
               1/2*(kron([u(3);u(4)],[u(5);u(6)])'*[-6*u(1);0;0;0]+kron([u(5);u(6)],[u(3);u(4)])'*[-6*u(1);0;0;0]) + ...
               1/6*kron(kron([u(3);u(4)],[u(3);u(4)]),[u(3);u(4)])'*[-6, 0, 0, 0, 0, 0, 0, 0]' + ...
               u(5);
               (1)*u(7)+(0)*u(8) + ...
               nu1*((1)*u(5)+(0)*u(6)           + nu1*((1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2))                       + nu2*(u(1))             + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[0;0;0;0]       + u(4)) + ...
               nu2*((1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2)) + ...
               nu3*(u(1)) + ...
               1/2*(kron([u(3);u(4)],[u(5);u(6)])'*[0;0;0;0]+kron([u(5);u(6)],[u(3);u(4)])'*[0;0;0;0]) + ...
               1/6*kron(kron([u(3);u(4)],[u(3);u(4)]),[u(3);u(4)])'*[0, 0, 0, 0, 0, 0, 0, 0]' + ...
               u(6);
               
               (1-3*u(1)^2)*u(9)+(-1)*u(10) + ...
               nu1*((1-3*u(1)^2)*u(7)+(-1)*u(8) + ...
               nu1*((1-3*u(1)^2)*u(5)+(-1)*u(6) + nu1*((1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1)) + nu2*(u(1)-u(1)^3-u(2)) + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[-6*u(1);0;0;0] + u(3)) + ...
               nu2*((1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1)) + ...
               nu3*(u(1)-u(1)^3-u(2)) + ...
               1/2*(kron([u(3);u(4)],[u(5);u(6)])'*[-6*u(1);0;0;0]+kron([u(5);u(6)],[u(3);u(4)])'*[-6*u(1);0;0;0]) + ...
               1/6*kron(kron([u(3);u(4)],[u(3);u(4)]),[u(3);u(4)])'*[-6, 0, 0, 0, 0, 0, 0, 0]' + ...
               u(5)) + ...
               nu2*((1-3*u(1)^2)*u(5)+(-1)*u(6) + nu1*((1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1)) + nu2*(u(1)-u(1)^3-u(2)) + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[-6*u(1);0;0;0] + u(3)) + ...
               nu3*((1-3*u(1)^2)*u(3)+(-1)*u(4) + nu1*(u(1)-u(1)^3-u(2)) + 1*u(1)) + ...
               nu4*(u(1)-u(1)^3-u(2)) + ...
               1/2*(kron([u(3);u(4)],[u(7);u(8)])'*[-6*u(1);0;0;0]+kron([u(7);u(8)],[u(3);u(4)])'*[-6*u(1);0;0;0] + kron([u(5);u(6)],[u(5);u(6)])'*[-6*u(1);0;0;0]) + ...
               1/6*(kron(kron([u(3);u(4)],[u(3);u(4)]),[u(5);u(6)])'*[-6, 0, 0, 0, 0, 0, 0, 0]' + kron(kron([u(3);u(4)],[u(5);u(6)]),[u(3);u(4)])'*[-6, 0, 0, 0, 0, 0, 0, 0]' + kron(kron([u(5);u(6)],[u(3);u(4)]),[u(3);u(4)])'*[-6, 0, 0, 0, 0, 0, 0, 0]') + ...
               u(7);
               (1)*u(9)+(0)*u(10) + ...
               nu1*((1)*u(7)+(0)*u(8) + ...
               nu1*((1)*u(5)+(0)*u(6)           + nu1*((1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2))                       + nu2*(u(1))             + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[0;0;0;0]       + u(4)) + ...
               nu2*((1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2)) + ...
               nu3*(u(1)) + ...
               1/2*(kron([u(3);u(4)],[u(5);u(6)])'*[0;0;0;0]+kron([u(5);u(6)],[u(3);u(4)])'*[0;0;0;0]) + ...
               1/6*kron(kron([u(3);u(4)],[u(3);u(4)]),[u(3);u(4)])'*[0, 0, 0, 0, 0, 0, 0, 0]' + ...
               u(6)) + ...
               nu2*((1)*u(5)+(0)*u(6)           + nu1*((1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2))                       + nu2*(u(1))             + 1/2*kron([u(3);u(4)],[u(3);u(4)])'*[0;0;0;0]       + u(4)) + ...
               nu3*((1)*u(3)+(0)*u(4) + nu1*(u(1)) + 1*u(2)) + ...
               nu4*(u(1)) + ...
               1/2*(kron([u(3);u(4)],[u(7);u(8)])'*[0;0;0;0]+kron([u(7);u(8)],[u(3);u(4)])'*[0;0;0;0] + kron([u(5);u(6)],[u(5);u(6)])'*[0;0;0;0]) + ...
               1/6*(kron(kron([u(3);u(4)],[u(3);u(4)]),[u(5);u(6)])'*[0, 0, 0, 0, 0, 0, 0, 0]' + kron(kron([u(3);u(4)],[u(5);u(6)]),[u(3);u(4)])'*[0, 0, 0, 0, 0, 0, 0, 0]' + kron(kron([u(5);u(6)],[u(3);u(4)]),[u(3);u(4)])'*[0, 0, 0, 0, 0, 0, 0, 0]') + ...
               u(8)];
           

%skip a few cycles
tran=1;

%initial cond
IC=[p1_lul' p2_lul' p3_lul' p4_lul'];

%solve
for k=1:tran
    
    %solve
    [T_src4,U_src4] = ode113(ISRC,0:dt:Tmax,[x0 IC],opts);
    
    %periodicity shooting thingy
    IC = U_src4(end,3:10);
    
end


%% 4st order LC

%perturbed LC
gamma_eps_v4 = v_u + U_src3(:,3)*eps + U_src3(:,5)*eps^2 + U_src3(:,7)*eps^3 + U_src4(:,9)*eps^4;
gamma_eps_n4 = n_u + U_src3(:,4)*eps + U_src3(:,6)*eps^2 + U_src3(:,8)*eps^3 + U_src4(:,10)*eps^4; 


%% plot

figure(1)
hold on
plot(t,v_p - gamma_eps_v1,'r','LineWidth',4)
plot(t,v_p - gamma_eps_v2,'g','LineWidth',4)
plot(t,v_p - gamma_eps_v3,'b','LineWidth',4)
plot(t,v_p - gamma_eps_v4,'c','LineWidth',4)
box on
axis square
xlim([0 Tmax])
xlabel('time t')
ylabel('Discrepancy \gamma(t;\epsilon) - iSRC')
legend('1','2','3','4')
title('x Error')
set(gca,'fontsize',15)

figure(2)
hold on
plot(t,n_p - gamma_eps_n1,'r','LineWidth',4)
plot(t,n_p - gamma_eps_n2,'g','LineWidth',4)
plot(t,n_p - gamma_eps_n3,'b','LineWidth',4)
plot(t,n_p - gamma_eps_n4,'c','LineWidth',4)
box on
axis square
xlim([0 Tmax])
xlabel('time t')
ylabel('Discrepancy \gamma(t;\epsilon) - iSRC')
legend('1','2','3','4')
title('y Error')
set(gca,'fontsize',15)

% xnorm=zeros(1,4);
% ynorm=zeros(1,4);
% 
% xnorm(1) = norm(abs(v_p - gamma_eps_v1));
% xnorm(2) = norm(abs(v_p - gamma_eps_v2)); 
% xnorm(3) = norm(abs(v_p - gamma_eps_v3)); 
% xnorm(4) = norm(abs(v_p - gamma_eps_v4)); 
% 
% ynorm(1) = norm(abs(n_p - gamma_eps_n1));
% ynorm(2) = norm(abs(n_p - gamma_eps_n2)); 
% ynorm(3) = norm(abs(n_p - gamma_eps_n3)); 
% ynorm(4) = norm(abs(n_p - gamma_eps_n4)); 
% 
% figure(71)
% hold on
% plot(1:4,xnorm,'k.','MarkerSize',40)
% plot(1:4,xnorm,'k-','LineWidth',3)
% plot(1:4,ynorm,'m.','MarkerSize',40)
% plot(1:4,ynorm,'m-','LineWidth',3)
% box on
% axis square
% xlabel('order')
% ylabel('abs(||discrepancy||)')
% title('L2 Error')
% set(gca,'fontsize',15)

figure(3)
hold on
plot(v_u,n_u,'color',[.7 .7 .7],'LineWidth',3)
plot(v_p,n_p,'k-','LineWidth',5)
plot(gamma_eps_v1,gamma_eps_n1,'r--','LineWidth',3)
plot(gamma_eps_v4,gamma_eps_n4,'c--','LineWidth',3)
box on
axis square
xlabel('x')
ylabel('y')
%legend('Unperturbed','Perturbed','1st order','4th order')
title('High Order LC Approximation')
set(gca,'fontsize',15)

