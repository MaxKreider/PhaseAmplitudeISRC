%% setup

%mr clean
clc
clf

%for the ODEs
format long
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

%time
Tmax = 6.663254495069035;
N = 2048+1;
dt = Tmax/(N-1);
t = 0:dt:Tmax;

Tp1 = 6.642127010587909;
T_s = 6.637406909544001;

%parameters
beta=.001;
delta=.1;
omega=2*pi/Tmax;

%choose 1st or 3rd order approximation (factor = 0 for 1st, factor = 1 for 3rd)
choice = 1;


%% base system

%unperturbed vector field
x0 = [-.575611 0.410697  -1.11779 0.672746];
F_p = @(t,u) [(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2); 
              (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
              u(3)-u(3)^3-u(4);
              u(3)];
          
%solve
[~,U]=ode113(F_p,0:dt:Tmax*1,x0,opts);

%ease of notation
v1 = U(:,1);
n1 = U(:,2);
v2 = U(:,3);
n2 = U(:,4);


%% perturbed system

%perturbed vector field
x0p = [-0.589511603659218   0.420915083744529  -1.089922687522716   0.663740923553011];
F_c = @(t,u)  Tp1/Tmax*[(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2) + delta*(u(3)-u(1)); 
              (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
              u(3)-u(3)^3-u(4) + delta*(u(1)-u(3));
              u(3)];

%solve
[~,U]=ode113(F_c,0:dt:Tmax*1,x0p,opts);

%ease of notation
v1p = U(:,1);
n1p = U(:,2);
v2p = U(:,3);
n2p = U(:,4);

% %compute unperturbed period T_p
% [~,loc]=findpeaks(v1p,'MinPeakHeight',0);
% timerz = T(loc);
% T_p1 = mean(diff(timerz))
% 
% [~,loc]=findpeaks(v2p,'MinPeakHeight',0);
% timerz = T(loc);
% T_p2 = mean(diff(timerz))


%% P section + nu + IC (hopf)

%define variational equation to find floquet multipliers and vectors
F_s = @(t,u) [(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2); 
              (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
               (.5-3*u(1)^2-u(2)^2)*u(3) + (-2*u(1)*u(2)-omega)*u(5);
               (.5-3*u(1)^2-u(2)^2)*u(4) + (-2*u(1)*u(2)-omega)*u(6);
               (-2*u(1)*u(2)+omega)*u(3) + (.5-u(1)^2-3*u(2)^2)*u(5);
               (-2*u(1)*u(2)+omega)*u(4) + (.5-u(1)^2-3*u(2)^2)*u(6)];

%solve variation equation
[~,U_s_hopf] = ode113(F_s,0:dt:Tmax,[x0(1:2) 1 0 0 1],opts);

%Monodromy matrix
M=zeros(2,2);
M(1,1) = U_s_hopf(end,3);
M(1,2) = U_s_hopf(end,4);
M(2,1) = U_s_hopf(end,5);
M(2,2) = U_s_hopf(end,6);
M_hopf = M;

%find floquet coordinates (eigenvectors)
[v,~] = eig(M);

%define new coordinate system, q2 is poincare section
q = v(:,1);

%choose several poincare sections for testing
theta = -.75;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
q = R*q;
n_hopf = null(q');

%plot
figure(1)
plot([x0(1)-q(1),x0(1)+q(1)],[x0(2)-q(2),x0(2)+q(2)],'color','g','LineWidth',1.5);

%Z0 
Z0_hopf_func = @(t,u) [(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2); 
                  (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
                  -(.5-3*u(1)^2-u(2)^2)*u(3) + -(-2*u(1)*u(2)+omega)*u(4);
                  -(-2*u(2)*u(1)-omega)*u(3) + -(.5-u(1)^2-3*u(2)^2)*u(4)];
                          
%solve it
la = F_p(0,x0);
IC = [x0(1:2) 0 2];
for i=1:3
    [~,U] = ode113(Z0_hopf_func,Tmax:-dt:0,IC,opts);
    IC = [x0(1:2) U(end,3:4)];
end
Z0_hopf = flip(U);
Z0_hopf = Z0_hopf(:,3:4);
Z0_hopf = Z0_hopf/(dot(Z0_hopf(1,:),la(1:2)));

%form integrand
integrand = zeros(1,length(0:dt:Tmax));
for i=1:length(0:dt:Tmax)
    integrand(i) = dot(Z0_hopf(i,:),[v2(i)-v1(i);0]);
end

%integrate
T1_hopf = cumtrapz(dt,integrand);
T1_hopf = -T1_hopf(end);

%find nu
nu_hopf = T1_hopf/Tmax;

%RHS
FM = [];
for j=1:length(U_s_hopf)
    
    %create df/deps
    dfdeps = [v2(j)-v1(j);0];
    
    %create fundamental matrix solution
    PHI = [U_s_hopf(j,3), U_s_hopf(j,4); U_s_hopf(j,5), U_s_hopf(j,6)];
    
    %create vector slot for FM
    la = F_p(0,[v1(j), n1(j), v2(j), n2(j)]);
    temp = PHI\(dfdeps+nu_hopf*la(1:2));
    
    %FM
    FM = [FM temp];
    
end

%integrate
FM = [cumtrapz(dt,FM(1,:));cumtrapz(dt,FM(2,:))];

%pls god pls
term3=[];
for i=1:length(U_s_hopf)
    
    %define PHI
    PHI = [U_s_hopf(j,3), U_s_hopf(j,4); U_s_hopf(j,5), U_s_hopf(j,6)];
    
    %sol
    term3(:,i) = PHI*FM(:,i);
    
end

%end point
FM = [term3(1,end); term3(2,end)];

%setup RHS
test=M-eye(2);
test=[test; n_hopf(1) n_hopf(2)];
FM = -[FM; 0];

%solve
p1_lul_hopf = test\FM;


%% P section + nu + IC (vdp)

%define variational equation to find floquet multipliers and vectors
F_s = @(t,u) [u(1)-u(1)^3-u(2); 
              u(1);
              (1-3*u(1)^2)*u(3) + (-1)*u(5);
              (1-3*u(1)^2)*u(4) + (-1)*u(6);
              (1)*u(3) + (0)*u(5);
              (1)*u(4) + (0)*u(6)];

%solve variation equation
[~,U_s_vdp] = ode113(F_s,0:dt:Tmax,[x0(3:4) 1 0 0 1],opts);

%Monodromy matrix
M=zeros(2,2);
M(1,1) = U_s_vdp(end,3);
M(1,2) = U_s_vdp(end,4);
M(2,1) = U_s_vdp(end,5);
M(2,2) = U_s_vdp(end,6);
M_vdp = M;

%find floquet coordinates (eigenvectors)
[v,~] = eig(M);

%define new coordinate system, q2 is poincare section
q = v(:,1);

%choose several poincare sections for testing
theta = .5;
% theta = 20;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
q = R*q;
n_vdp = null(q');

%Z0 
Z0_vdp_func = @(t,u) [u(1)-u(1)^3-u(2); 
                 u(1);
                 -(1-3*u(1)^2)*u(3) + -(1)*u(4);
                 -(-1)*u(3) + -(0)*u(4)];
                          
%solve it
la = F_p(0,x0);
IC = [x0(3:4) la(3:4)'];
for i=1:3
    [~,U] = ode113(Z0_vdp_func,Tmax:-dt:0,IC,opts);
    IC = [x0(3:4) U(end,3:4)];
end
Z0_vdp = flip(U);
Z0_vdp = Z0_vdp(:,3:4);
Z0_vdp = Z0_vdp/(dot(Z0_vdp(1,:),la(3:4)'));

%form integrand
integrand = zeros(1,length(0:dt:Tmax));
for i=1:length(0:dt:Tmax)
    integrand(i) = dot(Z0_vdp(i,:),[v1(i)-v2(i);0]);
end

%integrate
T1_vdp = cumtrapz(dt,integrand);
T1_vdp = -T1_vdp(end);

%find nu
nu_vdp = T1_vdp/Tmax;

%RHS
FM = [];
for j=1:length(U_s_vdp)
    
    %create df/deps
    dfdeps = [v1(j)-v2(j);0];
    
    %create fundamental matrix solution
    PHI = [U_s_vdp(j,3), U_s_vdp(j,4); U_s_vdp(j,5), U_s_vdp(j,6)];
    
    %create vector slot for FM
    la = F_p(0,[v1(j), n1(j), v2(j), n2(j)]);
    temp = PHI\(dfdeps+nu_vdp*la(3:4));
    
    %FM
    FM = [FM temp];
    
end

%integrate
FM = [cumtrapz(dt,FM(1,:));cumtrapz(dt,FM(2,:))];

%pls god pls
term3=[];
for i=1:length(U_s_vdp)
    
    %define PHI
    PHI = [U_s_vdp(j,3), U_s_vdp(j,4); U_s_vdp(j,5), U_s_vdp(j,6)];
    
    %sol
    term3(:,i) = PHI*FM(:,i);
    
end

%end point
FM = [term3(1,end); term3(2,end)];

%setup RHS
test=M-eye(2);
test=[test; n_vdp(1) n_vdp(2)];
FM = -[FM; 0];

%solve
p1_lul_vdp = test\FM;


%% iSRC

%iSRC equation
ISRC1 = @(t,u) [(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2); 
                (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
                u(3)-u(3)^3-u(4);
                u(3);
                (.5-3*u(1)^2-u(2)^2)*u(5) + (-2*u(2)*u(1)-omega)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2)) + u(3)-u(1);
                (-2*u(1)*u(2)+omega)*u(5) + (.5-u(1)^2-3*u(2)^2)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 0;
                (0)*u(5) + (0)*u(6) + (1-3*u(3)^2)*u(7) + (-1)*u(8) + nu_vdp*(u(3)-u(3)^3-u(4)) + u(1)-u(3);
                (0)*u(5) + (0)*u(6) + (1)*u(7) + (0)*u(8) + nu_vdp*(u(3))];
                    

%skip a few cycles
tran = 1;

%initial cond
IC1 = [p1_lul_hopf' p1_lul_vdp'];

%solve
for k=1:tran
    
    %solve
    [T_src,U_src] = ode113(ISRC1,0:dt:Tmax,[x0 IC1],opts);
    
    %periodicity shooting thingy
    IC1 = U_src(end,5:8);
    
end


%% nu2

%iSRC derivative
gamma_1_p_x = gradient(U_src(:,5),dt);
gamma_1_p_y = gradient(U_src(:,6),dt);
gamma_1_p_w = gradient(U_src(:,7),dt);
gamma_1_p_z = gradient(U_src(:,8),dt);

%H
H = zeros(4,length(t));
for i=1:length(t)
    H(1,i) = 1/2*kron(U_src(i,5:8),U_src(i,5:8))*[-6*v1(i);-2*n1(i);0;0;-2*n1(i);-2*v1(i);0;0;0;0;0;0;0;0;0;0];
    H(2,i) = 1/2*kron(U_src(i,5:8),U_src(i,5:8))*[-2*n1(i);-2*v1(i);0;0;-2*v1(i);-6*n1(i);0;0;0;0;0;0;0;0;0;0];
    H(3,i) = 1/2*kron(U_src(i,5:8),U_src(i,5:8))*[0;0;0;0;0;0;0;0;0;0;-6*v2(i);0;0;0;0;0];
    H(4,i) = 1/2*kron(U_src(i,5:8),U_src(i,5:8))*[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0];
end

%D
D = zeros(4,length(t));
for i=1:length(t)
    D(:,i) = [-1 0 1 0; 0 0 0 0; 1 0 -1 0; 0 0 0 0]*U_src(i,5:8)';
end

%change in timing
C = zeros(4,length(t));
for i=1:length(t)
   C(:,i) = [nu_hopf*gamma_1_p_x(i);nu_hopf*gamma_1_p_y(i);nu_vdp*gamma_1_p_w(i);nu_vdp*gamma_1_p_z(i)]; 
end

%form lol
integrand_1 = zeros(1,length(t));
integrand_2 = zeros(1,length(t));
RHS = zeros(4,length(t));
for i=1:length(t)
    RHS(:,i) = H(:,i)+D(:,i)+C(:,i);
    integrand_1(i) = dot([Z0_hopf(i,:)';0;0],RHS(:,i));
    integrand_2(i) = dot([0;0;Z0_vdp(i,:)'],RHS(:,i));
end

%integrate
T2_hopf = cumtrapz(t,integrand_1);
T2_hopf = -T2_hopf(end);
T2_vdp = cumtrapz(t,integrand_2);
T2_vdp = -T2_vdp(end);

%find nu
nu2_hopf = T2_hopf/Tmax;
nu2_vdp = T2_vdp/Tmax;


%% IC

%RHS
FM = [];
for j=1:length(U_s_hopf)
        
    %create fundamental matrix solution
    PHI = [U_s_hopf(j,3), U_s_hopf(j,4); U_s_hopf(j,5), U_s_hopf(j,6)];
    
    %create vector slot for FM
    la = F_p(0,[U_s_hopf(j,1), U_s_hopf(j,2), U_s_vdp(j,1), U_s_vdp(j,2)]);
    temp = PHI\(RHS(1:2,j)+nu2_hopf*la(1:2));
    
    %FM
    FM = [FM temp];
    
end

%integrate
FM = [cumtrapz(dt,FM(1,:));cumtrapz(dt,FM(2,:))];

%pls god pls
term3=[];
for i=1:length(U_s_hopf)
    
    %define PHI
    PHI = [U_s_hopf(j,3), U_s_hopf(j,4); U_s_hopf(j,5), U_s_hopf(j,6)];
    
    %sol
    term3(:,i) = PHI*FM(:,i);
    
end

%end point
FM = [term3(1,end); term3(2,end)];

%setup RHS
test=M_hopf-eye(2);
test=[test; n_hopf(1) n_hopf(2)];
FM = -[FM; 0];

%solve
p2_lul_hopf = test\FM;


%RHS
FM = [];
for j=1:length(U_s_vdp)
        
    %create fundamental matrix solution
    PHI = [U_s_vdp(j,3), U_s_vdp(j,4); U_s_vdp(j,5), U_s_vdp(j,6)];
    
    %create vector slot for FM
    la = F_p(0,[U_s_hopf(j,1), U_s_hopf(j,2), U_s_vdp(j,1), U_s_vdp(j,2)]);
    temp = PHI\(RHS(3:4,j)+nu2_vdp*la(3:4));
    
    %FM
    FM = [FM temp];
    
end

%integrate
FM = [cumtrapz(dt,FM(1,:));cumtrapz(dt,FM(2,:))];

%pls god pls
term3=[];
for i=1:length(U_s_vdp)
    
    %define PHI
    PHI = [U_s_vdp(j,3), U_s_vdp(j,4); U_s_vdp(j,5), U_s_vdp(j,6)];
    
    %sol
    term3(:,i) = PHI*FM(:,i);
    
end

%end point
FM = [term3(1,end); term3(2,end)];

%setup RHS
test=M_vdp-eye(2);
test=[test; n_vdp(1) n_vdp(2)];
FM = -[FM; 0];

%solve
p2_lul_vdp = test\FM;


%% iSRC

%iSRC equation
ISRC = @(t,u) [(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2); 
               (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
               u(3)-u(3)^3-u(4);
               u(3);
               
               (.5-3*u(1)^2-u(2)^2)*u(5) + (-2*u(2)*u(1)-omega)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2)) + u(3)-u(1);
                (-2*u(1)*u(2)+omega)*u(5) + (.5-u(1)^2-3*u(2)^2)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 0;
                (0)*u(5) + (0)*u(6) + (1-3*u(3)^2)*u(7) + (-1)*u(8) + nu_vdp*(u(3)-u(3)^3-u(4)) + u(1)-u(3);
                (0)*u(5) + (0)*u(6) + (1)*u(7) + (0)*u(8) + nu_vdp*(u(3));
               
               (.5-3*u(1)^2-u(2)^2)*u(9) + (-2*u(2)*u(1)-omega)*u(10) + (0)*u(11) + (0)*u(12) + nu_hopf*((.5-3*u(1)^2-u(2)^2)*u(5) + (-2*u(2)*u(1)-omega)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2)) + u(3)-u(1)) + nu2_hopf*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2)) + 1/2*(kron([u(5);u(6);u(7);u(8)],[u(5);u(6);u(7);u(8)])'*[-6*u(1);-2*u(2);0;0;-2*u(2);-2*u(1);0;0;0;0;0;0;0;0;0;0])-u(5)+u(7);
               (-2*u(1)*u(2)+omega)*u(9) + (.5-u(1)^2-3*u(2)^2)*u(10) + (0)*u(11) + (0)*u(12) + nu_hopf*((-2*u(1)*u(2)+omega)*u(5) + (.5-u(1)^2-3*u(2)^2)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 0) + nu2_hopf*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 1/2*(kron([u(5);u(6);u(7);u(8)],[u(5);u(6);u(7);u(8)])'*[-2*u(2);-2*u(1);0;0;-2*u(1);-6*u(2);0;0;0;0;0;0;0;0;0;0]) + 0;
               (0)*u(9) + (0)*u(10) + (1-3*u(3)^2)*u(11) + (-1)*u(12) + nu_vdp*((0)*u(5) + (0)*u(6) + (1-3*u(3)^2)*u(7) + (-1)*u(8) + nu_vdp*(u(3)-u(3)^3-u(4)) + u(1)-u(3)) + nu2_vdp*(u(3)-u(3)^3-u(4)) + 1/2*(kron([u(5);u(6);u(7);u(8)],[u(5);u(6);u(7);u(8)])'*[0;0;0;0;0;0;0;0;0;0;-6*u(3);0;0;0;0;0]) + u(5)-u(7);
               (0)*u(9) + (0)*u(10) + (1)*u(11) + (0)*u(12) + nu_vdp*((0)*u(5) + (0)*u(6) + (1)*u(7) + (0)*u(8) + nu_vdp*(u(3))) + nu2_vdp*(u(3)) + 0 + 0];

%skip a few cycles
tran=1;

%initial cond
IC=[p1_lul_hopf' p1_lul_vdp' p2_lul_hopf' p2_lul_vdp'];

%solve
for k=1:tran
    
    %solve
    [T_src2,U_src2] = ode113(ISRC,0:dt:Tmax,[x0 IC],opts);
    
    %periodicity shooting thingy
    IC = U_src2(end,5:12);
    
end


%% nu3

%iSRC derivative
gamma_2_p_x = gradient(U_src2(:,9),dt);
gamma_2_p_y = gradient(U_src2(:,10),dt);
gamma_2_p_w = gradient(U_src2(:,11),dt);
gamma_2_p_z = gradient(U_src2(:,12),dt);

%H
H = zeros(4,length(t));
for i=1:length(t)
    H(1,i) = 1/2*(kron(U_src2(i,5:8),U_src2(i,9:12))+kron(U_src2(i,9:12),U_src2(i,5:8)))*[-6*v1(i);-2*n1(i);0;0;-2*n1(i);-2*v1(i);0;0;0;0;0;0;0;0;0;0] + 1/6*kron(kron(U_src2(i,5:8),U_src2(i,5:8)),U_src2(i,5:8))*[-6, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
    H(2,i) = 1/2*(kron(U_src2(i,5:8),U_src2(i,9:12))+kron(U_src2(i,9:12),U_src2(i,5:8)))*[-2*n1(i);-2*v1(i);0;0;-2*v1(i);-6*n1(i);0;0;0;0;0;0;0;0;0;0] + 1/6*kron(kron(U_src2(i,5:8),U_src2(i,5:8)),U_src2(i,5:8))*[0, -2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
    H(3,i) = 1/2*(kron(U_src2(i,5:8),U_src2(i,9:12))+kron(U_src2(i,9:12),U_src2(i,5:8)))*[0;0;0;0;0;0;0;0;0;0;-6*v2(i);0;0;0;0;0] + 1/6*kron(kron(U_src2(i,5:8),U_src2(i,5:8)),U_src2(i,5:8))*[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
    H(4,i) = 1/2*(kron(U_src2(i,5:8),U_src2(i,9:12))+kron(U_src2(i,9:12),U_src2(i,5:8)))*[0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0] + 0;
end

%D
D = zeros(4,length(t));
for i=1:length(t)
    D(:,i) = [-1 0 1 0; 0 0 0 0; 1 0 -1 0; 0 0 0 0]*U_src2(i,9:12)';
end

%change in timing
C = zeros(4,length(t));
for i=1:length(t)
   C(:,i) = [nu_hopf*gamma_2_p_x(i);nu_hopf*gamma_2_p_y(i);nu_vdp*gamma_2_p_w(i);nu_vdp*gamma_2_p_z(i)] + ... 
            [nu2_hopf*gamma_1_p_x(i);nu2_hopf*gamma_1_p_y(i);nu2_vdp*gamma_1_p_w(i);nu2_vdp*gamma_1_p_z(i)]; 
end

%form lol
integrand_1 = zeros(1,length(t));
integrand_2 = zeros(1,length(t));
RHS = zeros(4,length(t));
for i=1:length(t)
    RHS(:,i) = H(:,i)+D(:,i)+C(:,i);
    integrand_1(i) = dot([Z0_hopf(i,:)';0;0],RHS(:,i));
    integrand_2(i) = dot([0;0;Z0_vdp(i,:)'],RHS(:,i));
end

%integrate
T3_hopf = cumtrapz(t,integrand_1);
T3_hopf = -T3_hopf(end);
T3_vdp = cumtrapz(t,integrand_2);
T3_vdp = -T3_vdp(end);

%find nu
nu3_hopf = T3_hopf/Tmax;
nu3_vdp = T3_vdp/Tmax;


%% IC

%RHS
FM = [];
for j=1:length(U_s_hopf)
        
    %create fundamental matrix solution
    PHI = [U_s_hopf(j,3), U_s_hopf(j,4); U_s_hopf(j,5), U_s_hopf(j,6)];
    
    %create vector slot for FM
    la = F_p(0,[U_s_hopf(j,1), U_s_hopf(j,2), U_s_vdp(j,1), U_s_vdp(j,2)]);
    temp = PHI\(RHS(1:2,j)+nu3_hopf*la(1:2));
    
    %FM
    FM = [FM temp];
    
end

%integrate
FM = [cumtrapz(dt,FM(1,:));cumtrapz(dt,FM(2,:))];

%pls god pls
term3=[];
for i=1:length(U_s_hopf)
    
    %define PHI
    PHI = [U_s_hopf(j,3), U_s_hopf(j,4); U_s_hopf(j,5), U_s_hopf(j,6)];
    
    %sol
    term3(:,i) = PHI*FM(:,i);
    
end

%end point
FM = [term3(1,end); term3(2,end)];

%setup RHS
test=M_hopf-eye(2);
test=[test; n_hopf(1) n_hopf(2)];
FM = -[FM; 0];

%solve
p3_lul_hopf = test\FM;

%RHS
FM = [];
for j=1:length(U_s_vdp)
        
    %create fundamental matrix solution
    PHI = [U_s_vdp(j,3), U_s_vdp(j,4); U_s_vdp(j,5), U_s_vdp(j,6)];
    
    %create vector slot for FM
    la = F_p(0,[U_s_hopf(j,1), U_s_hopf(j,2), U_s_vdp(j,1), U_s_vdp(j,2)]);
    temp = PHI\(RHS(3:4,j)+nu3_vdp*la(3:4));
    
    %FM
    FM = [FM temp];
    
end

%integrate
FM = [cumtrapz(dt,FM(1,:));cumtrapz(dt,FM(2,:))];

%pls god pls
term3=[];
for i=1:length(U_s_vdp)
    
    %define PHI
    PHI = [U_s_vdp(j,3), U_s_vdp(j,4); U_s_vdp(j,5), U_s_vdp(j,6)];
    
    %sol
    term3(:,i) = PHI*FM(:,i);
    
end

%end point
FM = [term3(1,end); term3(2,end)];

%setup RHS
test=M_vdp-eye(2);
test=[test; n_vdp(1) n_vdp(2)];
FM = -[FM; 0];

%solve
p3_lul_vdp = test\FM;


%% iSRC

%iSRC equation
ISRC = @(t,u) [(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2); 
               (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
               u(3)-u(3)^3-u(4);
               u(3);
               
               (.5-3*u(1)^2-u(2)^2)*u(5) + (-2*u(2)*u(1)-omega)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2)) + u(3)-u(1);
               (-2*u(1)*u(2)+omega)*u(5) + (.5-u(1)^2-3*u(2)^2)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 0;
               (0)*u(5) + (0)*u(6) + (1-3*u(3)^2)*u(7) + (-1)*u(8) + nu_vdp*(u(3)-u(3)^3-u(4)) + u(1)-u(3);
               (0)*u(5) + (0)*u(6) + (1)*u(7) + (0)*u(8) + nu_vdp*(u(3));
               
               (.5-3*u(1)^2-u(2)^2)*u(9) + (-2*u(2)*u(1)-omega)*u(10) + (0)*u(11) + (0)*u(12) + nu_hopf*((.5-3*u(1)^2-u(2)^2)*u(5) + (-2*u(2)*u(1)-omega)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2)) + u(3)-u(1)) + nu2_hopf*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2)) + 1/2*(kron([u(5);u(6);u(7);u(8)],[u(5);u(6);u(7);u(8)])'*[-6*u(1);-2*u(2);0;0;-2*u(2);-2*u(1);0;0;0;0;0;0;0;0;0;0]) + -u(5)+u(7);
               (-2*u(1)*u(2)+omega)*u(9) + (.5-u(1)^2-3*u(2)^2)*u(10) + (0)*u(11) + (0)*u(12) + nu_hopf*((-2*u(1)*u(2)+omega)*u(5) + (.5-u(1)^2-3*u(2)^2)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 0) + nu2_hopf*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 1/2*(kron([u(5);u(6);u(7);u(8)],[u(5);u(6);u(7);u(8)])'*[-2*u(2);-2*u(1);0;0;-2*u(1);-6*u(2);0;0;0;0;0;0;0;0;0;0]) + 0;
               (0)*u(9) + (0)*u(10) + (1-3*u(3)^2)*u(11) + (-1)*u(12) + nu_vdp*((0)*u(5) + (0)*u(6) + (1-3*u(3)^2)*u(7) + (-1)*u(8) + nu_vdp*(u(3)-u(3)^3-u(4)) + u(1)-u(3)) + nu2_vdp*(u(3)-u(3)^3-u(4)) + 1/2*(kron([u(5);u(6);u(7);u(8)],[u(5);u(6);u(7);u(8)])'*[0;0;0;0;0;0;0;0;0;0;-6*u(3);0;0;0;0;0]) + u(5)-u(7);
               (0)*u(9) + (0)*u(10) + (1)*u(11) + (0)*u(12) + nu_vdp*((0)*u(5) + (0)*u(6) + (1)*u(7) + (0)*u(8) + nu_vdp*(u(3))) + nu2_vdp*(u(3)) + 0 + 0;
               
               (.5-3*u(1)^2-u(2)^2)*u(13) + (-2*u(2)*u(1)-omega)*u(14) + (0)*u(15) + (0)*u(16) + nu_hopf*((.5-3*u(1)^2-u(2)^2)*u(9) + (-2*u(2)*u(1)-omega)*u(10) + (0)*u(11) + (0)*u(12) + nu_hopf*((.5-3*u(1)^2-u(2)^2)*u(5) + (-2*u(2)*u(1)-omega)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2)) + u(3)-u(1)) + nu2_hopf*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2)) + 1/2*(kron([u(5);u(6);u(7);u(8)],[u(5);u(6);u(7);u(8)])'*[-6*u(1);-2*u(2);0;0;-2*u(2);-2*u(1);0;0;0;0;0;0;0;0;0;0]) + -u(5)+u(7)) + nu2_hopf*((.5-3*u(1)^2-u(2)^2)*u(5) + (-2*u(2)*u(1)-omega)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2)) + u(3)-u(1)) + nu3_hopf*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2)) + 1/2*(kron([u(5);u(6);u(7);u(8)]',[u(9);u(10);u(11);u(12)]')+kron([u(9);u(10);u(11);u(12)]',[u(5);u(6);u(7);u(8)]'))*[-6*u(1);-2*u(2);0;0;-2*u(2);-2*u(1);0;0;0;0;0;0;0;0;0;0] + 1/6*kron(kron([u(5);u(6);u(7);u(8)]',[u(5);u(6);u(7);u(8)]'),[u(5);u(6);u(7);u(8)]')*[-6, 0, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]'-u(9)+u(11);
               (-2*u(1)*u(2)+omega)*u(13) + (.5-u(1)^2-3*u(2)^2)*u(14) + (0)*u(15) + (0)*u(16) + nu_hopf*((-2*u(1)*u(2)+omega)*u(9) + (.5-u(1)^2-3*u(2)^2)*u(10) + (0)*u(11) + (0)*u(12) + nu_hopf*((-2*u(1)*u(2)+omega)*u(5) + (.5-u(1)^2-3*u(2)^2)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 0) + nu2_hopf*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 1/2*(kron([u(5);u(6);u(7);u(8)],[u(5);u(6);u(7);u(8)])'*[-2*u(2);-2*u(1);0;0;-2*u(1);-6*u(2);0;0;0;0;0;0;0;0;0;0]) + 0) + nu2_hopf*((-2*u(1)*u(2)+omega)*u(5) + (.5-u(1)^2-3*u(2)^2)*u(6) + (0)*u(7) + (0)*u(8) + nu_hopf*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 0) + nu3_hopf*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 1/2*(kron([u(5);u(6);u(7);u(8)]',[u(9);u(10);u(11);u(12)]')+kron([u(9);u(10);u(11);u(12)]',[u(5);u(6);u(7);u(8)]'))*[-2*u(2);-2*u(1);0;0;-2*u(1);-6*u(2);0;0;0;0;0;0;0;0;0;0] + 1/6*kron(kron([u(5);u(6);u(7);u(8)]',[u(5);u(6);u(7);u(8)]'),[u(5);u(6);u(7);u(8)]')*[0, -2, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]';
               (0)*u(13) + (0)*u(14) + (1-3*u(3)^2)*u(15) + (-1)*u(16) + nu_vdp*((0)*u(9) + (0)*u(10) + (1-3*u(3)^2)*u(11) + (-1)*u(12) + nu_vdp*((0)*u(5) + (0)*u(6) + (1-3*u(3)^2)*u(7) + (-1)*u(8) + nu_vdp*(u(3)-u(3)^3-u(4)) + u(1)-u(3)) + nu2_vdp*(u(3)-u(3)^3-u(4)) + 1/2*(kron([u(5);u(6);u(7);u(8)],[u(5);u(6);u(7);u(8)])'*[0;0;0;0;0;0;0;0;0;0;-6*u(3);0;0;0;0;0]) + u(5)-u(7)) + nu2_vdp*((0)*u(5) + (0)*u(6) + (1-3*u(3)^2)*u(7) + (-1)*u(8) + nu_vdp*(u(3)-u(3)^3-u(4)) + u(1)-u(3)) + nu3_vdp*(u(3)-u(3)^3-u(4)) + 1/2*(kron([u(5);u(6);u(7);u(8)]',[u(9);u(10);u(11);u(12)]')+kron([u(9);u(10);u(11);u(12)]',[u(5);u(6);u(7);u(8)]'))*[0;0;0;0;0;0;0;0;0;0;-6*u(3);0;0;0;0;0] + 1/6*kron(kron([u(5);u(6);u(7);u(8)]',[u(5);u(6);u(7);u(8)]'),[u(5);u(6);u(7);u(8)]')*[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]' + u(9)-u(11);
               (0)*u(13) + (0)*u(14) + (1)*u(15) + (0)*u(16) + nu_vdp*((0)*u(9) + (0)*u(10) + (1)*u(11) + (0)*u(12) + nu_vdp*((0)*u(5) + (0)*u(6) + (1)*u(7) + (0)*u(8) + nu_vdp*(u(3))) + nu2_vdp*(u(3)) + 0 + 0) + nu2_vdp*((0)*u(5) + (0)*u(6) + (1)*u(7) + (0)*u(8) + nu_vdp*(u(3))) + nu3_vdp*(u(3)) + 0];

%skip a few cycles
tran=3;

%initial cond
IC=[p1_lul_hopf' p1_lul_vdp' p2_lul_hopf' p2_lul_vdp' p3_lul_hopf' p3_lul_vdp'];

%solve
for k=1:tran
    
    %solve
    [T_src3,U_src3] = ode113(ISRC,0:dt:Tmax,[x0 IC],opts);
    
    %periodicity shooting thingy
    IC = U_src3(end,5:16);
    
end


%% second iSRC (beta)

%unperturbed vector field
x0p = [-0.589511603659218   0.420915083744529  -1.089922687522716   0.663740923553011];
F_c = @(t,u)  [(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2) + delta*(u(3)-u(1)); 
              (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
              u(3)-u(3)^3-u(4) + delta*(u(1)-u(3));
              u(3)];

%solve
[~,U]=ode113(F_c,0:dt:Tp1*1,x0p,opts);

%ease of notation
v1c = U(:,1);
n1c = U(:,2);
v2c = U(:,3);
n2c = U(:,4);

% %compute unperturbed period T_p
% [~,loc]=findpeaks(v1c,'MinPeakHeight',0);
% timerz = T(loc);
% T_p1 = mean(diff(timerz))
% 
% [~,loc]=findpeaks(v2c,'MinPeakHeight',0);
% timerz = T(loc);
% T_p2 = mean(diff(timerz))


%% perturbed system

%unperturbed vector field
x0x = [-0.585919903526354   0.425903044085868  -1.089980141107930   0.663938942045873];
F_p1 = @(t,u) T_s/Tp1*[(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2) + delta*(u(3)-u(1)); 
              (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
              u(3)-u(3)^3-u(4) + delta*(u(1)-u(3)) + beta*(u(3)-u(3)^3-u(4));
              u(3) + beta*(u(3))];

%solve
[~,U]=ode113(F_p1,0:dt:Tp1*1,x0x,opts);

% %compute unperturbed period T_p
% [~,loc]=findpeaks(vp1c,'MinPeakHeight',0);
% timerz = T(loc);
% T_p1 = mean(diff(timerz))
% 
% [~,loc]=findpeaks(vp2c,'MinPeakHeight',0);
% timerz = T(loc);
% T_p2 = mean(diff(timerz))


%% P section (coupled system)

%define variational equation to find floquet multipliers and vectors
F_s = @(t,u) [(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2) + delta*(u(3)-u(1)); 
              (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
              u(3)-u(3)^3-u(4) + delta*(u(1)-u(3));
              u(3);
              (.5-3*u(1)^2-u(2)^2-delta)*u(5) + (-2*u(2)*u(1)-omega)*u(9) + (delta)*u(13) + (0)*u(17);
              (.5-3*u(1)^2-u(2)^2-delta)*u(6) + (-2*u(2)*u(1)-omega)*u(10) + (delta)*u(14) + (0)*u(18);
              (.5-3*u(1)^2-u(2)^2-delta)*u(7) + (-2*u(2)*u(1)-omega)*u(11) + (delta)*u(15) + (0)*u(19);
              (.5-3*u(1)^2-u(2)^2-delta)*u(8) + (-2*u(2)*u(1)-omega)*u(12) + (delta)*u(16) + (0)*u(20);
              (-2*u(1)*u(2)+omega)*u(5) + (.5-u(1)^2-3*u(2)^2)*u(9) + (0)*u(13) + (0)*u(17);
              (-2*u(1)*u(2)+omega)*u(6) + (.5-u(1)^2-3*u(2)^2)*u(10) + (0)*u(14) + (0)*u(18);
              (-2*u(1)*u(2)+omega)*u(7) + (.5-u(1)^2-3*u(2)^2)*u(11) + (0)*u(15) + (0)*u(19);
              (-2*u(1)*u(2)+omega)*u(8) + (.5-u(1)^2-3*u(2)^2)*u(12) + (0)*u(16) + (0)*u(20);
              (delta)*u(5) + (0)*u(9) + (1-3*u(3)^2-delta)*u(13) + (-1)*u(17);
              (delta)*u(6) + (0)*u(10) + (1-3*u(3)^2-delta)*u(14) + (-1)*u(18);
              (delta)*u(7) + (0)*u(11) + (1-3*u(3)^2-delta)*u(15) + (-1)*u(19);
              (delta)*u(8) + (0)*u(12) + (1-3*u(3)^2-delta)*u(16) + (-1)*u(20);
              (0)*u(5) + (0)*u(9) + (1)*u(13) + (0)*u(17);
              (0)*u(6) + (0)*u(10) + (1)*u(14) + (0)*u(18);
              (0)*u(7) + (0)*u(11) + (1)*u(15) + (0)*u(19);
              (0)*u(8) + (0)*u(12) + (1)*u(16) + (0)*u(20)];

%solve variation equation
[~,U_s] = ode113(F_s,0:dt:Tp1,[x0p 1 0 0 0 0 1 0 0 0 0 1 0 0 0 0 1],opts);

%Monodromy matrix
M=zeros(4,4);
M(1,1) = U_s(end,5);
M(1,2) = U_s(end,6);
M(1,3) = U_s(end,7);
M(1,4) = U_s(end,8);
M(2,1) = U_s(end,9);
M(2,2) = U_s(end,10);
M(2,3) = U_s(end,11);
M(2,4) = U_s(end,12);
M(3,1) = U_s(end,13);
M(3,2) = U_s(end,14);
M(3,3) = U_s(end,15);
M(3,4) = U_s(end,16);
M(4,1) = U_s(end,17);
M(4,2) = U_s(end,18);
M(4,3) = U_s(end,19);
M(4,4) = U_s(end,20);

%find floquet coordinates (eigenvectors)
[v,lambda] = eig(M);

%define new coordinate system, q2 is poincare section
q1 = v(:,2);
q2=null(q1');

%choose several poincare sections for testing
theta1 = 88.285;
theta2 = 88.285;
R1 = [cosd(theta1) -sind(theta1); sind(theta1) cosd(theta1)];
R2 = [cosd(theta2) -sind(theta2); sind(theta2) cosd(theta2)];
RR = [R1, zeros(2,2); zeros(2,2), R2];
q2 = RR*q2;
ns=null(q2');


%% nu

%Z0 
zzz = @(t,u) [(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2) + delta*(u(3)-u(1)); 
                    (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
                    u(3)-u(3)^3-u(4) + delta*(u(1)-u(3));
                    u(3);
                    -(.5-3*u(1)^2-u(2)^2-delta)*u(5) + -(-2*u(1)*u(2)+omega)*u(6) + -(delta)*u(7) + (0)*u(8);
                    -(-2*u(2)*u(1)-omega)*u(5) + -(.5-u(1)^2-3*u(2)^2)*u(6) + (0)*u(7) + (0)*u(8);
                    -(delta)*u(5) + -(0)*u(6) + -(1-3*u(3)^2-delta)*u(7) + -(1)*u(8);
                    -(0)*u(5) + (0)*u(6) + -(-1)*u(7) + (0)*u(8)];
                          
%solve it
IC = [x0p F_c(0,x0p)'];
for i=1:3
    [~,U] = ode113(zzz,Tp1:-dt:0,IC,opts);
    IC = [x0p U(end,5:8)];
end
Z0_s = flip(U);
Z0_s = Z0_s(:,5:8);
Z0_s = Z0_s/(dot(Z0_s(1,:),F_c(0,x0p)));

%form integrand
integrand = zeros(1,length(0:dt:T_s));
for i=1:length(0:dt:T_s)
    integrand(i) = dot(Z0_s(i,:),[0;0;U_s(i,3)-U_s(i,3)^3-U_s(i,4);U_s(i,3)]);
end

%integrate
T1s = cumtrapz(dt,integrand);
T1s = -T1s(end);

%find nu
nus = T1s/T_s;


%% IC

%RHS
FM = [];
for j=1:length(U_s)
    
    %create df/deps
    dfdeps = [0;0;U_s(j,3)-U_s(j,3)^3-U_s(j,4);U_s(j,3)];
    
    %create fundamental matrix solution
    PHI = [U_s(j,5), U_s(j,6), U_s(j,7), U_s(j,8);
           U_s(j,9), U_s(j,10), U_s(j,11), U_s(j,12);
           U_s(j,13), U_s(j,14), U_s(j,15), U_s(j,16);
           U_s(j,17), U_s(j,18), U_s(j,19), U_s(j,20)];

    %create vector slot for FM
    temp = PHI\(dfdeps+nus*F_c(0,[U_s(j,1) U_s(j,2) U_s(j,3) U_s(j,4)]));
    
    %FM
    FM = [FM temp];
    
end

%integrate
FM = [cumtrapz(dt,FM(1,:));cumtrapz(dt,FM(2,:));cumtrapz(dt,FM(3,:));cumtrapz(dt,FM(4,:))];

%pls god pls
term3=[];
for i=1:length(U_s)
    
    %create fundamental matrix solution
    PHI = [U_s(j,5), U_s(j,6), U_s(j,7), U_s(j,8);
           U_s(j,9), U_s(j,10), U_s(j,11), U_s(j,12);
           U_s(j,13), U_s(j,14), U_s(j,15), U_s(j,16);
           U_s(j,17), U_s(j,18), U_s(j,19), U_s(j,20)];
    
    %sol
    term3(:,i) = PHI*FM(:,i);
    
end

%end point
FM = [term3(1,end); term3(2,end); term3(3,end); term3(4,end)];

%setup RHS
testme=M-eye(4);
testme=[testme; ns(1) ns(2) ns(3) ns(4)];
FM = -[FM; 0];

%solve
p1_lul_s = testme\FM;
p_eps_lul = x0p'+beta*p1_lul_s;


%% iSRC (coupled)

%iSRC equation
ISRC1 = @(t,u) [(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2) + delta*(u(3)-u(1)); 
                (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
                u(3)-u(3)^3-u(4) + delta*(u(1)-u(3));
                u(3);
                (.5-3*u(1)^2-u(2)^2-delta)*u(5) + (-2*u(2)*u(1)-omega)*u(6) + (delta)*u(7) + (0)*u(8) + nus*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2) + delta*(u(3)-u(1))) + 0;
                (-2*u(1)*u(2)+omega)*u(5) + (.5-u(1)^2-3*u(2)^2)*u(6) + (0)*u(7) + (0)*u(8) + nus*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 0;
                (delta)*u(5) + (0)*u(6) + (1-3*u(3)^2-delta)*u(7) + (-1)*u(8) + nus*(u(3)-u(3)^3-u(4) + delta*(u(1)-u(3))) + u(3)-u(3)^3-u(4);
                (0)*u(5) + (0)*u(6) + (1)*u(7) + (0)*u(8) + nus*(u(3)) + u(3)];
                    

%skip a few cycles
tran=1;

%initial cond
IC1=p1_lul_s';

%solve
for k=1:tran
    
    %solve
    [T_src_s,U_src_s] = ode113(ISRC1,0:dt:Tp1,[x0p IC1],opts);
    
    %periodicity shooting thingy
    IC1 = U_src_s(end,5:8);
    
end


%% put it all together

%true coupled system (scaled to Tmax periodicity)
x0x = [-0.585919903526354   0.425903044085868  -1.089980141107930   0.663938942045873];
F_p1 = @(t,u) T_s/Tmax*[(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2) + delta*(u(3)-u(1)); 
              (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
              u(3)-u(3)^3-u(4) + delta*(u(1)-u(3)) + beta*(u(3)-u(3)^3-u(4));
              u(3) + beta*(u(3))];
[~,U]=ode113(F_p1,0:dt:Tmax*1,x0x,opts);
vp1c = U(:,1);
np1c = U(:,2);
vp2c = U(:,3);
np2c = U(:,4);

%true base system (scaled to Tmax periodicity)
x0 = [-.575611 0.410697  -1.11779 0.672746];
F_p = @(t,u) [(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2); 
              (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
              u(3)-u(3)^3-u(4);
              u(3)];
[T,U]=ode113(F_p,0:dt:Tmax*1,x0,opts);
v1 = U(:,1);
n1 = U(:,2);
v2 = U(:,3);
n2 = U(:,4);

%scale 2nd iSRC to be Tmax periodic
ISRC1 = @(t,u) Tp1/Tmax*[(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2) + delta*(u(3)-u(1)); 
                (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
                u(3)-u(3)^3-u(4) + delta*(u(1)-u(3));
                u(3);
                (.5-3*u(1)^2-u(2)^2-delta)*u(5) + (-2*u(2)*u(1)-omega)*u(6) + (delta)*u(7) + (0)*u(8) + nus*((.5-u(1)^2-u(2)^2)*u(1) - omega*u(2) + delta*(u(3)-u(1))) + 0;
                (-2*u(1)*u(2)+omega)*u(5) + (.5-u(1)^2-3*u(2)^2)*u(6) + (0)*u(7) + (0)*u(8) + nus*((.5-u(1)^2-u(2)^2)*u(2) + omega*u(1)) + 0;
                (delta)*u(5) + (0)*u(6) + (1-3*u(3)^2-delta)*u(7) + (-1)*u(8) + nus*(u(3)-u(3)^3-u(4) + delta*(u(1)-u(3))) + u(3)-u(3)^3-u(4);
                (0)*u(5) + (0)*u(6) + (1)*u(7) + (0)*u(8) + nus*(u(3)) + u(3)];
                    

%skip a few cycles
tran=1;

%initial cond
IC1=p1_lul_s';

%solve
for k=1:tran
    
    %solve
    [T_src_s,U_src_s] = ode113(ISRC1,0:dt:Tmax,[x0p IC1],opts);
    
    %periodicity shooting thingy
    IC1 = U_src_s(end,5:8);
    
end

figure(1)
subplot(2,2,1)
hold on
plot(T,vp1c-v1,'k-','LineWidth',10)
plot(T,U_src_s(:,5)*beta*1 + U_src3(:,5)*delta + (U_src3(:,9)*delta^2  + U_src3(:,13)*delta^3)*choice,'color',[.7 .7 .7],'LineWidth',5)
title('x')
set(gca,'fontsize',15)
box on

subplot(2,2,2)
hold on
plot(T,np1c-n1,'k-','LineWidth',10)
plot(T,U_src_s(:,6)*beta*1 + U_src3(:,6)*delta + (U_src3(:,10)*delta^2 + U_src3(:,14)*delta^3)*choice,'color',[.7 .7 .7],'LineWidth',5)
title('y')
set(gca,'fontsize',15)
box on

subplot(2,2,3)
hold on
plot(T,vp2c-v2,'k-','LineWidth',10)
plot(T,U_src_s(:,7)*beta*1 + U_src3(:,7)*delta + (U_src3(:,11)*delta^2 + U_src3(:,15)*delta^3)*choice,'color',[.7 .7 .7],'LineWidth',5)
title('w')
set(gca,'fontsize',15)
box on

subplot(2,2,4)
hold on
plot(T,np2c-n2,'k-','LineWidth',10)
plot(T,U_src_s(:,8)*beta*1 + U_src3(:,8)*delta + (U_src3(:,12)*delta^2 + U_src3(:,16)*delta^3)*choice,'color',[.7 .7 .7],'LineWidth',5)
title('z')
set(gca,'fontsize',15)
box on
