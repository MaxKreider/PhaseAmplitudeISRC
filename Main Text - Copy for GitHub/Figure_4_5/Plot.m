%% setup

%mr clean
clc
clf

%for the ODEs
format long
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

%time
N = 2048+1;
Tmax = 24.407647685044648;
dt = Tmax/(N-1);
t = 0:Tmax/(N-1):Tmax;

%parameter
mu = 0.13;
tau = 24.2;
B = .1;

%phase-amplitude parameters
omega = 1;
kappa = -0.060882514099964;

%Jacobian
Jac = @(t,u) pi/12*[mu-4*mu*u(1)^2, 1; -(24/tau)^2, B];


%% load data

%K's
K_0_1=load('K0_circadian.mat');
K_1_1=load('K1_circadian.mat');
K_2_1=load('K2_circadian.mat');
K_3_1=load('K3_circadian.mat');
K_4_1=load('K4_circadian.mat');
K_5_1=load('K5_circadian.mat');
K_6_1=load('K6_circadian.mat');
K_7_1=load('K7_circadian.mat');
K_8_1=load('K8_circadian.mat');
K_9_1=load('K9_circadian.mat');
K_10_1=load('K10_circadian.mat');
K_11_1=load('K11_circadian.mat');
K_12_1=load('K12_circadian.mat');
K_13_1=load('K13_circadian.mat');
K_14_1=load('K14_circadian.mat');
K_15_1=load('K15_circadian.mat');

%Z's
z0K_1=load('z0K_circadian.mat');

%I's
I0K_1=load('I0K_circadian.mat');

%put data in correct form
K_0_1 = cell2mat(struct2cell(K_0_1));
K_1_1 = cell2mat(struct2cell(K_1_1));
K_2_1 = cell2mat(struct2cell(K_2_1));
K_3_1 = cell2mat(struct2cell(K_3_1));
K_4_1 = cell2mat(struct2cell(K_4_1));
K_5_1 = cell2mat(struct2cell(K_5_1));
K_6_1 = cell2mat(struct2cell(K_6_1));
K_7_1 = cell2mat(struct2cell(K_7_1));
K_8_1 = cell2mat(struct2cell(K_8_1));
K_9_1 = cell2mat(struct2cell(K_9_1));
K_10_1 = cell2mat(struct2cell(K_10_1));
K_11_1 = cell2mat(struct2cell(K_11_1));
K_12_1 = cell2mat(struct2cell(K_12_1));
K_13_1 = cell2mat(struct2cell(K_13_1));
K_14_1 = cell2mat(struct2cell(K_14_1));
K_15_1 = cell2mat(struct2cell(K_15_1));


z0K_1 = cell2mat(struct2cell(z0K_1));
I0K_1 = cell2mat(struct2cell(I0K_1));


%% establish perturbation function thing

%welp
G = pi/12*[ones(1,length(K_0_1(2,:)));K_0_1(2,:)];


%% compute change in timing nu

%form integrand
integrand_1=zeros(1,length(t));
for i=1:length(t)
    integrand_1(i) = dot(z0K_1(:,i),G(:,i));
end

%spacing
spacing = 0:dt:Tmax;

%integrate
T1_1 = cumtrapz(spacing,integrand_1);
T1_1 = -T1_1(end);

%find nu
nu1 = T1_1/Tmax;


%% compute IC

%poincare section (isochron by default)
q2 = [0;1];

%rotate if you wish (and I do wish)
theta = 0;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
q2 = R*q2;

%find normal vector to poincare section
n=null(q2');

%define matrix A
A = [0 0; 0 1-exp(kappa*Tmax); n(1) n(2)];

%define b
integrand=zeros(1,length(t));
for i=1:length(t)
    integrand(i) = exp(-kappa*t(i))*dot(I0K_1(:,i),G(:,i));
end
entry = cumtrapz(t,integrand);
entry = entry(end);
entry = exp(kappa*Tmax)*entry;
b = [0;entry;0];

%solve
p1_lul_1 = A\b;


%% solve isrc

% figure(4)
% hold on
% xlabel('time t')
% ylabel('SRC phase 1')
% title('Phase SRC 1')
% xlim([0 Tmax])
% box on
% axis square
% set(gca,'FontSize',15)

%overlay iPRC relation for funsies
funsies = zeros(1,length(t));
for i=1:length(t)
    if i>1
        integrand = dot(z0K_1(:,1:i),G(:,1:i));
        integrand = cumtrapz(dt,integrand);
        integrand = integrand(end);
    else
        integrand = 0;
    end
    funsies(i) = nu1*omega*t(i) + integrand;
end
% plot(t,p1_lul_1(1)+funsies,'k-','LineWidth',2)
phase_boi = p1_lul_1(1)+funsies;

% figure(5)
% hold on
% xlabel('time t')
% ylabel('SRC amplitude 1')
% title('Amplitude SRC 1')
% xlim([0 Tmax])
% box on
% axis square
% set(gca,'FontSize',15)

%overlay iIRC relation for funsies
funsies = zeros(1,length(t));
for j=1:length(t)
    if j>1
    integrand = exp(-kappa*t(1:j)).*dot(I0K_1(:,1:j),G(:,1:j));
    integrand = cumtrapz(dt,integrand);
    integrand = integrand(end);
    else
        integrand=0;
    end
    funsies(j) = p1_lul_1(2)*exp(kappa*t(j)) + exp(kappa*t(j))*integrand;
end
% plot(t,funsies,'k-','LineWidth',2)
amp_boi = funsies;


%% start le loop

%eps_vec
evec = -0.04:0.005:0.04;

%vector of amplitudes
true_amp_x = zeros(1,length(evec));
iSRC_amp_x = zeros(1,length(evec));
approx_amp_x = zeros(1,length(evec));

%count
count = 1;

for eps=-0.04:0.005:0.04

    %% compute approximated LC in phase amplitude coordinates
    
    %perturbed LC
    gamma_eps_v_1 = mod(t + phase_boi*eps,Tmax);
    gamma_eps_n_1 = amp_boi*eps;
    
    
    %% use phase-amplitude to plot perturbed trajectory
    
    %solution
    sol_1 = zeros(2,length(t));
    for i=1:length(t)
        sol_1(:,i) = K_function(gamma_eps_v_1(i),gamma_eps_n_1(i),K_0_1,K_1_1,K_2_1,K_3_1,K_4_1,K_5_1,K_6_1,K_7_1,K_8_1,K_9_1,K_10_1,K_11_1,K_12_1,K_13_1,K_14_1,K_15_1);
    end
    
    %system
    F = @(t,u) [pi/12*(u(2)+mu*(u(1)-(4/3)*u(1)^3)+B); -pi/12*((24/tau)^2*u(1)-B*u(2))];
    F_p = @(t,u) [pi/12*(u(2)+mu*(u(1)-(4/3)*u(1)^3)+(B+eps)); -pi/12*((24/tau)^2*u(1)-(B+eps)*u(2))];
    
    %initial conditions
    x0 = [1.321694467420288   0.128378509658959];
    [T_u,U_u] = ode113(F,0:dt:Tmax,x0,opts);
    [T_p,U_p] = ode113(F_p,0:dt:150*Tmax,x0,opts);
    
    %perturbed period
    [~,loc]=findpeaks(U_p(:,1),'MinPeakHeight',0);
    timerz = T_p(loc);
    T_p = mean(diff(timerz));
    
    %chop it
    U_p = U_p(round(end-1*T_p/dt):end,:);
    
%     %visualize
%     figure(1)
%     hold on
%     plot(U_u(:,1),U_u(:,2),'k','LineWidth',3)
%     plot(U_p(:,1),U_p(:,2),'r','LineWidth',3)
%     title('LC')
%     xlabel('x')
%     ylabel('y')
%     set(gca,'fontsize',15)
%     box on
%     axis square
%     hold on
%     
%     %and again
%     figure(1)
%     plot(sol_1(1,:),sol_1(2,:),'c--','LineWidth',3)
    
    
    %% keep track of max
    
    %original system
    [~,by] = max(U_u(:,1));
    max0 = U_u(by,:);
    
    %perturbed system
    [~,my] = max(U_p(:,1));
    maxK = U_p(my,:);
    
    %iSRC
    [~,tie] = max(sol_1(1,:));
    maxQ = [sol_1(1,tie);sol_1(2,tie)];
    
%     %plot it
%     figure(1)
%     plot(maxK(1),maxK(2),'r.','MarkerSize',40)
%     plot(maxQ(1),maxQ(2),'c.','MarkerSize',40)
    
    %time at which max occurs
    tm = by;
    
    
    %% analytic method
    
    %first derivative of gamma_0
    g_0_p = zeros(2,length(t));
    for i=1:length(t)
        g_0_p(:,i) = F(t(i),U_u(i,:));
    end
    
    %second derivative of gamma_0
    g_0_pp = gradient(g_0_p,dt);
    
    %first derivative of gamma_1
    g_1_p = zeros(2,length(t));
    for i=1:length(t)
        g_1_p(:,i) = Jac(t(i),U_u(i,:))*sol_1(:,i) + nu1*F(t(i),U_u(i,:)) + pi/12*[1;U_u(i,2)];
    end
    
    %second derivative of gamma_1
    g_1_pp = gradient(g_1_p,dt);
    
    %solve for delta
    delta = -(g_0_p(1,1)+eps*g_1_p(1,1))/(g_0_pp(1,1)+eps*g_1_pp(1,1));
    
    %change it to an index
    delta_I = round(delta/dt);
    
    %write out new index for convenience
    new_index_max = mod(tm+delta_I,length(sol_1(1,:)));
    if new_index_max == 0
        new_index_max = length(sol_1(1,:))-1;
    end
    
%     %plot it
%     figure(1)
%     plot(sol_1(1,new_index_max),sol_1(2,new_index_max),'g.','MarkerSize',40)
    
    %keep track of amplitude 
    true_amp_x(count) = maxK(1);
    iSRC_amp_x(count) = maxQ(1);
    approx_amp_x(count) = sol_1(1,new_index_max);
    
    %update count
    count = count + 1;
    
end

%visualize results
figure(1)
hold on
plot(evec,true_amp_x,'k.','MarkerSize',30)
plot(evec,true_amp_x,'k-','LineWidth',3)
plot(evec,iSRC_amp_x,'c.','MarkerSize',30)
plot(evec,iSRC_amp_x,'c-','LineWidth',3)
plot(evec,approx_amp_x,'color',[255, 182, 193 256]/256,'Marker','.','MarkerSize',30)
plot(evec,approx_amp_x,'color',[255, 182, 193 256]/256,'LineWidth',3)
title('Circadian Temperature Oscillation')
xlabel('\epsilon')
ylabel('max(x)')
box on
axis square
set(gca,'fontsize',15)


%% functionals

%K
function[K] = K_function(theta,psi,K0,K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15)

%turn phase into index
index = mod(theta,24.407647685044648)/(24.407647685044648)*length(K1);
index = round(index);
if index == 0
    index = 1;
end

%compute K at given phase and amplitude
K = K0(:,index) + psi*K1(:,index) + psi^2*K2(:,index) + psi^3*K3(:,index) + psi^4*K4(:,index) + ...
    psi^5*K5(:,index) + psi^6*K6(:,index) + psi^7*K7(:,index) + psi^8*K8(:,index) + ...
    psi^9*K9(:,index) + psi^10*K10(:,index) + psi^11*K11(:,index) + psi^12*K12(:,index) + ...
    psi^13*K13(:,index) + psi^14*K14(:,index) + psi^15*K15(:,index);
end

