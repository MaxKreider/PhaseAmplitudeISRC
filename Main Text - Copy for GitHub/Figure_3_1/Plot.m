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

%rotate if you wish (I don't wish)
theta = 0;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
q2 = R*q2;

%find normal vector to poincare section
n = null(q2');

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


%% I copy/pasted this title :)

%perturbation
eps = 0.02;

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

%initial conditions
x0 = [1.321694467420288   0.128378509658959];
x01 = [1.32427 0.148761];
[T_u,U_u] = ode113(F,0:dt:Tmax,x0,opts);

%visualize
figure(1)
subplot(1,2,2)
hold on
plot(U_u(:,1),U_u(:,2),'k','LineWidth',3)
plot(sol_1(1,:),sol_1(2,:),'r','LineWidth',3)
title('LC')
xlabel('x')
ylabel('y')
set(gca,'fontsize',15)
box on
axis square
hold on


%% plot isochrons

%now, fix a value of theta (let's set theta = 0 first for simplicity)
theta_ind = 1;

%let's define a vector of sigma's to compute the isochron
sigma = -.5:.001:.5;

%store the Kn's nicely
Kx = [K_0_1(1,:); K_1_1(1,:); K_2_1(1,:); K_3_1(1,:); K_4_1(1,:); K_5_1(1,:); K_6_1(1,:); K_7_1(1,:); ...
    K_8_1(1,:); K_9_1(1,:); K_10_1(1,:); K_11_1(1,:); K_12_1(1,:); K_13_1(1,:); K_14_1(1,:); K_15_1(1,:)];
Ky = [K_0_1(2,:); K_1_1(2,:); K_2_1(2,:); K_3_1(2,:); K_4_1(2,:); K_5_1(2,:); K_6_1(2,:); K_7_1(2,:); ...
    K_8_1(2,:); K_9_1(2,:); K_10_1(2,:); K_11_1(2,:); K_12_1(2,:); K_13_1(2,:); K_14_1(2,:); K_15_1(2,:)];

%initialize an isochron vector
iso_0_x = zeros(1,length(sigma));
iso_0_y = zeros(1,length(sigma));

%solve
for k=0:128:N
    for i=1:length(sigma)
        for j=0:15
            iso_0_x(i) = iso_0_x(i) + Kx(j+1,theta_ind+k)*sigma(i)^j;
            iso_0_y(i) = iso_0_y(i) + Ky(j+1,theta_ind+k)*sigma(i)^j;
        end
    end

    %visualize
    plot(iso_0_x,iso_0_y,'color',[.7 .7 .7],'LineWidth',3)

    %reset
    iso_0_x = zeros(1,length(sigma));
    iso_0_y = zeros(1,length(sigma));
end

%poincare section
for i=1:length(sigma)
    for j=0:15
        iso_0_x(i) = iso_0_x(i) + Kx(j+1,theta_ind+0)*sigma(i)^j;
        iso_0_y(i) = iso_0_y(i) + Ky(j+1,theta_ind+0)*sigma(i)^j;
    end
end

%visualize
plot(iso_0_x,iso_0_y,'color','c','LineWidth',3)


%% amplitude curves

%let's define a sigma value
sigma = -.5:.25:.5;

%define a nice colormap for the phase
mymapboi = parula(length(sigma));

%initialize an isochron vector
isostable_x = zeros(1,length(Kx(1,:)));
isostable_y = zeros(1,length(Kx(1,:)));

%solve
for i=1:length(sigma)
    for j=0:15
        isostable_x = isostable_x + Kx(j+1,:)*sigma(i)^j;
        isostable_y = isostable_y + Ky(j+1,:)*sigma(i)^j;
    end
    
    %visualize
    plot(isostable_x,isostable_y,'color',[255, 182, 193 100]/256,'LineWidth',3)
    
    %reset
    isostable_x = zeros(1,length(Kx(1,:)));
    isostable_y = zeros(1,length(Kx(1,:)));
    
end

%visualize
figure(1)
hold on
plot(U_u(:,1),U_u(:,2),'k','LineWidth',3)
plot(sol_1(1,:),sol_1(2,:),'color',[0.9100 0.4100 0.1700],'LineWidth',3)
title('Cartesian')
xlabel('x')
ylabel('y')
set(gca,'fontsize',15)
box on
axis square
hold on
% xlim([x0(1)-.01 x0(1)+.01])
% ylim([x0(2)-.06 x0(2)+.06])

%% plot isochrons

%now, fix a value of theta (let's set theta = 0 first for simplicity)
theta_ind = 1;

%let's define a vector of sigma's to compute the isochron
sigma = -.5:.001:.5;

%store the Kn's nicely
Kx = [K_0_1(1,:); K_1_1(1,:); K_2_1(1,:); K_3_1(1,:); K_4_1(1,:); K_5_1(1,:); K_6_1(1,:); K_7_1(1,:); ...
    K_8_1(1,:); K_9_1(1,:); K_10_1(1,:); K_11_1(1,:); K_12_1(1,:); K_13_1(1,:); K_14_1(1,:); K_15_1(1,:)];
Ky = [K_0_1(2,:); K_1_1(2,:); K_2_1(2,:); K_3_1(2,:); K_4_1(2,:); K_5_1(2,:); K_6_1(2,:); K_7_1(2,:); ...
    K_8_1(2,:); K_9_1(2,:); K_10_1(2,:); K_11_1(2,:); K_12_1(2,:); K_13_1(2,:); K_14_1(2,:); K_15_1(2,:)];

%initialize an isochron vector
iso_0_x = zeros(1,length(sigma));
iso_0_y = zeros(1,length(sigma));

%solve
for k=0:128:N
    for i=1:length(sigma)
        for j=0:15
            iso_0_x(i) = iso_0_x(i) + Kx(j+1,theta_ind+k)*sigma(i)^j;
            iso_0_y(i) = iso_0_y(i) + Ky(j+1,theta_ind+k)*sigma(i)^j;
        end
    end

    %visualize
    plot(iso_0_x,iso_0_y,'color',[.7 .7 .7],'LineWidth',3)

    %reset
    iso_0_x = zeros(1,length(sigma));
    iso_0_y = zeros(1,length(sigma));
end

%poincare section
for i=1:length(sigma)
    for j=0:15
        iso_0_x(i) = iso_0_x(i) + Kx(j+1,theta_ind+0)*sigma(i)^j;
        iso_0_y(i) = iso_0_y(i) + Ky(j+1,theta_ind+0)*sigma(i)^j;
    end
end

%visualize
plot(iso_0_x,iso_0_y,'color','c','LineWidth',3)

%visualize
figure(1)
hold on
plot(U_u(:,1),U_u(:,2),'k','LineWidth',3)
plot(sol_1(1,:),sol_1(2,:),'color',[0.9100 0.4100 0.1700],'LineWidth',3)
title('Cartesian')
xlabel('x')
ylabel('y')
set(gca,'fontsize',15)
box on
axis square
hold on

%visualize
plot(iso_0_x,iso_0_y,'color','c','LineWidth',3)


%% phase amplitude figure

%visualize
% figure(2)
subplot(1,2,1)
hold on
plot(t,zeros(1,length(t)),'k','LineWidth',3)
plot(t+eps*gamma_eps_v_1,eps*circshift(gamma_eps_n_1,1028),'color',[0.9100 0.4100 0.1700],'LineWidth',3)
title('Phase-Amplitude')
xlabel('x')
ylabel('y')
set(gca,'fontsize',15)
box on
axis square
hold on
xticklabels({'-T/2','-T/4','0','T/4','T/2'})
xticks([1 24.407647685044648/4+1 24.407647685044648/2 3*24.407647685044648/4-1 24.407647685044648-1])
xlim([0 Tmax])


%% isochrons
for index=1:128:N
    xval = t(index);
    plot([xval xval],[-10^(-3) 1*3.7*10^(-3)],'color',[.7 .7 .7],'LineWidth',2)
    plot([xval xval],[-8.5*10^(-3) 8.5*10^(-3)],'color',[.7 .7 .7],'LineWidth',2)
end
ylim([-10^(-3) 3.7*10^(-3)])
ylim([-8.5*10^(-3) 8.5*10^(-3)])

%% amplitude
index = -8*10^(-3):8*10^(-3)/2:8*10^(-3);
for i=1:length(index)
    plot([0 Tmax],[index(i) index(i)],'color',[255, 182, 193 100]/256,'LineWidth',2)
end

%visualize
hold on
plot(t,zeros(1,length(t)),'k','LineWidth',3)
plot(t+eps*gamma_eps_v_1,eps*circshift(gamma_eps_n_1,1028),'color',[0.9100 0.4100 0.1700],'LineWidth',3)
title('Phase-Amplitude')
xlabel('\theta')
ylabel('\sigma')
set(gca,'fontsize',15)
box on
axis square
hold on
xticklabels({'-T/2','-T/4','0','T/4','T/2'})
xticks([1 24.407647685044648/4+1 24.407647685044648/2 3*24.407647685044648/4-1 24.407647685044648-1])
xlim([0 Tmax])

for index=1029
    xval = t(index);
    plot([xval xval],[-10^(-3) 1*3.7*10^(-3)],'color','c','LineWidth',2)
    plot([xval xval],[-8.5*10^(-3) 8.5*10^(-3)],'color','c','LineWidth',2)
end


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

