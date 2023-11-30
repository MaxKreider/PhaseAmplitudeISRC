%% setup

%mr clean
clc
clf

%ODE options
format long
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

%time
N = 2048*8+1;
Tmax = 6.663550284765624;
dt=Tmax/(N-1);
t = 0:dt:Tmax;

T_p = 6.710088675031566;
fac = Tmax/T_p;

%parameters
eps = .1;
delta = -.03;

%phase-amplitude parameters
omega_1 = 1;
omega_2 = 1;
kappa_1 = -1.059378788499514;
kappa_2 = -1.055266912086093;


%% load data

%K's
K_0_1=load('K0_vdp_1.mat');
K_1_1=load('K1_vdp_1.mat');
K_2_1=load('K2_vdp_1.mat');
K_3_1=load('K3_vdp_1.mat');
K_4_1=load('K4_vdp_1.mat');
K_5_1=load('K5_vdp_1.mat');
K_6_1=load('K6_vdp_1.mat');
K_7_1=load('K7_vdp_1.mat');
K_8_1=load('K8_vdp_1.mat');
K_9_1=load('K9_vdp_1.mat');
K_10_1=load('K10_vdp_1.mat');
K_11_1=load('K11_vdp_1.mat');
K_12_1=load('K12_vdp_1.mat');
K_13_1=load('K13_vdp_1.mat');
K_14_1=load('K14_vdp_1.mat');
K_15_1=load('K15_vdp_1.mat');

K_0_2=load('K0_vdp_2.mat');
K_1_2=load('K1_vdp_2.mat');
K_2_2=load('K2_vdp_2.mat');
K_3_2=load('K3_vdp_2.mat');
K_4_2=load('K4_vdp_2.mat');
K_5_2=load('K5_vdp_2.mat');
K_6_2=load('K6_vdp_2.mat');
K_7_2=load('K7_vdp_2.mat');
K_8_2=load('K8_vdp_2.mat');
K_9_2=load('K9_vdp_2.mat');
K_10_2=load('K10_vdp_2.mat');
K_11_2=load('K11_vdp_2.mat');
K_12_2=load('K12_vdp_2.mat');
K_13_2=load('K13_vdp_2.mat');
K_14_2=load('K14_vdp_2.mat');
K_15_2=load('K15_vdp_2.mat');

%Z's
z0K_1=load('z0K_vdp_1.mat');
z0K_2=load('z0K_vdp_2.mat');

%I's
I0K_1=load('I0K_vdp_1.mat');
I0K_2=load('I0K_vdp_2.mat');

%make the data have the correct structure
K_0_1 = cell2mat(struct2cell(K_0_1));
K_0_2 = cell2mat(struct2cell(K_0_2));
K_1_1 = cell2mat(struct2cell(K_1_1));
K_1_2 = cell2mat(struct2cell(K_1_2));
K_2_1 = cell2mat(struct2cell(K_2_1));
K_2_2 = cell2mat(struct2cell(K_2_2));
K_3_1 = cell2mat(struct2cell(K_3_1));
K_3_2 = cell2mat(struct2cell(K_3_2));
K_4_1 = cell2mat(struct2cell(K_4_1));
K_4_2 = cell2mat(struct2cell(K_4_2));
K_5_1 = cell2mat(struct2cell(K_5_1));
K_5_2 = cell2mat(struct2cell(K_5_2));
K_6_1 = cell2mat(struct2cell(K_6_1));
K_6_2 = cell2mat(struct2cell(K_6_2));
K_7_1 = cell2mat(struct2cell(K_7_1));
K_7_2 = cell2mat(struct2cell(K_7_2));
K_8_1 = cell2mat(struct2cell(K_8_1));
K_8_2 = cell2mat(struct2cell(K_8_2));
K_9_1 = cell2mat(struct2cell(K_9_1));
K_9_2 = cell2mat(struct2cell(K_9_2));
K_10_1 = cell2mat(struct2cell(K_10_1));
K_10_2 = cell2mat(struct2cell(K_10_2));
K_11_1 = cell2mat(struct2cell(K_11_1));
K_11_2 = cell2mat(struct2cell(K_11_2));
K_12_1 = cell2mat(struct2cell(K_12_1));
K_12_2 = cell2mat(struct2cell(K_12_2));
K_13_1 = cell2mat(struct2cell(K_13_1));
K_13_2 = cell2mat(struct2cell(K_13_2));
K_14_1 = cell2mat(struct2cell(K_14_1));
K_14_2 = cell2mat(struct2cell(K_14_2));
K_15_1 = cell2mat(struct2cell(K_15_1));
K_15_2 = cell2mat(struct2cell(K_15_2));

z0K_1 = cell2mat(struct2cell(z0K_1));
z0K_2 = cell2mat(struct2cell(z0K_2));
I0K_1 = cell2mat(struct2cell(I0K_1));
I0K_2 = cell2mat(struct2cell(I0K_2));


%% x-y dynamics

%function
F = @(t,u) [u(1)-u(1)^3-u(2); 
            u(1);
            1/fac*(u(3)-u(3)^3-u(4)+eps);
            1/fac*(u(3)+eps)];

%solve for uncoupled system
x0 = [-.146512 -1.2451 -0.785781525641688 1.301404835817183];
[T,U] = ode113(F,0:dt:Tmax*1,x0,opts);


%% coupled dynamics

%function
F_c = @(t,u) [u(1)-u(1)^3-u(2)+delta*(u(3)-u(1)); 
              u(1);
              1/fac*(u(3)-u(3)^3-u(4)+eps);
              1/fac*(u(3)+eps)];

%solve for system
x02 = [-0.13447  -1.28322  -0.785781525641688   1.301404835817183];
[Tc,Uc] = ode113(F_c,0:dt:Tmax,x02,opts);

%compute unperturbed period T_p
[~,loc]=findpeaks(Uc(:,1),'MinPeakHeight',0);
timerz = Tc(loc);
Tp = mean(diff(timerz));

%visualize unperturbed LC
figure(1)
hold on
plot(Tc,U(:,1),'k*','LineWidth',4)
plot(Tc,U(:,2),'m*','LineWidth',4)

% figure(2)
% hold on
% plot(Tc,Uc(:,3),'k*','LineWidth',4)
% plot(Tc,Uc(:,4),'m*','LineWidth',4)

%visualize perturbed LC
figure(1)
hold on
plot(T,Uc(:,1),'color',[.7 .7 .7],'LineWidth',4)
plot(T,Uc(:,2),'color',[255, 182, 193 256]/256,'LineWidth',4)
xlabel('time t')
ylabel('x(t), y(t)')
set(gca,'fontsize',15)
box on
axis square
xlim([0 Tmax])
title('Oscillator 1')

% figure(2)
% hold on
% plot(T,U(:,3),'color',[.7 .7 .7],'LineWidth',4)
% plot(T,U(:,4),'color',[255, 182, 193 256]/256,'LineWidth',4)
% xlabel('time t')
% ylabel('w(t), z(t)')
% set(gca,'fontsize',15)
% box on
% axis square
% xlim([0 Tmax])
% title('Oscillator 2')

% figure(4)
% hold on
% plot(T,U(:,3)-Uc(:,3),'color','k','LineWidth',4)
% plot(T,U(:,4)-Uc(:,4),'color','k','LineWidth',4)


%% establish coupling function thing

%welp
G1 = (K_0_2(1,:)-K_0_1(1,:));
G = [G1;zeros(1,length(G1))];


%% compute change in timing nu

%form integrand
integrand_1=zeros(1,length(t));
integrand_2=zeros(1,length(t));
for i=1:length(t)
    integrand_1(i) = dot(z0K_1(:,i),G(:,i));
    integrand_2(i) = dot(-z0K_2(:,i),[0;0]);
end

%spacing
spacing = T_p/(N-1);

%integrate
T1_1 = cumtrapz(spacing,integrand_1);
T1_1 = -T1_1(end);
T1_2 = cumtrapz(spacing,integrand_2);
T1_2 = -T1_2(end);

%find nu
nu1 = T1_1/T_p;
nu2 = T1_2/T_p;

%check to make sure it works
T_p_1 = T_p+delta*T1_1;
T_p_2 = T_p+delta*T1_2;


%% compute IC 1

%poincare section (isochron by default)
q2 = [0;1];

%rotate if you wish
theta = 0;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
q2 = R*q2;

%find normal vector to poincare section
n=null(q2');

%define matrix A
A = [0 0; 0 1-exp(kappa_1*Tmax); n(1) n(2)];

%define b
integrand=zeros(1,length(t));
for i=1:length(t)
    integrand(i) = exp(-kappa_1*t(i))*dot(I0K_1(:,i),G(:,i));
end
entry = cumtrapz(t,integrand);
entry = entry(end);
entry = exp(kappa_1*Tmax)*entry;
b = [0;entry;0];

%solve
p1_lul_1 = A\b;


%% compute IC 2

%poincare section (isochron by default)
q2 = [0;1];

%rotate if you wish
theta = 0;
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
q2 = R*q2;

%find normal vector to poincare section
n=null(q2');

%define matrix A
A = [0 0; 0 1-exp(kappa_2*Tmax); n(1) n(2)];

%define b
integrand=zeros(1,length(t));
for i=1:length(t)
    integrand(i) = exp(-kappa_2*t(i))*dot(I0K_2(:,i),[0;0]);
end
entry = cumtrapz(t,integrand);
entry = entry(end);
entry = exp(kappa_2*Tmax)*entry;
b = [0;entry;0];

%solve
p1_lul_2 = A\b;


%% solve isrc

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
    funsies(i) = nu1*omega_1*t(i) + integrand;
end
%plot(t,p1_lul_1(1)+funsies,'m--','LineWidth',2)
phase1 = funsies;

%overlay iIRC relation for funsies
funsies = zeros(1,length(t));
for j=1:length(t)
    if j>1
    integrand = exp(-kappa_1*t(1:j)).*dot(I0K_1(:,1:j),G(:,1:j));
    integrand = cumtrapz(dt,integrand);
    integrand = integrand(end);
    else
        integrand=0;
    end
    funsies(j) = p1_lul_1(2)*exp(kappa_1*t(j)) + exp(kappa_1*t(j))*integrand;
end
%plot(t,funsies,'m--','LineWidth',2)
amp1 = funsies;

%overlay iPRC relation for funsies
funsies = zeros(1,length(t));
for i=1:length(t)
    if i>1
        integrand = dot(z0K_2(:,1:i),zeros(2,i));
        integrand = cumtrapz(dt,integrand);
        integrand = integrand(end);
    else
        integrand = 0;
    end
    funsies(i) = nu2*omega_2*t(i) + integrand;
end
%plot(t,p1_lul_2(1)+funsies,'m--','LineWidth',2)
phase2 = funsies;

%overlay iIRC relation for funsies
funsies = zeros(1,length(t));
for j=1:length(t)
    if j>1
    integrand = exp(-kappa_2*t(1:j)).*dot(I0K_2(:,1:j),zeros(2,j));
    integrand = cumtrapz(dt,integrand);
    integrand = integrand(end);
    else
        integrand=0;
    end
    funsies(j) = p1_lul_2(2)*exp(kappa_2*t(j)) + exp(kappa_2*t(j))*integrand;
end
%plot(t,funsies,'m--','LineWidth',2)
amp2 = funsies;


%% compute approximated LC in phase amplitude coordinates

gamma_eps_v_1 = mod(omega_1*t + phase1*delta,Tmax);
gamma_eps_n_1 = amp1*delta;
gamma_eps_v_2 = mod(omega_2*t + phase2*delta,Tmax);
gamma_eps_n_2 = amp2*delta;

gamma_eps_v_11 = mod(omega_1*t,Tmax);
gamma_eps_n_11 = amp1*delta*0;
gamma_eps_v_22 = mod(omega_2*t,Tmax);
gamma_eps_n_22 = amp2*delta*0;


%% use phase-amplitude to plot perturbed trajectory

%solution
sol_1 = zeros(2,length(t));
sol_2 = zeros(2,length(t));
sol_11 = zeros(2,length(t));
sol_22 = zeros(2,length(t));
for i=1:length(t)
    sol_1(:,i) = K_function(1*gamma_eps_v_1(i),1*gamma_eps_n_1(i),K_0_1,K_1_1,K_2_1,K_3_1,K_4_1,K_5_1,K_6_1,K_7_1,K_8_1,K_9_1,K_10_1,K_11_1,K_12_1,K_13_1,K_14_1,K_15_1,T_p);
    sol_2(:,i) = K_function(1*gamma_eps_v_2(i),1*gamma_eps_n_2(i),K_0_2,K_1_2,K_2_2,K_3_2,K_4_2,K_5_2,K_6_2,K_7_2,K_8_2,K_9_2,K_10_2,K_11_2,K_12_2,K_13_2,K_14_2,K_15_2,T_p);
    
    sol_11(:,i) = K_function(1*gamma_eps_v_11(i),1*gamma_eps_n_11(i),K_0_1,K_1_1,K_2_1,K_3_1,K_4_1,K_5_1,K_6_1,K_7_1,K_8_1,K_9_1,K_10_1,K_11_1,K_12_1,K_13_1,K_14_1,K_15_1,T_p);
    sol_22(:,i) = K_function(1*gamma_eps_v_22(i),1*gamma_eps_n_22(i),K_0_2,K_1_2,K_2_2,K_3_2,K_4_2,K_5_2,K_6_2,K_7_2,K_8_2,K_9_2,K_10_2,K_11_2,K_12_2,K_13_2,K_14_2,K_15_2,T_p);
   
end

sol_1(1,:) = sol_1(1,:)-sol_11(1,:);
sol_1(2,:) = sol_1(2,:)-sol_11(2,:);
sol_2(1,:) = sol_2(1,:)-sol_22(1,:);
sol_2(2,:) = sol_2(2,:)-sol_22(2,:);

%now plot difference
figure(3)
hold on
plot(T,Uc(:,1)-U(:,1),'k*','LineWidth',12)
plot(T,Uc(:,2)-U(:,2),'m*','LineWidth',12)

% figure(4)
% hold on
% plot(T,U(:,3)-Uc(:,3),'k*','LineWidth',12)
% plot(T,U(:,4)-Uc(:,4),'m*','LineWidth',12)

%plot it
figure(3)
plot(t,sol_1(1,:),'color',[.7 .7 .7],'marker','*','LineWidth',6)
plot(t,sol_1(2,:),'color',[255, 182, 193 256]/256,'marker','*','LineWidth',6)
xlim([0 Tmax])
box on
axis square
set(gca,'fontsize',15)
xlabel('time t')
ylabel('x(t), y(t)')
title('iSRC')

% figure(4)
% plot(t,sol_2(1,:),'color',[.7 .7 .7],'marker','*','LineWidth',6)
% plot(t,sol_2(2,:),'color',[255, 182, 193 256]/256,'marker','*','LineWidth',6)
% xlim([0 Tmax])
% box on
% axis square


%% compute quantitative comparison

%numerator
numer = cumtrapz(dt,abs(Uc(:,1) - U(:,1) - sol_1(1,:)').^2 + abs(Uc(:,2) - U(:,2) - sol_1(2,:)').^2);
numer = sqrt(numer(end));

%mean value
meaner1 = cumtrapz(dt,Uc(:,1))/T_p;
meaner1 = meaner1(end);
meaner2 = cumtrapz(dt,Uc(:,2))/T_p;
meaner2 = meaner2(end);

%denom
denomer = cumtrapz(dt,abs(Uc(:,1) - meaner1).^2 + abs(Uc(:,2) - meaner2).^2);
denomer = sqrt(denomer(end));

%compute E_delta
E_delta = 1/delta*numer/denomer


%% functionals

%K
function[K] = K_function(theta,psi,K0,K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,T)

%turn phase into index
index = mod(theta,T)/(T)*length(K1);
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

