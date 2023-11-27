%% setup

%keep it clean
clc
clf

%for the ODEs
format long
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

%time
N = 2048+1;
Tmax = 6.663550284765624;
dt=Tmax/(N-1);
t = 0:Tmax/(N-1):Tmax;

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

%put data in correct form
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
F = @(t,u) [fac*(u(1)-u(1)^3-u(2)); 
            fac*u(1);
            (u(3)-u(3)^3-u(4)+eps);
            (u(3)+eps)];

%solve for uncoupled system
x0 = [-.146512 -1.2451   -0.785781525641688 1.301404835817183];
[T,U] = ode113(F,0:dt:T_p*1,x0,opts);


%% coupled dynamics

%function
F_c = @(t,u) [fac*(u(1)-u(1)^3-u(2))+delta*(u(3)-u(1)); 
              fac*u(1);
              (u(3)-u(3)^3-u(4)+eps);
              (u(3)+eps)];

%solve for system
x02 = [-0.134496906250048  -1.283462958304520  -0.785781525641688   1.301404835817183];
[Tc,Uc] = ode113(F_c,0:dt:1*T_p,x02,opts);

%compute unperturbed period T_p
[~,loc]=findpeaks(Uc(:,1),'MinPeakHeight',0);
timerz = Tc(loc);
Tp = mean(diff(timerz));


%% establish coupling function thing

for j=1:length(K_0_2)

%welp
G1 = (circshift(K_0_2(1,:),j)-K_0_1(1,:));
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
nu1(j) = T1_1/T_p;
nu2(j) = T1_2/T_p;

end

figure(1)
plot(t,nu2-nu1,'LineWidth',2)
xlabel('phase offset')
ylabel('frequency matching curve')
xlim([0 Tmax])
box on
axis square
yline(0)
set(gca,'fontsize',12)
title('Frequency Matching Curve: Entrainment')


