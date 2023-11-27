%% setup

%mr clean
clc
clf

%ODE options
format long
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

%time
N = 2048*8+1;
Tmax = 27.58055;
dt=Tmax/(N-1);
t = 0:dt:Tmax;

T_p = 27.631037974211299;
fac = Tmax/T_p;

%system parameters
tau_m = 10;
Delta = 0.3;
J = 21;
Theta = 4;
tau_d = 5;

%parameters
eps = 0.25;

%phase-amplitude parameters
omega_1 = 1;
kappa_2 = -0.059997791426366;
kappa_1 = -0.407952188760009;


%% load data

%K's
K_00 = load('K00.mat');
K_10 = load('K10.mat');
K_01 = load('K01.mat');
K_20 = load('K20.mat');
K_11 = load('K11.mat');
K_02 = load('K02.mat');
K_30 = load('K30.mat');
K_21 = load('K21.mat');
K_12 = load('K12.mat');
K_03 = load('K03.mat');
K_40 = load('K40.mat');
K_31 = load('K31.mat');
K_22 = load('K22.mat');
K_13 = load('K13.mat');
K_04 = load('K04.mat');
K_50 = load('K50.mat');
K_41 = load('K41.mat');
K_32 = load('K32.mat');
K_23 = load('K23.mat');
K_14 = load('K14.mat');
K_05 = load('K05.mat');
K_60 = load('K60.mat');
K_51 = load('K51.mat');
K_42 = load('K42.mat');
K_33 = load('K33.mat');
K_24 = load('K24.mat');
K_15 = load('K15.mat');
K_06 = load('K06.mat');
K_70 = load('K70.mat');
K_61 = load('K61.mat');
K_52 = load('K52.mat');
K_43 = load('K43.mat');
K_34 = load('K34.mat');
K_25 = load('K25.mat');
K_16 = load('K16.mat');
K_07 = load('K07.mat');
K_80 = load('K80.mat');
K_71 = load('K71.mat');
K_62 = load('K62.mat');
K_53 = load('K53.mat');
K_44 = load('K44.mat');
K_35 = load('K35.mat');
K_26 = load('K26.mat');
K_17 = load('K17.mat');
K_08 = load('K08.mat');

%Z's
Z0 = load('Z0.mat');

%I's
I0_1=load('I0_2.mat');
I0_2=load('I0_1.mat');

%make the data have the correct structure
K_00 = cell2mat(struct2cell(K_00));
K_10 = cell2mat(struct2cell(K_10));
K_01 = cell2mat(struct2cell(K_01));
K_20 = cell2mat(struct2cell(K_20));
K_11 = cell2mat(struct2cell(K_11));
K_02 = cell2mat(struct2cell(K_02));
K_30 = cell2mat(struct2cell(K_30));
K_21 = cell2mat(struct2cell(K_21));
K_12 = cell2mat(struct2cell(K_12));
K_03 = cell2mat(struct2cell(K_03));
K_40 = cell2mat(struct2cell(K_40));
K_31 = cell2mat(struct2cell(K_31));
K_22 = cell2mat(struct2cell(K_22));
K_13 = cell2mat(struct2cell(K_13));
K_04 = cell2mat(struct2cell(K_04));
K_50 = cell2mat(struct2cell(K_50));
K_41 = cell2mat(struct2cell(K_41));
K_32 = cell2mat(struct2cell(K_32));
K_23 = cell2mat(struct2cell(K_23));
K_14 = cell2mat(struct2cell(K_14));
K_05 = cell2mat(struct2cell(K_05));
K_60 = cell2mat(struct2cell(K_60));
K_51 = cell2mat(struct2cell(K_51));
K_42 = cell2mat(struct2cell(K_42));
K_33 = cell2mat(struct2cell(K_33));
K_24 = cell2mat(struct2cell(K_24));
K_15 = cell2mat(struct2cell(K_15));
K_06 = cell2mat(struct2cell(K_06));
K_70 = cell2mat(struct2cell(K_70));
K_61 = cell2mat(struct2cell(K_61));
K_52 = cell2mat(struct2cell(K_52));
K_43 = cell2mat(struct2cell(K_43));
K_34 = cell2mat(struct2cell(K_34));
K_25 = cell2mat(struct2cell(K_25));
K_16 = cell2mat(struct2cell(K_16));
K_07 = cell2mat(struct2cell(K_07));
K_80 = cell2mat(struct2cell(K_80));
K_71 = cell2mat(struct2cell(K_71));
K_62 = cell2mat(struct2cell(K_62));
K_53 = cell2mat(struct2cell(K_53));
K_44 = cell2mat(struct2cell(K_44));
K_35 = cell2mat(struct2cell(K_35));
K_26 = cell2mat(struct2cell(K_26));
K_17 = cell2mat(struct2cell(K_17));
K_08 = cell2mat(struct2cell(K_08));

Z0 = cell2mat(struct2cell(Z0));
I0_1 = cell2mat(struct2cell(I0_1));
I0_2 = cell2mat(struct2cell(I0_2));


%% unperturbed dynamics

%function
F = @(t,u) [1/tau_m*(u(1)^2 - (pi*tau_m*u(2))^2 - J*tau_m*u(3) + Theta); ...
            1/tau_m*(Delta/(pi*tau_m)+2*u(2)*u(1)); ...
            1/tau_d*(-u(3)+u(2))];

%solve for uncoupled system
x0 = [2.28564   0.0688145   0.0239971];
[T,U] = ode113(F,0:dt:Tmax,x0,opts);


%% perturbed dynamics

%function
F_c = @(t,u) 1/fac*[1/tau_m*(u(1)^2 - (pi*tau_m*u(2))^2 - (J+eps)*tau_m*u(3) + Theta); ...
                    1/tau_m*(Delta/(pi*tau_m)+2*u(2)*u(1)); ...
                    1/tau_d*(-u(3)+u(2))];

%solve for system
x02 = [2.25266 0.0681603 0.0239081]';
[Tc,Uc] = ode113(F_c,0:dt:1*Tmax*1,x02,opts);

% %compute perturbed period T_p
% [~,loc]=findpeaks(Uc(:,1),'MinPeakHeight',0);
% timerz = Tc(loc);
% Tp = mean(diff(timerz))

% %visualize perturbed LC
% figure(1)
% hold on
% plot3(U(:,1),U(:,2),U(:,3),'k','LineWidth',4)
% plot3(Uc(:,1),Uc(:,2),Uc(:,3),'r.','MarkerSize',20)
% plot3(x0(1),x0(2),x0(3),'k.','MarkerSize',50)
% plot3(x02(1),x02(2),x02(3),'r.','MarkerSize',50)


%% establish coupling function thing

%welp
G1 = -K_00(3,:);
G = [G1;zeros(1,length(G1));zeros(1,length(G1))];


%% compute change in timing nu

%form integrand
integrand=zeros(1,length(t));
for i=1:length(t)
    integrand(i) = dot(Z0(i,:),G(:,i));
end

%spacing
spacing = T_p/(N-1);

%integrate
T1 = cumtrapz(spacing,integrand);
T1 = -T1(end);

%find nu
nu1 = T1/Tmax;


%% compute IC 1

%poincare section (isochron by default)
q2 = [0;1;1];

%find normal vector to poincare section
n=null(q2');

%define matrix A
A = [0 0 0; 0 1-exp(kappa_1*Tmax) 0; 0 0 1-exp(kappa_2*Tmax); n(1) n(2) n(3)];

%define b
integrand_1=zeros(1,length(t));
integrand_2=zeros(1,length(t));
for i=1:length(t)
    integrand_1(i) = exp(kappa_1*(Tmax-t(i)))*dot(I0_1(i,:),G(:,i));
    integrand_2(i) = exp(kappa_2*(Tmax-t(i)))*dot(I0_2(i,:),G(:,i));
end
entry_1 = cumtrapz(t,integrand_1);
entry_2 = cumtrapz(t,integrand_2);
entry_1 = entry_1(end);
entry_2 = entry_2(end);
b = [0;entry_1;entry_2;0];

%solve
p1_lul_1 = A\b;


%% solve isrc

%overlay iPRC relation for funsies
funsies = zeros(1,length(t));
for i=1:length(t)
    if i>1
        integrand = dot(Z0(1:i,:)',G(:,1:i));
        integrand = cumtrapz(dt,integrand);
        integrand = integrand(end);
    else
        integrand = 0;
    end
    funsies(i) = nu1*omega_1*t(i) + integrand;
end
phase1 = funsies; %%%%% maybe add p1_lul_1(1)

%overlay iIRC relation for funsies
funsies = zeros(1,length(t));
for j=1:length(t)
    if j>1
    integrand = exp(-kappa_1*t(1:j)).*dot(I0_1(1:j,:)',G(:,1:j));
    integrand = cumtrapz(dt,integrand);
    integrand = integrand(end);
    else
        integrand=0;
    end
    funsies(j) = p1_lul_1(2)*exp(kappa_1*t(j)) + exp(kappa_1*t(j))*integrand;
end
amp1 = funsies;

%overlay iIRC relation for funsies
funsies = zeros(1,length(t));
for j=1:length(t)
    if j>1
    integrand = exp(-kappa_2*t(1:j)).*dot(I0_2(1:j,:)',G(:,1:j));
    integrand = cumtrapz(dt,integrand);
    integrand = integrand(end);
    else
        integrand=0;
    end
    funsies(j) = p1_lul_1(3)*exp(kappa_2*t(j)) + exp(kappa_2*t(j))*integrand;
end
amp2 = funsies;


%% compute approximated LC in phase amplitude coordinates

gamma_eps_x_1 = mod(omega_1*t + phase1*eps,Tmax+dt);
gamma_eps_y_1 = amp1*eps;
gamma_eps_z_1 = amp2*eps;

gamma_eps_x_11 = mod(omega_1*t,Tmax+dt);
gamma_eps_y_11 = amp1*eps*0;
gamma_eps_z_11 = amp2*eps*0;


%% use phase-amplitude to plot perturbed trajectory

%solution
sol_1 = zeros(3,length(t));
sol_11 = zeros(3,length(t));
for i=1:length(t)
    sol_1(:,i) = K_function(gamma_eps_x_1(i),gamma_eps_y_1(i),gamma_eps_z_1(i),K_00,K_01,K_02,K_03,K_04,K_05,K_06,K_07,K_08,K_10,K_11,K_12,...
                            K_13,K_14,K_15,K_16,K_17,K_20,K_21,K_22,K_23,K_24,K_25,K_26,...
                            K_30,K_31,K_32,K_33,K_34,K_35,K_40,K_41,K_42,K_43,K_44,K_50,K_51,K_52,K_53,K_60,K_61,K_62,K_70,K_71,K_80,Tmax);
    
    sol_11(:,i) = K_function(gamma_eps_x_11(i),gamma_eps_y_11(i),gamma_eps_z_11(i),K_00,K_01,K_02,K_03,K_04,K_05,K_06,K_07,K_08,K_10,K_11,K_12,...
                            K_13,K_14,K_15,K_16,K_17,K_20,K_21,K_22,K_23,K_24,K_25,K_26,...
                            K_30,K_31,K_32,K_33,K_34,K_35,K_40,K_41,K_42,K_43,K_44,K_50,K_51,K_52,K_53,K_60,K_61,K_62,K_70,K_71,K_80,Tmax);
end
sol(1,:) = sol_1(1,:)-sol_11(1,:);
sol(2,:) = sol_1(2,:)-sol_11(2,:);
sol(3,:) = sol_1(3,:)-sol_11(3,:);
                        
% figure(1)
% plot3(U(:,1)'+sol_1(1,:),U(:,2)'+sol_1(2,:),U(:,3)'+sol_1(3,:),'g.','MarkerSize',30)
% plot3(U(1,1)'+sol_1(1,1),U(1,2)'+sol_1(2,1),U(1,3)'+sol_1(3,1),'b.','MarkerSize',40)
                        
%now plot difference
figure(1)
hold on
plot(T,Uc(:,1)-U(:,1),'k*','LineWidth',12)

figure(2)
hold on
plot(T,Uc(:,2)-U(:,2),'k*','LineWidth',12)

figure(3)
hold on
plot(T,Uc(:,3)-U(:,3),'k*','LineWidth',12)

%plot it
figure(1)
plot(t,sol(1,:),'color','r','marker','*','LineWidth',6)
xlim([0 Tmax])
box on
axis square
set(gca,'fontsize',15)
xlabel('time t')
ylabel('x(t)')
title('x-iSRC')

figure(2)
plot(t,sol(2,:),'color','g','marker','*','LineWidth',6)
xlim([0 Tmax])
box on
axis square
set(gca,'fontsize',15)
xlabel('time t')
ylabel('y(t)')
title('y-iSRC')

figure(3)
plot(t,sol(3,:),'color','c','marker','*','LineWidth',6)
xlim([0 Tmax])
box on
axis square
set(gca,'fontsize',15)
xlabel('time t')
ylabel('z(t)')
title('z-iSRC')

figure(4)
hold on
plot(t,U(:,1),'k-','Linewidth',6)
plot(t,Uc(:,1),'color',[.7 .7 .7],'Linewidth',4)
plot(t,sol_1(1,:),'--','color','r','Linewidth',4)
xlim([0 Tmax])
box on
axis square
set(gca,'fontsize',15)
xlabel('time t')
ylabel('x(t)')
title('x-component')

figure(5)
hold on
plot(t,U(:,2),'k-','Linewidth',6)
plot(t,Uc(:,2),'color',[.7 .7 .7],'Linewidth',4)
plot(t,sol_1(2,:),'--','color','g','Linewidth',4)
xlim([0 Tmax])
box on
axis square
set(gca,'fontsize',15)
xlabel('time t')
ylabel('y(t)')
title('y-component')

figure(6)
hold on
plot(t,U(:,3),'k-','Linewidth',6)
plot(t,Uc(:,3),'color',[.7 .7 .7],'Linewidth',4)
plot(t,sol_1(3,:),'--','color','c','Linewidth',4)
xlim([0 Tmax])
box on
axis square
set(gca,'fontsize',15)
xlabel('time t')
ylabel('z(t)')
title('z-component')


%% compute quantitative comparison

%numerator
numer = cumtrapz(dt,abs(Uc(:,1) - U(:,1) - sol_1(1,:)).^2 + abs(Uc(:,2) - U(:,2) - sol_1(2,:)).^2 + abs(Uc(:,3) - U(:,3) - sol_1(3,:)).^2);
numer = sqrt(numer(end));

%mean value
meaner1 = cumtrapz(dt,Uc(:,1))/T_p;
meaner1 = meaner1(end);
meaner2 = cumtrapz(dt,Uc(:,2))/T_p;
meaner2 = meaner2(end);
meaner3 = cumtrapz(dt,Uc(:,3))/T_p;
meaner3 = meaner3(end);

%denom
denomer = cumtrapz(dt,abs(Uc(:,1) - meaner1).^2 + abs(Uc(:,2) - meaner2).^2 + abs(Uc(:,3) - meaner3).^2);
denomer = denomer(end);

%compute E_delta
E_delta = 1/eps*numer/denomer


%% functionals

%K
function[K] = K_function(theta,psi_1,psi_2,K00,K01,K02,K03,K04,K05,K06,K07,K08,K10,K11,K12,K13,K14,K15,K16,K17,K20,K21,K22,K23,K24,K25,K26,...
                         K30,K31,K32,K33,K34,K35,K40,K41,K42,K43,K44,K50,K51,K52,K53,K60,K61,K62,K70,K71,K80,T)

%turn phase into index
index = mod(theta,T)/(T)*length(K00);
index = round(index);
if index == 0
    index = 1;
end

%compute K at given phase and amplitude
K = K00(:,index) + ...
    psi_1*K10(:,index) + psi_2*K01(:,index) + ...
    psi_2^2*K02(:,index) + psi_1*psi_2*K11(:,index) + psi_1^2*K20(:,index) + ...
    psi_2^3*K03(:,index) + psi_1*psi_2^2*K12(:,index) + psi_1^2*psi_2*K21(:,index) + psi_1^3*K30(:,index) + ...
    psi_2^4*K04(:,index) + psi_1*psi_2^3*K13(:,index) + psi_1^2*psi_2^2*K22(:,index) + psi_1^3*psi_2*K31(:,index) + psi_1^4*K40(:,index) + ...
    psi_2^5*K05(:,index) + psi_1*psi_2^4*K14(:,index) + psi_1^2*psi_2^3*K23(:,index) + psi_1^3*psi_2^2*K32(:,index) + psi_1^4*psi_2*K41(:,index) + psi_1^5*K50(:,index) + ...
    psi_2^6*K06(:,index) + psi_1*psi_2^5*K15(:,index) + psi_1^2*psi_2^4*K24(:,index) + psi_1^3*psi_2^3*K33(:,index) + psi_1^4*psi_2^2*K42(:,index) + ...
            + psi_1^5*psi_2*K51(:,index) + psi_1^6*K60(:,index) + ...
    psi_2^7*K07(:,index) + psi_1*psi_2^6*K16(:,index) + psi_1^2*psi_2^5*K25(:,index) + psi_1^3*psi_2^4*K34(:,index) + psi_1^4*psi_2^3*K43(:,index) + ...
            + psi_1^5*psi_2^2*K52(:,index) + psi_1^6*psi_2*K61(:,index) + psi_1^7*K70(:,index) + ...     
    psi_2^8*K08(:,index) + psi_1*psi_2^7*K17(:,index) + psi_1^2*psi_2^6*K26(:,index) + psi_1^3*psi_2^5*K35(:,index) + psi_1^4*psi_2^4*K44(:,index) + ...
            + psi_1^5*psi_2^3*K53(:,index) + psi_1^6*psi_2^2*K62(:,index) + psi_1^7*psi_2*K71(:,index) + psi_1^8*K80(:,index);  
end

