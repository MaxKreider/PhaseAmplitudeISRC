%% setup

%keep it clean
clc
clf

%for the ODEs
format long
opts = odeset('RelTol',1e-12,'AbsTol',1e-14);

%time
N = 2048*4+1;
Tmax = 6.663550284765624;
dt=Tmax/(N-1);
t = 0:Tmax/(N-1):Tmax;

%make sure oscillators have identical periods for eps = 0.1
T_p = 6.710088675031566;
fac = Tmax/T_p;

%fixed parameters
eps = .1;

%coupling parameter (to vary)
delta = linspace(-0.51,0.15,21);

%period vectors
T1 = zeros(1,length(delta));
T2 = T1;

%% compute period as a function of delta

for j = 1:length(delta)
    
    %define vector field for each delta
    F = @(t,u) [fac*(u(1)-u(1)^3-u(2))+delta(j)*(u(3)-u(1));
        fac*u(1);
        u(3)-u(3)^3-u(4)+eps;
        u(3)+eps];
    
    %integration time
    if j > 13 && j < 20
        factor = 500;
    else
        factor = 120;
    end
    
    %solve for system
    x0 = [0.787155062518509  -1.065028290519620  -1.023916642131385   1.173324317167155];
    [T,U] = ode113(F,0:dt:1*T_p*factor,x0,opts);
    
    %ease of notation
    v1 = U(:,1);
    n1 = U(:,2);
    v2 = U(:,3);
    n2 = U(:,4);
    
    %chop off
    v1 = v1(round(end-25*Tmax/dt):end);
    v2 = v2(round(end-25*Tmax/dt):end);
        
    %find phases
    [~,phi1]=findpeaks(v1,'MinPeakHeight',0);
    [~,phi2]=findpeaks(v2,'MinPeakHeight',0);
    
    %take last 10 entries
    phi1 = phi1(end-5:end);
    phi2 = phi2(end-5:end);
    
    %find difference
    phase_lag = abs(phi2-phi1);
    M = max(phase_lag);
    phase_change = abs(phase_lag - M);
    
    for i=1:length(phase_change)
        if  phase_change(i) > 2
            fprintf('I have failed')
            delta(j)
            break
        end
    end
    
    %compute period 1st
    [~,loc]=findpeaks(v1,'MinPeakHeight',0);
    timerz = T(loc);
    T1(j) = mean(diff(timerz));
    
    %compute period 2nd
    [~,loc]=findpeaks(v2,'MinPeakHeight',0);
    timerz = T(loc);
    T2(j) = mean(diff(timerz));
    
    j
end


%% visualize

figure(1)
hold on
plot(delta,T2,'b-','LineWidth',2)
plot(delta,T2,'b.','MarkerSize',30)
plot(delta,T1,'k-','LineWidth',2)
plot(delta,T1,'k.','MarkerSize',30)
title('Period vs \delta')
xlabel('Coupling Strength \delta')
ylabel('Period')
box on
axis square
set(gca,'fontsize',12)
ylim([T_p-0.01 T_p+0.01])



