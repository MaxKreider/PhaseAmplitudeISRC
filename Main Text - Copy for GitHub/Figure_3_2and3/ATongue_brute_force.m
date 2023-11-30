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

%parameters
beta = linspace(-.4,.4,200);
delta = logspace(log(0.11),log(.5),200);
omega=2*pi/Tmax;

%fail indicator
fail = 1;

%solution matrix
sol_matrix = ones(length(delta),length(beta));


%% system

for p=1:length(delta)
    tic
    for q=1:length(beta)
        p
        q
        %unperturbed vector field
        x0 = [0.365302629541246  -0.606783646555696   0.731381671068370  -1.071908705393243];
        F_p = @(t,u) [(.5-u(1)^2-u(2)^2)*u(1) - omega*u(2) + delta(p)*(u(3)-u(1));
            (.5-u(1)^2-u(2)^2)*u(2) + omega*u(1);
            u(3)-u(3)^3-u(4) + delta(p)*(u(1)-u(3)) + beta(q)*(u(3)-u(3)^3-u(4));
            u(3) + beta(q)*(u(3))];
        
        %solve
        [T,U]=ode113(F_p,0:dt:Tmax*50,x0,opts);
        
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
                fail = 0;
                break
            end
        end
        
        %update
        sol_matrix(p,q) = fail;
        fail = 1;
        
    end
    toc
end

figure(1)
imagesc(sol_matrix)
colorbar

writematrix(sol_matrix,'AtongueMatrix2.txt')

