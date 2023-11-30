%% setup

%mr clean
clc
clf

%for the ODEs
format long
opts = odeset('RelTol',1e-10,'AbsTol',1e-12);

%time
N = 2048+1;
Tmax = 2*pi;
dt = Tmax/(N-1);
t = 0:dt:Tmax;

%parameter
mu = 1;
omega = 1;

%Jacobian
Jac = @(t,u) [-3*u(1)^2-u(2)^2+mu, -2*u(1)*u(2)-omega;
    -2*u(1)*u(2)+omega, -3*u(2)^2-u(1)^2+mu];


%% dynamics

%eps_vec
evec = -.95:.05:1.95;

%vector of amplitudes
true_amp_x = zeros(1,length(evec));
iSRC_amp_x = zeros(1,length(evec));
approx_amp_x = zeros(1,length(evec));

%count
count = 1;

%loop it
for eps=-.95:.05:1.95
        
    %system
    F = @(t,u) [-u(1)^3-u(1)*u(2)^2-omega*u(2)+mu*u(1); -u(2)^3-u(1)^2*u(2)+omega*u(1)+mu*u(2)];
    F_p = @(t,u) [-u(1)^3-u(1)*u(2)^2-omega*u(2)+(mu+eps)*u(1); -u(2)^3-u(1)^2*u(2)+omega*u(1)+(mu+eps)*u(2)];
    
    %initial conditions
    x0 = [sqrt(mu) 0];
    x01 = [sqrt(mu+eps) 0];
    [T_u,U_u] = ode113(F,0:dt:Tmax,x0,opts);
    [T_p,U_p] = ode113(F_p,0:dt:Tmax,x01,opts);
    
%     %visualize
%     figure(69)
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
    
    
    %% compute change in timing nu
    
    %normal forms are good for something I suppose
    nu = 0;
    
    
    %% poincare section
    
    %define variational equation to find floquet multipliers and vectors
    F_v = @(t,u) [-u(1)^3-u(1)*u(2)^2-omega*u(2)+mu*u(1); ...
        -u(2)^3-u(1)^2*u(2)+omega*u(1)+mu*u(2); ...
        (-3*u(1)^2-u(2)^2+mu)*u(3)+(-2*u(1)*u(2)-omega)*u(5); ...
        (-3*u(1)^2-u(2)^2+mu)*u(4)+(-2*u(1)*u(2)-omega)*u(6); ...
        (-2*u(1)*u(2)+omega)*u(3)+(-3*u(2)^2-u(1)^2+mu)*u(5); ...
        (-2*u(1)*u(2)+omega)*u(4)+(-3*u(2)^2-u(1)^2+mu)*u(6)];
    
    %solve variation equation
    [T,U] = ode113(F_v,0:dt:Tmax,[x0 1 0 0 1],opts);
    
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
    theta = 0;
    R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
    q2 = R*q2;
    n=null(q2');
    
%     %plot
%     figure(69)
%     plot([x0(1)-q2(1),x0(1)+q2(1)],[x0(2)-q2(2),x0(2)+q2(2)],'color','g','LineWidth',1.5);
    
    
    %% IC for iSRC
    
    %RHS
    FM = [];
    for j=1:length(U)
        
        %create df/deps
        dfdeps = [U(j,1);U(j,2)];
        
        %create fundamental matrix solution
        PHI = [U(j,3), U(j,4); U(j,5), U(j,6)];
        
        %create vector slot for FM
        b = F_v(0,[U(j,1), U(j,2), U(j,3), U(j,4), U(j,5), U(j,6)]);
        temp = PHI\(dfdeps+nu*b(1:2));
        
        %FM
        FM = [FM temp];
        
    end
    
    %integrate
    FM = [cumtrapz(dt,FM(1,:));cumtrapz(dt,FM(2,:))];
    
    %pls god pls
    term3=[];
    for i=1:length(T)
        
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
    ISRC = @(t,u) [-u(1)^3-u(1)*u(2)^2-omega*u(2)+mu*u(1); ...
        -u(2)^3-u(1)^2*u(2)+omega*u(1)+mu*u(2);
        (-3*u(1)^2-u(2)^2+mu)*u(3)+(-2*u(1)*u(2)-omega)*u(4) + nu*(-u(1)^3-u(1)*u(2)^2-omega*u(2)+mu*u(1)) + u(1);
        (-2*u(1)*u(2)+omega)*u(3)+(-3*u(2)^2-u(1)^2+mu)*u(4) + nu*(-u(2)^3-u(1)^2*u(2)+omega*u(1)+mu*u(2)) + u(2)];
    
    %skip a few cycles
    tran=2;
    
    %initial cond
    IC=p1_lul';
    
    %solve
    for k=1:tran
        
        %solve
        [T_src,U_src] = ode113(ISRC,0:dt:Tmax,[x0 IC],opts);
        
        %periodicity shooting thingy
        IC = U_src(end,3:4);
        
    end
    
    %perturbed LC
    gamma_eps_v = U_src(:,1) + U_src(:,3)*eps;
    gamma_eps_n = U_src(:,2) + U_src(:,4)*eps;
    
%     %plot
%     figure(69)
%     plot(gamma_eps_v,gamma_eps_n,'--','color','c','LineWidth',2)
    
    
    %% max :)
    
    %original system
    [~,by] = max(U_u(:,1));
    max0 = U_u(by,:);
    
    %perturbed system
    [~,my] = max(U_p(:,1));
    maxK = U_p(my,:);
    
    %iSRC
    [~,tie] = max(gamma_eps_v);
    maxQ = [gamma_eps_v(tie);gamma_eps_n(tie)];
    
%     %plot it
%     figure(69)
%     plot(max0(1),max0(2),'k.','MarkerSize',40)
%     plot(maxK(1),maxK(2),'r.','MarkerSize',40)
%     plot(maxQ(1),maxQ(2),'c.','MarkerSize',40)
    
    %time at which max occurs
    tm = by;
    
    
    %% min :(
    
    %original system
    [~,bym] = min(U_u(:,1));
    max0m = U_u(bym,:);
    
    %perturbed system
    [~,mym] = min(U_p(:,1));
    maxKm = U_p(mym,:);
    
    %iSRC
    [~,tiem] = min(gamma_eps_v);
    maxQm = [gamma_eps_v(tiem);gamma_eps_n(tiem)];
    
%     %plot it
%     figure(69)
%     plot(max0m(1),max0m(2),'k.','MarkerSize',40)
%     plot(maxKm(1),maxKm(2),'r.','MarkerSize',40)
%     plot(maxQm(1),maxQm(2),'c.','MarkerSize',40)
    
    %time at which min occurs
    tm2 = bym;
    
    
    %% my boi Newton
    
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
        g_1_p(:,i) = Jac(t(i),U_u(i,:))*U_src(i,3:4)' + nu*F(t(i),U_u(i,:)) + [U_u(i,1);U_u(i,2)];
    end
    
    %second derivative of gamma_1
    g_1_pp = gradient(g_1_p,dt);
    
    %solve for delta
    delta = -(g_0_p(1,1)+eps*g_1_p(1,1))/(g_0_pp(1,1)+eps*g_1_pp(1,1));
    
    %change it to an index
    delta_I = round(delta/dt);
    
    %write out new index bc this is annoying af
    new_index_max = mod(tm+delta_I,length(U_src(:,3)));
    if new_index_max == 0
        new_index_max = length(U_src(:,3))-1;
    end
    new_index_min = mod(tm2+delta_I,length(U_src(:,3)));
    if new_index_min == 0
        new_index_min = length(U_src(:,3))-1;
    end
    
%     %plot it
%     figure(69)
%     plot(gamma_eps_v(new_index_max),gamma_eps_n(new_index_max),'g.','MarkerSize',40)
%     plot(gamma_eps_v(new_index_min),gamma_eps_n(new_index_min),'g.','MarkerSize',40)
        
    %keep track of amplitude (max-min)
    true_amp_x(count) = maxK(1);
    iSRC_amp_x(count) = maxQ(1);
    approx_amp_x(count) = gamma_eps_v(new_index_max);
    
    %update count
    count = count + 1;
    
end

%visualize results
figure(1)
hold on
plot(evec,true_amp_x,'k.','MarkerSize',40)
plot(evec,true_amp_x,'k-','LineWidth',4)
plot(evec,iSRC_amp_x,'c.','MarkerSize',40)
plot(evec,iSRC_amp_x,'c-','LineWidth',4)
plot(evec,approx_amp_x,'color',[255, 182, 193 256]/256,'Marker','.','MarkerSize',40)
plot(evec,approx_amp_x,'color',[255, 182, 193 256]/256,'LineWidth',4)
plot(evec,1/2*(evec)+1,'color',[.5 .5 .5],'LineWidth',4)
title('Amplitude vs. Perturbation')
xlabel('\epsilon')
ylabel('Amplitude')
box on
axis square
set(gca,'fontsize',15)

