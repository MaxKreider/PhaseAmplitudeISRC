%% setup

%mr clean
clc
clf

%for the ODEs
format long


%% load data

load('K0_vdp_1.mat','K_0')
load('K1_vdp_1.mat','K_1')
load('K2_vdp_1.mat','K_2')
load('K3_vdp_1.mat','K_3')
load('K4_vdp_1.mat','K_4')
load('K5_vdp_1.mat','K_5')
load('K6_vdp_1.mat','K_6')
load('K7_vdp_1.mat','K_7')
load('K8_vdp_1.mat','K_8')
load('K9_vdp_1.mat','K_9')
load('K10_vdp_1.mat','K_10')
load('K11_vdp_1.mat','K_11')
load('K12_vdp_1.mat','K_12')
load('K13_vdp_1.mat','K_13')
load('K14_vdp_1.mat','K_14')
load('K15_vdp_1.mat','K_15')


%% organize

%time
N = 2048*8+1;
Tmax = 6.663550284765624;
t = 0:Tmax/(N-1):Tmax;

%store the Kn's nicely
Kx = [K_0(1,:); K_1(1,:); K_2(1,:); K_3(1,:); K_4(1,:); K_5(1,:); K_6(1,:); K_7(1,:); ...
      K_8(1,:); K_9(1,:); K_10(1,:); K_11(1,:); K_12(1,:); K_13(1,:); K_14(1,:); K_15(1,:)];
Ky = [K_0(2,:); K_1(2,:); K_2(2,:); K_3(2,:); K_4(2,:); K_5(2,:); K_6(2,:); K_7(2,:); ...
      K_8(2,:); K_9(2,:); K_10(2,:); K_11(2,:); K_12(2,:); K_13(2,:); K_14(2,:); K_15(2,:)];

%store their derivatives nicely
dKx = gradient(Kx,Tmax/(N-1));
dKy = gradient(Ky,Tmax/(N-1));


%% compute A_i

%initialize
A_0 = zeros(2,2,length(Kx));
A_1 = zeros(2,2,length(Kx));
A_2 = zeros(2,2,length(Kx));
A_3 = zeros(2,2,length(Kx));
A_4 = zeros(2,2,length(Kx));
A_5 = zeros(2,2,length(Kx));
A_6 = zeros(2,2,length(Kx));
A_7 = zeros(2,2,length(Kx));
A_8 = zeros(2,2,length(Kx));
A_9 = zeros(2,2,length(Kx));
A_10 = zeros(2,2,length(Kx));
A_11 = zeros(2,2,length(Kx));
A_12 = zeros(2,2,length(Kx));
A_13 = zeros(2,2,length(Kx));

%loop
for i=1:length(Kx)
    A_0(:,:,i) = [[dKx(1,i);dKy(1,i)],1*[Kx(2,i);Ky(2,i)]];
    A_1(:,:,i) = [[dKx(2,i);dKy(2,i)],2*[Kx(3,i);Ky(3,i)]];
    A_2(:,:,i) = [[dKx(3,i);dKy(3,i)],3*[Kx(4,i);Ky(4,i)]];
    A_3(:,:,i) = [[dKx(4,i);dKy(4,i)],4*[Kx(5,i);Ky(5,i)]];
    A_4(:,:,i) = [[dKx(5,i);dKy(5,i)],5*[Kx(6,i);Ky(6,i)]];
    A_5(:,:,i) = [[dKx(6,i);dKy(6,i)],6*[Kx(7,i);Ky(7,i)]];
    A_6(:,:,i) = [[dKx(7,i);dKy(7,i)],7*[Kx(8,i);Ky(8,i)]];
    A_7(:,:,i) = [[dKx(8,i);dKy(8,i)],8*[Kx(9,i);Ky(9,i)]];
    A_8(:,:,i) = [[dKx(9,i);dKy(9,i)],9*[Kx(10,i);Ky(10,i)]];
    A_9(:,:,i) = [[dKx(10,i);dKy(10,i)],10*[Kx(11,i);Ky(11,i)]];
    A_10(:,:,i) = [[dKx(11,i);dKy(11,i)],11*[Kx(12,i);Ky(12,i)]];
    A_11(:,:,i) = [[dKx(12,i);dKy(12,i)],12*[Kx(13,i);Ky(13,i)]];
    A_12(:,:,i) = [[dKx(13,i);dKy(13,i)],13*[Kx(14,i);Ky(14,i)]];
    A_13(:,:,i) = [[dKx(14,i);dKy(14,i)],14*[Kx(15,i);Ky(15,i)]];
end


%% Z0, I0

%compute
z0K = zeros(2,length(Kx));
I0K = zeros(2,length(Kx));
for i=1:length(Kx)
    
    %det
    big_D = det(A_0(:,:,i));
    
    %z0
    z0K(:,i) = [Ky(2,i);-Kx(2,i)]/big_D;
    
    %I0
    I0K(:,i) = [-dKy(1,i);dKx(1,i)]/big_D;
end

% %visualize
% figure(100)
% hold on
% plot(t,I0K(1,:)/max(abs(I0K(1,:))),'k-','LineWidth',5)
% plot(t,I0K(2,:)/max(abs(I0K(2,:))),'k-','LineWidth',5)
% plot(t_lo,I0(1,:)/max(abs(I0(1,:))),'r--','LineWidth',3)
% plot(t_lo,I0(2,:)/max(abs(I0(2,:))),'c--','LineWidth',3)
% title('I0 Comparison')
% xlabel('phase \phi')
% ylabel('I')
% box on 
% axis square
% set(gca,'fontsize',15)
% xlim([0 6.664])


%% Z1, I1

%compute
soll = zeros(2,2,length(Kx));
for i=1:length(Kx)
    %solve
    sol1(:,:,i) = -[z0K(:,i)';I0K(:,i)']*A_1(:,:,i)*[z0K(:,i)';I0K(:,i)'];
end

%store
z1K = squeeze(sol1(1,:,:));
I1K = squeeze(sol1(2,:,:));

% figure(101)
% hold on
% plot(t,I1K(1,:)/max(abs(I1K(1,:))),'k-','LineWidth',5)
% plot(t,I1K(2,:)/max(abs(I1K(2,:))),'k-','LineWidth',5)
% plot(t_lo,I1(1,:)/max(abs(I1(1,:))),'r--','LineWidth',3)
% plot(t_lo,I1(2,:)/max(abs(I1(2,:))),'c--','LineWidth',3)
% title('I1 Comparison')
% xlabel('phase \phi')
% ylabel('I')
% box on 
% axis square
% set(gca,'fontsize',15)
% xlim([0 6.664])


%% Z2, I2

%compute
soll = zeros(2,2,length(Kx));
for i=1:length(Kx)
    %solve
    sol1(:,:,i) = -[z0K(:,i)';I0K(:,i)']*(A_2(:,:,i)*[z0K(:,i)';I0K(:,i)'] + A_1(:,:,i)*[z1K(:,i)';I1K(:,i)']);
end

%store
z2K = squeeze(sol1(1,:,:));
I2K = squeeze(sol1(2,:,:));

% figure(102)
% hold on
% plot(t,I2K(1,:)/max(abs(I2K(1,:))),'k-','LineWidth',5)
% plot(t,I2K(2,:)/max(abs(I2K(2,:))),'k-','LineWidth',5)
% plot(t_lo,I2(1,:)/max(abs(I2(1,:))),'r--','LineWidth',3)
% plot(t_lo,I2(2,:)/max(abs(I2(2,:))),'c--','LineWidth',3)
% title('I2 Comparison')
% xlabel('phase \phi')
% ylabel('z')
% box on 
% axis square
% set(gca,'fontsize',15)
% xlim([0 6.664])


%% Z3, I3

%compute
soll = zeros(2,2,length(Kx));
for i=1:length(Kx)
    %solve
    sol1(:,:,i) = -[z0K(:,i)';I0K(:,i)']*(A_3(:,:,i)*[z0K(:,i)';I0K(:,i)'] + A_2(:,:,i)*[z1K(:,i)';I1K(:,i)'] + A_1(:,:,i)*[z2K(:,i)';I2K(:,i)']);
end

%store
z3K = squeeze(sol1(1,:,:));
I3K = squeeze(sol1(2,:,:));

% figure(103)
% hold on
% plot(t,I3K(1,:)/max(abs(I3K(1,:))),'k-','LineWidth',5)
% plot(t,I3K(2,:)/max(abs(I3K(2,:))),'k-','LineWidth',5)
% plot(t_lo,I3(1,:)/max(abs(I3(1,:))),'r--','LineWidth',3)
% plot(t_lo,I3(2,:)/max(abs(I3(2,:))),'c--','LineWidth',3)
% title('I3 Comparison')
% xlabel('phase \phi')
% ylabel('I')
% box on 
% axis square
% set(gca,'fontsize',15)
% xlim([0 6.664])


%% Z4, I4

%compute
soll = zeros(2,2,length(Kx));
for i=1:length(Kx)
    %solve
    sol1(:,:,i) = -[z0K(:,i)';I0K(:,i)']*(A_4(:,:,i)*[z0K(:,i)';I0K(:,i)'] + A_3(:,:,i)*[z1K(:,i)';I1K(:,i)'] + ...
                    A_2(:,:,i)*[z2K(:,i)';I2K(:,i)'] + A_1(:,:,i)*[z3K(:,i)';I3K(:,i)']);
end

%store
z4K = squeeze(sol1(1,:,:));
I4K = squeeze(sol1(2,:,:));

% figure(104)
% hold on
% plot(t,I4K(1,:)/max(abs(I4K(1,:))),'k-','LineWidth',5)
% plot(t,I4K(2,:)/max(abs(I4K(2,:))),'k-','LineWidth',5)
% plot(t_lo,I4(1,:)/max(abs(I4(1,:))),'r--','LineWidth',3)
% plot(t_lo,I4(2,:)/max(abs(I4(2,:))),'c--','LineWidth',3)
% title('I4 Comparison')
% xlabel('phase \phi')
% ylabel('I')
% box on 
% axis square
% set(gca,'fontsize',15)
% xlim([0 6.664])


%% Z5, I5

%compute
soll = zeros(2,2,length(Kx));
for i=1:length(Kx)
    %solve
    sol1(:,:,i) = -[z0K(:,i)';I0K(:,i)']*(A_5(:,:,i)*[z0K(:,i)';I0K(:,i)'] + A_4(:,:,i)*[z1K(:,i)';I1K(:,i)'] + ...
                    A_3(:,:,i)*[z2K(:,i)';I2K(:,i)'] + A_2(:,:,i)*[z3K(:,i)';I3K(:,i)'] + A_1(:,:,i)*[z4K(:,i)';I4K(:,i)']);
end

%store
z5K = squeeze(sol1(1,:,:));
I5K = squeeze(sol1(2,:,:));

% figure(105)
% hold on
% plot(t,I5K(1,:)/max(abs(I5K(1,:))),'k-','LineWidth',5)
% plot(t,I5K(2,:)/max(abs(I5K(2,:))),'k-','LineWidth',5)
% plot(t_lo,I5(1,:)/max(abs(I5(1,:))),'r--','LineWidth',3)
% plot(t_lo,I5(2,:)/max(abs(I5(2,:))),'c--','LineWidth',3)
% title('I5 Comparison')
% xlabel('phase \phi')
% ylabel('I')
% box on 
% axis square
% set(gca,'fontsize',15)
% xlim([0 6.664])


%% Z6, I6

%compute
soll = zeros(2,2,length(Kx));
for i=1:length(Kx)
    %solve
    sol1(:,:,i) = -[z0K(:,i)';I0K(:,i)']*(A_6(:,:,i)*[z0K(:,i)';I0K(:,i)'] + A_5(:,:,i)*[z1K(:,i)';I1K(:,i)'] + ...
                    A_4(:,:,i)*[z2K(:,i)';I2K(:,i)'] + A_3(:,:,i)*[z3K(:,i)';I3K(:,i)'] + A_2(:,:,i)*[z4K(:,i)';I4K(:,i)'] + ...
                    A_1(:,:,i)*[z5K(:,i)';I5K(:,i)']);
end

%store
z6K = squeeze(sol1(1,:,:));
I6K = squeeze(sol1(2,:,:));

% figure(106)
% hold on
% plot(t,I6K(1,:)/max(abs(I6K(1,:))),'k-','LineWidth',5)
% plot(t,I6K(2,:)/max(abs(I6K(2,:))),'k-','LineWidth',5)
% plot(t_lo,I6(1,:)/max(abs(I6(1,:))),'r--','LineWidth',3)
% plot(t_lo,I6(2,:)/max(abs(I6(2,:))),'c--','LineWidth',3)
% title('I6 Comparison')
% xlabel('phase \phi')
% ylabel('I')
% box on 
% axis square
% set(gca,'fontsize',15)
% xlim([0 6.664])


%% Z7, I7

%compute
soll = zeros(2,2,length(Kx));
for i=1:length(Kx)
    %solve
    sol1(:,:,i) = -[z0K(:,i)';I0K(:,i)']*(A_7(:,:,i)*[z0K(:,i)';I0K(:,i)'] + A_6(:,:,i)*[z1K(:,i)';I1K(:,i)'] + ...
                    A_5(:,:,i)*[z2K(:,i)';I2K(:,i)'] + A_4(:,:,i)*[z3K(:,i)';I3K(:,i)'] + A_3(:,:,i)*[z4K(:,i)';I4K(:,i)'] + ...
                    A_2(:,:,i)*[z5K(:,i)';I5K(:,i)'] + A_1(:,:,i)*[z6K(:,i)';I6K(:,i)']);
end

%store
z7K = squeeze(sol1(1,:,:));
I7K = squeeze(sol1(2,:,:));

% figure(107)
% hold on
% plot(t,I7K(1,:)/max(abs(I7K(1,:))),'k-','LineWidth',5)
% plot(t,I7K(2,:)/max(abs(I7K(2,:))),'k-','LineWidth',5)
% plot(t_lo,I7(1,:)/max(abs(I7(1,:))),'r--','LineWidth',3)
% plot(t_lo,I7(2,:)/max(abs(I7(2,:))),'c--','LineWidth',3)
% title('I7 Comparison')
% xlabel('phase \phi')
% ylabel('I')
% box on 
% axis square
% set(gca,'fontsize',15)
% xlim([0 6.664])


%% Z8, I8

%compute
soll = zeros(2,2,length(Kx));
for i=1:length(Kx)
    %solve
    sol1(:,:,i) = -[z0K(:,i)';I0K(:,i)']*(A_8(:,:,i)*[z0K(:,i)';I0K(:,i)'] + A_7(:,:,i)*[z1K(:,i)';I1K(:,i)'] + ...
                    A_6(:,:,i)*[z2K(:,i)';I2K(:,i)'] + A_5(:,:,i)*[z3K(:,i)';I3K(:,i)'] + A_4(:,:,i)*[z4K(:,i)';I4K(:,i)'] + ...
                    A_3(:,:,i)*[z5K(:,i)';I5K(:,i)'] + A_2(:,:,i)*[z6K(:,i)';I6K(:,i)'] + A_1(:,:,i)*[z7K(:,i)';I7K(:,i)']);
end

%store
z8K = squeeze(sol1(1,:,:));
I8K = squeeze(sol1(2,:,:));

% figure(108)
% hold on
% plot(t,I8K(1,:)/max(abs(I8K(1,:))),'k-','LineWidth',5)
% plot(t,I8K(2,:)/max(abs(I8K(2,:))),'k-','LineWidth',5)
% plot(t_lo,I8(1,:)/max(abs(I8(1,:))),'r--','LineWidth',3)
% plot(t_lo,I8(2,:)/max(abs(I8(2,:))),'c--','LineWidth',3)
% title('I8 Comparison')
% xlabel('phase \phi')
% ylabel('I')
% box on 
% axis square
% set(gca,'fontsize',15)
% xlim([0 6.664])


%% Z8, I8

%compute
soll = zeros(2,2,length(Kx));
for i=1:length(Kx)
    %solve
    sol1(:,:,i) = -[z0K(:,i)';I0K(:,i)']*(A_8(:,:,i)*[z0K(:,i)';I0K(:,i)'] + A_7(:,:,i)*[z1K(:,i)';I1K(:,i)'] + ...
                    A_6(:,:,i)*[z2K(:,i)';I2K(:,i)'] + A_5(:,:,i)*[z3K(:,i)';I3K(:,i)'] + A_4(:,:,i)*[z4K(:,i)';I4K(:,i)'] + ...
                    A_3(:,:,i)*[z5K(:,i)';I5K(:,i)'] + A_2(:,:,i)*[z6K(:,i)';I6K(:,i)'] + A_1(:,:,i)*[z7K(:,i)';I7K(:,i)']);
end

%store
z8K = squeeze(sol1(1,:,:));
I8K = squeeze(sol1(2,:,:));

% figure(108)
% hold on
% plot(t,I8K(1,:)/max(abs(I8K(1,:))),'k-','LineWidth',5)
% plot(t,I8K(2,:)/max(abs(I8K(2,:))),'k-','LineWidth',5)
% plot(t_lo,I8(1,:)/max(abs(I8(1,:))),'r--','LineWidth',3)
% plot(t_lo,I8(2,:)/max(abs(I8(2,:))),'c--','LineWidth',3)
% title('I8 Comparison')
% xlabel('phase \phi')
% ylabel('I')
% box on 
% axis square
% set(gca,'fontsize',15)
% xlim([0 6.664])


%% Z9, I9

%compute
soll = zeros(2,2,length(Kx));
for i=1:length(Kx)
    %solve
    sol1(:,:,i) = -[z0K(:,i)';I0K(:,i)']*(A_9(:,:,i)*[z0K(:,i)';I0K(:,i)'] + A_8(:,:,i)*[z1K(:,i)';I1K(:,i)'] + ...
                    A_7(:,:,i)*[z2K(:,i)';I2K(:,i)'] + A_6(:,:,i)*[z3K(:,i)';I3K(:,i)'] + A_5(:,:,i)*[z4K(:,i)';I4K(:,i)'] + ...
                    A_4(:,:,i)*[z5K(:,i)';I5K(:,i)'] + A_3(:,:,i)*[z6K(:,i)';I6K(:,i)'] + A_2(:,:,i)*[z7K(:,i)';I7K(:,i)'] + ...
                    A_1(:,:,i)*[z8K(:,i)';I8K(:,i)']);
end

%store
z9K = squeeze(sol1(1,:,:));
I9K = squeeze(sol1(2,:,:));

% figure(109)
% hold on
% plot(t,I9K(1,:)/max(abs(I9K(1,:))),'k-','LineWidth',5)
% plot(t,I9K(2,:)/max(abs(I9K(2,:))),'k-','LineWidth',5)
% title('I9 Comparison')
% xlabel('phase \phi')
% ylabel('I')
% box on 
% axis square
% set(gca,'fontsize',15)
% xlim([0 6.664])


%% Z10, I10

%compute
soll = zeros(2,2,length(Kx));
for i=1:length(Kx)
    %solve
    sol1(:,:,i) = -[z0K(:,i)';I0K(:,i)']*(A_10(:,:,i)*[z0K(:,i)';I0K(:,i)'] + A_9(:,:,i)*[z1K(:,i)';I1K(:,i)'] + ...
                    A_8(:,:,i)*[z2K(:,i)';I2K(:,i)'] + A_7(:,:,i)*[z3K(:,i)';I3K(:,i)'] + A_6(:,:,i)*[z4K(:,i)';I4K(:,i)'] + ...
                    A_5(:,:,i)*[z5K(:,i)';I5K(:,i)'] + A_4(:,:,i)*[z6K(:,i)';I6K(:,i)'] + A_3(:,:,i)*[z7K(:,i)';I7K(:,i)'] + ...
                    A_2(:,:,i)*[z8K(:,i)';I8K(:,i)'] + A_1(:,:,i)*[z9K(:,i)';I9K(:,i)']);
end

%store
z10K = squeeze(sol1(1,:,:));
I10K = squeeze(sol1(2,:,:));

% figure(109)
% hold on
% plot(t,I10K(1,:)/max(abs(I10K(1,:))),'k-','LineWidth',5)
% plot(t,I10K(2,:)/max(abs(I10K(2,:))),'k-','LineWidth',5)
% title('I10 Comparison')
% xlabel('phase \phi')
% ylabel('I')
% box on 
% axis square
% set(gca,'fontsize',15)
% xlim([0 6.664])


%% save data

save('z0K_vdp_1.mat','z0K')
save('z1K_vdp_1.mat','z1K')
save('z2K_vdp_1.mat','z2K')
save('z3K_vdp_1.mat','z3K')
save('z4K_vdp_1.mat','z4K')
save('z5K_vdp_1.mat','z5K')
save('z6K_vdp_1.mat','z6K')
save('z7K_vdp_1.mat','z7K')
save('z8K_vdp_1.mat','z8K')

save('I0K_vdp_1.mat','I0K')
save('I1K_vdp_1.mat','I1K')
save('I2K_vdp_1.mat','I2K')
save('I3K_vdp_1.mat','I3K')
save('I4K_vdp_1.mat','I4K')
save('I5K_vdp_1.mat','I5K')
save('I6K_vdp_1.mat','I6K')
save('I7K_vdp_1.mat','I7K')
save('I8K_vdp_1.mat','I8K')


%% 3D plots (Alberto)

%angle and amplitude vectors
theta = flip(t);
sigma = -.75:.01:.75;

%angle and amplitude meshes
[THETA,SIGMA] = meshgrid(theta,sigma);

%8th order PRC (Alberto)
Zx8 = z0K(1,:) + SIGMA.*z1K(1,:) + SIGMA.^2.*z2K(1,:) + SIGMA.^3.*z3K(1,:) + ...
      SIGMA.^4.*z4K(1,:) + SIGMA.^5.*z5K(1,:) + SIGMA.^6.*z6K(1,:) + SIGMA.^7.*z7K(1,:) + SIGMA.^8.*z8K(1,:);
Zy8 = z0K(2,:) + SIGMA.*z1K(2,:) + SIGMA.^2.*z2K(2,:) + SIGMA.^3.*z3K(2,:) + ...
      SIGMA.^4.*z4K(2,:) + SIGMA.^5.*z5K(2,:) + SIGMA.^6.*z6K(2,:) + SIGMA.^7.*z7K(2,:) + SIGMA.^8.*z8K(2,:); 
%   
%8th order
Ix8 = I0K(1,:) + SIGMA.*I1K(1,:) + SIGMA.^2.*I2K(1,:) + SIGMA.^3.*I3K(1,:) + ...
      SIGMA.^4.*I4K(1,:) + SIGMA.^5.*I5K(1,:) + SIGMA.^6.*I6K(1,:) + SIGMA.^7.*I7K(1,:) + SIGMA.^8.*I8K(1,:);
Iy8 = I0K(2,:) + SIGMA.*I1K(2,:) + SIGMA.^2.*I2K(2,:) + SIGMA.^3.*I3K(2,:) + ...
      SIGMA.^4.*I4K(2,:) + SIGMA.^5.*I5K(2,:) + SIGMA.^6.*I6K(2,:) + SIGMA.^7.*I7K(2,:) + SIGMA.^8.*I8K(2,:);   

figure(1)
surf(SIGMA,THETA,Zx8,'EdgeColor','none')
title('8th Order PRC (x) - P')
ylabel('phase \theta')
xlabel('amplitude \psi')
box on 
axis square

figure(2)
surf(SIGMA,THETA,Zy8,'EdgeColor','none')
title('8th Order PRC (y) - P')
ylabel('phase \theta')
xlabel('amplitude \psi')
box on 
axis square  
  
figure(3)
surf(SIGMA,THETA,Ix8,'EdgeColor','none')
title('8th Order IRC (x) - P')
ylabel('phase \theta')
xlabel('amplitude \psi')
box on 
axis square

figure(4)
surf(SIGMA,THETA,Iy8,'EdgeColor','none')
title('8th Order IRC (y) - P')
ylabel('phase \theta')
xlabel('amplitude \psi')
box on 
axis square
