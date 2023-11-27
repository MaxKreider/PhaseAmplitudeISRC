%% make the plots

%load in data
ATongue = readmatrix('ATongueMatrix2.txt');

%extend the matrix for plotting purposes
ATongue2 = ATongue*30;
ATongue2 = [ones(50,200)*80; ATongue2];

%syntax
ATongue = flip(ATongue);
ATongue2 = flip(ATongue2);

%Figure 3.2
figure(1)
colormap parula
imagesc(ATongue)
xlabel('T_f / T_g')
xticks(100)
xticklabels({'1'})
ylabel('\delta')
yticks(200)
yticklabels({'0'})
set(gca,'fontsize',15)
box on
axis square
hold on
plot(100,15,'g.','MarkerSize',40)
plot([100 100],[15 200],'g-','LineWidth',1.5)
plot(100,100,'gv','MarkerFaceColor','g','MarkerSize',10)

%Figure 3.2
figure(2)
colormap parula
imagesc(ATongue)
xlabel('T_f / T_g')
xticks(100)
xticklabels({'1'})
ylabel('\delta')
yticks(200)
yticklabels({'0'})
set(gca,'fontsize',15)
box on
axis square
hold on
plot(130,15,'r.','MarkerSize',40)
plot([130 130],[15 200],'r-','LineWidth',1.5)
plot(130,100,'rv','MarkerFaceColor','r','MarkerSize',10)

%figure 3.3
figure(3)
colormap parula
imagesc(ATongue2)
xlabel('\beta')
xticks(100)
xticklabels({'0'})
ylabel('\delta')
yticks(200)
yticklabels({'0'})
set(gca,'fontsize',15)
axis square
hold on
plot(130,15,'color',[.7 .7 .7],'marker','.','MarkerSize',40)
text(128,25,'4','color',[.7 .7 .7])
plot(100,15,'color',[.7 .7 .7],'marker','.','MarkerSize',40)
text(103,25,'3','color',[.7 .7 .7])
plot(100,200,'color',[.7 .7 .7],'marker','.','MarkerSize',40)
text(97.5,212,'2','color',[.7 .7 .7])
plot(130,200,'color',[.7 .7 .7],'marker','.','MarkerSize',40)
text(127.5,212,'1','color',[.7 .7 .7])
plot([130 100],[15 15],'color',[.7 .7 .7],'LineWidth',1.5)
plot(115,15,'color',[.7 .7 .7],'marker','>','MarkerFaceColor',[.7 .7 .7],'MarkerSize',7)
plot([130 100],[200 200],'color',[.7 .7 .7],'LineWidth',1.5)
plot(115,200,'color',[.7 .7 .7],'marker','<','MarkerFaceColor',[.7 .7 .7],'MarkerSize',7)
plot([100 100],[15 200],'color',[.7 .7 .7],'LineWidth',1.5)
plot(100,107,'color',[.7 .7 .7],'marker','^','MarkerFaceColor',[.7 .7 .7],'MarkerSize',7)

c = [0.2422    0.1504    0.6603;
     0.9769    0.9839    0.0805;
     1 1 1];
 colormap(c)

