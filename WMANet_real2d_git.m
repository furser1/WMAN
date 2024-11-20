%% load data 

load real2d.mat

%% Patching data
w1 =32;
w2 =32;
s1z =2;
s2z =2;
dn_patch = yc_patch(DataNoisy,1,w1,w2,s1z,s2z);
save('Input_Patches_2Dreal1','dn_patch','-V7.3');
filename = 'Input_Patches_2Dreal1.csv';
writematrix(dn_patch, filename, 'Delimiter', ',');


%% process with mssa \svmf \patch 

real1demssa=fxymssa(DataNoisy,5,120,0.002,3,0);

outc = csvread('output_2dreal1patch.csv');
[n1,n2]=size(Dn);
w1 =32;
w2 =32;
s1z =2;
s2z =2;
real1depatch=yc_patch_inv(outc',1,n1,n2,w1,w2,s1z,s2z);


figure;yc_imagesc(DataNoisy)
title('Noisy')
clim([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:50:260]);
set(gca,'XTickLabel',[0:50:260],'FontSize',12);
set(gca,'YTick',[0:50:513]);
set(gca,'YTickLabel',[0:50*dt:513*dt],'FontSize',12);
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  
print(img, '-dpng', '-r600', './real1datan.png')    

