%% load data 

load real3d.mat

%% Patching data
w1 =12;
w2 =12;
w3 =12;
s1z =4;
s2z =2;
s3z =2;
dn_patch = yc_patch3d(DataNoisey,1,w1,w2,w3,s1z,s2z,s3z);


% It is better to save .mat as -V 7.3 or later because the size of the generated patch is large.
save('Input_Patches_3Dreal1','dn_patch','-V7.3');
filename = 'Input_Patches_3Dreal1.csv';% 你的CSV数据文件名
writematrix(dn_patch, filename, 'Delimiter', ',');




c1 = csvread('output_3dreal1.csv');
w1 =12;
w2 =12;
w3 =12;
s1z =4;
s2z =2;
s3z = 2;
out1=yc_patch3d_inv(c1',1,n1,n2,n3,w1,w2,w3,s1z,s2z,s3z);


c1 = csvread('output_3dreal1.csv');

% Loading the Data and the Output Patches of the DDUL
% Synthetic or Field Example, zero for synthetic and one for field example.

% UnPatching
[n1,n2,n3] = size(DataClean);
w1 =12;
w2 =12;
w3 =12;
s1z =4;
s2z =2;
s3z = 2;
out1=yc_patch3d_inv(c1',1,n1,n2,n3,w1,w2,w3,s1z,s2z,s3z);

%% process with mssa \svmf \patch 

real1demssa=fxymssa(DataNoisy,5,120,0.002,3,0);

outc = csvread('output_2dreal1patch.csv');
[n1,n2]=size(Dn);
w1 =32;
w2 =32;
s1z =2;
s2z =2;
real1depatch=yc_patch_inv(outc',1,n1,n2,w1,w2,s1z,s2z);

%% image

figure;yc_imagesc(squeeze(real2_ls_mssa(:,:,18)))
title('MSSA')
clim([0,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:10:61]);
set(gca,'XTickLabel',[0:10:61],'FontSize',12);
set(gca,'YTick',[0:20:181]);
set(gca,'YTickLabel',[0:20*dt:181*dt],'FontSize',12);
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
colormap('jet')
img =gcf; 
print(img, '-dpng', '-r600', './real2_ls_mssa.png')  
