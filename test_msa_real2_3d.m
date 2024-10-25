 seismic_data1= readsegyfile('Kerry3D.segy');
  dt=0.0004;
dat = reshape(seismic_data1,1252,735,287);
for i =1:6
slx = 30*i+50;
sly = 200;
Sx_dat(:,:,i) = (dat(:,slx,:));
Sy_dat(:,:,i) = (dat(:,:,sly));
end

figure;
for i=1:6
    subplot(2,3,i)
    imagesc(squeeze(Sx_dat(:,:,i)))
    subtitle(sprintf("x slice is %d",30*i+50))
end
colormap(seis(3))

figure;
for i=1:6
    subplot(2,3,i)
    imagesc(squeeze(Sy_dat(:,:,i)))
    subtitle(sprintf("y slice is %d",30*i+50))
end
colormap(seis(3))


figure;
for i=1:6
    subplot(2,3,i)
    imagesc(squeeze(dat1(:,i,:)))
    subtitle(sprintf("x slice is %d",30*i+50))
end
colormap(seis(3))

dat1=dat(420:540,220:280,120:180);


figure;
for i=1:6
    subplot(2,3,i)
    imagesc(squeeze(dat1(:,i,:)))
    subtitle(sprintf("y slice is %d",30*i+50))
end
colormap(seis(3))

save datreal.mat dat1
csvwrite('datreal.csv',dat1);

clc;
clear;
close all;
cd D:\Desktop\MSAnet\realdata
load("datreal.mat");
DataClean=dat1;
cd ..
% preparing the patches where w1, w2, and w3 are the patch size, while the s1z, s2z, and s3z are the number of shift samples between neighbor windows. 
% The default values are w1,w2,w3 =15,15,15 and s1z,s2z,s3z=1,1,1.
rng('default')
rng(2023055);
[Nt,Nx,Ny] = size(DataClean);

Dc=reshape(DataClean,Nt,Nx*Ny);
Dc0=scalecyk(Dc,1);
Dc=Dc0;
figure;
imagesc(Dc)
Input = zeros(size(Dc));
mask=rand(1,Nx*Ny);
mask(logical(mask<0.95))=0;
mask(logical(mask>=0.95))=1;
Input0 = zeros(size(Input));
for i=1:Nt
    Input0(i,:)=Dc(i,:)+randn(1,Nx*Ny).*mask*1.3;
end

Input=Input0+0.2*randn(Nt,Nx*Ny);


figure;
imagesc(reshape(Dc,Nt,Nx*Ny))
caxis([-1,1])
colormap(seis(1))


Dn= Input;
figure;
imagesc(reshape(Dn,Nt,Nx*Ny))
caxis([-1,1])
colormap(seis(1))

figure;
imagesc(reshape(Dn-Dc,Nt,Nx*Ny))
caxis([-1,1])
colormap(seis(1))



fid = fopen('real2datanoisy.bin','wb');
fwrite(fid,Dn,'float');
fclose(fid);
figure;
plot(Dn(:,15));
hold on;
plot(Dc(:,15));
hold on;

DataNoisey=reshape(Dn,Nt,Nx,Ny);
w1 =12;
w2 =12;
w3 =12;
s1z =4;
s2z =2;
s3z =2;
dn_patch = yc_patch3d(DataNoisey,1,w1,w2,w3,s1z,s2z,s3z);



% It is better to save .mat as -V 7.3 or later because the size of the generated patch is large.
save('Input_Patches_3Dreal1','dn_patch','-V7.3');
filename = 'Input_Patches_3Dreal2.csv';% 你的CSV数据文件名
writematrix(dn_patch, filename, 'Delimiter', ',');



outa = csvread('output_3dreal1.csv');

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
real2dewma=yc_patch3d_inv(outa',1,n1,n2,n3,w1,w2,w3,s1z,s2z,s3z);
%% fmssa  \ patchunet
real2demssa=fxymssa(DataNoisey,5,120,0.002,3,0);


outb = csvread('output_3dreal2patch.csv');


% Loading the Data and the Output Patches of the DDUL
% Synthetic or Field Example, zero for synthetic and one for field example.
[n1,n2,n3]=size(DataNoisey);
w1 =12;
w2 =12;
w3 =12;
s1z =4;
s2z =2;
s3z = 2;
real2depatch=yc_patch3d_inv(outb',1,n1,n2,n3,w1,w2,w3,s1z,s2z,s3z);

real2dataclean=Dc;
real2datan=DataNoisey;

figure;
imagesc(reshape(real2depatch,Nt,Nx*Ny))
colormap(seis(1))
caxis([-1,1])



fid = fopen('real2dataclean.bin','wb');
fwrite(fid,real2dataclean,'float');
fclose(fid);

fid = fopen('real2datan.bin','wb');
fwrite(fid,real2datan,'float');
fclose(fid);

fid = fopen('real2demssa.bin','wb');
fwrite(fid,real2demssa,'float');
fclose(fid);


fid = fopen('real2dewma.bin','wb');
fwrite(fid,real2dewma,'float');
fclose(fid);

fid = fopen('real2depatch.bin','wb');
fwrite(fid,real2depatch,'float');
fclose(fid);

 fid1 = fopen('real2_ls_mssa.bin','rb'); 
[A1] = fread(fid1,'float');
real2_ls_mssa=reshape(A1,181,61,31);

fid1 = fopen('real2_ls_svmf.bin','rb'); 
[A1] = fread(fid1,'float');
real2_ls_svmf=reshape(A1,181,61,31);

fid1 = fopen('real2_ls_patch.bin','rb'); 
[A1] = fread(fid1,'float');
real2_ls_patch=reshape(A1,181,61,31);

fid1 = fopen('real2_ls_wma.bin','rb'); 
[A1] = fread(fid1,'float');
real2_ls_wma=reshape(A1,181,61,31);
%% mean var
mean_mssa=mean(reshape(real2_ls_mssa,1,181*61*31));
mean_svmf=mean(reshape(real2_ls_svmf,1,181*61*31));
mean_patch=mean(reshape(real2_ls_patch,1,181*61*31));
mean_wma=mean(reshape(real2_ls_wma,1,181*61*31));

var_mssa=var(reshape(real2_ls_mssa,1,181*61*31));
var_svmf=var(reshape(real2_ls_svmf,1,181*61*31));
var_patch=var(reshape(real2_ls_patch,1,181*61*31));
var_wma=var(reshape(real2_ls_wma,1,181*61*31));
%% images

figure;yc_imagesc(squeeze(real2_ls_mssa(:,:,18)))
title('MSSA')
caxis([0,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:10:61]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:10:61],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:20:181]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:20*dt:181*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%grid on;% 显示网格
colorbar;
colormap('jet')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real2_ls_mssa.png')  


figure;yc_imagesc(squeeze(real2_ls_svmf(:,:,18)))
title('SVMF')
caxis([0,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:10:61]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:10:61],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:20:181]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:20*dt:181*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%grid on;% 显示网格
colorbar;
colormap('jet')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real2_ls_svmf.png')  

figure;yc_imagesc(squeeze(real2_ls_patch(:,:,18)))
title('PATCHUNET')
caxis([0,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:10:61]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:10:61],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:20:181]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:20*dt:181*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%grid on;% 显示网格
colorbar;
colormap('jet')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real2_ls_patch.png')  

figure;yc_imagesc(squeeze(real2_ls_wma(:,:,18)))
title('WMA')
caxis([0,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:10:61]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:10:61],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:20:181]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:20*dt:181*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%grid on;% 显示网格
colorbar;
colormap('jet')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real2_ls_wma.png')  