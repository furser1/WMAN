clc;
clear;

cd D:\Desktop\MSAnet\realdata

[c1]=textread('segy2d.txt','%f','headerlines',1);

c2=reshape(c1,1501,188,345);
c11=reshape(c1,1501,188*345);
figure;imagesc(c11(:,188*200:188*220))
sl=115;
%dat = squeeze(c2(225:225+511,91:180,sl));
dat = squeeze(c2(225:225+511,sl,50:310));
dt=0.002;
figure;
imagesc(dat)
colormap(seis(1))
dat1=bandpass_filter(fliplr(dat),dt,1,50,85,110);

dat2=scalecyk(dat1,1);
figure;
imagesc(dat2)
colormap(seis(1))
dt=0.002;
rng('default')
rng(20224);

[Nt,Nx] = size(dat1);
mask=rand(1,Nx);
mask(logical(mask<0.95))=0;
mask(logical(mask>=0.95))=1;
Dn = zeros(size(dat));
for i=1:Nt
    Dn(i,:)=dat2(i,:)+randn(1,Nx).*mask*1;
end

DataNoisy=Dn+0.2*randn(Nt,Nx);
figure;
imagesc(DataNoisy)
colormap(seis(2))
caxis([-1,1])
%%%%%%%
cd ..
fid1 = fopen('real1noisy.bin','rb'); 
[A1] = fread(fid1,'float');
real1noisy=reshape(A1,512,261);


%% 数据patching
w1 =32;
w2 =32;
s1z =2;
s2z =2;
dn_patch = yc_patch(DataNoisy,1,w1,w2,s1z,s2z);
% It is better to save .mat as -V 7.3 or later because the size of the generated patch is large.
save('Input_Patches_2Dreal1','dn_patch','-V7.3');
filename = 'Input_Patches_2Dreal1.csv';% 你的CSV数据文件名
writematrix(dn_patch, filename, 'Delimiter', ',');
%% 去噪数据结果显示

fid1 = fopen('orreal_wma.bin','rb'); 
[A1] = fread(fid1,'float');
real1dewma=reshape(A1,512,261);

fid1 = fopen('orereal_wma.bin','rb'); 
[A1] = fread(fid1,'float');
real1errwma=reshape(A1,512,261);

%% fmssa \ PatchNet
real1demssa=fxymssa(DataNoisy,5,120,0.002,3,0);


fid1 = fopen('real1desvmf.bin','rb'); 
[A1] = fread(fid1,'float');
real1desvmf=reshape(A1,512,261);



outc = csvread('output_2dreal1patch.csv');


% Loading the Data and the Output Patches of the DDUL
% Synthetic or Field Example, zero for synthetic and one for field example.

% UnPatching
[n1,n2]=size(Dn);
w1 =32;
w2 =32;
s1z =2;
s2z =2;
real1depatch=yc_patch_inv(outc',1,n1,n2,w1,w2,s1z,s2z);
%% save data

save real1noisy.mat DataNoisy

save real1demssa.mat real1demssa
save real1depatch.mat real1depatch
save real1dewma.mat real1dewma



fid = fopen('real1noisy.bin','wb');
fwrite(fid,DataNoisy,'float');
fclose(fid);


fid = fopen('real1dewma.bin','wb');
fwrite(fid,real1dewma,'float');
fclose(fid);

%% image


figure;yc_imagesc(DataNoisy)
title('Noisy')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:50:260]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:50:260],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:513]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:513*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real1datan.png')    



figure;
imagesc(real1dewma)
title('Proposed')
caxis([-1,1])
colormap(seis(1))
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:50:260]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:50:260],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:513]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:513*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real1dewma.png')    


figure;
imagesc(real1demssa)
title('MSSA')
caxis([-1,1])
dt=2/1000;
colormap(seis(1))
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:50:260]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:50:260],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:513]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:513*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real1demssa.png')    


figure;
imagesc(real1desvmf)
title('SVMF')
caxis([-1,1])
dt=2/1000;
colormap(seis(1))
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:50:260]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:50:260],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:513]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:513*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real1desvmf.png')    

figure;
imagesc(real1depatch)
title('PATCHUNET')
caxis([-1,1])
dt=2/1000;
colormap(seis(1))
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:50:260]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:50:260],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:513]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:513*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real1depatch.png')    


real1errwma=DataNoisy-real1dewma;
figure;yc_imagesc(real1errwma)
title('Proposed')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:50:260]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:50:260],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:513]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:513*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real1derrwma.png')    


real1errmssa=DataNoisy-real1demssa;
figure;yc_imagesc(real1errmssa)
title('MSSA')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:50:260]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:50:260],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:513]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:513*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real1errmssa.png')  


real1errpatch=DataNoisy-real1depatch;
figure;yc_imagesc(real1errpatch)
title('PATCHUNET')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:50:260]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:50:260],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:513]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:513*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real1errpatch.png')  

real1errsvmf=DataNoisy-real1desvmf;
figure;yc_imagesc(real1errsvmf)
title('SVMF')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:50:260]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:50:260],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:513]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:513*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real1errsvmf.png')  

reeal1zmsvmf=real1desvmf(170:270,50:126);
reeal1zmmssa=real1demssa(170:270,50:126);
reeal1zmpatch=real1depatch(170:270,50:126);
reeal1zmwma=real1dewma(170:270,50:126);


figure;yc_imagesc(reeal1zmsvmf)
title('SVMF')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,300,400])
set(gca,'XTick',[0:20:77]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[50:20:126],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:25:100]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[170*dt:25*dt:270*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real1zmsvmf.png')  


figure;yc_imagesc(reeal1zmmssa)
title('MSSA')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,300,400])
set(gca,'XTick',[0:20:77]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[50:20:126],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:25:100]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[170*dt:25*dt:270*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './real1zmmssa.png')  

figure;yc_imagesc(reeal1zmpatch)
title('PATCHUNET')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,300,400])
set(gca,'XTick',[0:20:77]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[50:20:126],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:25:100]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[170*dt:25*dt:270*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './reeal1zmpatch.png') 

figure;yc_imagesc(reeal1zmwma)
title('Proposed')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,300,400])
set(gca,'XTick',[0:20:77]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[50:20:126],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:25:100]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[170*dt:25*dt:270*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './reeal1zmwma.png') 



















