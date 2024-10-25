%%%%%%%%%%%%%%%%%%%%test _ok%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
clear;
close all;

cd D:\Desktop\test_skipNet\DDUL-main\Datasets
    load 3Dsyn_example
cd D:\Desktop\MSAnet

% preparing the patches where w1, w2, and w3 are the patch size, while the s1z, s2z, and s3z are the number of shift samples between neighbor windows. 
% The default values are w1,w2,w3 =15,15,15 and s1z,s2z,s3z=1,1,1.
rng('default')
rng(202305);
[Nt,Nx,Ny] = size(DataClean);
Dc=reshape(DataClean,Nt,Nx*Ny);
Input = zeros(size(Dc));
mask=rand(1,Nx*Ny);
mask(logical(mask<0.95))=0;
mask(logical(mask>=0.95))=1;
Input0 = zeros(size(Input));
for i=1:Nt
    Input0(i,:)=Dc(i,:)+randn(1,Nx*Ny).*mask*0.8;
end

Input=Input0+0.2*randn(Nt,Nx*Ny);


figure;
imagesc(reshape(DataClean,Nt,Nx*Ny))
caxis([-1,1])
colormap(seis(1))

figure;
imagesc(reshape(Input,Nt,Nx*Ny))
caxis([-1,1])
colormap(seis(1))
DataNoisey=reshape(Input,Nt,Nx,Ny);

Syn2dataclean=DataClean;
Syn2datan=DataNoisey;
dn = DataNoisey;
w1 =12;
w2 =12;
w3 =12;
s1z =1;
s2z =1;
s3z =1;
dn_patch = yc_patch3d(dn,1,w1,w2,w3,s1z,s2z,s3z);
% It is better to save .mat as -V 7.3 or later because the size of the generated patch is large.
save('Input_Patches_3Dsyn1','dn_patch','-V7.3');
filename = 'Input_Patches_3Dsyn2.csv';% 你的CSV数据文件名
writematrix(dn_patch, filename, 'Delimiter', ',');

%% WMANet

outa = csvread('output_3dsyn2.csv');

% Loading the Data and the Output Patches of the DDUL
% Synthetic or Field Example, zero for synthetic and one for field example.

% UnPatching
[n1,n2,n3] = size(DataClean);
n1=126;
n2=32;
n3=32;
[n1,n2,n3]=size(DataNoisy);
w1 =12;
w2 =12;
w3 =12;
s1z =1;
s2z =1;
s3z = 1;
Syn2dewma=yc_patch3d_inv(outa',1,n1,n2,n3,w1,w2,w3,s1z,s2z,s3z);
figure;yc_imagesc(squeeze(DataNoisy(:,14,:)));
figure;imagesc(Input)
title('Denoised Signal')
%% fmssa \ PatchNet
Syn2demssa=fxymssa(DataNoisey,5,120,0.002,3,0);

outb = csvread('output_3dsyn2patch.csv');


% Loading the Data and the Output Patches of the DDUL
% Synthetic or Field Example, zero for synthetic and one for field example.

% UnPatching
[n1,n2,n3] = size(DataClean);
n1=126;
n2=32;
n3=32;
[n1,n2,n3]=size(DataNoisy);
w1 =12;
w2 =12;
w3 =12;
s1z =1;
s2z =1;
s3z = 1;
Syn2depatch=yc_patch3d_inv(outb',1,n1,n2,n3,w1,w2,w3,s1z,s2z,s3z);
figure;yc_imagesc(squeeze(Syn2depatch(:,18,:)));

fid1 = fopen('syn2desvmf.bin','rb'); 
[A1] = fread(fid1,'float');
syn2desvmf=reshape(A1,126,32,32);

fid = fopen('Syn2dataclean.bin','wb');
fwrite(fid,Syn2dataclean,'float');
fclose(fid);

fid = fopen('Syn2datan.bin','wb');
fwrite(fid,Syn2datan,'float');
fclose(fid);

fid = fopen('Syn2demssa.bin','wb');
fwrite(fid,Syn2demssa,'float');
fclose(fid);


fid = fopen('Syn2dewma.bin','wb');
fwrite(fid,Syn2dewma,'float');
fclose(fid);

fid = fopen('Syn2depatch.bin','wb');
fwrite(fid,Syn2depatch,'float');
fclose(fid);
%%
fid1 = fopen('Syn2dataclean.bin','rb'); 
[A1] = fread(fid1,'float');
Syn2dataclean=reshape(A1,126,32,32);

fid1 = fopen('Syn2datan.bin','rb'); 
[A1] = fread(fid1,'float');
Syn2datan=reshape(A1,126,32,32);


fid1 = fopen('Syn2demssa.bin','rb'); 
[A1] = fread(fid1,'float');
Syn2demssa=reshape(A1,126,32,32);

fid1 = fopen('Syn2desvmf.bin','rb'); 
[A1] = fread(fid1,'float');
Syn2desvmf=reshape(A1,126,32,32);

fid1 = fopen('Syn2depatch.bin','rb'); 
[A1] = fread(fid1,'float');
Syn2depatch=reshape(A1,126,32,32);

fid1 = fopen('Syn2dewma.bin','rb'); 
[A1] = fread(fid1,'float');
Syn2dewma=reshape(A1,126,32,32);


save Syn2dataclean.mat Syn2dataclean
save Syn2datan.mat Syn2datan

save Syn2demssa.mat Syn2demssa
save Syn2desvmf.mat Syn2desvmf
save Syn2depatch.mat Syn2depatch
save Syn2dewma.mat Syn2dewma
%% SNR
SNR2d=snr(reshape(DataClean,1,126*32*32),reshape(DataNoisey-DataClean,1,126*32*32));
SNR2mssa=snr(reshape(syn2demssa,1,126*32*32),reshape(syn2demssa-DataClean,1,126*32*32));
SNR2svmf=snr(reshape(syn2desvmf,1,126*32*32),reshape(syn2desvmf-DataClean,1,126*32*32));
SNR2patch=snr(reshape(Syn2depatch,1,126*32*32),reshape(Syn2depatch-DataClean,1,126*32*32));
SNR2wma=snr(reshape(Syn2dewma,1,126*32*32),reshape(Syn2dewma-DataClean,1,126*32*32));

pSNR2mssa=psnr(reshape(Syn2demssa,1,126*32*32),reshape(Syn2dataclean,1,126*32*32));
pSNR2svmf=psnr(reshape(Syn2desvmf,1,126*32*32),reshape(Syn2dataclean,1,126*32*32));
pSNR2patch=psnr(reshape(Syn2depatch,1,126*32*32),reshape(Syn2dataclean,1,126*32*32));
pSNR2wma=psnr(reshape(Syn2dewma,1,126*32*32),reshape(Syn2dataclean,1,126*32*32));




%% 2D image %%



dx=5;
[nt,nx]=size(squeeze(syn2demssa(:,14,:)));
nk=4*(2^nextpow2(nx));%隧道中道数少，这样可以增大数据量、使成图效果更细腻，16可调
nf=4*(2^nextpow2(nt));%隧道中道数少，这样可以增大数据量、使成图效果更细腻，16可调
f=(-nf/2+1:nf/2)/nf/dt;
k=(-nk/2+1:nk/2)/nk/dx;

S=fftshift(fft2(squeeze(Syn2demssa(:,14,:)),nf,nk));%注意此处用了fftshift
Syn2fkclean=fliplr(S);  % 矩阵左右翻转

S=fftshift(fft2(squeeze(Syn2demssa(:,14,:)),nf,nk));%注意此处用了fftshift
Syn2fkmssa=fliplr(S);  % 矩阵左右翻转

S=fftshift(fft2(squeeze(Syn2desvmf(:,14,:)),nf,nk));%注意此处用了fftshift
Syn2fksvmf=fliplr(S);  % 矩阵左右翻转

S=fftshift(fft2(squeeze(Syn2depatch(:,14,:)),nf,nk));%注意此处用了fftshift
Syn2fkpatch=fliplr(S);  % 矩阵左右翻转

S=fftshift(fft2(squeeze(Syn2dewma(:,14,:)),nf,nk));%注意此处用了fftshift
Syn2fkwma=fliplr(S);  % 矩阵左右翻转




figure;imagesc(k,f,abs(Syn2fkclean(:,:)))
set(gcf,"Position",[150,120,450,620])
title('Clean')
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Wavenumber[c/m]','FontSize',15,'linewidth',1.5);
ylabel('f[Hz]','FontSize',15,'linewidth',1.5);
ylim([0 120])
xlim([-0.1 0.1])
%grid on;% 显示网格
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2fkclean.png')  




figure;imagesc(k,f,abs(Syn2fkmssa(:,:)))
set(gcf,"Position",[150,120,450,620])
title('MSSA')
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Wavenumber[c/m]','FontSize',15,'linewidth',1.5);
ylabel('f[Hz]','FontSize',15,'linewidth',1.5);
ylim([0 120])
xlim([-0.1 0.1])
%grid on;% 显示网格
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2fkmssa.png')  


figure;imagesc(k,f,abs(Syn2fksvmf(:,:)))
set(gcf,"Position",[150,120,450,620])
title('SVMF')
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Wavenumber[c/m]','FontSize',15,'linewidth',1.5);
ylabel('f[Hz]','FontSize',15,'linewidth',1.5);
ylim([0 120])
xlim([-0.1 0.1])
%grid on;% 显示网格
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './syn2fksvmf.png') 


figure;imagesc(k,f,abs(Syn2fkpatch(:,:)))
set(gcf,"Position",[150,120,450,620])
title('PATCHUNET')
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Wavenumber[c/m]','FontSize',15,'linewidth',1.5);
ylabel('f[Hz]','FontSize',15,'linewidth',1.5);
ylim([0 120])
xlim([-0.1 0.1])
%grid on;% 显示网格
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './syn2fkpatch.png') 

figure;imagesc(k,f,abs(Syn2fkwma(:,:)))
set(gcf,"Position",[150,120,450,620])
title('Porposed')
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Wavenumber[c/m]','FontSize',15,'linewidth',1.5);
ylabel('f[Hz]','FontSize',15,'linewidth',1.5);
ylim([0 120])
xlim([-0.1 0.1])
%grid on;% 显示网格
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './syn2fkwma.png') 




figure;
set(gcf,"Position",[150,120,650,320])
plot(DataClean(:,20,15),'LineWidth',1.5,'Color','b')
hold on;
plot(Syn2demssa(:,20,15),'LineWidth',1.5,'Color','r')
xlim([0,126]);                                   
set(gca,'XTick',(0:20:126));%设置要多少个刻度要从1/0开始 纵轴
set(gca,'XTickLabel',(0:20*dt:126*dt),'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Time (s)','FontSize',15,'linewidth',1.5);
ylabel('Amplitude','FontSize',15,'linewidth',1.5);
legend('Clean','MSSA')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2plotmssa.png')  

figure;
set(gcf,"Position",[150,120,650,320])
plot(DataClean(:,20,15),'LineWidth',1.5,'Color','b')
hold on;
plot(Syn2desvmf(:,20,15),'LineWidth',1.5,'Color','r')
xlim([0,126]);                                   
set(gca,'XTick',(0:20:126));%设置要多少个刻度要从1/0开始 纵轴
set(gca,'XTickLabel',(0:20*dt:126*dt),'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Time (s)','FontSize',15,'linewidth',1.5);
ylabel('Amplitude','FontSize',15,'linewidth',1.5);
legend('Clean','SVMF')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2plotsvmf.png')  

figure;
set(gcf,"Position",[150,120,650,320])
plot(DataClean(:,20,15),'LineWidth',1.5,'Color','b')
hold on;
plot(Syn2depatch(:,20,15),'LineWidth',1.5,'Color','r')
xlim([0,126]);                                   
set(gca,'XTick',(0:20:126));%设置要多少个刻度要从1/0开始 纵轴
set(gca,'XTickLabel',(0:20*dt:126*dt),'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Time (s)','FontSize',15,'linewidth',1.5);
ylabel('Amplitude','FontSize',15,'linewidth',1.5);
legend('Clean','PATCHUNET')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2plotpatch.png')  

figure;
set(gcf,"Position",[150,120,650,320])
plot(DataClean(:,20,15),'LineWidth',1.5,'Color','b')
hold on;
plot(Syn2dewma(:,20,15),'LineWidth',1.5,'Color','r')
xlim([0,126]);                                   
set(gca,'XTick',(0:20:126));%设置要多少个刻度要从1/0开始 纵轴
set(gca,'XTickLabel',(0:20*dt:126*dt),'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴

set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Time (s)','FontSize',15,'linewidth',1.5);
ylabel('Amplitude','FontSize',15,'linewidth',1.5);
legend('Clean','Proposed')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2plotwma.png')  
























figure('units','normalized','Position',[0.4 0.1 0.4, 0.3]);
plot([1:126]*dt,DataClean(:,20,15),'color','b','linewidth',3);
hold on;
plot([1:126]*dt,DataNoisey(:,20,15),'color','g','linewidth',3);
hold on;
plot([1:126]*dt,syn2dewma(:,20,15),'color','r','linewidth',3);
legend('Clean data','Noisy data','Denoised data');
xlim([1 126]*dt)
set(gca,'FontName','Times New Roman','FontSize',14);
title('Inline = 18','FontName','Times New Roman','Fontsize',16,'Fontweight','Bold');
xlabel('{{Times} (s)}','FontName','Times New Roman','FontSize',14);
ylabel('{{{Amplitude}}}','FontName','Times New Roman','FontSize',14);
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2Amp.png')    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;yc_imagesc(squeeze(DataClean(:,18,:)))
title('Clean data')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:10:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:10:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%grid on;% 显示网格
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2dataclean.png')  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;yc_imagesc(squeeze(DataNoisey(:,18,:)))
title('Noisy Data')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:10:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:10:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%grid on;% 显示网格
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2datan.png')    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[q1,f1]=fftplot2d(squeeze(DataClean(:,18,:)),dt,1);
[q2,f2]=fftplot2d(squeeze(out1(:,18,:)),dt,1);
maxf=200;
nfft=(2^nextpow2(n1));
df=1/dt/nfft;

figure('units','normalized','Position',[0.4 0.1 0.4, 0.3]);
nn1=floor(maxf/df);
plot(f1(1:nn1,1),q1(1:nn1,1),'color','b','linewidth',3);
hold on;
plot(f2(1:nn1,1),q2(1:nn1,1),'color','r','linewidth',3);
legend('Clean data','Denoised data');
xlim([0 125])
set(gca,'FontName','Times New Roman','FontSize',14);
title('Spectrum','FontName','Times New Roman','Fontsize',16,'Fontweight','Bold');
xlabel('{{{Frequence}} (Hz)}','FontName','Times New Roman','FontSize',14);
ylabel('{{{Amplitude}}}','FontName','Times New Roman','FontSize',14);
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2Amps.png')    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;yc_imagesc(squeeze(out1(:,18,:)));
title('Denoised Signal')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:20:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:20:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%grid on;% 显示网格
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2ded.png')    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;yc_imagesc(squeeze(DataNoisey(:,18,:)-out1(:,18,:)))
title('Removed Noise Section')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:10:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:10:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%grid on;% 显示网格
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2derr.png')    


% 
% 
% 
c1 = csvread('loss_3dsyn2.csv',1,0);
figure;
set(gcf,"Position",[150,120,650,320])
plot(c1(:,1),'LineWidth',1.5,'color','r')
hold on;
plot(c1(:,2),'LineWidth',1.5,'color','b')
xlim([0,81]);
%ylim([0,0.15])
set(gca,'XTick',(0:20:100));%设置要多少个刻度要从1/0开始 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Epoch','FontSize',15,'linewidth',1.5);
ylabel('Loss','FontSize',15,'linewidth',1.5);
legend('Training data','Validation data')
 saveas(gcf,'loss3dsyn2.png')
