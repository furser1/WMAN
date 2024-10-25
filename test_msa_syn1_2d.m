%%%%%%%%%%%%%%%%%%%%test _ok%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
close all;
cd D://Desktop/MSAnet

[DataClean,h,t] = levents();
DataClean=yc_scale(DataClean);
figwigb(DataClean,0.9,'Clean Data')
figure;
plot(DataClean(:,30))
  dt = 2./1000;
  [Nt,Nx]=size(DataClean);
%% % 加载奇异单道噪音以及高斯随机噪音

rng('default')
rng(2023120);
DataNoisey=zeros(size(DataClean));
erra = zeros(size(DataClean));
mask=rand(1,Nx);
mask(logical(mask<0.9))=0;
mask(logical(mask>=0.9))=1;
for i=1:Nt
    erra(i,:)=0.4*randn(1,Nx).*mask;
end


DataNoisey = DataClean + erra;
nrra = 0.15*randn(Nt,Nx);
%nrra=bandpass_filter(nrra,dt,2,20,40,100);

DataNoisey=DataNoisey+nrra;
figimag(DataClean)
colormap(seis(1))
figimag(DataNoisey)
colormap(seis(1))
[p,q]=fftplot2d(DataNoisey,dt,1);
figure
plot(q(1:end/2),p(1:end/2))


figure;

wigbcyk(DataNoisey,0.9)

set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:30:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'FontName','Times New Roman','Linewidth',1.5,'Fontsize',15,'Fontweight','bold') ;
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);

img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn2derr.png')    
load Syn1datacl







figwigb(erra,0.4,'erra noise')



figimag(DataNoisey)
colormap(seis(1))
colorbar


 cd D:/Desktop/MSAnet
%% 数据patching
dn = DataNoisey;
w1 =48;
w2 =48;
s1z =1;
s2z =1;
dn_patch = yc_patch(dn,1,w1,w2,s1z,s2z);
% It is better to save .mat as -V 7.3 or later because the size of the generated patch is large.
save('Input_Patches_2Dsyn1','dn_patch','-V7.3');
filename = 'Input_Patches_2Dsyn1.csv';% 你的CSV数据文件名
writematrix(dn_patch, filename, 'Delimiter', ',');
%% 去噪数据结果显示
outa = csvread('output_2dsyn1.csv');
% Loading the Data and the Output Patches of the DDUL
% Synthetic or Field Example, zero for synthetic and one for field example.

% UnPatching
[n1,n2]=size(DataClean);
w1 =48;
w2 =48;
s1z =1;
s2z =1;
syn1dewma=yc_patch_inv(outa',1,n1,n2,w1,w2,s1z,s2z);

syn1errwma=DataNoisey-syn1dewma;

%% process with mssa \svmf \patch 
syn1demssa=fxmssa(DataNoisey,5,120,dt,3,0);

fid1 = fopen('syn1desvmf.bin','rb'); 
[A1] = fread(fid1,'float');
syn1desvmf=reshape(A1,501,120);

outc = csvread('output_2dsyn1patch.csv');
[n1,n2]=size(DataClean);
w1 =48;
w2 =48;
s1z =1;
s2z =1;
syn1depatch=yc_patch_inv(outc',1,n1,n2,w1,w2,s1z,s2z);

close all
%% save data
save Syn1dataclean.mat DataClean
save Syn1datan.mat DataNoisey
save Syn1demssa.mat syn1demssa
save Syn1dewma.mat syn1dewma
save syn1desvmf.mat syn1desvmf
save Syn1depatch.mat syn1depatch
save Syn1errwma.mat syn1errwma
%% fk 

[nt,nx]=size(syn1demssa);
nk=4*(2^nextpow2(nx));%隧道中道数少，这样可以增大数据量、使成图效果更细腻，16可调
nf=4*(2^nextpow2(nt));%隧道中道数少，这样可以增大数据量、使成图效果更细腻，16可调
f=(-nf/2+1:nf/2)/nf/dt;
k=(-nk/2+1:nk/2)/nk/dx;

S=fftshift(fft2(DataClean,nf,nk));%注意此处用了fftshift
Syn1fkclean=fliplr(S);  % 矩阵左右翻转

S=fftshift(fft2(syn1demssa,nf,nk));%注意此处用了fftshift
Syn1fkmssa=fliplr(S);  % 矩阵左右翻转

S=fftshift(fft2(syn1desvmf,nf,nk));%注意此处用了fftshift
Syn1fksvmf=fliplr(S);  % 矩阵左右翻转

S=fftshift(fft2(syn1depatch,nf,nk));%注意此处用了fftshift
Syn1fkpatch=fliplr(S);  % 矩阵左右翻转

S=fftshift(fft2(syn1dewma,nf,nk));%注意此处用了fftshift
Syn1fkwma=fliplr(S);  % 矩阵左右翻转

save Syn1fkclean.mat Syn1fkclean
save Syn1fkmssa.mat Syn1fkmssa
save Syn1fksvmf.mat Syn1fksvmf
save Syn1fkpatch.mat Syn1fkpatch
save Syn1fkwma.mat Syn1fkwma

figure;imagesc(k,f,abs(Syn1fkclean(:,:)))
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
print(img, '-dpng', '-r600', './Syn1fkclean.png')  


figure;imagesc(k,f,abs(Syn1fkmssa(:,:)))
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
print(img, '-dpng', '-r600', './Syn1fkmssa.png')  


figure;imagesc(k,f,abs(Syn1fksvmf(:,:)))
set(gcf,"Position",[150,120,450,620])
title('svmf')
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Wavenumber[c/m]','FontSize',15,'linewidth',1.5);
ylabel('f[Hz]','FontSize',15,'linewidth',1.5);
ylim([0 120])
xlim([-0.1 0.1])
%grid on;% 显示网格
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './syn1fksvmf.png') 


figure;imagesc(k,f,abs(Syn1fkpatch(:,:)))
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
print(img, '-dpng', '-r600', './syn1fkpatch.png') 

figure;imagesc(k,f,abs(Syn1fkwma(:,:)))
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
print(img, '-dpng', '-r600', './syn1fkwma.png') 

 figure;
   imagesc(k,f,abs(S(:,:))); 
%  imagesc(k,f(nf/2:nf),abs(S(nf/2:nf,:)));  % 这里为了只显示f>0的fk图像，因此f、S值只能取一半
   xlabel('Wavenumber[c/m]'); ylabel('f[Hz]');title('syn1fkmssa');
   set(gca,'YDir','normal'); % 让纵轴正值向上

figure;
imagesc(real(S))
%% SNR
SNR1d=snr(reshape(DataClean,1,501*120),reshape(DataNoisey-DataClean,1,501*120));
SNR1mssa=snr(reshape(syn1demssa,1,501*120),reshape(syn1demssa-DataClean,1,501*120));
SNR1svmf=snr(reshape(syn1desvmf,1,501*120),reshape(syn1desvmf-DataClean,1,501*120));
SNR1patch=snr(reshape(syn1depatch,1,501*120),reshape(syn1depatch-DataClean,1,501*120));
SNR1wma=snr(reshape(syn1dewma,1,501*120),reshape(syn1dewma-DataClean,1,501*120));


pSNR1d=psnr(reshape(DataClean,1,501*120),reshape(DataClean,1,501*120));
pSNR1mssa=psnr(reshape(syn1demssa,1,501*120),reshape(DataClean,1,501*120));
pSNR1svmf=psnr(reshape(syn1desvmf,1,501*120),reshape(DataClean,1,501*120));
pSNR1patch=psnr(reshape(syn1depatch,1,501*120),reshape(DataClean,1,501*120));
pSNR1wma=psnr(reshape(syn1dewma,1,501*120),reshape(DataClean,1,501*120));




%% 2D image %%

figure;yc_imagesc(DataClean)

title('Clean')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:30:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
%grid on;% 显示网格
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1dataclean.png')    


figure;yc_imagesc(DataNoisey-DataClean)
title('Noise')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:30:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1n.png')  

figure;yc_imagesc(DataNoisey)
title('Noisy')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:30:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1datan.png')    



figure;
imagesc(syn1dewma)
title('Proposed')
caxis([-0.9,0.9])
colormap(seis(1))
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:30:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1dewma.png')    


figure;
imagesc(syn1demssa)
title('MSSA')
caxis([-1,1])
dt=2/1000;
colormap(seis(1))
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:30:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1demssa.png')    


figure;
imagesc(syn1desvmf)
title('SVMF')
caxis([-1,1])
dt=2/1000;
colormap(seis(1))
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:30:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1desvmf.png')    

figure;
imagesc(syn1depatch)
title('PATCHUNET')
caxis([-1,1])
dt=2/1000;
colormap(seis(1))
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:30:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1depatch.png')    


syn1errwma=DataNoisey-syn1dewma;
figure;yc_imagesc(syn1errwma)
title('Proposed')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:30:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1derrwma.png')    


syn1errmssa=DataNoisey-syn1demssa;
figure;yc_imagesc(syn1errmssa)
title('MSSA')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:30:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1errmssa.png')  


syn1errpatch=DataNoisey-syn1depatch;
figure;yc_imagesc(syn1errpatch)
title('PATCHUNET')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:30:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1errpatch.png')  

syn1errsvmf=DataNoisey-syn1desvmf;
figure;yc_imagesc(syn1errsvmf)
title('SVMF')
caxis([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);%设置要多少个刻度要从1/0开始 横轴
set(gca,'XTickLabel',[0:30:120],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 横轴
set(gca,'YTick',[0:50:501]);%设置要多少个刻度要从1/0开始 纵轴
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1errsvmf.png')  



figure;
set(gcf,"Position",[150,120,650,320])
line55=[DataClean(:,55),syn1dewma(:,55)];
plot(DataClean(:,55),'LineWidth',1.5,'Color','b')
hold on;
plot(syn1demssa(:,55),'LineWidth',1.5,'Color','r')
xlim([0,501]);                                   
set(gca,'XTick',(0:100:501));%设置要多少个刻度要从1/0开始 纵轴
set(gca,'XTickLabel',(0:100*dt:501*dt),'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Time (s)','FontSize',15,'linewidth',1.5);
ylabel('Amplitude','FontSize',15,'linewidth',1.5);
legend('Clean','MSSA')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1plotmssa.png')  

figure;
set(gcf,"Position",[150,120,650,320])
plot(DataClean(:,55),'LineWidth',1.5,'Color','b')
hold on;
plot(syn1desvmf(:,55),'LineWidth',1.5,'Color','r')
xlim([0,501]);                                   
set(gca,'XTick',(0:100:501));%设置要多少个刻度要从1/0开始 纵轴
set(gca,'XTickLabel',(0:100*dt:501*dt),'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Time (s)','FontSize',15,'linewidth',1.5);
ylabel('Amplitude','FontSize',15,'linewidth',1.5);
legend('Clean','SVMF')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1plotsvmf.png')  

figure;
set(gcf,"Position",[150,120,650,320])
plot(DataClean(:,55),'LineWidth',1.5,'Color','b')
hold on;
plot(syn1depatch(:,55),'LineWidth',1.5,'Color','r')
xlim([0,501]);                                   
set(gca,'XTick',(0:100:501));%设置要多少个刻度要从1/0开始 纵轴
set(gca,'XTickLabel',(0:100*dt:501*dt),'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Time (s)','FontSize',15,'linewidth',1.5);
ylabel('Amplitude','FontSize',15,'linewidth',1.5);
legend('Clean','PATCHUNET')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1plotpatch.png')  

figure;
set(gcf,"Position",[150,120,650,320])
plot(DataClean(:,55),'LineWidth',1.5,'Color','b')
hold on;
plot(syn1dewma(:,55),'LineWidth',1.5,'Color','r')
xlim([0,501]);                                   
set(gca,'XTick',(0:100:501));%设置要多少个刻度要从1/0开始 纵轴
set(gca,'XTickLabel',(0:100*dt:501*dt),'FontSize',12);%给坐标刻度加上内容，这个值随意更改 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Time (s)','FontSize',15,'linewidth',1.5);
ylabel('Amplitude','FontSize',15,'linewidth',1.5);
legend('Clean','Proposed')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1plotwma.png')  


c1 = csvread('loss_2dsyn1.csv',1,0);
figure;
set(gcf,"Position",[150,120,650,320])
plot(c1(:,1),'LineWidth',1.5,'Color','b')
hold on;
plot(c1(:,2),'LineWidth',1.5,'Color','r')

xlim([0,48]);
%ylim([-1,1])
set(gca,'XTick',(0:5:100));%设置要多少个刻度要从1/0开始 纵轴
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Epoch','FontSize',15,'linewidth',1.5);
ylabel('Loss','FontSize',15,'linewidth',1.5);
legend('Training data','Validation data')
img =gcf;  %获取当前画图的句柄
print(img, '-dpng', '-r600', './Syn1loss.png')  
save Syn1loss.mat c1
