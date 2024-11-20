%% Unpatching csv data
outa = csvread('output_2dsyn1.csv');
n1=501;
n2=120;

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

%% fk transform 
[nt,nx]=size(syn1demssa);
nk=4*(2^nextpow2(nx));
nf=4*(2^nextpow2(nt));
f=(-nf/2+1:nf/2)/nf/dt;
k=(-nk/2+1:nk/2)/nk/dx;

S=fftshift(fft2(DataClean,nf,nk));
Syn1fkclean=fliplr(S);  
S=fftshift(fft2(syn1demssa,nf,nk));
Syn1fkmssa=fliplr(S);  

S=fftshift(fft2(syn1desvmf,nf,nk));
Syn1fksvmf=fliplr(S);  

S=fftshift(fft2(syn1depatch,nf,nk));
Syn1fkpatch=fliplr(S);  

S=fftshift(fft2(syn1dewma,nf,nk));
Syn1fkwma=fliplr(S);  


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




%% 2D image data and F-K %%

figure;yc_imagesc(DataClean)
title('Clean')
clim([-1,1])
dt=2/1000;
set(gcf,"Position",[150,120,450,620])
set(gca,'XTick',[0:30:120]);
set(gca,'XTickLabel',[0:30:120],'FontSize',12);
set(gca,'YTick',[0:50:501]);
set(gca,'YTickLabel',[0:50*dt:501*dt],'FontSize',12);
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Tracenumber','FontSize',15,'linewidth',1.5);
ylabel('Time (s)','FontSize',15,'linewidth',1.5);
colorbar;
img =gcf;  
print(img, '-dpng', '-r600', './Syn1dataclean.png')    


figure;imagesc(k,f,abs(Syn1fkwma(:,:)))
set(gcf,"Position",[150,120,450,620])
title('Porposed')
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Wavenumber[c/m]','FontSize',15,'linewidth',1.5);
ylabel('f[Hz]','FontSize',15,'linewidth',1.5);
ylim([0 120])
xlim([-0.1 0.1])
colorbar;
img =gcf;  
print(img, '-dpng', '-r600', './syn1fkwma.png') 
