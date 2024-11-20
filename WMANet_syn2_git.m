%% load data 
load syn3d.mat
%% Patching data
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


%% Unpatching csv data

outa = csvread('output_3dsyn2.csv');

n1=126;
n2=32;
n3=32;
w1 =12;
w2 =12;
w3 =12;
s1z =1;
s2z =1;
s3z = 1;
Syn2dewma=yc_patch3d_inv(outa',1,n1,n2,n3,w1,w2,w3,s1z,s2z,s3z);

%% fmssa \ PatchNet
Syn2demssa=fxymssa(DataNoisey,5,120,0.002,3,0);

outb = csvread('output_3dsyn2patch.csv');


% Loading the Data and the Output Patches of the DDUL
% Synthetic or Field Example, zero for synthetic and one for field example.

% UnPatching
n1=126;
n2=32;
n3=32;
w1 =12;
w2 =12;
w3 =12;
s1z =1;
s2z =1;
s3z = 1;
Syn2depatch=yc_patch3d_inv(outb',1,n1,n2,n3,w1,w2,w3,s1z,s2z,s3z);
%% fk transform



dx=5;
[nt,nx]=size(squeeze(syn2demssa(:,14,:)));
nk=4*(2^nextpow2(nx));
nf=4*(2^nextpow2(nt));
f=(-nf/2+1:nf/2)/nf/dt;
k=(-nk/2+1:nk/2)/nk/dx;

S=fftshift(fft2(squeeze(Syn2demssa(:,14,:)),nf,nk));
Syn2fkclean=fliplr(S);  

S=fftshift(fft2(squeeze(Syn2demssa(:,14,:)),nf,nk));
Syn2fkmssa=fliplr(S);  

S=fftshift(fft2(squeeze(Syn2desvmf(:,14,:)),nf,nk));
Syn2fksvmf=fliplr(S);  

S=fftshift(fft2(squeeze(Syn2depatch(:,14,:)),nf,nk));
Syn2fkpatch=fliplr(S); 

S=fftshift(fft2(squeeze(Syn2dewma(:,14,:)),nf,nk));
Syn2fkwma=fliplr(S); 
%% 2D image data and F-K %%


figure;
set(gcf,"Position",[150,120,650,320])
plot(DataClean(:,20,15),'LineWidth',1.5,'Color','b')
hold on;
plot(Syn2depatch(:,20,15),'LineWidth',1.5,'Color','r')
xlim([0,126]);                                   
set(gca,'XTick',(0:20:126));
set(gca,'XTickLabel',(0:20*dt:126*dt),'FontSize',12);
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Time (s)','FontSize',15,'linewidth',1.5);
ylabel('Amplitude','FontSize',15,'linewidth',1.5);
legend('Clean','PATCHUNET')
img =gcf;  
print(img, '-dpng', '-r600', './Syn2plotpatch.png')  



figure;imagesc(k,f,abs(Syn2fkclean(:,:)))
set(gcf,"Position",[150,120,450,620])
title('Clean')
set(gca,'Linewidth',1.5,'Fontsize',15,'Fontweight','bold');
set(gcf,'Color','w')
xlabel('Wavenumber[c/m]','FontSize',15,'linewidth',1.5);
ylabel('f[Hz]','FontSize',15,'linewidth',1.5);
ylim([0 120])
xlim([-0.1 0.1])

colorbar;
img =gcf;  
print(img, '-dpng', '-r600', './Syn2fkclean.png')  


