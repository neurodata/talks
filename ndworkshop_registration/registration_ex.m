% for a registration example
% we need a few different brain modalities
% in the middle MNI space
clear all;
close all;

addpath /cis/home/dtward/Functions/avwQuiet
rng(3);
colors = rand(256,3);
colors(1,:) = 0;


atlas = '/cis/home/dtward/Documents/mricloud_atlases/DTIatlas/Adult/labels168/s02a_168labels.img';
atlas = '/cis/home/dtward/Documents/mricloud_atlases/Adult27-55/Adt27-55_01_Adt27-55_01_FullLabels.img';
avw = avw_img_read(atlas);
I = avw.img;
I = padarray(I,([217 217 217] - size(I))/2);


Ishow = flip(squeeze(I(:,round(size(I,2)/2),:))',1);
Im = mod(Ishow,256)+1;
IRGB = reshape( [colors(Im,1),colors(Im,2),colors(Im,3)] , [size(Im,1),size(Im,2),3]);
figure('position',[1 1 512 512],'paperpositionmode','auto','inverthardcopy','off');
axes('position',[0,0,1,1]);
axis off;
imagesc(IRGB);
axis image
fontsize = 30;
axes('position',[0,0,1,1],'visible','off')
text(0.5,0.9,'MNI space','color','w','fontsize',fontsize,'horizontalalignment','center')

saveas(gcf,'reg_mni.png')
% xI = 1 : size(I,2);
% yI = 1 : size(I,1);
% zI = 1 : size(I,3);
% xI = xI - mean(xI);
% yI = yI - mean(yI);
% zI = zI - mean(zI);
xI = linspace(-1,1,size(I,2));
yI = linspace(-1,1,size(I,1));
zI = linspace(-1,1,size(I,3));
[XI,YI,ZI] = meshgrid(xI,yI,zI);

%%
% load some other data
image1 = '/cis/home/dtward/Documents/mricloud_atlases/DTIatlas/Adult/b0/s02a_b0.img';
avw = avw_img_read(image1);
I1 = avw.img;
x = linspace(-1,1,size(I1,2));
y = linspace(-1,1,size(I1,1));
z = linspace(-1,1,size(I1,3));
[X,Y,Z] = meshgrid(x,y,z);

% a good transform to bring it to  MNI
A = diag([0.9,0.85,0.9,1]);
B = inv(A);
Xs = B(1,1)*XI + B(1,2)*YI + B(1,3)*ZI + B(1,4);
Ys = B(2,1)*XI + B(2,2)*YI + B(2,3)*ZI + B(2,4);
Zs = B(3,1)*XI + B(3,2)*YI + B(3,3)*ZI + B(3,4);

F = griddedInterpolant({y,x,z},I1,'linear','nearest');
I1D = F(Ys,Xs,Zs);

Ishow = flip(squeeze(I1D(:,round(size(I,2)/2),:))',1);

figure('position',[1 1 512 512],'paperpositionmode','auto','inverthardcopy','off');
colormap gray
axes('position',[0,0,1,1]);
axis off;
imagesc(Ishow);
axis image
axes('position',[0,0,1,1],'visible','off')
text(0.5,0.9,'Registered 1','color','w','fontsize',fontsize,'horizontalalignment','center')

saveas(gcf,'reg_1b.png')


% a bad transform that brings it away
linear = 0.05;
quadratic = 0.05;
Xs = XI + linear*(randn*Xs + randn*Ys + randn*Zs) + quadratic*(randn*Xs.^2 + randn.*Xs.*Ys + randn.*Xs.*Zs + randn.*Ys.^2 + randn.*Ys.*Zs + randn.*Zs.^2);
Ys = YI + linear*(randn*Xs + randn*Ys + randn*Zs) + quadratic*(randn*Xs.^2 + randn.*Xs.*Ys + randn.*Xs.*Zs + randn.*Ys.^2 + randn.*Ys.*Zs + randn.*Zs.^2);
Zs = ZI + linear*(randn*Xs + randn*Ys + randn*Zs) + quadratic*(randn*Xs.^2 + randn.*Xs.*Ys + randn.*Xs.*Zs + randn.*Ys.^2 + randn.*Ys.*Zs + randn.*Zs.^2);

F = griddedInterpolant({y,x,z},I1,'linear','nearest');
I1D = F(Ys,Xs,Zs);

Ishow = flip(squeeze(I1D(:,round(size(I,2)/2),:))',1);

figure('position',[1 1 512 512],'paperpositionmode','auto','inverthardcopy','off');
colormap gray
axes('position',[0,0,1,1]);
axis off;
imagesc(Ishow);
axis image
axes('position',[0,0,1,1],'visible','off')
text(0.5,0.9,'Image 1 (b0)','color','w','fontsize',fontsize,'horizontalalignment','center')

saveas(gcf,'reg_1a.png')



%%

image1 = '/cis/home/dtward/Documents/mricloud_atlases/DTIatlas/Adult/dwi/s02a_dwi.img';
avw = avw_img_read(image1);
I1 = avw.img;
x = linspace(-1,1,size(I1,2));
y = linspace(-1,1,size(I1,1));
z = linspace(-1,1,size(I1,3));
[X,Y,Z] = meshgrid(x,y,z);

% a good transform to bring it to  MNI
A = diag([0.9,0.85,0.9,1]);
B = inv(A);
Xs = B(1,1)*XI + B(1,2)*YI + B(1,3)*ZI + B(1,4);
Ys = B(2,1)*XI + B(2,2)*YI + B(2,3)*ZI + B(2,4);
Zs = B(3,1)*XI + B(3,2)*YI + B(3,3)*ZI + B(3,4);

F = griddedInterpolant({y,x,z},I1,'linear','nearest');
I1D = F(Ys,Xs,Zs);

Ishow = flip(squeeze(I1D(:,round(size(I,2)/2),:))',1);

figure('position',[1 1 512 512],'paperpositionmode','auto','inverthardcopy','off');
colormap gray
axes('position',[0,0,1,1]);
axis off;
imagesc(Ishow);
axis image
axes('position',[0,0,1,1],'visible','off')
text(0.5,0.9,'Registered 2','color','w','fontsize',fontsize,'horizontalalignment','center')

saveas(gcf,'reg_2b.png')


% a bad transform that brings it away
Xs = XI + linear*(randn*Xs + randn*Ys + randn*Zs) + quadratic*(randn*Xs.^2 + randn.*Xs.*Ys + randn.*Xs.*Zs + randn.*Ys.^2 + randn.*Ys.*Zs + randn.*Zs.^2);
Ys = YI + linear*(randn*Xs + randn*Ys + randn*Zs) + quadratic*(randn*Xs.^2 + randn.*Xs.*Ys + randn.*Xs.*Zs + randn.*Ys.^2 + randn.*Ys.*Zs + randn.*Zs.^2);
Zs = ZI + linear*(randn*Xs + randn*Ys + randn*Zs) + quadratic*(randn*Xs.^2 + randn.*Xs.*Ys + randn.*Xs.*Zs + randn.*Ys.^2 + randn.*Ys.*Zs + randn.*Zs.^2);

F = griddedInterpolant({y,x,z},I1,'linear','nearest');
I1D = F(Ys,Xs,Zs);

Ishow = flip(squeeze(I1D(:,round(size(I,2)/2),:))',1);

figure('position',[1 1 512 512],'paperpositionmode','auto','inverthardcopy','off');
colormap gray
axes('position',[0,0,1,1]);
axis off;
imagesc(Ishow);
axis image
axes('position',[0,0,1,1],'visible','off')
text(0.5,0.9,'Image 2 (DWI)','color','w','fontsize',fontsize,'horizontalalignment','center')

saveas(gcf,'reg_2a.png')


%%

image1 = '/cis/home/dtward/Documents/mricloud_atlases/DTIatlas/Adult/fa/s02a_fa.img';
avw = avw_img_read(image1);
I1 = avw.img;
I1 = padarray(I1,[1,1,1]);
x = linspace(-1,1,size(I1,2));
y = linspace(-1,1,size(I1,1));
z = linspace(-1,1,size(I1,3));
[X,Y,Z] = meshgrid(x,y,z);

% a good transform to bring it to  MNI
A = diag([0.9,0.85,0.9,1]);
B = inv(A);
Xs = B(1,1)*XI + B(1,2)*YI + B(1,3)*ZI + B(1,4);
Ys = B(2,1)*XI + B(2,2)*YI + B(2,3)*ZI + B(2,4);
Zs = B(3,1)*XI + B(3,2)*YI + B(3,3)*ZI + B(3,4);

F = griddedInterpolant({y,x,z},I1,'linear','nearest');
I1D = F(Ys,Xs,Zs);

Ishow = flip(squeeze(I1D(:,round(size(I,2)/2),:))',1);

figure('position',[1 1 512 512],'paperpositionmode','auto','inverthardcopy','off');
colormap gray
axes('position',[0,0,1,1]);
axis off;
imagesc(Ishow);
axis image
axes('position',[0,0,1,1],'visible','off')
text(0.5,0.9,'Registered 3','color','w','fontsize',fontsize,'horizontalalignment','center')

saveas(gcf,'reg_3b.png')


% a bad transform that brings it away
Xs = XI + linear*(randn*Xs + randn*Ys + randn*Zs) + quadratic*(randn*Xs.^2 + randn.*Xs.*Ys + randn.*Xs.*Zs + randn.*Ys.^2 + randn.*Ys.*Zs + randn.*Zs.^2);
Ys = YI + linear*(randn*Xs + randn*Ys + randn*Zs) + quadratic*(randn*Xs.^2 + randn.*Xs.*Ys + randn.*Xs.*Zs + randn.*Ys.^2 + randn.*Ys.*Zs + randn.*Zs.^2);
Zs = ZI + linear*(randn*Xs + randn*Ys + randn*Zs) + quadratic*(randn*Xs.^2 + randn.*Xs.*Ys + randn.*Xs.*Zs + randn.*Ys.^2 + randn.*Ys.*Zs + randn.*Zs.^2);

F = griddedInterpolant({y,x,z},I1,'linear','nearest');
I1D = F(Ys,Xs,Zs);

Ishow = flip(squeeze(I1D(:,round(size(I,2)/2),:))',1);

figure('position',[1 1 512 512],'paperpositionmode','auto','inverthardcopy','off');
colormap gray
axes('position',[0,0,1,1]);
axis off;
imagesc(Ishow);
axis image
axes('position',[0,0,1,1],'visible','off')
text(0.5,0.9,'Image 3 (DWI)','color','w','fontsize',fontsize,'horizontalalignment','center')

saveas(gcf,'reg_3a.png')


%%

image1 = '/cis/home/dtward/Documents/mricloud_atlases/DTIatlas/Adult/trace/s02a_trace.img';
avw = avw_img_read(image1);
I1 = avw.img;
I1 = padarray(I1,[1,1,1]);
x = linspace(-1,1,size(I1,2));
y = linspace(-1,1,size(I1,1));
z = linspace(-1,1,size(I1,3));
[X,Y,Z] = meshgrid(x,y,z);

% a good transform to bring it to  MNI
A = diag([0.9,0.85,0.9,1]);
B = inv(A);
Xs = B(1,1)*XI + B(1,2)*YI + B(1,3)*ZI + B(1,4);
Ys = B(2,1)*XI + B(2,2)*YI + B(2,3)*ZI + B(2,4);
Zs = B(3,1)*XI + B(3,2)*YI + B(3,3)*ZI + B(3,4);

F = griddedInterpolant({y,x,z},I1,'linear','nearest');
I1D = F(Ys,Xs,Zs);

Ishow = flip(squeeze(I1D(:,round(size(I,2)/2),:))',1);

figure('position',[1 1 512 512],'paperpositionmode','auto','inverthardcopy','off');
colormap gray
axes('position',[0,0,1,1]);
axis off;
imagesc(Ishow);
axis image
axes('position',[0,0,1,1],'visible','off')
text(0.5,0.9,'Registered 4','color','w','fontsize',fontsize,'horizontalalignment','center')

saveas(gcf,'reg_4b.png')


% a bad transform that brings it away
Xs = XI + linear*(randn*Xs + randn*Ys + randn*Zs) + quadratic*(randn*Xs.^2 + randn.*Xs.*Ys + randn.*Xs.*Zs + randn.*Ys.^2 + randn.*Ys.*Zs + randn.*Zs.^2);
Ys = YI + linear*(randn*Xs + randn*Ys + randn*Zs) + quadratic*(randn*Xs.^2 + randn.*Xs.*Ys + randn.*Xs.*Zs + randn.*Ys.^2 + randn.*Ys.*Zs + randn.*Zs.^2);
Zs = ZI + linear*(randn*Xs + randn*Ys + randn*Zs) + quadratic*(randn*Xs.^2 + randn.*Xs.*Ys + randn.*Xs.*Zs + randn.*Ys.^2 + randn.*Ys.*Zs + randn.*Zs.^2);

F = griddedInterpolant({y,x,z},I1,'linear','nearest');
I1D = F(Ys,Xs,Zs);

Ishow = flip(squeeze(I1D(:,round(size(I,2)/2),:))',1);

figure('position',[1 1 512 512],'paperpositionmode','auto','inverthardcopy','off');
colormap gray
axes('position',[0,0,1,1]);
axis off;
imagesc(Ishow);
axis image
axes('position',[0,0,1,1],'visible','off')
text(0.5,0.9,'Image 4 (trace)','color','w','fontsize',fontsize,'horizontalalignment','center')

saveas(gcf,'reg_4a.png')
