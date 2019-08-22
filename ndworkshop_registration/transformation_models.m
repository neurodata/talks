clear all;
close all;
fclose all;

rng(1);
% here the goal will be
% 1. demonstrate some transformation groups
% 2. demonstrate similar versus dissimilar
% 3. demonstrate priors
fontsize = 30;
% load up an mricloud atlas
addpath /cis/home/dtward/Functions/avwQuiet
avw = avw_img_read(['/cis/home/dtward/Documents/mricloud_atlases/Adult27-55/Adt27-55_01_Adt27-55_01_MNI.img']);


I = flip(squeeze(avw.img(round(end/2)+5,:,:))',1);
I = padarray(I,[1 1],'both');
x = 1 : size(I,2);
y = 1 : size(I,1);
xI = x - mean(x);
yI = y - mean(y);

xI = linspace(-1,1,size(I,2));
yI = linspace(-size(I,1)/size(I,2),size(I,1)/size(I,2),size(I,1));



figure;
imagesc(xI,yI,I)
axis image
colormap gray;

[XI,YI] = meshgrid(xI,yI);


% now we need a sample space, make it a bit bigger, make it square
xJ = 1 : round(size(I,2)*1.1);
xJ = xJ - mean(xJ);
yJ = xJ;

xJ = linspace(-1.1,1.1,length(xJ));
yJ = linspace(-1.1,1.1,length(yJ));

[XJ,YJ] = meshgrid(xJ,yJ);

% show the first image
Xs = XJ;
Ys = YJ;
F = griddedInterpolant({yI,xI},I,'linear','nearest');
ID = F(Ys,Xs);

figure('position',[1 1 512 512],'paperpositionmode','auto');
axes('position',[0,0,1,1]);
imagesc(xJ,yJ,ID)
axis image
colormap gray
IDeye = ID;
contours = linspace(-2,2,21);

hold on;
contour(xJ,yJ,Xs,contours,'r','linewidth',2)
contour(xJ,yJ,Ys,contours,'r','linewidth',2)
contour(xJ,yJ,Xs,[0,0],'b','linewidth',2)
contour(xJ,yJ,Ys,[0,0],'b','linewidth',2)
hold off;

text(0,yJ(1),'Identity','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')

saveas(gcf,'tform_models_id.png')


%%
% get a rigid transform
close all;
for n = 1 : 6
t = randn(2,1)*0.1;
theta = randn*15*pi/180;
A = [cos(theta),-sin(theta),t(1);
    sin(theta),cos(theta),t(2);
    0,0,1];
B = inv(A);
Xs = B(1,1)*XJ + B(1,2)*YJ + B(1,3);
Ys = B(2,1)*XJ + B(2,2)*YJ + B(2,3);

F = griddedInterpolant({yI,xI},I,'linear','nearest');
ID = F(Ys,Xs);

figure('position',[1 1 512 512],'paperpositionmode','auto');
axes('position',[0,0,1,1]);
imagesc(xJ,yJ,ID)
axis image
colormap gray


hold on;
contour(xJ,yJ,Xs,contours,'r','linewidth',2)
contour(xJ,yJ,Ys,contours,'r','linewidth',2)
contour(xJ,yJ,Xs,[0,0],'b','linewidth',2)
contour(xJ,yJ,Ys,[0,0],'b','linewidth',2)
hold off;

text(0,yJ(1),'Rigid','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
saveas(gcf,['tform_models_rigid_' num2str(n) '.png'])


end

%%
% get a similitude transform
close all;
for n = 1 : 6
t = randn(2,1)*0.1;
theta = randn*15*pi/180;
s = exp(randn*0.15);

A = [cos(theta),-sin(theta),t(1);
    sin(theta),cos(theta),t(2);
    0,0,1]*diag([s,s,1]);
B = inv(A);
Xs = B(1,1)*XJ + B(1,2)*YJ + B(1,3);
Ys = B(2,1)*XJ + B(2,2)*YJ + B(2,3);

F = griddedInterpolant({yI,xI},I,'linear','nearest');
ID = F(Ys,Xs);

figure('position',[1 1 512 512],'paperpositionmode','auto');
axes('position',[0,0,1,1]);
imagesc(xJ,yJ,ID)
axis image
colormap gray


hold on;
contour(xJ,yJ,Xs,contours,'r','linewidth',2)
contour(xJ,yJ,Ys,contours,'r','linewidth',2)
contour(xJ,yJ,Xs,[0,0],'b','linewidth',2)
contour(xJ,yJ,Ys,[0,0],'b','linewidth',2)
hold off;

text(0,yJ(1),'Similarity','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
saveas(gcf,['tform_models_similarity_' num2str(n) '.png'])


end

%%
% get a affine transform
close all;
for n = 1 : 6
t = randn(2,1)*0.1;


A = [expm(randn(2,2)*0.15),t;0,0,1];
B = inv(A);
Xs = B(1,1)*XJ + B(1,2)*YJ + B(1,3);
Ys = B(2,1)*XJ + B(2,2)*YJ + B(2,3);

F = griddedInterpolant({yI,xI},I,'linear','nearest');
ID = F(Ys,Xs);

figure('position',[1 1 512 512],'paperpositionmode','auto');
axes('position',[0,0,1,1]);
imagesc(xJ,yJ,ID)
axis image
colormap gray


hold on;
contour(xJ,yJ,Xs,contours,'r','linewidth',2)
contour(xJ,yJ,Ys,contours,'r','linewidth',2)
contour(xJ,yJ,Xs,[0,0],'b','linewidth',2)
contour(xJ,yJ,Ys,[0,0],'b','linewidth',2)
hold off;

text(0,yJ(1),'Affine','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
saveas(gcf,['tform_models_affine_' num2str(n) '.png'])


end


%%
close all;
for n = 1 : 6
    close all;
    for n2 = 1 : 6

vx = randn(size(XJ));
vy = randn(size(XJ));

rs = [10,20,30,40,50,60];
r = rs(n);
mags = [1:6]*0.005;
mag = mags(n2);

[X_,Y_] = meshgrid(-r:r);
K = exp(-(X_.^2 + Y_.^2)/2/(r/3).^2);
vx = convn(vx,K,'same');
vy = convn(vy,K,'same');

std_ = std([vx(:);vy(:)]);
vx = vx/std_*mag;
vy = vy/std_*mag;

Xs = XJ + vx;
Ys = YJ + vy;

F = griddedInterpolant({yI,xI},I,'linear','nearest');
ID = F(Ys,Xs);

figure('position',[1 1 512 512],'paperpositionmode','auto');
axes('position',[0,0,1,1]);
imagesc(xJ,yJ,ID)
axis image
colormap gray


hold on;
contour(xJ,yJ,Xs,contours,'r','linewidth',2)
contour(xJ,yJ,Ys,contours,'r','linewidth',2)
contour(xJ,yJ,Xs,[0,0],'b','linewidth',2)
contour(xJ,yJ,Ys,[0,0],'b','linewidth',2)
hold off;

text(0,yJ(1),'Displacement','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
saveas(gcf,['tform_models_disp_' num2str(n) '_' num2str(n2) '.png'])


end
end


%%
close all;
for n = 1 : 6
    
    W = exp(-(XJ.^2 + YJ.^2)/2/0.2^2)*n*0.15;
    
    Xs = XJ + W;
    Ys = YJ + W;
    
    F = griddedInterpolant({yI,xI},I,'linear','nearest');
    ID = F(Ys,Xs);
    
    figure('position',[1 1 512 512],'paperpositionmode','auto');
    axes('position',[0,0,1,1]);
    imagesc(xJ,yJ,ID)
    axis image
    colormap gray
    
    
    hold on;
    contour(xJ,yJ,Xs,contours,'r','linewidth',2)
    contour(xJ,yJ,Ys,contours,'r','linewidth',2)
    contour(xJ,yJ,Xs,[0,0],'b','linewidth',2)
    contour(xJ,yJ,Ys,[0,0],'b','linewidth',2)
    hold off;
    
    text(0,yJ(1),'Displacement','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
    saveas(gcf,['tform_models_disp_' num2str(n)  '.png'])
    
    
end

%%
close all;
nt = 5;
for n = 1 : 6
    
    W = exp(-(XJ.^2 + YJ.^2)/2/0.2^2)*n*0.15;    
    vx = W;
    vy = W;
    Xs = XJ;
    Ys = YJ;
    for t = 1 : nt
        F = griddedInterpolant({yJ,xJ},Xs-XJ,'linear','nearest');
        Xs = F(YJ+vy/nt,XJ+vx/nt) + XJ+vx/nt;
        F = griddedInterpolant({yJ,xJ},Ys-YJ,'linear','nearest');
        Ys = F(YJ+vy/nt,XJ+vx/nt) + YJ+vy/nt;
    end
    
    
    
    
    F = griddedInterpolant({yI,xI},I,'linear','nearest');
    ID = F(Ys,Xs);
    
    figure('position',[1 1 512 512],'paperpositionmode','auto');
    axes('position',[0,0,1,1]);
    imagesc(xJ,yJ,ID)
    axis image
    colormap gray
    
    
    hold on;
    contour(xJ,yJ,Xs,contours,'r','linewidth',2)
    contour(xJ,yJ,Ys,contours,'r','linewidth',2)
    contour(xJ,yJ,Xs,[0,0],'b','linewidth',2)
    contour(xJ,yJ,Ys,[0,0],'b','linewidth',2)
    hold off;
    
    text(0,yJ(1),'Flow','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
    saveas(gcf,['tform_models_flow_' num2str(n)  '.png'])
    
    
end