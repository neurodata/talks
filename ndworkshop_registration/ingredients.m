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



figure;
imagesc(xI,yI,I)
axis image
colormap gray;

[XI,YI] = meshgrid(x,y);


% now we need a sample space, make it a bit bigger, make it square
xJ = 1 : round(size(I,2)*1.1);
xJ = xJ - mean(xJ);
yJ = xJ;
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
text(0,yJ(1),'Identity','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
saveas(gcf,'ex2_id.png')
IDeye = ID;
%%
% get a rigid transform
t = randn(2,1)*5;
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
text(0,yJ(1),'Rigid','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
saveas(gcf,'ex2_rigid.png')


%%
% affine transform
A = expm([randn(2,2)*15*pi/180,randn(2,1)*5;0,0,0]);

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
text(0,yJ(1),'Affine','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
saveas(gcf,'ex2_affine.png')


%%
% spline
rng(2)
% need some displacement
vx = randn(size(XJ));
vy = randn(size(XJ));
r = 75;
[X_,Y_] = meshgrid(-r:r);
K = exp(-(X_.^2 + Y_.^2)/2/(r/3).^2);
K = K/sum(K(:))*750;
vx = convn(vx,K,'same');
vy = convn(vy,K,'same');
W = exp(-(XJ.^2 + YJ.^2)/2/20^2);
% vx = vx.*W;
% vy = vy.*W;
% vy = vy*0;
vx = vx + W*100;

Xs = XJ - vx;
Ys = YJ - vy;

F = griddedInterpolant({yI,xI},I,'linear','nearest');
ID = F(Ys,Xs);

figure('position',[1 1 512 512],'paperpositionmode','auto');
axes('position',[0,0,1,1]);
imagesc(xJ,yJ,ID)
axis image
colormap gray
text(0,yJ(1),'Smooth displacement','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
saveas(gcf,'ex2_disp_bad.png')

%%
% diffeomorphism
nt = 10;
dt = 1.0/nt;
Xs = XJ;
Ys = YJ;
for t = 1 : nt
    F = griddedInterpolant({yJ,xJ},Xs-XJ,'linear','nearest');
    Xs = F(YJ-vy*dt,XJ-vx*dt) + (XJ-vx*dt);
    
    F = griddedInterpolant({yJ,xJ},Ys-YJ,'linear','nearest');
    Ys = F(YJ-vy*dt,XJ-vx*dt) + (YJ-vy*dt);
end

F = griddedInterpolant({yI,xI},I,'linear','nearest');
ID = F(Ys,Xs);

figure('position',[1 1 512 512],'paperpositionmode','auto');
axes('position',[0,0,1,1]);
imagesc(xJ,yJ,ID)
axis image
colormap gray
text(0,yJ(1),'Smooth flow','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
saveas(gcf,'ex2_flow_bad.png')


%%
% spline
rng(2)
% need some displacement
vx = randn(size(XJ));
vy = randn(size(XJ));
vx = convn(vx,K,'same');
vy = convn(vy,K,'same');
W = exp(-(XJ.^2 + YJ.^2)/2/20^2);
% vx = vx.*W;
% vy = vy.*W;
% vy = vy*0;
% vx = vx + W*100;

Xs = XJ - vx;
Ys = YJ - vy;

F = griddedInterpolant({yI,xI},I,'linear','nearest');
ID = F(Ys,Xs);

figure('position',[1 1 512 512],'paperpositionmode','auto');
axes('position',[0,0,1,1]);
imagesc(xJ,yJ,ID)
axis image
colormap gray
text(0,yJ(1),'Smooth displacement','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
saveas(gcf,'ex2_disp.png')
% diffeomorphism
nt = 10;
dt = 1.0/nt;
Xs = XJ;
Ys = YJ;
for t = 1 : nt
    F = griddedInterpolant({yJ,xJ},Xs-XJ,'linear','nearest');
    Xs = F(YJ-vy*dt,XJ-vx*dt) + (XJ-vx*dt);
    
    F = griddedInterpolant({yJ,xJ},Ys-YJ,'linear','nearest');
    Ys = F(YJ-vy*dt,XJ-vx*dt) + (YJ-vy*dt);
end

F = griddedInterpolant({yI,xI},I,'linear','nearest');
ID = F(Ys,Xs);

figure('position',[1 1 512 512],'paperpositionmode','auto');
axes('position',[0,0,1,1]);
imagesc(xJ,yJ,ID)
axis image
colormap gray
text(0,yJ(1),'Smooth flow','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
saveas(gcf,'ex2_flow.png')


%%
% now how likely



%%
for n = 1 : 5
% spline
rng(2)
% need some displacement
vx = randn(size(XJ))*n/2;
vy = randn(size(XJ))*n/2;
vx = convn(vx,K,'same');
vy = convn(vy,K,'same');
W = exp(-(XJ.^2 + YJ.^2)/2/20^2);
% vx = vx.*W;
% vy = vy.*W;
% vy = vy*0;
% vx = vx + W*100;

Xs = XJ - vx;
Ys = YJ - vy;

F = griddedInterpolant({yI,xI},I,'linear','nearest');
ID = F(Ys,Xs);

figure('position',[1 1 512 512],'paperpositionmode','auto');
axes('position',[0,0,1,1]);
imagesc(xJ,yJ,ID)
axis image
colormap gray
text(0,yJ(1),'Smooth displacement','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
% saveas(gcf,'ex2_disp.png')
% diffeomorphism
nt = 10;
dt = 1.0/nt;
Xs = XJ;
Ys = YJ;
for t = 1 : nt
    F = griddedInterpolant({yJ,xJ},Xs-XJ,'linear','nearest');
    Xs = F(YJ-vy*dt,XJ-vx*dt) + (XJ-vx*dt);
    
    F = griddedInterpolant({yJ,xJ},Ys-YJ,'linear','nearest');
    Ys = F(YJ-vy*dt,XJ-vx*dt) + (YJ-vy*dt);
end

F = griddedInterpolant({yI,xI},I,'linear','nearest');
ID = F(Ys,Xs);

figure('position',[1 1 512 512],'paperpositionmode','auto');
axes('position',[0,0,1,1]);
imagesc(xJ,yJ,ID)
axis image
colormap gray
if n == 1
text(0,yJ(1),'Likely','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
elseif n == 5
    text(0,yJ(1),'Unlikely','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
end
saveas(gcf,['ex2_flow_n_'  num2str(n) '.png'])
end


%%
% just simple translation
for n = 1 : 5
    
    
    A = expm([randn(2,2)*15*pi/180,randn(2,1)*5;0,0,0]);
    
    A = [eye(2)*1.01,[(n - 3)*10;0];[0,0,1]];
    
    t = [(n - 3)*10;0];
    theta = 15*pi/180*(n-3)/2;
    A = [cos(theta),-sin(theta),t(1);
        sin(theta),cos(theta),t(2);
        0,0,1];
    A = A*diag([1.02,1.02,1]);


    
    B = inv(A);
    Xs = B(1,1)*XJ + B(1,2)*YJ + B(1,3);
    Ys = B(2,1)*XJ + B(2,2)*YJ + B(2,3);
    
    
    
    F = griddedInterpolant({yI,xI},I,'linear','nearest');
    ID = F(Ys,Xs);
    
    figure('position',[1 1 512 512],'paperpositionmode','auto');
    axes('position',[0,0,1,1]);
    Ishow = cat(3,IDeye,ID,IDeye);
    Ishow = (Ishow - min(Ishow(:)))/(max(Ishow(:)) - min(Ishow(:)));
    imagesc(xJ,yJ,Ishow);
    axis image
    colormap gray
    if n == 1
    text(0,yJ(1),'Bad','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
    elseif n == 3
        text(0,yJ(1),'Good','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
    elseif n == 5
        text(0,yJ(1),'Bad','color','w','verticalalignment','top','fontsize',fontsize,'horizontalalignment','center')
    end
    saveas(gcf,['ex2_affine_n_' num2str(n) '.png'])
    



end