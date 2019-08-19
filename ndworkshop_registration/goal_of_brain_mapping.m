% make overview figs
% make figures for my what is the goal of brain mapping slides

clear all
close all;
fclose all;
addpath /cis/home/dtward/Functions/plotting
addpath /cis/home/dtward/Functions/nrrd/
addpath /cis/home/dtward/Functions/avwQuiet/
addpath /cis/home/dtward/Functions/downsample/

% input arguments in this cell
template_name = '/cis/home/dtward/Documents/ARA/Mouse_CCF/average_template_50.nrrd';
label_name = '/cis/home/dtward/Documents/ARA/Mouse_CCF/annotation_50.nrrd';
target_name = '/cis/home/dtward/Documents/ailey-mouse-May-2019/gad2cre_9um_9um_5um.tif';
% pixel size is required here as the tif data structure does not store it
dxJ0 = [9 9 5];





%% 
% allen atlas
[I,meta] = nrrdread(template_name);
I = double(I);
dxI = diag(sscanf(meta.spacedirections,'(%d,%d,%d) (%d,%d,%d) (%d,%d,%d)',[3,3]))';



% want padding of 1mm
npad = round(1000/dxI(1))*0;
I = padarray(I,[1,1,1]*npad,0,'both');


% scale it for numerical stability, since its scale doesn't matter
I = I - mean(I(:));
I = I/std(I(:));
nxI = [size(I,2),size(I,1),size(I,3)];
xI = (0:nxI(1)-1)*dxI(1);
yI = (0:nxI(2)-1)*dxI(2);
zI = (0:nxI(3)-1)*dxI(3);
xI = xI - mean(xI);
yI = yI - mean(yI);
zI = zI - mean(zI);
danfigure(1);
sliceView(xI,yI,zI,I);

[XI,YI,ZI] = meshgrid(xI,yI,zI);
fxI = (0:nxI(1)-1)/nxI(1)/dxI(1);
fyI = (0:nxI(2)-1)/nxI(2)/dxI(2);
fzI = (0:nxI(3)-1)/nxI(3)/dxI(3);
[FXI,FYI,FZI] = meshgrid(fxI,fyI,fzI);

[L,meta] = nrrdread(label_name);
L = padarray(L,[1,1,1]*npad,0,'both');



% vikram mouse
info = imfinfo(target_name);
%%
% downsample to about same res as atlas
down = round(dxI./dxJ0);
for f = 1 : length(info)
    disp(['File ' num2str(f) ' of ' num2str(length(info))])
    J_ = double(imread(target_name,f));
    if f == 1
        nxJ0 = [size(J_,2),size(J_,1),length(info)];
        nxJ = floor(nxJ0./down);
        J = zeros(nxJ(2),nxJ(1),nxJ(3));
    end
    % downsample J_
    Jd = zeros(nxJ(2),nxJ(1));
    for i = 1 : down(1)
        for j = 1 : down(2)
            Jd = Jd + J_(i:down(2):down(2)*nxJ(2), j:down(1):down(1)*nxJ(1))/down(1)/down(2);
        end
    end
    
    slice = floor( (f-1)/down(3) ) + 1;
    if slice > nxJ(3)
        break;
    end
    J(:,:,slice) = J(:,:,slice) + Jd/down(3);
    
    if ~mod(f-1,10)
    danfigure(1234);
    imagesc(J(:,:,slice));
    axis image
    drawnow;
    end
end
dxJ = dxJ0.*down;
xJ = (0:nxJ(1)-1)*dxJ(1);
yJ = (0:nxJ(2)-1)*dxJ(2);
zJ = (0:nxJ(3)-1)*dxJ(3);

xJ = xJ - mean(xJ);
yJ = yJ - mean(yJ);
zJ = zJ - mean(zJ);

% J = avw.img;
% J(isnan(J)) = 0;
% danfigure(2);
% sliceView(xJ,yJ,zJ,J)
% [XJ,YJ,ZJ] = meshgrid(xJ,yJ,zJ);
J0 = J; % save it
J0_orig = J0;

danfigure(2);
sliceView(xJ,yJ,zJ,J0_orig);



%%
% make my two figures
rng(1);
colors = rand(256,3);
% colors = 0.5*colors + 0.5;
colors(1,:) = 0;
Lm = mod(L,256)+1;
LRGB = reshape([colors(Lm,1),colors(Lm,2),colors(Lm,3)],[size(Lm),3]);

Ishow = I;

figure;
imagesc(zI,yI,squeeze(Ishow(:,round(size(I,2)/2),:)))
axis image
colormap gray
axis off;


qlim = [0.01,0.99];
clim = quantile(I(:),qlim);

Ishow = (I - clim(1))/diff(clim);

label_alpha = 0.5;
Ishow = bsxfun(@plus, Ishow*(1-label_alpha) , LRGB*label_alpha);



figure;
imagesc(zI,yI,squeeze(Ishow(:,round(size(I,2)/2),:,:)))
axis image
colormap gray
axis off;


out = squeeze(Ishow(:,round(size(I,2)/2),:,:));
imwrite(out,'ex1_allen_slice.png')

%%
out = squeeze(Ishow(round(size(I,1)/2),:,:,:));
figure;
imagesc(zI,xI,out)
axis image
colormap gray
axis off;
imwrite(out,'ex1_allen_slice.png')


%%
Ishow = J;
out = squeeze(Ishow(:,:,round(size(J,2)/2),:));

qlim = [0.45,0.99];
clim = quantile(out(:),qlim);

out = (out - clim(1))/diff(clim);
out(out < 0) = 0;
out(out > 1) = 1;

% out = out.^1.5;
out = out(30:end-20,:);
figure;
imagesc(xJ, yJ(30:end-20),out)
axis image
colormap gray;
imwrite(out,'ex1_target_slice.png')



%%
% now a deforming grid
rng(2);
nx = 136;
x = 1:nx;
[X,Y] = meshgrid(x);

phiX = X + 0.1*Y + 0.001*X.^2 + 0.002*X.*Y;
phiY = Y + 0.1*X + 0.001*X.^2 - 0.001*Y.^2 - 0.001*X.*Y;


% need some displacement
vx = randn(size(X));
vy = randn(size(X));
r = 70;
[X_,Y_] = meshgrid(-r:r);
K = exp(-(X_.^2 + Y_.^2)/2/(r/3).^2);
K = K/sum(K(:))*500;
vx = convn(vx,K,'same');
vy = convn(vy,K,'same');
phiX = X + 0.002*X.*Y + vx;
phiY = Y + vy;


down = 15;
figure;
hold on;
for i = 1 : down : nx
    plot(phiX(i,:),phiY(i,:),'k','linewidth',2)
end
for i = 1 : down : nx
    plot(phiX(:,i),phiY(:,i),'k','linewidth',2)
end

axis tight
axis off
set(gcf,'paperpositionmode','auto')
saveas(gcf,['ex1_warp.png'])
