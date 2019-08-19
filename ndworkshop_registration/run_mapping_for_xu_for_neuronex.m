% in this code we jointly estimate an affine and deformable mapping from
% atlas to target
% intensity differences are accounted for using a cubic polynomial
% transform
% missing tissue and artifacts are accounted for using a three component
% mixture model

clear all
close all;
fclose all;

tic

addpath /cis/home/dtward/Functions/plotting
addpath /cis/home/dtward/Functions/avwQuiet
addpath /cis/home/dtward/Functions/frame2Gif
addpath /cis/home/dtward/Functions/nrrd


%%
% input parameters
% open images
atlas_name = 'FluoroAtlas_Downsample.img';
atlas_name = '/cis/home/dtward/Documents/ARA/Mouse_CCF/average_template_25.nrrd';
atlas_name = '/cis/home/dtward/Documents/ARA/Mouse_CCF/average_template_50.nrrd';

target_name = '180517_ch1_Downsample.img';
target_name = '/cis/home/hliu/data/Data_public/Daniel/Josh_mice/180517_Downsample/180517_Downsample.img';

label_name = '/cis/home/dtward/Documents/ARA/Mouse_CCF/annotation/ccf_2017/annotation_25.nrrd';
label_name = '/cis/home/dtward/Documents/ARA/Mouse_CCF/annotation/ccf_2017/annotation_50.nrrd';

output_prefix = 'test_01_'; % all outputs should be prepended with this prefix, it can contain directories too

% save frames
save_all_iters_until = 0;
save_every_n_iters = 10;
write_gifs_every_n_iters = 10;

% other parameters
nT = 5; % timesteps in numerical integration of flow
sigmaM = 1; % we will standardize J so that std(J) is 1
sigmaR = 1e0;
sigmaA = sigmaM*10; % artifact, large variance
CA = 3; % estimate of artifact inensity, +3 sigma
sigmaB = sigmaM; % background
CB = -1; % estimate of background intensity, -1 sigma
niter = 500;
naffine = 100; % affine only
nM = 5; % number of m steps per e step
nMaffine = 1; % number of m steps per e step durring affine only
% a = dxI(1)*5;
a = 0.1250;
p = 2;
eV = 5e-3;
eL = 2e-4;
eT = 1e-3;
post_affine_reduce = 0.1; % after deformation starts, reduce step size for affine
eL = 1e-4; 

% for neuronex, let's get this done fast
naffine = 50;
niter = 200;
nM = 1;


%%
% load images
if ~isempty(strfind(atlas_name,'.img'))
    avw = avw_img_read(atlas_name);
    I = avw.img;
elseif ~isempty(strfind(atlas_name,'.nrrd'))
    [vol,meta] = nrrdread(atlas_name);
    I = double(vol);
end

% let's zero bad the atlas so there are no unfortunate boundary conditions pad it
I = padarray(I,[1,1,1],0,'both');
nxI = [size(I,2),size(I,1),size(I,3)];
% dxI = double( avw.hdr.dime.pixdim([3,2,4]) );
% dxI = [0.025,0.025,0.025];
dxI = [1,1,1]*0.05; % we are using the 50 micron atlas    
% for numerical stability we can standardize I, it doesn't change anything
% since we will estimate intensity transforms anyway
I = I - mean(I(:));
I = I / std(I(:));

% set up atlas domain
xI = (0 : nxI(1)-1)*dxI(1);
yI = (0 : nxI(2)-1)*dxI(2);
zI = (0 : nxI(3)-1)*dxI(3);
xI = xI - mean(xI);
yI = yI - mean(yI);
zI = zI - mean(zI);
[XI,YI,ZI] = meshgrid(xI,yI,zI);

fxI = (0:nxI(1)-1)/nxI(1)/dxI(1);
fyI = (0:nxI(2)-1)/nxI(2)/dxI(2);
fzI = (0:nxI(3)-1)/nxI(3)/dxI(3);
[FXI,FYI,FZI] = meshgrid(fxI,fyI,fzI);

% loat target
avw = avw_img_read(target_name);
J = avw.img;
nxJ = double( avw.hdr.dime.dim([3,2,4]) );
% dxJ = double( avw.hdr.dime.pixdim([3,2,4]) );
% dxJ = [0.025,0.025,0.025];
% header information for this example image is not correct
dxJ = [1,1,1]*0.05;
% let's standardize J as well to make things simpler
J = J - mean(J(:));
J = J / std(J(:));
xJ = (0 : nxJ(1)-1)*dxJ(1);
yJ = (0 : nxJ(2)-1)*dxJ(2);
zJ = (0 : nxJ(3)-1)*dxJ(3);
xJ = xJ - mean(xJ);
yJ = yJ - mean(yJ);
zJ = zJ - mean(zJ);
[XJ,YJ,ZJ] = meshgrid(xJ,yJ,zJ);


danfigure(1);
sliceView(xI,yI,zI,I)
danfigure(2);
sliceView(xJ,yJ,zJ,J)


%%
% set up differential operator to enforce smoothness
LL = (1 - 2 * a^2 * ( (cos(2*pi*dxI(1)*FXI) - 1)/dxI(1)^2 + (cos(2*pi*dxI(2)*FYI) - 1)/dxI(2)^2 + (cos(2*pi*dxI(3)*FZI) - 1)/dxI(3)^2 )).^(2*p);
Khat = 1.0./LL;
dt = 1/nT;
%%
% initialize
A = eye(4);
% we will need to encode a permutation
A = [0,0,1,0;1,0,0,0;0,1,0,0;0,0,0,1]'*A;
A = diag([-1,1,1,1])*A;

vtx = zeros([size(I),nT]);
vty = zeros([size(I),nT]);
vtz = zeros([size(I),nT]);
It = zeros([size(I),nT]);
It(:,:,:,1) = I;

% we need an initial linear transformation to compute our first weight
Jq = quantile(J(:),[0.1 0.9]);
Iq = quantile(I(:),[0.1,0.9]);
coeffs = [mean(Jq)-mean(Iq)*diff(Jq)/diff(Iq); diff(Jq)/diff(Iq)];
coeffs = [coeffs;0;0];





%%
% start
Esave = [];
EMsave = [];
ERsave = [];
EAsave = [];
EBsave = [];
Asave = [];
frame_errW = [];
frame_W = [];
frame_I = [];
frame_errRGB = [];
for it = 1 : niter
    
    % deform image
    phiinvx = XI;
    phiinvy = YI;
    phiinvz = ZI;
    for t = 1 : nT*(it > naffine)
        % sample image
        if t > 1
            F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
            It(:,:,:,t) = F(phiinvy,phiinvx,phiinvz);
        end
        % update diffeo, add and subtract identity for better boundary conditions
        Xs = XI - vtx(:,:,:,t)*dt;
        Ys = YI - vty(:,:,:,t)*dt;
        Zs = ZI - vtz(:,:,:,t)*dt;
        F = griddedInterpolant({yI,xI,zI},phiinvx-XI,'linear','nearest');
        phiinvx = F(Ys,Xs,Zs) + Xs;
        F = griddedInterpolant({yI,xI,zI},phiinvy-YI,'linear','nearest');
        phiinvy = F(Ys,Xs,Zs) + Ys;
        F = griddedInterpolant({yI,xI,zI},phiinvz-ZI,'linear','nearest');
        phiinvz = F(Ys,Xs,Zs) + Zs;
        
    end
    % get final deformed image 
    % note that this image is not needed for the code so it is commented
    % out
    % if you need it for an output, copy this part to the end before
    % writing outputs
    %     F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
    %     phiI = F(phiinvy,phiinvx,phiinvz);
    
    % now apply affine, go to sampling of J
    % ideally I should fix the double interpolation here
    % for now just leave it
    B = inv(A);
    Xs = B(1,1)*XJ + B(1,2)*YJ + B(1,3)*ZJ + B(1,4);
    Ys = B(2,1)*XJ + B(2,2)*YJ + B(2,3)*ZJ + B(2,4);
    Zs = B(3,1)*XJ + B(3,2)*YJ + B(3,3)*ZJ + B(3,4);
    % we could have computed phiI and then transformed it again
    % this would be double interpolation and create blurring
    % the below section is commented out, and replaced with applying all
    % the transforms (diffeo and affine) at once
    %     F = griddedInterpolant({yI,xI,zI},phiI,'linear','nearest');
    %     AphiI = F(Ys,Xs,Zs);
    F = griddedInterpolant({yI,xI,zI},phiinvx-XI,'linear','nearest');
    phiinvBx = F(Ys,Xs,Zs) + Xs;
    F = griddedInterpolant({yI,xI,zI},phiinvy-YI,'linear','nearest');
    phiinvBy = F(Ys,Xs,Zs) + Ys;
    F = griddedInterpolant({yI,xI,zI},phiinvz-ZI,'linear','nearest');
    phiinvBz = F(Ys,Xs,Zs) + Zs;
    F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
    AphiI = F(phiinvBy,phiinvBx,phiinvBz);
    
    
    % now apply the cubic intensity transformation
    D = [ones(size(AphiI(:))),AphiI(:),AphiI(:).^2,AphiI(:).^3];
    fAphiI = reshape(D*coeffs,size(J));
    
    
    danfigure(3);
    sliceView(xJ,yJ,zJ,fAphiI);
    
    
    err = fAphiI - J;
    
    danfigure(6666);
    sliceView(xJ,yJ,zJ,cat(4,J,fAphiI,J));
    
    
    
    % now a weight
    doENumber = nMaffine;
    if it > naffine; doENumber = nM; end
    if ~mod(it-1,doENumber)
        WM = 1/sqrt(2*pi*sigmaM^2)*exp(-1.0/2.0/sigmaM^2*err.^2);
        WA = 1/sqrt(2*pi*sigmaA^2)*exp(-1.0/2.0/sigmaA^2*(CA - J).^2);
        WB = 1/sqrt(2*pi*sigmaB^2)*exp(-1.0/2.0/sigmaB^2*(CB - J).^2);
        
        Wsum = WM + WA + WB;
        WM = WM./Wsum;
        WA = WA./Wsum;
        WB = WB./Wsum;
        
        % there is probably a numerically better way to do this
        % if I took logs and then subtractedthe max and then took
        % exponentials maybe
        % not worth it, the numerics are okay
        danfigure(45);
%         sliceView(xJ,yJ,zJ,WM);
        sliceView(xJ,yJ,zJ,cat(4,WM,WA,WB));
        
    end
    errW = err.*WM;
    danfigure(4);
    sliceView(xJ,yJ,zJ,errW);
    
    
    % cost
    vtxhat = fft(fft(fft(vtx,[],1),[],2),[],3);
    vtyhat = fft(fft(fft(vty,[],1),[],2),[],3);
    vtzhat = fft(fft(fft(vtz,[],1),[],2),[],3);
    ER = sum(sum(sum(LL.*sum(abs(vtxhat).^2 + abs(vtyhat).^2 + abs(vtzhat).^2,4))))/2/sigmaR^2*dt*prod(dxI)/(size(I,1)*size(I,2)*size(I,3));
    EM = sum(sum(sum((fAphiI - J).^2.*WM)))*prod(dxJ)/2/sigmaM^2;
    EA = sum(sum(sum((CA - J).^2.*WA)))*prod(dxJ)/2/sigmaA^2;
    EB = sum(sum(sum((CB - J).^2.*WB)))*prod(dxJ)/2/sigmaB^2;
    E = ER + EM + EA + EB;
    % note that this energy will not necessarily decrease
    % the marginal likelihood will increase, but we don't know how to
    % calculate that easily.  This is the whole point of the EM algorithm
    fprintf(1,'Energy %g, reg %g, match %g, artifact %g, background %g\n',E,ER,EM,EA,EB);
    Esave = [Esave,E];
    ERsave = [ERsave,ER];
    EMsave = [EMsave,EM];
    EAsave = [EAsave,EA];
    EBsave = [EBsave,EB];
    
    % at the end, we'll update affine and coeffs
    % first let's do affine gradient
    % use a right perturbation A \mapsto A expm(e dA) and take derivative
    % with respect to e at e = 0
    [fAphiI_x,fAphiI_y,fAphiI_z] = gradient(fAphiI,dxJ(1),dxJ(2),dxJ(3));
    grad = zeros(4,4);
    for r = 1 : 3
        for c = 1 : 4
            dA = (double((1:4)'==r)) * double(((1:4)==c));
            
            AdAB = A * dA * B;
            AdABX = AdAB(1,1)*XJ + AdAB(1,2)*YJ + AdAB(1,3)*ZJ + AdAB(1,4);
            AdABY = AdAB(2,1)*XJ + AdAB(2,2)*YJ + AdAB(2,3)*ZJ + AdAB(2,4);
            AdABZ = AdAB(3,1)*XJ + AdAB(3,2)*YJ + AdAB(3,3)*ZJ + AdAB(3,4);
            
            grad(r,c) = -sum(sum(sum(errW.*(fAphiI_x.*AdABX + fAphiI_y.*AdABY + fAphiI_z.*AdABZ))))*prod(dxJ)/sigmaM^2;
            
        end
    end
    % if your affine needs to be rigid, uncomment the line below
    %     grad(1:2,1:2) = grad(1:2,1:2) - grad(1:2,1:2)';
    
    % now pull back the error,
    % again we will avoid double interpolation
    %     F = griddedInterpolant({yJ,xJ,zJ},errW,'linear','nearest'); % this should really be zero padded
    %     Xs = A(1,1)*XI + A(1,2)*YI + A(1,3)*ZI + A(1,4);
    %     Ys = A(2,1)*XI + A(2,2)*YI + A(2,3)*ZI + A(2,4);
    %     Zs = A(3,1)*XI + A(3,2)*YI + A(3,3)*ZI + A(3,4);
    %     lambda1 = -F(Ys,Xs,Zs)/abs(det(A))/sigmaM^2;
    % I should zero pad this so interpolation works properly
    %     lambda1pad = padarray(lambda1,[1,1,1],0,'both');
    %     xIpad = [xI(1)-dxI(1),xI,xI(end)+dxI(1)];
    %     yIpad = [yI(1)-dxI(2),yI,yI(end)+dxI(2)];
    %     zIpad = [zI(1)-dxI(3),zI,zI(end)+dxI(3)];
    % and flow back
    phi1tinvx = XI;
    phi1tinvy = YI;
    phi1tinvz = ZI;
    for t = nT*(it>naffine) : -1 : 1
        % update diffeo (note plus instead of minus)
        Xs = XI + vtx(:,:,:,t)*dt;
        Ys = YI + vty(:,:,:,t)*dt;
        Zs = ZI + vtz(:,:,:,t)*dt;
        F = griddedInterpolant({yI,xI,zI},phi1tinvx-XI,'linear','nearest');
        phi1tinvx = F(Ys,Xs,Zs) + Xs;
        F = griddedInterpolant({yI,xI,zI},phi1tinvy-YI,'linear','nearest');
        phi1tinvy = F(Ys,Xs,Zs) + Ys;
        F = griddedInterpolant({yI,xI,zI},phi1tinvz-ZI,'linear','nearest');
        phi1tinvz = F(Ys,Xs,Zs) + Zs;
        % determinant of jacobian
        [phi1tinvx_x,phi1tinvx_y,phi1tinvx_z] = gradient(phi1tinvx,dxI(1),dxI(2),dxI(3));
        [phi1tinvy_x,phi1tinvy_y,phi1tinvy_z] = gradient(phi1tinvy,dxI(1),dxI(2),dxI(3));
        [phi1tinvz_x,phi1tinvz_y,phi1tinvz_z] = gradient(phi1tinvz,dxI(1),dxI(2),dxI(3));
        detjac = phi1tinvx_x.*(phi1tinvy_y.*phi1tinvz_z - phi1tinvy_z.*phi1tinvz_y) ...
            - phi1tinvx_y.*(phi1tinvy_x.*phi1tinvz_z - phi1tinvy_z.*phi1tinvz_x) ...
            + phi1tinvx_z.*(phi1tinvy_x.*phi1tinvz_y - phi1tinvy_y.*phi1tinvz_x);
        % again, we are avoiding double interpolation
        %         F = griddedInterpolant({yIpad,xIpad,zIpad},lambda1pad,'linear','nearest'); % should be zero padded so that data out of FOV doesn't contribute
        %         lambda = F(phi1tinvy,phi1tinvx,phi1tinvz).*detjac;
                
        % deform the error once
        Aphi1tinvx = A(1,1)*phi1tinvx + A(1,2)*phi1tinvy + A(1,3)*phi1tinvz + A(1,4);
        Aphi1tinvy = A(2,1)*phi1tinvx + A(2,2)*phi1tinvy + A(2,3)*phi1tinvz + A(2,4);
        Aphi1tinvz = A(3,1)*phi1tinvx + A(3,2)*phi1tinvy + A(3,3)*phi1tinvz + A(3,4);
        F = griddedInterpolant({yJ,xJ,zJ},(-errW/sigmaM^2),'linear','nearest');
        lambda = F(Aphi1tinvy,Aphi1tinvx,Aphi1tinvz).*detjac.*abs(det(A));
        
        % get the gradient of the image
        fI = It(:,:,:,t);        
        D = [ones(size(fI(:))),fI(:),fI(:).^2,fI(:).^3];        
        fI = reshape(D*coeffs,size(I));
        [fI_x,fI_y,fI_z] = gradient(fI,dxI(1),dxI(2),dxI(3));
        
        % set up the gradient
        gradx = fI_x.*lambda;
        grady = fI_y.*lambda;
        gradz = fI_z.*lambda;
        
        % kernel and reg
        % note I could add more smoothness here if I wanted
        gradx = ifftn(fftn(gradx).*Khat + vtxhat(:,:,:,t)/sigmaR^2,'symmetric');
        grady = ifftn(fftn(grady).*Khat + vtyhat(:,:,:,t)/sigmaR^2,'symmetric');
        gradz = ifftn(fftn(gradz).*Khat + vtzhat(:,:,:,t)/sigmaR^2,'symmetric');
        
        % now update
        vtx(:,:,:,t) = vtx(:,:,:,t) - gradx*eV;
        vty(:,:,:,t) = vty(:,:,:,t) - grady*eV;
        vtz(:,:,:,t) = vtz(:,:,:,t) - gradz*eV;
    end
    
    % now update coeffs (make sure to include weight properly)
    D = [ones(size(AphiI(:))),AphiI(:),AphiI(:).^2,AphiI(:).^3];
    coeffs = ( D'*bsxfun(@times,D,WM(:)) ) \ (D'*bsxfun(@times,J(:),WM(:)));
    % I can also update my constants
    CB = sum(WB(:).*J(:))/sum(WB(:));
    CA = sum(WA(:).*J(:))/sum(WA(:));
    
    % update A
    e = [ones(3)*eL,ones(3,1)*eT;0,0,0,0];
    if it > naffine
        e = e * post_affine_reduce;
    end
    A = A * expm(-e.*grad);
    
    
    danfigure(8);
    Asave = [Asave,A(:)];
    subplot(1,3,1)
    plot(Asave([1,2,3,5,6,7,9,10,11],:)')
    title('linear part')
    subplot(1,3,2)
    plot(Asave([13,14,15],:)')
    ylabel mm
    title('translation part')
    legend('x','y','z','location','best')
    subplot(1,3,3)
    plot([Esave;ERsave;EMsave;EAsave;EBsave]')
    legend('tot','reg','match','artifact','background')
    title('Energy')
    drawnow;
    
   
    % write out some data
    if it <= save_all_iters_until || ~mod(it-1,save_every_n_iters)
        frame_errRGB = [frame_errRGB,getframe(6666)];
        frame_errW = [frame_errW,getframe(4)];
        frame_I = [frame_I,getframe(3)];
        frame_W = [frame_W,getframe(45)];
    end
    if it <= save_all_iters_until || ~mod(it-1,write_gifs_every_n_iters)
        frame2Gif(frame_errRGB,[output_prefix 'err_RGB.gif'])
        frame2Gif(frame_errW,[output_prefix 'err_W.gif'])
        frame2Gif(frame_I,[output_prefix 'AphiI.gif'])
        frame2Gif(frame_W,[output_prefix 'W.gif'])
        saveas(8, [output_prefix 'plots.png'])
    end
    

    
end

toc

%%
% now the algorithm is done
% we want to define any outputs here
% right now there are no outputs

% this image is not used in the algorithm, but you may want it as an output
F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
phiI = F(phiinvy,phiinvx,phiinvz);


%%
% map the labels (TO DO)
[L,meta] = nrrdread(label_name);
sliceView((L>0)+randn(size(L))*0.1)
labels = unique(L(:));

% build a sparse matrix of colors
colors = sparse(1,3);
rng(1);
for i = 1 : length(labels)
    colors(labels(i)+1,:) = rand(1,3); % this takes a surprisingly long time
end
colors(1,:) = 0;
% 
R = reshape(full(colors(L+1,1)),size(L));
G = reshape(full(colors(L+1,2)),size(L));
B = reshape(full(colors(L+1,3)),size(L));
figure
sliceView(xI,yI,zI,cat(4,R,G,B));

% deform the RGB
for c = 1 : 3    
    F = griddedInterpolant({yI,xI,zI},double(L(:,:,c)),'linear','nearest');
    AphiL(:,:,c) = F(phiinvBy,phiinvBx,phiinvBz);
end
% an alpha blend
alpha = 0.5;



%%
% I need to make a picture showing the deformation as an isosurface
% I can show the atlas before and after transform, and the overlayd
qlim = [0.01,0.99];
climI = quantile(I(:),qlim);
climJ = quantile(J(:),qlim);

danfigure(99);
Ishow = I;
imagesc(zI,yI,squeeze(Ishow(:,round(end/2),:)),climI)
axis image
colormap gray

hold on;
contour(zI,yI,squeeze(XI(:,round(end/2),:)),-10:10,'r','linewidth',2)
contour(zI,yI,squeeze(YI(:,round(end/2),:)),-10:10,'r','linewidth',2)
contour(zI,yI,squeeze(ZI(:,round(end/2),:)),-10:10,'r','linewidth',2)
hold off



F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
phiI = F(phiinvy,phiinvx,phiinvz);
Ishow = phiI;
danfigure(100);
imagesc(zI,yI,squeeze(Ishow(:,round(end/2),:)),climI)
axis image
colormap gray


hold on;
contour(zI,yI,squeeze(phi1tinvx(:,round(end/2),:)),-10:10,'r','linewidth',2)
contour(zI,yI,squeeze(phi1tinvy(:,round(end/2),:)),-10:10,'r','linewidth',2)
contour(zI,yI,squeeze(phi1tinvz(:,round(end/2),:)),-10:10,'r','linewidth',2)
hold off



Ishow = I-phiI;
danfigure(101);
imagesc(zI,yI,squeeze(Ishow(:,round(end/2),:)))
axis image
colormap gray

hold on;
contour(zI,yI,squeeze(phi1tinvx(:,round(end/2),:)),-10:10,'r','linewidth',2)
contour(zI,yI,squeeze(phi1tinvy(:,round(end/2),:)),-10:10,'r','linewidth',2)
contour(zI,yI,squeeze(phi1tinvz(:,round(end/2),:)),-10:10,'r','linewidth',2)
hold off


% apply Ainv to J
F = griddedInterpolant({yJ,xJ,zJ},J,'linear','nearest');
Xs = A(1,1)*XI + A(1,2)*YI + A(1,3)*ZI + A(1,4);
Ys = A(2,1)*XI + A(2,2)*YI + A(2,3)*ZI + A(2,4);
Zs = A(3,1)*XI + A(3,2)*YI + A(3,3)*ZI + A(3,4);
AiJ = F(Ys,Xs,Zs);


Ishow = AiJ;
danfigure(102);
imagesc(zI,yI,squeeze(Ishow(:,round(end/2),:)),climJ)
axis image
colormap gray


%%
close all
% okay
% need initial A
A_ = eye(4);
% we will need to encode a permutation
A_ = [0,0,1,0;1,0,0,0;0,1,0,0;0,0,0,1]'*A_;
A_ = diag([-1,1,1,1])*A_;
B = inv(A_);

Xs = B(1,1)*XJ + B(1,2)*YJ + B(1,3)*ZJ + B(1,4);
Ys = B(2,1)*XJ + B(2,2)*YJ + B(2,3)*ZJ + B(2,4);
Zs = B(3,1)*XJ + B(3,2)*YJ + B(3,3)*ZJ + B(3,4);
F = griddedInterpolant({yI,xI,zI},I,'linear','nearest');
I0 = F(Ys,Xs,Zs);

I0s = squeeze(I0(:,round(end/2),:));
I0s = (I0s - climI(1))/diff(climI);
Js = squeeze(J(:,round(end/2),:));
Js = (Js - climJ(1))/diff(climJ);


danfigure(1)
set(1,'paperpositionmode','auto','inverthardcopy','off','color','w')

Ishow = cat(3,I0s,Js,I0s);
imagesc(zJ,yJ,Ishow)
hold on;
contour(zJ,yJ,squeeze(Xs(:,round(end/2),:)),-10:10,'r','linewidth',2)
contour(zJ,yJ,squeeze(Ys(:,round(end/2),:)),-10:10,'r','linewidth',2)
contour(zJ,yJ,squeeze(Zs(:,round(end/2),:)),-10:10,'r','linewidth',2)
hold off
axis off
title('Before registration','fontsize',20)

saveas(gcf,'ardent0.png')

%
% now the deformed atlas
Is = squeeze(AphiI(:,round(end/2),:));
Is = (Is - climI(1))/diff(climI);
danfigure(3)
set(3,'paperpositionmode','auto','inverthardcopy','off','color','w')
Ishow = cat(3,Is,Js,Is);
imagesc(zJ,yJ,Ishow)
hold on;
contour(zJ,yJ,squeeze(phiinvBx(:,round(end/2),:)),-10:10,'r','linewidth',2)
% contour(zJ,yJ,squeeze(phiinvBy(:,round(end/2),:)),-20:20,'r','linewidth',2)
contour(zJ,yJ,squeeze(phiinvBz(:,round(end/2),:)),-10:10,'r','linewidth',2)
hold off
axis off
title('After registration','fontsize',20)

saveas(gcf,'ardent1.png')

danfigure(5);
set(5,'paperpositionmode','auto','inverthardcopy','off','color','w')
imagesc(zJ,yJ,I0s,[0,1])
axis image
colormap gray
axis off
title('Allen CCF Atlas','fontsize',20)
saveas(gcf,'ardent2.png')

danfigure(6);
set(6,'paperpositionmode','auto','inverthardcopy','off','color','w')
imagesc(zJ,yJ,Js,[0,1])
axis image
colormap gray
axis off
title('Serial Two Photon Experiment','fontsize',20)
saveas(gcf,'ardent3.png')