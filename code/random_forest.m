tfs=20;
lfs=15;

%%

p=0.2;
q=0.02;

A11=round(rand(n)<p);
A12=round(rand(n)<q);

A=[A11 A12; A12 A11];

spy(A)

%%

X0=[rand(n,2); rand(n,2)-[1 1]];
X1=[rand(n,2)-[0, 1]; rand(n,2)-[1 0]];

figure(1), clf
subplot(121)
cla, hold all
plot(X0(:,1),X0(:,2),'bo')
plot(X1(:,1),X1(:,2),'rx')
axis('square')
set(gca,'XTick',[],'YTick',[])
title('Data Matrix','fontsize',tfs)
xlabel('dimension 1','fontsize',lfs)
ylabel('dimension 2','fontsize',lfs)

subplot(122)
cla, hold all
A=[0 1; 1 0];
imagesc(A), colormap('default'), colorbar
axis('square')
set(gca,'XTick',[],'YTick',[])
title('Probability Matrix','fontsize',tfs)
xlabel('dimension 1','fontsize',lfs)
ylabel('dimension 2','fontsize',lfs)


% export_fig(['../images/.png']);
print('../images/sample_XOR.png','-dpng')


%%
p=10;
d=5;

A=zeros(p,p);
for i=1:d
   j=randperm(p);
   A(i,j(1))=1;
end


%%
figure(1), clf, 
colormap('gray')

subplot(121)
imagesc(A), axis('equal')
set(gca,'XTick',[],'YTick',[])
title('RF Projection Matrix','fontsize',tfs)
ylabel('features','fontsize',lfs)
xlabel('dimensions','fontsize',lfs)

n=50;
x=randn(p,n);
subplot(122)
imagesc(x), %axis('equal')
set(gca,'XTick',[],'YTick',[])
title('Data Matrix','fontsize',tfs)
xlabel('samples','fontsize',lfs)
ylabel('dimensions','fontsize',lfs)

% if print_fig,
    export_fig(['../images/RF1.png']);
% end

%%

A=zeros(p,p);
j=randperm(p*d);
A(j(1:d))=1;
A=A';

%%
figure(2), clf, 

subplot(121)
imagesc(A), axis('equal')
set(gca,'XTick',[],'YTick',[])
title('RerF Projection Matrix','fontsize',tfs)
set(gca,'XColor','none','YColor','none')
ylabel('features','fontsize',lfs)
xlabel('dimensions','fontsize',lfs)


n=50;
x=randn(p,n);
subplot(122)
imagesc(x), %axis('equal')
set(gca,'XTick',[],'YTick',[])
title('Data Matrix','fontsize',tfs)
xlabel('samples','fontsize',lfs)
ylabel('dimensions','fontsize',lfs)

% export_fig(['../images/RerF1.png']);
print('../images/RerF1.png','-dpng')

%%

left=0.03;
bottom=0.1;
width=0.3;
height=0.8;
eps=0.03;

pos = [left bottom width height];

figure(3), clf, colormap('gray') 
subplot('position',pos)
imagesc(A), 
axis('equal')
set(gca,'XTick',[],'YTick',[])
title('RerF Projection Matrix','fontsize',tfs)
ylabel('features','fontsize',lfs)
xlabel('dimensions','fontsize',lfs)
set(gca,'XColor','none','YColor','none')

pos(1)=pos(1)+width+eps;
subplot('position',pos)
imagesc(eye(p))
axis('equal')
set(gca,'XTick',[],'YTick',[])
title('Dictionary','fontsize',tfs)
set(gca,'XColor','none','YColor','none')


pos(1)=pos(1)+width+eps;
subplot('position', pos)
imagesc(x), %axis('equal')
set(gca,'XTick',[],'YTick',[])
title('Data Matrix','fontsize',tfs)
xlabel('samples','fontsize',lfs)
ylabel('dimensions','fontsize',lfs)

print('../images/RF_Dict1.png','-dpng')

%% random dictionary

pos = [left bottom width height];

figure(3), clf, colormap('gray') 
subplot('position',pos)
imagesc(A), 
axis('equal')
set(gca,'XTick',[],'YTick',[])
title('RerF Projection Matrix','fontsize',tfs)
ylabel('features','fontsize',lfs)
xlabel('dimensions','fontsize',lfs)
set(gca,'XColor','none','YColor','none')

pos(1)=pos(1)+width+eps;
subplot('position',pos)
imagesc(rand(p))
axis('equal')
set(gca,'XTick',[],'YTick',[])
title('Dictionary','fontsize',tfs)
set(gca,'XColor','none','YColor','none')


pos(1)=pos(1)+width+eps;
subplot('position', pos)
imagesc(x), %axis('equal')
set(gca,'XTick',[],'YTick',[])
title('Data Matrix','fontsize',tfs)
xlabel('samples','fontsize',lfs)
ylabel('dimensions','fontsize',lfs)

print('../images/RF_Dict2.png','-dpng')


%% structured dictionary

pos = [left bottom width height];

figure(3), clf, colormap('gray') 
subplot('position',pos)
imagesc(A), 
axis('equal')
set(gca,'XTick',[],'YTick',[])
title('RerF Projection Matrix','fontsize',tfs)
ylabel('features','fontsize',lfs)
xlabel('dimensions','fontsize',lfs)
set(gca,'XColor','none','YColor','none')

pos(1)=pos(1)+width+eps;
subplot('position',pos)

z=[1 1 zeros(1,8)];
ZZ=zeros(10);
for i=1:10
    ZZ(i,:)=circshift(z,i-2);
end
imagesc(ZZ)
axis('equal')
set(gca,'XTick',[],'YTick',[])
title('Dictionary','fontsize',tfs)
set(gca,'XColor','none','YColor','none')


pos(1)=pos(1)+width+eps;
subplot('position', pos)
imagesc(x), %axis('equal')
set(gca,'XTick',[],'YTick',[])
title('Data Matrix','fontsize',tfs)
xlabel('samples','fontsize',lfs)
ylabel('dimensions','fontsize',lfs)

print('../images/RF_Dict3.png','-dpng')
