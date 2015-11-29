
clear all;
close all;
base1 = 'C:\Users\vhasfckefays\Projects\MFStrength\ProVelMats\';
load ([base1 'VX601007TPad']);
bcVXmask=VX;
load('maskshftShrink7T');
mask=maskshft;
base2 = 'C:\Users\vhasfckefays\Projects\MFStrength\ProVelMatsPadUnWrp\';
load([base2 'VX4040247TPadUnWr']);
bcVX=VX;
slbg=27;
slsml=slbg;
slbg2=slbg;
bg= squeeze(mask(1,slbg,:,:));
bg2= squeeze(bcVXmask(1,slbg2,:,:));
bg (bg<0)=5;
bg2 (bg2<0)=2;
im = squeeze(bcVX(1,slsml,:,:));
%im2 = squeeze(bcVX(1,slsml2,:,:));
im (im<0)=2;
%im2 (im2<0)=2;

hf = figure('units','normalized','position',[.2 .2 .6 .6]);
ax1 = subplot(2,3,1);
ibg = imagesc(bg);
axis off
title('Big')
ax2 = subplot(2,3,4);
iim = imagesc(im);
axis off
title('Small')

ax3 = subplot(2,3,[2 5]);%[2:3, 5:6]
ibg = imagesc(bg);
axis off
hold on
%   Overlay the image, and set the transparency
wdth= size(im,2);
shft= (size(bg,2)-size(im,2))/2;
full= size(bg,1);
shftPlus= shft+0;

iim = imagesc(im,'XData', [1+shftPlus wdth+shftPlus],'YData',[1 full]); %[shft wdth+shft]
set(iim,'Alphadata',0.5);
title(sprintf('Using transparency while overlaying images'));

ax4 = subplot(2,3,[3 6]);%[2:3, 5:6]
ibg = imagesc(bg);
axis off
hold on
%   Overlay the image, and set the transparency
iim = imagesc(bg2,'XData',[1 full],'YData',[1 full]); 
set(iim,'Alphadata',0.5);
title(sprintf('Using transparency while overlaying images'));