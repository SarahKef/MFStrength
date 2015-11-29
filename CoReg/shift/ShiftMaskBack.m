clear all;
close all;
F = 'Shrink'; % 'Shrink';
f = '7T';
load (['maskshft' F f]);
path = ['C:\Users\vhasfckefays\Projects\MF Strength 4DFlow\OrgVelMats\' f '\VX\'];
sk = 'VX';
SL = 60; 
ExtShft3T = 1;
ExtShft7T = 0;
cs = '404024';
V = load([path sk cs]);
V = V.bcVX;
V (V<0) = -1;
V (V>0) = 0;
FOV1 = size(maskshft,4);
FOV2 = size(V,4);
T=1;
mask=zeros(T, size(V,2), size(V,3), size(V,4));

for sl = ((SL-size(V,2))/2)+1:size(V,2)+((SL-size(V,2))/2)
    d = squeeze (maskshft(:,sl,:,:));
    d = shiftl(d,0,((FOV1 - FOV2)/2)+ eval(['ExtShft' f]),1);
    d = d(:, 1:size(V,4));
    mask(T,sl - ((SL-size(V,2))/2) ,:,:)= d;
end

save (['mask' sk cs F f], 'mask');

sk = 15;
shft = 0;
vel = imagesc(squeeze(V(1,sk,:,:)));
msk = squeeze(mask(1,sk,:,:));
axis off
hold on
%   Overlay the image, and set the transparency
iim = imagesc(msk,'XData',[1+shft size(V,4)+shft],'YData',[1 size(V,3)]); 
set(iim,'Alphadata',0.6);