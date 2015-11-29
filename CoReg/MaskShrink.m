clear all;
close all;
mask = load ('maskshft7T');
mask = mask.maskshft;
maskshft = zeros(size(mask));
for i=1:60
b=squeeze(mask(1,i,:,:));
maskshft(1,i,:,:) = bwmorph(b,'thin',1);
end
save maskshftShrink7T maskshft