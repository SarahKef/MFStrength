clear all;
close all;
load ('maskMRACut7T');
T=1;
SL=60;
d=cell(T,SL);
y=cell(T,SL);
maskshft=zeros(1,60,304,304);
for t=1:T
    for sl= 1:SL
        d{t,sl}= squeeze (mask(t,sl,:,:));
        y{t,sl}= shiftu(d{t,sl},0,5,1);
        y{t,sl}= shiftl(y{t,sl},0,0,1);
        maskshft(t,sl,:,:)=(y{t,sl});
    end
end
save maskshft7T maskshft