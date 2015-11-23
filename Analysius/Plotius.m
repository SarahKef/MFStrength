clear all;
close all;
Config{1,1} = 'C:\Users\vhasfckefays\Projects\MFStrength\AnalysisVelMatsPad\Stat';
Config{1,2} = 'C:\Users\vhasfckefays\Projects\MFStrength\AnalysisVelMatsPad';
Config{1,3} = [0 1.5 0 1.1];
Config{2,1} = 'C:\Users\vhasfckefays\Projects\MFStrength\AnalysisVelMatsPadUnWrp\Stat';
Config{2,2} = 'C:\Users\vhasfckefays\Projects\MFStrength\AnalysisVelMatsPadUnWrp';
Config{2,3} = [-0.2 2 -0.2 2.5];
Config{3,1} = '/Users/sarah/code/MFStrength/StatData/Stat';
Config{3,2} = '/Users/sarah/code/MFStrength/StatData';
Config{3,3} = [-0.2 2 -0.2 2.5];
configIdx = 2;

load(Config{configIdx,1}); %PadUnWrp
% dim:1 3T,7T dim:2 Vx,Vy,Vz,Mag dim:3 all cases dim:4 mean,std,RMS dim:5 all over,ROI
Cases = {'60x100x24', '60x71x24','60x40x24','40x40x24', '40x40x16', '40x40x8', '40x40x4'};
BaseFolder = Config{configIdx,2};
%% <std>
% 3T vs. 7T over cases


fh = figure(1);
res = 1; %1:all over, 2:ROI
metric = 4;

x=1:7;
std(:,1) = squeeze(Stat(1,metric,:,2,res)) ./ squeeze(Stat(1,metric,1,1,res));
std(:,2) = squeeze(Stat(2,metric,:,2,res)) ./ squeeze(Stat(2,metric,1,1,res));
mycolor=[0 0 1;1 0 0];
b = bar(std, 'grouped');
colormap(mycolor);
legend('3T', '7T', 'Location', 'NorthWest');
hold on;
plot(x-0.15,std(:,1),'bo-','linewidth',2);
plot(x+0.15,std(:,2),'ro-','linewidth',2);
set(gca, 'XTickLabel', Cases);
ylabel('Relative <\sigma_{Vmag}>');
xlabel('Cases');
title('Relative Mean of St. dev. of V_{mag} over all voxels');
set(fh, 'Position', [520   378   790   420]);

%% Distribution Histograms 
FileSuf = {'60100', '6071','6040','404024', '404016', '40408', '40404'};
FilePre = {'VX','VY','VZ','Mag'};
FileFolder = {'3T','7T'}; 
Metric = {'mean', '\sigma', 'RMS'};

Case = [4, 1, 2]; %dim:1 FilePre, dim:2 FileSuf, dim:3 Metric

figure(2);
CustAxis = [0 0.6 0 1.1];
col = {'b', 'r'};
minThreshold = 0.0001;
for cs = 1:4
    Case(2) = cs + 3;
    for i = 1:2
        FileName = [FilePre{Case(1)} FileSuf{Case(2)} FileFolder{i} 'Stat'];
        FilePath = [BaseFolder filesep FileName]; 
        Val = load (FilePath);
        Val = Val.VoxStat;
        Val = reshape(Val(Case(3),:,:,:), 1, numel(Val(Case(3),:,:,:)));

        subplot(2,4,4*(i-1)+cs);
        [c, x] = hist(Val(Val > minThreshold), 2000);
        c = c ./ max(c);
        widthIdx = (c >= 0.5);
        minX = min(x(widthIdx));
        maxX = max(x(widthIdx));
        width = maxX - minX;
        if configIdx == 3
            skw = skewness(c);
        else
            skw = 0;
        end
        bar(x, c, col{i});
        axis(CustAxis);
        line([minX maxX], [0.5 0.5], 'Color', 'k', 'LineWidth', 2);
        xlabel('\sigma_{|V|}');
        ylabel('Normalized Frequency');
        title(sprintf('Histogram of %s of |V| of %s: %s\nWidth = %f Skewness = %f', Metric{Case(3)}, FileFolder{i}, FileSuf{Case(2)}, width, skw));
    end
end

%% TI vs. <V>mag

FileSuf = {'60100', '6071','6040','404024', '404016', '40408', '40404'};
FilePre = {'VX','VY','VZ','Mag','TI'};
FileFolder = {'3T','7T'}; 
Metric = {'mean', '\sigma', 'RMS'}; % TI60407TMag  Mag60407TStat

FileCfg {1,1} = 'Mag'; 
FileCfg {1,2} = 'Stat'; 
FileCfg {2,1} = 5;
FileCfg {2,2} = 4; 
FileCfg {3,1} = 'TI'; 
FileCfg {3,2} = 'VoxStat'; 
TIperMag = cell (size(FileSuf, 2), size(FileFolder, 2));


for j = 1:size(FileSuf, 2)
    for k = 1:size(FileFolder, 2)
        for i = 1:2
            FileName = [FilePre{FileCfg{2,i}} FileSuf{j} FileFolder{k}];
            FilePath = [BaseFolder filesep FileName FileCfg{1,i}]; 
            Val = load(FilePath);
            eval(['X = squeeze(Val.' FileCfg{3,i} '(1,:,:,:)' ')' ';']);
            A  = (reshape(squeeze(X), [1, numel(X)]))';
            if (i == 1)
                TI = A;
            else
                idx = (A > 0);
                TIperMag{j,k}(:,1) = A(idx);
                TIperMag{j,k}(:,2) = TI(idx);
            end
        end
    end
end

figure(3);
count = 0;
binsNo = 100;
SNR = zeros([size(TIperMag) 5]);
Lines = cell(size(TIperMag));
for k =1:size(FileFolder, 2)
    for j = 1:size(FileSuf, 2)
        count = count + 1;
        subplot(size(FileFolder, 2), size(FileSuf, 2), count);
        axis(Config{configIdx,3});
        hold on;
        
        [SortedData, Idx] = sort(TIperMag{j,k}(:,1), 1);
        SortedData(:,2) = TIperMag{j,k}(Idx,2);
        binX = round(linspace(1, size(TIperMag{j,k}, 1), binsNo+1))';
        binYmin = zeros(binsNo, 2);
        maxLimit = 0.9;
        binYmax = zeros(round(binsNo*maxLimit), 2);
        for i=1:binsNo
            binYmin(i, :) = min(SortedData(binX(i):binX(i+1), :));
            if (i > round(binsNo*maxLimit)) 
                continue;
            end
            binYmax(i, :) = max(SortedData(binX(i):binX(i+1), :));
        end
        
        pLow = polyfit(binYmin(:,1), binYmin(:,2), 1);
        pHigh = polyfit(binYmax(:,1), binYmax(:,2), 1);
        pMid(1) = tan(0.5 .* (atan(pLow(1)) + atan(pHigh(1))));
        pMid(2) = pLow(2) + (pHigh(2) - pLow(2)).*(pMid(1) - pLow(1))./(pHigh(1) - pLow(1));
        
        predictedTI = polyval(pMid, TIperMag{j,k}(:,1));
        noiseIdx = predictedTI <= TIperMag{j,k}(:,2);
        signalIdx = predictedTI > TIperMag{j,k}(:,2);
        plot(TIperMag{j,k}(noiseIdx,1), TIperMag{j,k}(noiseIdx,2), 'r.');
        plot(TIperMag{j,k}(signalIdx,1), TIperMag{j,k}(signalIdx,2), 'g.');
%        plot(TIperMag{j,k}(:,1), TIperMag{j,k}(:,2), '.');
%        plot(binYmin(:,1), binYmin(:,2), 'm.');
%        plot(binYmax(:,1), binYmax(:,2), 'b.');

%        line(Config{configIdx,3}(1:2), polyval(pLow, Config{configIdx,3}(1:2)), 'Color', 'r');
%        line(Config{configIdx,3}(1:2), polyval(pHigh, Config{configIdx,3}(1:2)), 'Color', 'g');
        line(Config{configIdx,3}(1:2), polyval(pMid, Config{configIdx,3}(1:2)), 'Color', 'k', 'LineWidth' , 2);
        
        noiseNo = sum(noiseIdx);
        signalNo = sum(signalIdx);
        totalError = sum(TIperMag{j,k}(:,2));
        noiseRatio = noiseNo./(noiseNo + signalNo);
        signalRatio = signalNo./(noiseNo + signalNo);
        
        SNR(j,k,:) = [noiseNo signalNo noiseRatio signalRatio totalError];
        Lines{j,k} = [pMid; pLow; pHigh];
        
        %title(sprintf('Noise: %d%%, Signal: %d%%\nTotal Error=%.2f', round(100*noiseRatio), round(100*signalRatio), totalError));
        title(sprintf('Noise: %0.2f%%, Signal: %0.2f%%\nTotal Error=%.2f', 100*noiseRatio, 100*signalRatio, totalError));
        xlabel('<V_{mag}>');
        ylabel('$$\sqrt{\Sigma(\delta_i^2)}$$','Interpreter','latex');
    end
end
%%
figure(4);
hold on;
plot(SNR(:,1,2) ./ SNR(:,1,1),'bo-', 'LineWidth', 2);
plot(SNR(:,2,2) ./ SNR(:,2,1),'ro-', 'LineWidth', 2);
legend('3T SNR', '7T SNR');
set(gca, 'XTickLabel', Cases);
xlabel('Cases');

%%
figure(5);
hold on;
plot(SNR(:,1,3),'rx-');
plot(SNR(:,1,4),'b.-');
plot(SNR(:,2,3),'kx-');
plot(SNR(:,2,4),'g.-');
legend('3T Noise Ratio', '3T Signal Ratio', '7T Noise Ratio', '7T Signal Ratio', 'Location', 'Best');
set(gca, 'XTickLabel', Cases);
xlabel('Cases');


%save([BaseFolder filesep 'results'], 'SNR', 'Lines', 'TIperMag');


