clear all;
close all;
Config{1,1} = 'C:\Users\vhasfckefays\Projects\MFStrength\AnalysisVelMatsPad\Stat';
Config{1,2} = 'C:\Users\vhasfckefays\Projects\MFStrength\AnalysisVelMatsPad';
Config{1,3} = [0 1.5 0 1.1];
Config{2,1} = 'C:\Users\vhasfckefays\Projects\MFStrength\AnalysisVelMatsPadUnWrp\Stat';
Config{2,2} = 'C:\Users\vhasfckefays\Projects\MFStrength\AnalysisVelMatsPadUnWrp';
Config{2,3} = [-0.2 2.2 -0.2 2.5];
Config{3,1} = '/Users/sarah/code/MFStrength/StatData/Stat';
Config{3,2} = '/Users/sarah/code/MFStrength/StatData';
Config{3,3} = [-0.2 2 -0.2 2.5];
configIdx = 2;
ROI{1} = 1:40;
ROI{2} = 1:304;
ROI{3} = 1:304;

load(Config{configIdx,1}); %PadUnWrp
% dim:1 3T,7T dim:2 Vx,Vy,Vz,Mag dim:3 all cases dim:4 mean,std,RMS dim:5 all over,ROI
%Cases = {'60x100x24', '60x71x24','60x40x24','40x40x24', '40x40x16', '40x40x8', '40x40x4'};
Cases = {'60x100-24°', '60x71-24°','60x40-24°','40x40-24°', '40x40-16°', '40x40-8°', '40x40-4°'};
BaseFolder = Config{configIdx,2};
%% <std>
% 3T vs. 7T over cases


fh = figure(1);
res = 2; %1:all over, 2:ROI
metric = 4;

x=1:7;
std(:,1) = squeeze(Stat(1,metric,:,2,res)) ./ squeeze(Stat(1,metric,1,1,res));
std(:,2) = squeeze(Stat(2,metric,:,2,res)) ./ squeeze(Stat(2,metric,1,1,res));
mycolor=[0 0 1;1 0 0];
b = bar(std, 'grouped');
colormap(mycolor);
hold on;
h1 = plot(x-0.15,std(:,1),'bo-.','linewidth',2);
h2 = plot(x+0.15,std(:,2),'rd-','linewidth',2);
set(gca, 'XTickLabel', Cases);
ylabel('$$\frac{<\sigma_{V_{mag}}>}{<\overline{V}_{mag}>}$$','Interpreter','latex');
xlabel('Cases');
title('Relative mean of st. dev. of $$\overline{V}_{mag}$$ over all voxels','Interpreter','latex');
set(fh, 'Position', [520   378   790   420]);
legend([h1 h2], {'3T', '7T'}, 'Location', 'NorthWest');

%% Distribution Histograms 
FileSuf = {'60100', '6071','6040','404024', '404016', '40408', '40404'};
FilePre = {'VX','VY','VZ','Mag'};
FileFolder = {'3T','7T'}; 
Metric = {'mean', '\sigma', 'RMS'};

Case = [4, 1, 2]; %dim:1 FilePre, dim:2 FileSuf, dim:3 Metric

figure(2);
CustAxis = [0 0.6 0 100];
col = {'b', 'r'};
minThreshold = 0.0001;
for cs = 1:4
    Case(2) = cs + 3;
    for i = 1:2
        FileName = [FilePre{Case(1)} FileSuf{Case(2)} FileFolder{i} 'Stat'];
        FilePath = [BaseFolder filesep FileName]; 
        Val = load (FilePath);
        Val = Val.VoxStat(:,ROI{1},ROI{2},ROI{3});
        Val = reshape(Val(Case(3),:,:,:), 1, numel(Val(Case(3),:,:,:)));

        subplot(2,4,4*(i-1)+cs);
        [c, x] = hist(Val(Val > minThreshold), 2000);
        c = c ;%./ max(c);
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
        xlabel('\sigma_{V_{mag}}');
        ylabel('Normalized Frequency');
         if configIdx == 3
             title(sprintf('Histogram of %s of V_{mag} of %s: %s\nWidth = %f Skewness = %f', Metric{Case(3)}, FileFolder{i}, Cases{Case(2)}, width, skw));
         else
             title(sprintf('Histogram of %s of V_{mag} of %s: %s\nWidth = %f', Metric{Case(3)}, FileFolder{i}, Cases{Case(2)}, width));
         end
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
snapshotCase = [4, 2, 27 0.3 0.35]; %dim:1 FileSuf, dim:2 FileFolder, dim:3 Z slice, dim:4 Vmag
snapshot = cell(1,2);
baseSnapshot = cell(size(FileFolder, 2), 2);

for j = 1:size(FileSuf, 2)
    for k = 1:size(FileFolder, 2)
        for i = 1:2
            FileName = [FilePre{FileCfg{2,i}} FileSuf{j} FileFolder{k}];
            FilePath = [BaseFolder filesep FileName FileCfg{1,i}]; 
            Val = load(FilePath);
            eval(['X = squeeze(Val.' FileCfg{3,i} '(1,ROI{1},ROI{2},ROI{3})' ')' ';']);
            if (j == 1)
                baseSnapshot{k,i} = squeeze(X(snapshotCase(3),ROI{2},ROI{3}));
            end
            if (j == snapshotCase(1) && k == snapshotCase(2))
                snapshot{i} = squeeze(X(snapshotCase(3),ROI{2},ROI{3}));
            end
            
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
set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
count = 0;
binsNo = 100;
SNR = zeros([size(TIperMag) 7]);
Lines = cell(size(TIperMag));
baseIdx = cell(size(FileFolder));
baseNoise = cell(size(FileFolder));
for k =1:size(FileFolder, 2)
    for j = 1:size(FileSuf, 2)
        count = count + 1;
        subplot(size(FileFolder, 2), size(FileSuf, 2), count);
        axis(Config{configIdx,3});
        hold on;
        
        [SortedData, Idx] = sort(TIperMag{j,k}(:,1), 1);
        SortedData(:,2) = TIperMag{j,k}(Idx,2);
        binX = round(linspace(1, size(TIperMag{j,k}, 1), binsNo+1))';
        binYmin = zeros(0, 2);
        maxLimit = 0.9;
        minLimit = 0.5;
        binYmax = zeros(0, 2);
        for i=1:binsNo
            if (j ~= 7 || k ~= 1 || i < binsNo*minLimit)
            [~, I] = min(SortedData(binX(i):binX(i+1), 2));
            binYmin = [binYmin; SortedData(binX(i) + I - 1, :)]; %#ok<*AGROW>
            end
            if (i > round(binsNo*maxLimit)) 
                continue;
            end
            [~, I] = max(SortedData(binX(i):binX(i+1), 2));
            newPoint = SortedData(binX(i) + I - 1, :);
            if (newPoint(2) < max(binYmax(:,2)))
                continue;
            end
            binYmax = [binYmax; SortedData(binX(i) + I - 1, :)];
        end
        
        pLow = polyfit(binYmin(:,1), binYmin(:,2), 1);
        pHigh = polyfit(binYmax(:,1), binYmax(:,2), 1);
        pMid(1) = tan(0.5 .* (atan(pLow(1)) + atan(pHigh(1))));
        pMid(2) = pLow(2) + (pHigh(2) - pLow(2)).*(pMid(1) - pLow(1))./(pHigh(1) - pLow(1));
        
        %pMid = [0 0.1];
        
        predictedTI = polyval(pMid, TIperMag{j,k}(:,1));
        noiseIdx = predictedTI <= TIperMag{j,k}(:,2);
        signalIdx = predictedTI > TIperMag{j,k}(:,2);
        
        if (j == 1)
            baseIdx{k} = signalIdx;
            predictedTI = polyval(pMid, baseSnapshot{k, 2});
            baseNoise{k} = (baseSnapshot{k, 1} > 0) & (predictedTI < baseSnapshot{k, 1});
        end
        
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
        noiseMean = mean(TIperMag{j,k}(noiseIdx,2));
        signalMean = mean(TIperMag{j,k}(signalIdx,2));
        noiseRatio = noiseNo./(noiseNo + signalNo);
        signalRatio = signalNo./(noiseNo + signalNo);
        
        SNR(j,k,:) = [noiseNo signalNo noiseRatio signalRatio totalError noiseMean signalMean];
        Lines{j,k} = [pMid; pLow; pHigh];
        
        %title(sprintf('Noise: %d%%, Signal: %d%%\nTotal Error=%.2f', round(100*noiseRatio), round(100*signalRatio), totalError));
        
        if (count == j) 
            title(sprintf('%s\n\nQuant: %0.2f%%\nUnquant: %0.2f%%\nMeanQuant = %.3f', Cases{count}, 100*signalRatio, 100*noiseRatio, signalMean), 'FontSize', 12,'fontweight','bold');
        else
            title(sprintf('Quant: %0.2f%%\nUnquant: %0.2f%%\nMeanQuant = %.3f', 100*signalRatio, 100*noiseRatio, signalMean),'FontSize', 12,'fontweight','bold');
            xlabel('$$\overline{V}_{mag}$$','Interpreter','latex','FontSize', 12,'fontweight','bold');
        end
        if (j == 1)
            ylabel(sprintf('%s\n\n\n%s', FileFolder{k}, '$$\sqrt{\Sigma(\sigma_i^2)}$$'),'Interpreter','latex','FontSize', 12,'fontweight','bold');
        end
        
        if (j == snapshotCase(1) && k == snapshotCase(2))
            predictedTI = polyval(pMid, snapshot{2});
            idxSnapshot = (snapshot{2} >= snapshotCase(4)) & (snapshot{2} <= snapshotCase(5));
            idxNoise1 = (snapshot{1} > 0) & (predictedTI < snapshot{1});
            idxNoise2 = (snapshot{1} > 0) & (predictedTI < snapshot{1}) & (baseNoise{k} ~= 1);
            idxNoise3 = (snapshot{1} > 0) & (predictedTI < snapshot{1}) & (baseNoise{k} ~= 1) & idxSnapshot;
            idxNoise4 = (snapshot{1} > 0) & (predictedTI < snapshot{1}) & idxSnapshot;
            idxSignal = (snapshot{1} > 0) & (predictedTI > snapshot{1});
            idxSignal2 = (snapshot{1} > 0) & (predictedTI > snapshot{1}) & idxSnapshot;
            slice1 = idxNoise1 - idxSignal;
            slice2 = idxNoise2 - idxSignal;
            slice3 = idxNoise3 - idxSignal;
            slice4 = idxNoise3 + idxNoise2 - idxSignal - idxSignal2;
            slice5 = idxNoise4 + idxNoise1 - idxSignal - idxSignal2;
        end
    end
end
%%
figure(4);
hold on;
plot(SNR(:,1,2) ./ SNR(:,1,1),'bo-.', 'LineWidth', 2);
plot(SNR(:,2,2) ./ SNR(:,2,1),'rd-', 'LineWidth', 2);
legend('3T FNR', '7T FNR');
set(gca, 'XTickLabel', Cases);
xlabel('Cases');
ylabel('FNR');

%%
figure(5);
hold on;
plot(SNR(:,1,3),'bo-.', 'LineWidth', 2, 'MarkerSize', 10);
plot(SNR(:,1,4),'bs--', 'LineWidth', 2, 'MarkerSize', 10);
plot(SNR(:,2,3),'rd-', 'LineWidth', 2, 'MarkerSize', 10);
plot(SNR(:,2,4),'rx:', 'LineWidth', 2, 'MarkerSize', 10);
legend('3T Noise Ratio', '3T Flow Ratio', '7T Noise Ratio', '7T Flow Ratio', 'Location', 'Best');
set(gca, 'XTickLabel', Cases);
xlabel('Cases');
ylabel('Ratio');


%%
figure(6);
imagesc(slice1);
set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);

%%
figure(7);
set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
count = 0;
TIperMagUnbiased = cell(TIperMag);
for k =1:size(FileFolder, 2)
    for j = 1:size(FileSuf, 2)
        count = count + 1;
        subplot(size(FileFolder, 2), size(FileSuf, 2), count);
        axis(Config{configIdx,3});
        hold on;
        
        TIperMagUnbiased{j,k} = TIperMag{j,k}(baseIdx{k},:);
        predictedTI = polyval(Lines{j,k}(1, :), TIperMagUnbiased{j,k}(:,1));
        noiseIdx = predictedTI <= TIperMagUnbiased{j,k}(:,2);
        signalIdx = predictedTI > TIperMagUnbiased{j,k}(:,2);

        plot(TIperMagUnbiased{j,k}(noiseIdx,1), TIperMagUnbiased{j,k}(noiseIdx,2), 'r.');
        plot(TIperMagUnbiased{j,k}(signalIdx,1), TIperMagUnbiased{j,k}(signalIdx,2), 'g.');
%        plot(TIperMagUnbiased{j,k}(:,1), TIperMagUnbiased{j,k}(:,2), '.');

        %line(Config{configIdx,3}(1:2), polyval(Lines{j,k}(1, :), Config{configIdx,3}(1:2)), 'Color', 'k', 'LineWidth' , 2);
        
        noiseNo = sum(noiseIdx);
        signalNo = sum(signalIdx);
        totalError = sum(TIperMagUnbiased{j,k}(:,2));
        noiseRatio = noiseNo./(noiseNo + signalNo);
        signalRatio = signalNo./(noiseNo + signalNo);
        

        %title(sprintf('Noise: %d%%, Signal: %d%%\nTotal Error=%.2f', round(100*noiseRatio), round(100*signalRatio), totalError));
        title(sprintf('Signal: %0.2f%%\nNoise: %0.2f%%\nTotal Error=%.2f', 100*signalRatio, 100*noiseRatio, totalError));
        xlabel('<V_{mag}>');
        ylabel('$$\sqrt{\Sigma(\delta_i^2)}$$','Interpreter','latex');
        
    end
end


%%
figure(8);
imagesc(slice2);
set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);

%%
figure(9);
imagesc(slice3);
set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);

%%
figure(10);
imagesc(slice4);
set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);

%%
figure(11);
imagesc(slice5);
set (gcf, 'Position', [0,0,1596,826]);
map = [0, 1.0, 0
    0, 0.4, 0
    1, 1, 1
    1, 0.6, 0
    1, 0, 0];
colormap(map); 
caxis([-2 2]);

%%
figure(12);
set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
count = 0;
vLines = [0.1 0.15 0.3 0.35 0.5 0.55 1.3 1.35];
for k = 2
    for j = [1 7]
        count = count + 1;
        subplot(1, 2, count);
        axis(Config{configIdx,3});
        hold on;
        
        predictedTI = polyval(Lines{j,k}(1, :), TIperMag{j,k}(:,1));
        noiseIdx = predictedTI <= TIperMag{j,k}(:,2);
        signalIdx = predictedTI > TIperMag{j,k}(:,2);

        plot(TIperMag{j,k}(noiseIdx,1), TIperMag{j,k}(noiseIdx,2), 'r.');
        plot(TIperMag{j,k}(signalIdx,1), TIperMag{j,k}(signalIdx,2), 'g.');
%        plot(TIperMagUnbiased{j,k}(:,1), TIperMagUnbiased{j,k}(:,2), '.');

        line(Config{configIdx,3}(1:2), polyval(Lines{j,k}(1, :), Config{configIdx,3}(1:2)), 'Color', 'k', 'LineWidth' , 2);
        
        for vl = vLines
            line([vl vl], Config{configIdx,3}(1:2), 'Color', 'k', 'LineWidth' , 1);
        end
        
        noiseNo = sum(noiseIdx);
        signalNo = sum(signalIdx);
        totalError = sum(TIperMag{j,k}(:,2));
        noiseRatio = noiseNo./(noiseNo + signalNo);
        signalRatio = signalNo./(noiseNo + signalNo);
        

        %title(sprintf('Noise: %d%%, Signal: %d%%\nTotal Error=%.2f', round(100*noiseRatio), round(100*signalRatio), totalError));
        title(sprintf('Signal: %0.2f%%\nNoise: %0.2f%%\nTotal Error=%.2f', 100*signalRatio, 100*noiseRatio, totalError));
        xlabel('<V_{mag}>');
        ylabel('$$\sqrt{\Sigma(\delta_i^2)}$$','Interpreter','latex');
        
    end
end


%% <Vmag>
% 3T vs. 7T over cases


fh = figure(13);
res = 2; %1:all over, 2:ROI
metric = 4;

x=1:7;
vmag(:,1) = squeeze(Stat(1,metric,:,1,res));
vmag(:,2) = squeeze(Stat(2,metric,:,1,res));
mycolor=[0 0 1;1 0 0];
b = bar(vmag, 'grouped');
colormap(mycolor);
hold on;
h1 = plot(x-0.15,vmag(:,1),'bo-.','linewidth',2);
h2 = plot(x+0.15,vmag(:,2),'rd-','linewidth',2);
set(gca, 'XTickLabel', Cases);
ylabel('$$<\overline{V}_{mag}>$$','Interpreter','latex');
xlabel('Cases');
title('Spatial mean of temporal mean of $$\overline{V}_{mag}$$ over all voxels','Interpreter','latex');
set(fh, 'Position', [520   378   790   420]);
legend([h1 h2], {'3T', '7T'}, 'Location', 'NorthWest');

%%
figure(14);
hold on;
plot(SNR(:,1,7) ,'bo-.', 'LineWidth', 2);
plot(SNR(:,2,7) ,'rd-', 'LineWidth', 2);
legend('3T', '7T');
set(gca, 'XTickLabel', Cases);
xlabel('Cases');
ylabel('Variation Mean of Quantifiable Voxels (m/s)');


%%
figure(15);
set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
count = 0;
vLines = [0.3 0.35];
for k = [1 2]
    for j = 1
        count = count + 1;
        subplot(1, 2, count);
        axis([0 1.5 0 1]);
        hold on;
        
        predictedTI = polyval(Lines{j,k}(1, :), TIperMag{j,k}(:,1));
        noiseIdx1 = (predictedTI <= TIperMag{j,k}(:,2)) & (TIperMag{j,k}(:,1) >= vLines(1)) & (TIperMag{j,k}(:,1) <= vLines(2));
        noiseIdx2 = (predictedTI <= TIperMag{j,k}(:,2)) & ((TIperMag{j,k}(:,1) < vLines(1)) | (TIperMag{j,k}(:,1) > vLines(2)));
        signalIdx1 = (predictedTI > TIperMag{j,k}(:,2)) & (TIperMag{j,k}(:,1) >= vLines(1)) & (TIperMag{j,k}(:,1) <= vLines(2));
        signalIdx2 = (predictedTI > TIperMag{j,k}(:,2)) & ((TIperMag{j,k}(:,1) < vLines(1)) | (TIperMag{j,k}(:,1) > vLines(2)));

        plot(TIperMag{j,k}(noiseIdx2,1), TIperMag{j,k}(noiseIdx2,2), '.', 'Color', [1 0.6 0]);
        plot(TIperMag{j,k}(signalIdx2,1), TIperMag{j,k}(signalIdx2,2), '.', 'Color', [0 0.4 0]);
        plot(TIperMag{j,k}(noiseIdx1,1), TIperMag{j,k}(noiseIdx1,2), 'r.', 'MarkerSize', 15);
        plot(TIperMag{j,k}(signalIdx1,1), TIperMag{j,k}(signalIdx1,2), 'g.');

        line(Config{configIdx,3}(1:2), polyval(Lines{j,k}(1, :), Config{configIdx,3}(1:2)), 'Color', 'k', 'LineWidth' , 2);
        
        for vl = vLines
            line([vl vl], Config{configIdx,3}(1:2), 'Color', 'k', 'LineWidth' , 1);
        end

        
    end
end
