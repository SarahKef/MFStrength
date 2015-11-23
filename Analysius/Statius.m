clear all;
close all;
FileSuf = {'60100', '6071','6040','404024', '404016', '40408', '40404'};
FilePre = {'VX','VY','VZ','Mag'};
FileFolder = {'3T','7T'}; 
% PadUnWrp
InputBaseFolder = 'C:\Users\vhasfckefays\Projects\MF Strength 4DFlow\ProVelMats';
OutputBaseFolder = 'C:\Users\vhasfckefays\Projects\MF Strength 4DFlow\AnalysisVelMats';
Stat = zeros(max(size(FileFolder)), max(size(FilePre)), max(size(FileSuf)), 3, 2);
ROI{1} = 30;
ROI{2} = 1:304;
ROI{3} = 1:304;

for i = 1:size(FileFolder,2)
    load ([InputBaseFolder '\' 'maskshft' FileFolder{i}]);
    for k = 1:size(FileSuf,2)
        for j = 1:size(FilePre,2)
            FileName = [FilePre{j} FileSuf{k} FileFolder{i} 'Pad']; %'UnWr'
            FilePath = [InputBaseFolder '\' FileName]; 
            
            Val = load (FilePath);
            eval(['Val = Val.' FilePre{j} ';']);
            Mask=repmat(maskshft,size(Val,1),1);
            Val = Val.*Mask;
            VoxStat = calcVoxStat(Val, 1); 
            for r = 1:3
                Stat(i, j, k, r, 1) = mean(mean(mean(VoxStat(r,:,:,:))));
                Stat(i, j, k, r, 2) = mean(mean(mean(VoxStat(r,ROI{1},ROI{2},ROI{3}))));
            end
            
            if exist('TI', 'var')
                TI = TI + VoxStat(2,:,:,:).^2;
            else
                TI = VoxStat(2,:,:,:).^2;
            end
            
            save([OutputBaseFolder '\' FilePre{j} FileSuf{k} FileFolder{i} 'Stat'], 'VoxStat');
            fprintf('i=%d, j=%d, k=%d\n', i, j, k);
        end

        TI = sqrt(TI);
        save([OutputBaseFolder '\' 'TI' FileSuf{k} FileFolder{i} 'Mag'], 'TI');
        
        clear TI;
    end
end
save([OutputBaseFolder '\' 'Stat'], 'Stat');
