clear all;
close all;
FOV=[304 304];
FileSuf = {'60100', '6071','6040','404024', '404016', '40408', '40404'};
FilePre = {'VX','VY','VZ'};
FileFolder = {'3T'}; %'7T',
F = ''; %''; 'Shrink'
InputBaseFolder = 'C:\Users\vhasfckefays\Projects\MF Strength 4DFlow\ProVelMatsPadUnWrp3T';
OutputBaseFolder = 'C:\Users\vhasfckefays\Projects\MF Strength 4DFlow\ProVelMatsPadUnWrp3T';


for i = 1:size(FileFolder,2)
    for k = 1:size(FileSuf,2)
        for j = 1:size(FilePre,2)
            FileName = [FilePre{j} FileSuf{k} F FileFolder{i} 'PadUnWr'];
            load([OutputBaseFolder '\' FileName], FilePre{j});
        end
        Mag = sqrt(VX.^2 + VY.^2 + VZ.^2);
        save([OutputBaseFolder '\' 'Mag' FileSuf{k} F FileFolder{i} 'PadUnWr'], 'Mag');
    end
end
