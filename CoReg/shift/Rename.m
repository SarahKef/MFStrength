clear all;
close all;

FileSuf = {'60100', '6071','6040','404024', '404016', '40408', '40404'};
FilePre = {'VX','VY','VZ'};
FileFolder = {'7T'};

InputBaseFolder = 'C:\Users\vhasfckefays\Projects\MF Strength 4DFlow\ProVelMatsPadUnWrpShrink7T';
OutputBaseFolder = 'C:\Users\vhasfckefays\Projects\MF Strength 4DFlow\ProVelMatsPadUnWrp';

for i = 1:size(FileFolder,2)
    for k = 1:size(FileSuf,2)
        for j = 1:size(FilePre,2)
            FileName = [FilePre{j} FileSuf{k} 'Shrink7TPadUnWr'];
            load([InputBaseFolder '\' FileName], FilePre{j});
            save([OutputBaseFolder '\' FilePre{j} FileSuf{k} FileFolder{i} 'PadUnWr'], FilePre{j});
        end
    end
end