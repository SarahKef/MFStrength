clear all;
close all;
FOV=[304 304];
FileSuf = {'60100', '6071','6040','404024', '404016', '40408', '40404'};
FilePre = {'VX','VY','VZ'};
FileFolder = {'3T','7T'}; 
InputBaseFolder = 'C:\Users\vhasfckefays\Projects\MF Strength 4DFlow\OrgVelMats';
OutputBaseFolder = 'C:\Users\vhasfckefays\Projects\MF Strength 4DFlow\ProVelMats';

LShift = [0 1 1 1 1 1 1; 0 0 0 0 0 0 0];
SlShft = [0 0 0 10 10 10 10];
for i = 1:size(FileFolder,2)
    for j = 1:size(FilePre,2)
        for k = 1:size(FileSuf,2)
            FileName = [FilePre{j} FileSuf{k}];
            FilePath = [InputBaseFolder '\' FileFolder{i} '\' FilePre{j} '\' FileName]; 
            
            Val=load (FilePath);
            eval(['Val = Val.bc' FilePre{j} ';']);
            T=size(Val,1);
            SL=size(Val,2);
            eval([FilePre{j} '=single(zeros(T, 60, FOV(1), FOV(2)));']);

            for t=1:T
                for sl= 1:SL
                    d = squeeze(Val(t, sl, :, :));
                    d = padarray(d, (FOV - size(d))/2);
                    %d= shiftd(d,0,1,0);
                    d= shiftr(d, 0, LShift(i,k), 0);
                    eval([FilePre{j} '(t, sl + SlShft(k), :, :) = d;']);
                end
            end
            save([OutputBaseFolder '\' FileName FileFolder{i} 'Pad'], FilePre{j});
            clear (FilePre{j});
            fprintf('i=%d, j=%d, k=%d\n', i, j, k);
            
        end
    end
end

for i = 1:size(FileFolder,2)
    for k = 1:size(FileSuf,2)
        for j = 1:size(FilePre,2)
            FileName = [FilePre{j} FileSuf{k}];
            load([OutputBaseFolder '\' FileName FileFolder{i} 'Pad'], FilePre{j});
        end
        Mag = sqrt(VX.^2 + VY.^2 + VZ.^2);
        save([OutputBaseFolder '\' 'Mag' FileSuf{k} FileFolder{i} 'Pad'], 'Mag');
    end
end
