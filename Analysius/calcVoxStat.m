function VoxStat = calcVoxStat(Val, dim)
VoxStat = mean(Val,dim);
VoxStat = cat(dim, VoxStat, std(Val,0,dim));
VoxStat = cat(dim, VoxStat, sqrt(mean(Val.^2,dim)));