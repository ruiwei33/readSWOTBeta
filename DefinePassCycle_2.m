NODATA=NaN;
listing=dir(DataDir);
idx=[];
for i = 1:length(listing)
    file = listing(i).name;
    if contains(file, 'Node_Truth')
      idx = [idx,i];
    end
end

for j = 1: length(idx)
    CyclePass=cellfun(@str2num,strrep(extractBetween(listing(idx(j)).name, 'SP_', ['_', River]),'_',' '),'un',0);
    Dates(j)=cellfun(@str2num,extractBetween(listing(idx(j)).name, [River, '_'], 'T'),'un',0);
    Cycles(j) = CyclePass{1,1}(1);
    Passes(j) = CyclePass{1,1}(2);
end

Pass=unique(Passes);
nPass=length(Pass);
nObs=length(Passes);
date_str=cellfun(@mat2str,Dates,'UniformOutput',false);
dates=datetime(date_str, 'InputFormat','yyyyMMdd');