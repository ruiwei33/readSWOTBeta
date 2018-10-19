
%reach
Err.RchHeight.All=CalcErrorStats(Hr,Hrt);
Err.RchSlope.All=CalcErrorStats(Sr,Srt);
Err.RchSlopeEn.All=CalcErrorStats(Sre,Srte);
Err.RchWidth.All=CalcErrorStats(Wr,Wrt);

for i = 1:nPass
    ind.(['Pass' num2str(Pass(i))])=Passes==Pass(i);
end

for  j = 1:nPass
    Err.RchHeight.(['Pass' num2str(Pass(j))])=CalcErrorStats(Hr(:,ind.(...
        ['Pass' num2str(Pass(j))])),Hrt(:,ind.(['Pass' num2str(Pass(j))])));
    
    Err.RchSlope.(['Pass' num2str(Pass(j))])=CalcErrorStats(Sr(:,ind.(...
        ['Pass' num2str(Pass(j))])),Srt(:,ind.(['Pass' num2str(Pass(j))])));
    
    Err.RchSlopeEn.(['Pass' num2str(Pass(j))])=CalcErrorStats(Sre(:,ind.(...
        ['Pass' num2str(Pass(j))])),Srte(:,ind.(['Pass' num2str(Pass(j))])));
    
    Err.RchWidth.(['Pass' num2str(Pass(j))])=CalcErrorStats(Wr(:,ind.(...
        ['Pass' num2str(Pass(j))])),Wrt(:,ind.(['Pass' num2str(Pass(j))])));
end
%node
Err.NodeHeight=CalcErrorStats(Hn,Hnt);
Err.NodeWidth=CalcErrorStats(Wn,Wnt);