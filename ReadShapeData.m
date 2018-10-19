function [SWOTReaches,TrueReaches,SWOTNodes,TrueNodes] = ReadShapeData(Cycles, Passes,Dates,DataDir,River)

for i=1:length(Passes)
    BaseNameReaches=['SWOT_L2_HR_River_SP_' num2str(Cycles(i),'%03.f') '_' ...
        num2str(Passes(i),'%03.f') '_' River '_' num2str(Dates{1,i}) 'T000000_Reach' ];
    SWOTReaches(i).fname=[DataDir  BaseNameReaches '/' BaseNameReaches '.shp'];   
    [SWOTReaches(i).S,SWOTReaches(i).A] = shaperead(SWOTReaches(i).fname);

    BaseNameReaches=['SWOT_L2_HR_River_SP_' num2str(Cycles(i),'%03.f') '_' ...
        num2str(Passes(i),'%03.f') '_' River '_' num2str(Dates{1,i}) 'T000000_Reach_Truth' ];    
    TrueReaches(i).fname=[DataDir  BaseNameReaches '/' BaseNameReaches '.shp'];   
    [TrueReaches(i).S,TrueReaches(i).A] = shaperead(TrueReaches(i).fname);
    
    BaseNameNodes=['SWOT_L2_HR_River_SP_' num2str(Cycles(i),'%03.f') '_' ...
        num2str(Passes(i),'%03.f') '_' River '_' num2str(Dates{1,i}) 'T000000_Node' ];
    SWOTNodes(i).fname=[DataDir  BaseNameNodes '/' BaseNameNodes '.shp'];    
    [SWOTNodes(i).S,SWOTNodes(i).A] = shaperead(SWOTNodes(i).fname);


    BaseNameNodes=['SWOT_L2_HR_River_SP_' num2str(Cycles(i),'%03.f') '_' ...
        num2str(Passes(i),'%03.f') '_' River '_' num2str(Dates{1,i}) 'T000000_Node_Truth' ];
    TrueNodes(i).fname=[DataDir  BaseNameNodes '/' BaseNameNodes '.shp'];    
    [TrueNodes(i).S,TrueNodes(i).A] = shaperead(TrueNodes(i).fname);
    
end