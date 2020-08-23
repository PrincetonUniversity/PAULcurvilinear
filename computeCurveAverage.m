function [FinalX, FinalY_central, FinalY_UB, FinalY_LB]=computeCurveAverage(principalCurve_cell_T,xmin_T,xmax_T,gridPtUnit)
% principalCurve_cell_T is dim x n column vectors

FinalX = xmin_T:gridPtUnit:xmax_T;
curvePosX = cellfun(@(x) x(1,:),principalCurve_cell_T,'UniformOutput',false);
curvePosY = cellfun(@(x) x(2,:),principalCurve_cell_T,'UniformOutput',false);
[uniqueX,ia,~] = cellfun(@(x) unique(x(~isnan(x))),curvePosX,...
    'UniformOutput',false);
Select_Defined_unique = cellfun(@(x,y,ind) [x; y(ind)],uniqueX,curvePosY,ia,...
    'UniformOutput',false);
Select_Defined_unique_noNaN = cellfun(@(x) x(:,~isnan(x(2,:))),...
    Select_Defined_unique,'UniformOutput',false);% remove points were y is NaN
Select_Defined_unique_noNaN = Select_Defined_unique_noNaN(cellfun(@(x) size(x,2)>=2,...
    Select_Defined_unique_noNaN)); % prevent single point error

FinalY_cell = cellfun(@(x) interp1(x(1,:),x(2,:),FinalX,'linear','extrap'),...
    Select_Defined_unique_noNaN,'UniformOutput',false); 
FinalY = cell2mat(FinalY_cell);
FinalY_central = median(FinalY,1);

FinalY_LB = prctile(FinalY,2.5);
FinalY_UB = prctile(FinalY,97.5);
