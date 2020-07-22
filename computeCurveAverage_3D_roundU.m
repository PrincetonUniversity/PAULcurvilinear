function [FinalX, FinalY_central, FinalY_B95, FinalZ_central, FinalZ_B95]=computeCurveAverage_3D_roundU(principalCurve_cell_T,xmin_T,xmax_T,gridPtUnit)
% principleCurve_cell_T is dim x n column vectors

FinalX = xmin_T:gridPtUnit:xmax_T;
curvePosX = cellfun(@(x) x(1,:),principalCurve_cell_T,'UniformOutput',false);
curvePosY = cellfun(@(x) x(2,:),principalCurve_cell_T,'UniformOutput',false);
curvePosZ = cellfun(@(x) x(3,:),principalCurve_cell_T,'UniformOutput',false);
[uniqueX,ia,~] = cellfun(@(x) unique(x(~isnan(x))),curvePosX,...
    'UniformOutput',false);
Select_Defined_unique = cellfun(@(x,y,z,ind) [x; y(ind); z(ind)],uniqueX,curvePosY,curvePosZ,ia,...
    'UniformOutput',false);
Select_Defined_unique_noNaN = cellfun(@(x) x(:,~isnan(x(2,:))),...
    Select_Defined_unique,'UniformOutput',false);% remove points where y is NaN
Select_Defined_unique_noNaN = cellfun(@(x) x(:,~isnan(x(3,:))),...
    Select_Defined_unique_noNaN,'UniformOutput',false);% remove points where z is NaN
Select_Defined_unique_noNaN = Select_Defined_unique_noNaN(cellfun(@(x) size(x,2)>=2,...
    Select_Defined_unique_noNaN)); % prevent single point error

FinalY_cell = cellfun(@(x) interp1(x(1,:),x(2,:),FinalX,'linear','extrap'),...
    Select_Defined_unique_noNaN,'UniformOutput',false); 
FinalY = cell2mat(FinalY_cell);
FinalY_central = median(FinalY,1);
FinalZ_cell = cellfun(@(x) interp1(x(1,:),x(3,:),FinalX,'linear','extrap'),...
    Select_Defined_unique_noNaN,'UniformOutput',false); 
FinalZ = cell2mat(FinalZ_cell);
FinalZ_central = median(FinalZ,1);

% go thru a half circle, get uncertainties out
% for any point (y_i,z_i) the projection is (y_i*cosA, z_i*sinA), A \in
% [0,pi]
% Create large 3d matrix --> will not work for feature length>=1000 pixels
A = linspace(-pi,pi,360);
yranges = prctile(FinalY-FinalY_central,97.5,1); % may contain 0
zranges = prctile(FinalZ-FinalZ_central,97.5,1); % may contian 0
FinalY_3D = repmat((FinalY-FinalY_central)./yranges,1,1,length(A));
FinalZ_3D = repmat((FinalZ-FinalZ_central)./zranges,1,1,length(A));
FinalY_3D(isnan(FinalY_3D)) = 0;
FinalZ_3D(isnan(FinalZ_3D)) = 0;
A_3D = repmat(reshape(A,1,1,length(A)),size(FinalY_3D,1),size(FinalY_3D,2));
AllData = (FinalY_3D.*cos(A_3D)+FinalZ_3D.*sin(A_3D)); % projection to the desired dir
AllBounds = squeeze(prctile(AllData,97.5,1));
% % use both 97.5% and 2.5 -- similar
% AllBoundsL = squeeze(prctile(AllData,2.5,1)); 
% AllBounds(:,1:length(A)/2) = -AllBoundsL(:,length(A)/2+1:end);

% transfer the length along (cosA, sinA) back to (y,z) coordinate
A_3D_2 = repmat(A,size(FinalY_3D,2),1);
FinalY_B95 = AllBounds.* cos(A_3D_2).*yranges'+FinalY_central';
FinalZ_B95 = AllBounds.* sin(A_3D_2).*zranges'+FinalZ_central';
