function [AllPt1x,AllPt2x,AllPt1y,AllPt2y,AllAngles]=Division2LineSegs(FeatureEnhancedProj,imDimxy,indent,GaussFilterSigma)
%% Shifted divisions
numSubImxy=size(FeatureEnhancedProj,1)/imDimxy;
if indent == 0
    IMraw_split2D = mat2cell(FeatureEnhancedProj, repmat(imDimxy,1,numSubImxy),...
    repmat(imDimxy,1,numSubImxy));
else
    IMraw_split2D = mat2cell(FeatureEnhancedProj, ...
    [indent repmat(imDimxy,1,numSubImxy-1) imDimxy-indent],...
    [indent repmat(imDimxy,1,numSubImxy-1) imDimxy-indent]);
end
%% Smoothing (gaussfilt on each subimage)
IM_smoothed=cellfun(@imgaussfilt,IMraw_split2D,...
    repmat({GaussFilterSigma},size(IMraw_split2D)),'UniformOutput',false);
%% Binarize again on each sub (to ensure proper skeletonization)
IM_smoothed_bw=cellfun(@imbinarize,IM_smoothed,'UniformOutput',false);
IM_smoothed_bw_skel=cellfun(@bwskel,IM_smoothed_bw,'UniformOutput',false);
%% Extract Straight Lines using Hough Transform
AllLinesfromHough=cellfun(@HoughTransform2DLine_v2,IM_smoothed_bw_skel,'UniformOutput',false);

AllPt1s=cellfun(@(x) transpose({x(:).point1}), AllLinesfromHough,'UniformOutput',false);
AllPt1sMat=cellfun(@cell2mat,AllPt1s,'UniformOutput',false);
AllPt2s=cellfun(@(x) transpose({x(:).point2}), AllLinesfromHough,'UniformOutput',false);
AllPt2sMat=cellfun(@cell2mat,AllPt2s,'UniformOutput',false);
AllThetas=cellfun(@(x) transpose({x(:).theta}), AllLinesfromHough,'UniformOutput',false);
AllThetasMat=cellfun(@cell2mat,AllThetas,'UniformOutput',false);

if indent == 0
    [X,Y]=meshgrid(0:imDimxy:size(FeatureEnhancedProj,1)-1,0:imDimxy:size(FeatureEnhancedProj,2)-1);%create 2D indices
else
    [X,Y]=meshgrid([0 indent:imDimxy:(size(FeatureEnhancedProj,1))],...
    [0 indent:imDimxy:(size(FeatureEnhancedProj,2))]);%create 2D indices
end

IndexCell=repmat(JoinIndices(X,Y),1,1,size(AllPt1s,3));
[size1,size2]=cellfun(@size,AllPt1s,'UniformOutput',false);
IndexCell_rep=cellfun(@repmat,IndexCell,size1,size2,'UniformOutput',false);
AllPt1s_adjusted=cellfun(@plus, AllPt1sMat, IndexCell_rep,'UniformOutput',false);% Add to 2D coords
AllPt2s_adjusted=cellfun(@plus, AllPt2sMat, IndexCell_rep,'UniformOutput',false);

AllPt1x=cellfun(@(x) x(:,1),AllPt1s_adjusted(~cellfun(@isempty,AllPt1s_adjusted)),...
    'UniformOutput',false);AllPt1x=cell2mat(AllPt1x);
AllPt1y=cellfun(@(x) x(:,2),AllPt1s_adjusted(~cellfun(@isempty,AllPt1s_adjusted)),...
    'UniformOutput',false);AllPt1y=cell2mat(AllPt1y);
AllPt2x=cellfun(@(x) x(:,1),AllPt2s_adjusted(~cellfun(@isempty,AllPt2s_adjusted)),...
    'UniformOutput',false);AllPt2x=cell2mat(AllPt2x);
AllPt2y=cellfun(@(x) x(:,2),AllPt2s_adjusted(~cellfun(@isempty,AllPt2s_adjusted)),...
    'UniformOutput',false);AllPt2y=cell2mat(AllPt2y);
AllAngles=AllThetasMat(~cellfun(@isempty,AllThetasMat));AllAngles=cell2mat(AllAngles);
