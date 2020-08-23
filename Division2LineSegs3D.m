function [AllPt1x,AllPt2x,AllPt1y,AllPt2y,AllPt1z,AllPt2z]=Division2LineSegs3D(FeatureEnhancedIM,GridSize,indent,GaussFilterSigma)
%% Extend parameters to 3D nonequivalent
if isscalar(GridSize)
    GridSize = [GridSize GridSize size(FeatureEnhancedIM,3)];
elseif length(GridSize)==2
    GridSize = [GridSize(1) GridSize(1) GridSize(2)];
end
if isscalar(GaussFilterSigma)
    GaussFilterSigma = [GaussFilterSigma GaussFilterSigma GaussFilterSigma];
elseif length(GaussFilterSigma)==2
    GaussFilterSigma = [GaussFilterSigma(1) GaussFilterSigma(1) GaussFilterSigma(2)];
end

numSubImxy=size(FeatureEnhancedIM)./GridSize;
%% Shifted divisions (8 cases for 3D)
indent_code = 0;
if indent(1) == 0
    indent_code = indent_code + 1;
end
if indent(2) == 0
    indent_code = indent_code + 2;
end
if indent(3) == 0
    indent_code = indent_code + 5;
end
switch indent_code
    case 8
        dim1Dist=repmat(GridSize(1),1,numSubImxy(1));
        dim2Dist=repmat(GridSize(2),1,numSubImxy(2));
        dim3Dist=repmat(GridSize(3),1,numSubImxy(3));
    case 3
        dim1Dist=repmat(GridSize(1),1,numSubImxy(1));
        dim2Dist=repmat(GridSize(2),1,numSubImxy(2));
        dim3Dist=[indent(3) repmat(GridSize(3),1,numSubImxy(3)-1) GridSize(3)-indent(3)];
    case 6
        dim1Dist=repmat(GridSize(1),1,numSubImxy(1));
        dim2Dist=[indent(2) repmat(GridSize(2),1,numSubImxy(2)-1) GridSize(2)-indent(2)];
        dim3Dist=repmat(GridSize(3),1,numSubImxy(3));
    case 1
        dim1Dist=repmat(GridSize(1),1,numSubImxy(1));
        dim2Dist=[indent(2) repmat(GridSize(2),1,numSubImxy(2)-1) GridSize(2)-indent(2)];
        dim3Dist=[indent(3) repmat(GridSize(3),1,numSubImxy(3)-1) GridSize(3)-indent(3)];
    case 7
        dim1Dist=[indent(1) repmat(GridSize(1),1,numSubImxy(1)-1) GridSize(1)-indent(1)];
        dim2Dist=repmat(GridSize(2),1,numSubImxy(2));
        dim3Dist=repmat(GridSize(3),1,numSubImxy(3));
    case 2
        dim1Dist=[indent(1) repmat(GridSize(1),1,numSubImxy(1)-1) GridSize(1)-indent(1)];
        dim2Dist=repmat(GridSize(2),1,numSubImxy(2));
        dim3Dist=[indent(3) repmat(GridSize(3),1,numSubImxy(3)-1) GridSize(3)-indent(3)];
    case 5
        dim1Dist=[indent(1) repmat(GridSize(1),1,numSubImxy(1)-1) GridSize(1)-indent(1)];
        dim2Dist=[indent(2) repmat(GridSize(2),1,numSubImxy(2)-1) GridSize(2)-indent(2)];
        dim3Dist=repmat(GridSize(3),1,numSubImxy(3));
    otherwise
        dim1Dist=[indent(1) repmat(GridSize(1),1,numSubImxy(1)-1) GridSize(1)-indent(1)];
        dim2Dist=[indent(2) repmat(GridSize(2),1,numSubImxy(2)-1) GridSize(2)-indent(2)];
        dim3Dist=[indent(3) repmat(GridSize(3),1,numSubImxy(3)-1) GridSize(3)-indent(3)];
end
IMraw_split3D = mat2cell(FeatureEnhancedIM,dim1Dist,dim2Dist,dim3Dist);
%% Pre-processing B: Smoothing (gaussfilt on each subimage)
%********** this function imgaussfilt3 might be too slow; worth rewriting
%later...
IM_smoothed=cellfun(@imgaussfilt3,IMraw_split3D,...
    repmat({GaussFilterSigma},size(IMraw_split3D)),'UniformOutput',false);
maxSmoothedVox = max(max(max(cell2mat(IM_smoothed))));
%% Binarize again on each sub (to ensure proper skeletonization)
IM_smoothed = cellfun(@(x) x./maxSmoothedVox,IM_smoothed,'UniformOutput',false);
IM_smoothed_bw=cellfun(@imbinarize2,IM_smoothed,'UniformOutput',false);

volThresh = pi*(max(GaussFilterSigma.*sqrt(2)*1.5))^2*norm(GridSize);
volnum=cellfun(@nnz, IM_smoothed_bw);
noisyind = volnum>=volThresh;
IM_smoothed_bw(noisyind)={false(size(IM_smoothed_bw{1}))};

IM_smoothed_bw_skel=cellfun(@bwskel,IM_smoothed_bw,'UniformOutput',false);
IM_smoothed_bw_skel_cln=cellfun(@(x) bwareaopen(x,7),IM_smoothed_bw_skel,'UniformOutput',false);
%% Extract Straight Lines using Hough Transform
% [~,~,~,~,~,~,AllLinesfromHough]=cellfun(@hough3DLine,IM_smoothed_bw_skel_cln,'UniformOutput',false);

% HTspec = struct()

%%%%%%%%%%%%%%%%% change cellfun to for loop %%%%%%%%%%%%%%%%%%%%%
AllLinesfromHough = cell(size(IM_smoothed_bw_skel_cln));
nonemptyInd = find(cellfun(@(x) any(x,'all'),IM_smoothed_bw_skel_cln));
if ~isempty(nonemptyInd)
    AllLinesfromHough_nonempty = cell(size(nonemptyInd));
    IM_smoothed_bw_skel_cln_nonempty = IM_smoothed_bw_skel_cln(nonemptyInd);
    parfor sm=1:length(nonemptyInd)
        [~,~,~,~,~,~,lines_temp]=hough3DLine(IM_smoothed_bw_skel_cln_nonempty{sm});
        AllLinesfromHough_nonempty(sm) = {lines_temp};
    end
    AllLinesfromHough(nonemptyInd) = AllLinesfromHough_nonempty;
end
%%%%%%%%%%%%%%%%% change cellfun to for loop %%%%%%%%%%%%%%%%%%%%%


nonempty = ~cellfun(@isempty,AllLinesfromHough);
AllPt1s = cell(size(AllLinesfromHough)); 
AllPt1sMat = cell(size(AllLinesfromHough));
AllPt1s_adjusted = cell(size(AllLinesfromHough));
AllPt2s = cell(size(AllLinesfromHough)); 
AllPt2sMat = cell(size(AllLinesfromHough));
AllPt2s_adjusted = cell(size(AllLinesfromHough));

AllPt1s(nonempty) = cellfun(@(x) transpose({x(:).point1}), AllLinesfromHough(nonempty),'UniformOutput',false);
AllPt1sMat(nonempty)=cellfun(@cell2mat,AllPt1s(nonempty),'UniformOutput',false);
AllPt2s(nonempty) = cellfun(@(x) transpose({x(:).point2}), AllLinesfromHough(nonempty),'UniformOutput',false);
AllPt2sMat(nonempty)=cellfun(@cell2mat,AllPt2s(nonempty),'UniformOutput',false);

% move 3D indices to whole-image frame
switch indent_code
    case 8
        xx=0:GridSize(1):size(FeatureEnhancedIM,1)-1;
        yy=0:GridSize(2):size(FeatureEnhancedIM,2)-1;
        zz=0:GridSize(3):size(FeatureEnhancedIM,3)-1;
    case 3
        xx=0:GridSize(1):size(FeatureEnhancedIM,1)-1;
        yy=repmat(GridSize(2),1,numSubImxy(2));
        zz=[indent(3) repmat(GridSize(3),1,numSubImxy(3)-1) GridSize(3)-indent(3)];
    case 6
        xx=0:GridSize(1):size(FeatureEnhancedIM,1)-1;
        yy=[0 indent(2):GridSize(2):(size(FeatureEnhancedIM,2))];
        zz=0:GridSize(3):size(FeatureEnhancedIM,3)-1;
    case 1
        xx=0:GridSize(1):size(FeatureEnhancedIM,1)-1;
        yy=[0 indent(2):GridSize(2):(size(FeatureEnhancedIM,2))];
        zz=[indent(3) repmat(GridSize(3),1,numSubImxy(3)-1) GridSize(3)-indent(3)];
    case 7
        xx=[0 indent(1):GridSize(1):(size(FeatureEnhancedIM,1))];
        yy=0:GridSize(2):size(FeatureEnhancedIM,2)-1;
        zz=0:GridSize(3):size(FeatureEnhancedIM,3)-1;
    case 2
        xx=[0 indent(1):GridSize(1):(size(FeatureEnhancedIM,1))];
        yy=0:GridSize(2):size(FeatureEnhancedIM,2)-1;
        zz=[indent(3) repmat(GridSize(3),1,numSubImxy(3)-1) GridSize(3)-indent(3)];
    case 5
        xx=[0 indent(1):GridSize(1):(size(FeatureEnhancedIM,1))];
        yy=[0 indent(2):GridSize(2):(size(FeatureEnhancedIM,2))];
        zz=0:GridSize(3):size(FeatureEnhancedIM,3)-1;
    otherwise
        xx=[0 indent(1):GridSize(1):(size(FeatureEnhancedIM,1))];
        yy=[0 indent(2):GridSize(2):(size(FeatureEnhancedIM,2))];
        zz=[indent(3) repmat(GridSize(3),1,numSubImxy(3)-1) GridSize(3)-indent(3)];
end

[Y,X,Z] = meshgrid(yy,xx,zz);

IndexCell=JoinIndices(X,Y,Z);
[size1,size2,size3]=cellfun(@size,AllPt1s,'UniformOutput',false);
IndexCell_rep=cellfun(@repmat,IndexCell,size1,size2,size3,'UniformOutput',false);
AllPt1s_adjusted(nonempty)=cellfun(@plus, AllPt1sMat(nonempty), IndexCell_rep(nonempty),'UniformOutput',false);% Add to 2D coords
AllPt2s_adjusted(nonempty)=cellfun(@plus, AllPt2sMat(nonempty), IndexCell_rep(nonempty),'UniformOutput',false);

AllPt1x=cellfun(@(x) x(:,1),AllPt1s_adjusted(nonempty),'UniformOutput',false);
AllPt1x=cell2mat(AllPt1x(:));
AllPt1y=cellfun(@(x) x(:,2),AllPt1s_adjusted(nonempty),'UniformOutput',false);
AllPt1y=cell2mat(AllPt1y(:));
AllPt1z=cellfun(@(x) x(:,3),AllPt1s_adjusted(nonempty),'UniformOutput',false);
AllPt1z=cell2mat(AllPt1z(:));
AllPt2x=cellfun(@(x) x(:,1),AllPt2s_adjusted(nonempty),'UniformOutput',false);
AllPt2x=cell2mat(AllPt2x(:));
AllPt2y=cellfun(@(x) x(:,2),AllPt2s_adjusted(nonempty),'UniformOutput',false);
AllPt2y=cell2mat(AllPt2y(:));
AllPt2z=cellfun(@(x) x(:,3),AllPt2s_adjusted(nonempty),'UniformOutput',false);
AllPt2z=cell2mat(AllPt2z(:));

