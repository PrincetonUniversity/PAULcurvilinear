function p2 = plotUncertaintyTube(FinalX_B95,FinalY_B95,FinalZ_B95,color,facealpha,reduceFactor)
% plot uncertainty tube surface using tetramesh
% FinalX_B95 is of #points_along_curve x #curves, organized in order

[MeshIndex1_col,MeshIndex1_row] = meshgrid(1:size(FinalX_B95,2),1:size(FinalX_B95,1)-1);
[MeshIndex2_col,MeshIndex2_row] = meshgrid(1:size(FinalX_B95,2),2:size(FinalX_B95,1));
[MeshIndex3_col,MeshIndex3_row] = meshgrid([2:size(FinalX_B95,2) 1],1:size(FinalX_B95,1)-1);
[MeshIndex4_col,MeshIndex4_row] = meshgrid([2:size(FinalX_B95,2) 1],2:size(FinalX_B95,1));

ind1 = sub2ind(size(FinalX_B95),MeshIndex1_row,MeshIndex1_col);
ind2 = sub2ind(size(FinalX_B95),MeshIndex2_row,MeshIndex2_col);
ind3 = sub2ind(size(FinalX_B95),MeshIndex3_row,MeshIndex3_col);
ind4 = sub2ind(size(FinalX_B95),MeshIndex4_row,MeshIndex4_col);

p = patch('Faces',[ind1(:) ind2(:) ind4(:) ind3(:)],...
     'Vertices',[FinalX_B95(:),FinalY_B95(:),FinalZ_B95(:)],'visible','off');
 
% to save space and time, patches are reduced before plotting
nfv = reducepatch(p,reduceFactor);
p2 = patch('Faces',nfv.faces, 'Vertices',nfv.vertices,...
    'FaceColor',color,'FaceAlpha',facealpha,'EdgeAlpha',0);

