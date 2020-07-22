function dists=DistanceBtPoints(Pts1,Pts2)
% distance between each pair of points
% Pt1 and Pt2 are nx2 or nx3 matrices    
if size(Pts1,2) == 2
    dists = sqrt((Pts1(:,1)-Pts2(:,1)).^2+(Pts1(:,2)-Pts2(:,2)).^2);
elseif size(Pts1,2) == 3
    dists = sqrt((Pts1(:,1)-Pts2(:,1)).^2+(Pts1(:,2)-Pts2(:,2)).^2+(Pts1(:,3)-Pts2(:,3)).^2);
end
        