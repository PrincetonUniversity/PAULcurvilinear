function Area=calculatePolygonArea(Seg1_1x,Seg1_1y,Seg1_2x,Seg1_2y,Seg2_1x,Seg2_1y,Seg2_2x,Seg2_2y)

x=[Seg1_1x;Seg2_1x;Seg1_2x;Seg2_2x];
y=[Seg1_1y;Seg2_1y;Seg1_2y;Seg2_2y];
P=[x y];
CollBool = rank(bsxfun(@minus, P, P(1,:))) < 2; %sort points in ccw

if CollBool
    V=0;
else
    [~,V] = convhull(x,y);
end
Area=V;