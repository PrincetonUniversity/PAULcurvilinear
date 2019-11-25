function lines=HoughTransform2DLine_v2(BW)

[H,T,R] = hough(BW);
P  = houghpeaks(H,3,'threshold',max(H(:)));
if isempty(P)
    lines=[];
    return
else
    lines = houghlines(BW,T,R,P,'FillGap',4,'MinLength',7);
end
for i=1:length(lines)
   lines(i).midpt=0.5*(lines(i).point1+lines(i).point2);
end
