function ResultCell=JoinIndices(varargin)
% put grid indices (usually from meshgrid) into cells

X=varargin{1};Xcell=num2cell(X);
Y=varargin{2};Ycell=num2cell(Y);

if length(varargin)==3
    Z=varargin{3};Zcell=num2cell(Z);
    ResultCell=cellfun(@horzcat,Xcell,Ycell,Zcell,'UniformOutput',false);
else
    ResultCell=cellfun(@horzcat,Xcell,Ycell,'UniformOutput',false);
end