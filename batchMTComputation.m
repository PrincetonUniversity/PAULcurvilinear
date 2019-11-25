function batchMTComputation(OutputFileName,image2D,GridSize,NumOfShifts,PSFsigma,pixelSize,varargin)
% batch PAUL calculation for one or more images

beta = 0.5;
cRatio = 2;
l_KR = 20.4;
gridPtUnit = 0.1;

if ~isempty(varargin)
    if length(varargin) >=1 && ~isempty(varargin{1})
        beta = varargin{1};
    end
    if length(varargin) >=2 && ~isempty(varargin{2})
        cRatio = varargin{2};
    end
    if length(varargin) >=3 && ~isempty(varargin{3})
        l_KR = 10^(0.366*log10(varargin{3})-0.0503); 
    end
    if length(varargin) >=4 && ~isempty(varargin{4})
        gridPtUnit = varargin{4};
    end
end

if isnumeric(image2D)
    [FinalX_central_allGps,FinalY_central_allGps,...
        FinalX_UB95_allGps,FinalY_UB95_allGps,FinalX_LB95_allGps,FinalY_LB95_allGps,medium95PercError] = ...
        MTcomputation(image2D,NumOfShifts,GridSize,PSFsigma,pixelSize,beta,cRatio,l_KR,gridPtUnit);
    save([OutputFileName '_Division_' sprintf('%d',GridSize) 'x' sprintf('%d',NumOfShifts)],...
        'FinalX_central_allGps','FinalY_central_allGps','FinalX_UB95_allGps','FinalY_UB95_allGps',...
        'FinalX_LB95_allGps','FinalY_LB95_allGps','medium95PercError');
    disp([OutputFileName ' Division ' sprintf('%d',GridSize) ' x ' sprintf('%d',NumOfShifts) ' PAUL Results Saved'] );   
else
    NumOfImages=length(image2D);
    FinalX_central_allGps_allIMs=cell(NumOfImages,1);
    FinalY_central_allGps_allIMs=cell(NumOfImages,1);
    FinalX_UB95_allGps_allIMs=cell(NumOfImages,1);
    FinalY_UB95_allGps_allIMs=cell(NumOfImages,1);
    FinalX_LB95_allGps_allIMs=cell(NumOfImages,1);
    FinalY_LB95_allGps_allIMs=cell(NumOfImages,1);
    medium95PercError_allIMs=cell(NumOfImages,1);
    
    disp('Image Data Loaded');
    
    display([' Division ' sprintf('%d',GridSize) ' x ' sprintf('%d',NumOfShifts)] );
    
    parfor i=1:NumOfImages
    [FinalX_central_allGps,FinalY_central_allGps,...
        FinalX_UB95_allGps,FinalY_UB95_allGps,FinalX_LB95_allGps,FinalY_LB95_allGps,medium95PercError] = ...
        MTcomputation(image2D{i},NumOfShifts,GridSize,PSFsigma,pixelSize,beta,cRatio,l_KR,gridPtUnit);
        
        FinalX_central_allGps_allIMs(i)={FinalX_central_allGps};
        FinalY_central_allGps_allIMs(i)={FinalY_central_allGps};
        FinalX_UB95_allGps_allIMs(i)={FinalX_UB95_allGps};
        FinalY_UB95_allGps_allIMs(i)={FinalY_UB95_allGps};
        FinalX_LB95_allGps_allIMs(i)={FinalX_LB95_allGps};
        FinalY_LB95_allGps_allIMs(i)={FinalY_LB95_allGps};
        medium95PercError_allIMs(i)={medium95PercError};
        disp(i)
    end
    save([OutputFileName '_Division_' sprintf('%d',GridSize) 'x' sprintf('%d',NumOfShifts)],...
        'FinalX_central_allGps_allIMs','FinalY_central_allGps_allIMs','FinalX_UB95_allGps_allIMs',...
        'FinalY_UB95_allGps_allIMs','FinalX_LB95_allGps_allIMs','FinalY_LB95_allGps_allIMs',...
        'medium95PercError_allIMs');
    disp([OutputFileName ' Division ' sprintf('%d',GridSize) ' x ' sprintf('%d',NumOfShifts) ' PAUL Results Saved'] );
    
end