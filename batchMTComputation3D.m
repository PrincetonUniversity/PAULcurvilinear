function batchMTComputation3D(OutputFileName,image3D,GridSize,NumOfShifts,PSFsigma,pixelSize,varargin)
% batch PAUL calculation for one or more images

alpha = 0.5;
beta = 0.5;
cRatio = 5;
l_KR = 20.4;
gridPtUnit = 0.1;
num_seg_thresh = 5;

if ~isempty(varargin)
    if length(varargin) >=1 && ~isempty(varargin{1})
        alpha = varargin{1};
    end
    if length(varargin) >=2 && ~isempty(varargin{2})
        beta = varargin{2};
    end
    if length(varargin) >=3 && ~isempty(varargin{3})
        cRatio = varargin{3};
    end
    if length(varargin) >=4 && ~isempty(varargin{4})
        l_KR = 10^(0.366*log10(varargin{4})-0.0503); 
    end
    if length(varargin) >=5 && ~isempty(varargin{5})
        gridPtUnit = varargin{5};
    end
    if length(varargin) >=6 && ~isempty(varargin{6})
        num_seg_thresh = varargin{6};
    end
end

if isnumeric(image3D)
    [FinalX_central_allGps,FinalY_central_allGps,FinalZ_central_allGps,...
        FinalX_B95_allGps, FinalY_B95_allGps, FinalZ_B95_allGps,medium95PercError]=...
        MTcomputation3D(image3D,NumOfShifts,GridSize,PSFsigma,pixelSize,alpha,beta,cRatio,l_KR,gridPtUnit,num_seg_thresh);
    save([OutputFileName '_Division_' sprintf('%d',GridSize) 'x' sprintf('%d',NumOfShifts)],...
        'FinalX_central_allGps','FinalY_central_allGps','FinalZ_central_allGps',...
        'FinalX_B95_allGps','FinalY_B95_allGps','FinalZ_B95_allGps','medium95PercError');
    disp([OutputFileName ' Division ' sprintf('%d',GridSize) ' x ' sprintf('%d',NumOfShifts) ' PAUL Results Saved'] );
else
    NumOfImages=length(image3D);
    FinalX_central_allGps_allIMs=cell(NumOfImages,1);
    FinalY_central_allGps_allIMs=cell(NumOfImages,1);
    FinalZ_central_allGps_allIMs=cell(NumOfImages,1);
    FinalX_B95_allGps_allIMs=cell(NumOfImages,1);
    FinalY_B95_allGps_allIMs=cell(NumOfImages,1);
    FinalZ_B95_allGps_allIMs=cell(NumOfImages,1);
    medium95PercError_allIMs=cell(NumOfImages,1);
    
    disp('Image Data Loaded');
    
    display([' Division ' sprintf('%d',GridSize) ' x ' sprintf('%d',NumOfShifts)] );
    
    for i=1:NumOfImages
    [FinalX_central_allGps,FinalY_central_allGps,FinalZ_central_allGps,...
        FinalX_B95_allGps, FinalY_B95_allGps, FinalZ_B95_allGps,medium95PercError]=...
        MTcomputation3D(image3D{i},NumOfShifts,GridSize,PSFsigma,pixelSize,alpha,beta,cRatio,l_KR,gridPtUnit,num_seg_thresh);
        
        FinalX_central_allGps_allIMs(i)={FinalX_central_allGps};
        FinalY_central_allGps_allIMs(i)={FinalY_central_allGps};
        FinalZ_central_allGps_allIMs(i)={FinalZ_central_allGps};
        FinalX_B95_allGps_allIMs(i)={FinalX_B95_allGps};
        FinalY_B95_allGps_allIMs(i)={FinalY_B95_allGps};
        FinalZ_B95_allGps_allIMs(i)={FinalZ_B95_allGps};
        medium95PercError_allIMs(i)={medium95PercError};
        disp(i)
    end
    save([OutputFileName '_Division_' sprintf('%d',GridSize) 'x' sprintf('%d',NumOfShifts)],...
        'FinalX_central_allGps_allIMs','FinalY_central_allGps_allIMs','FinalZ_central_allGps_allIMs',...
        'FinalX_B95_allGps_allIMs','FinalY_B95_allGps_allIMs','FinalZ_B95_allGps_allIMs',...
        'medium95PercError_allIMs');
    disp([OutputFileName ' Division ' sprintf('%d',GridSize) ' x ' sprintf('%d',NumOfShifts) ' PAUL Results Saved'] );
    
end