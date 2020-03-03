# Prior-Apprised Unsupervised Learning (PAUL) of curvilinear features
Implementation of PAUL framwork on curvilinear features.

## Requirements
The code uses the following toolboxes and versions at the stage of development; earlier versions might also work:
* 'MATLAB' '9.6'
* 'Image Processing Toolbox' '10.4' 
* 'Statistics and Machine Learning Toolbox' '11.5' 
* 'Parallel Computing Toolbox' '7.0' 
* 'MATLAB Parallel Server' '7.0'

## Syntax
```
batchMTComputation(OutputFileName, image2D, GridSize, NumOfShifts, PSFsigma, pixelSize);
batchMTComputation(OutputFileName, image2D, GridSize, NumOfShifts, PSFsigma, pixelSize, beta);
batchMTComputation(OutputFileName, image2D, GridSize, NumOfShifts, PSFsigma, pixelSize, beta, cRatio);
batchMTComputation(OutputFileName, image2D, GridSize, NumOfShifts, PSFsigma, pixelSize, beta, cRatio, estimatedLp);
batchMTComputation(OutputFileName, image2D, GridSize, NumOfShifts, PSFsigma, pixelSize, beta, cRatio, estimatedLp, gridPtUnit);
```

## Input arguments
* `OutputFileName` — the string that contains the output file name.
* `image2D` — numeric array of a single input image or cells of numeric arrays of multiple images. The input image needs to be **square** (e.g. 512 x 512 pixels). The pixel values can be 8-bit, 16-bit or double. The image can be grayscale or RGB (will be converted into grayscale in the code).
* `GridSize` — the side of the sub image used in stage two. It should divide the image side (e.g. for 512 x 512 images, `GridSize` cannot fall out of {1, 2, 4, 8, 16, 32, 64, 128, 256, 512}).
* `NumOfShifts` — the number of shifts applied in stage two.
* `PSFsigma` — the standard deviation of the approximated Gaussian point-spread function (PSF), in the unit of micrometers.
* `pixelSize` — the side of the pixel, in the unit of micrometers.
* `beta` — the beta parameter in stage one. Default is 0.5. Put `[]` if using the default value.
* `cRatio` — the ratio between 2c^2 and maximum S^2 in stage one. Default is 2 (c=max(S)). Put `[]` if using the default value.
* `estimatedLp` — the estimated persistence length of the features to be detected, in the unit of micrometers. Default is 5200. Put `[]` if using the default value.
* `gridPtUnit` — the grid point density in the final PAUL principal curve, in the unit of pixels. Default is 0.1.


## Output
The output is saved in a .mat file containing the following variables:
* `FinalX_central_allGps` (or `FinalX_central_allGps_allIMs`), `FinalY_central_allGps` (or `FinalY_central_allGps_allIMs`) — x and y coordinates of the final PAUL principal curve, in the pixel space.
* `FinalX_UB95_allGps` (or `FinalX_UB95_allGps_allIMs`), `FinalY_UB95_allGps` (or `FinalY_UB95_allGps_allIMs`) — x and y coordinates of the upper bound of the final PAUL principal curve, in the pixel space.
* `FinalX_LB95_allGps` (or `FinalX_LB95_allGps_allIMs`), `FinalY_LB95_allGps` (or `FinalY_LB95_allGps_allIMs`) — x and y coordinates of the lower bound of the final PAUL principal curve, in the pixel space.
* `medium95PercError` (or `medium95PercError_allIMs`) — medians of all 95% bootstrap bounds, in the unit of pixels.

## Example

Run PAUL procedure on an example image (/examples/example_image.mat)
```matlab
>> batchMTComputation('example_result',example_image,32,8,0.1376,0.076,[],5,1000);
```

Plot the result
```matlab
load('/examples/example_image.mat');
load('/examples/example_result_Division_32x8.mat');
imagesc(example_image)
hold on
for i = 1:length(FinalX_central_allGps)
    plot(FinalX_central_allGps{i},FinalY_central_allGps{i},'Color',[0.3010, 0.7450, 0.9330],'LineWidth',0.6);
    plot(FinalX_UB95_allGps{i},FinalY_UB95_allGps{i},'--','Color',[0.9290, 0.6940, 0.1250],'LineWidth',0.5);
    plot(FinalX_LB95_allGps{i},FinalY_LB95_allGps{i},'--','Color',[0.9290, 0.6940, 0.1250],'LineWidth',0.5);
end
```
![example_image](example/example_image.png)
![example_result](example/example_image_result.png)


