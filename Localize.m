function localizationAll = Localize(imageFrame,threshold,rbox,rollingBallRadius,minAreaPixels,pixelSize)
% Localize the single molecule in the SINGLE frame "imageFrame" according to a threshold and a define box dimension
% Generally the input values are:
% rollingBallRadius = 5; % Rolling ball background removal
% threshold = 75;
% rbox = 4;
% minAreaPixels = 1; %9;
% pixelSize = 0.065;
% it returns the localization both in pixel and in um and the localization
% precsion (the veridicity of this last parametrs needs to be checked more...)

% initializing the rolling ball...
rollingBallDimension = 2 * rollingBallRadius - 1;
rollingBall = zeros (rollingBallDimension, rollingBallDimension);
for y = 1 : rollingBallDimension
    for x = 1 : rollingBallDimension
        if ( ( ( x - rollingBallRadius )^2 + ( y - rollingBallRadius )^2 ) < rollingBallRadius^2 )
            rollingBall(y,x) = 1;
        end
    end
end
rollingBall = rollingBall / sum(sum(rollingBall));

% inizializing variables
nLoc = 0;
errors = 0;
discarded = 0;

backgroundEstimate = conv2(imageFrame, rollingBall, 'same');
filteredImage = imageFrame - backgroundEstimate;
filteredImage = filteredImage .* (filteredImage>0);
% imtool(filteredImage);
% imwrite(uint16(filteredImage),'threshold.tif')
thresholdedImage = filteredImage .* (filteredImage >= threshold);
% imtool(thresholdedImage);

[height, width] = size (filteredImage);

while (1)
    %         clear('rowMax');
    %         clear('yI');
    %         clear('y');
    %         clear('x');
    %         clear('maxValue');
    [rowMax, yI] = max(thresholdedImage);
    [maxValue, x] = max(rowMax);
    y=yI(x);
    
    if (maxValue<=0)
        break;
    end
    
    if ( x <= rbox || x >= (width-rbox+1) || ...
            y <= rbox || y >= (height-rbox+1) )
        thresholdedImage(y,x) = 0;
        continue;
    end
    assert(x > rbox );
    assert(y > rbox );
    assert(x <= width-rbox);
    assert(y <= height-rbox);
    
    dataBox = filteredImage( (y-rbox):(y+rbox), (x-rbox):(x+rbox) );
    
    maxBox = thresholdedImage( (y-rbox):(y+rbox), (x-rbox):(x+rbox) );
    area = sum ( sum ( maxBox > 0 ) );
    
    thresholdedImage(y,x) = 0;
    
    if (area < minAreaPixels)
        continue;
    end
    
    modelParameters = IIT_G2DFit_gaussian2DFittingSupervised ( dataBox, 20, pixelSize, ...
        rbox*pixelSize, rbox*pixelSize, 1*pixelSize );
    modelParameters.kValue = pi*modelParameters.sx*modelParameters.sy*modelParameters.peak/(pixelSize.^2);
    
    %         newDataBox = imageFrame( (y-rbox):(y+rbox), (x-rbox):(x+rbox) );
    %         modelParameters = IIT_G2DFit_gaussian2DFittingSupervised ( newDataBox, 20, pixelSize, ...
    %                                                                    rbox*pixelSize, rbox*pixelSize, 1*pixelSize );
    %         modelParameters.kValue = 1500;
    
    if (modelParameters.sx == -1)
        errors = errors + 1;
        %ppp(:,:,errors) = dataBox;
        continue;
    end
    
    %         if (modelParameters.ux > 2*rbox+1 || modelParameters.uy > 2*rbox+1)
    %             discarded = discarded + 1;
    %             continue;
    %         end
    if ( modelParameters.ux <= 0 || modelParameters.ux >= (2*rbox+1)*pixelSize || ....
            modelParameters.uy <= 0 || modelParameters.uy >= (2*rbox+1)*pixelSize )
        discarded = discarded + 1;
        continue;
    end
    
    xM = ( x - rbox - 1 ) * pixelSize + modelParameters.ux;
    yM = ( y - rbox - 1 ) * pixelSize + modelParameters.uy;
    
    % 'bleaching' of the area where fluorophore was found and characterized
    sxPixels = ceil ( modelParameters.fwhmX / pixelSize );
    syPixels = ceil ( modelParameters.fwhmY / pixelSize );
    xMPixels = xM / pixelSize + 1;
    yMPixels = yM / pixelSize + 1;
    startingX = floor( max( ( xMPixels - sxPixels ), 1 ) );
    endingX = ceil( min( ( xMPixels + sxPixels ), width ) );
    startingY = floor( max( ( yMPixels - syPixels ), 1 ) );
    endingY = ceil( min( ( yMPixels + syPixels ), height ) );
    for yd = startingY : endingY
        for xd = startingX : endingX
            if ( ( ( ( xd - xMPixels )^2) / (sxPixels^2) + ( ( yd - yMPixels )^2) / (syPixels^2) ) <= 1 )
                thresholdedImage( yd , xd ) = 0;
            end
        end
    end
    
    %         imtool(thresholdedImage);
    
    bg = 12;
    lpX = 1.3 * sqrt( ( ((modelParameters.sx*pixelSize/2.0)^2 + pixelSize^2/12.0) / modelParameters.kValue) + ....
        ( 8 * pi * ( modelParameters.sx*pixelSize/2.0 )^4 * bg^2 / (pixelSize * modelParameters.kValue)^2 )  );
    lpY = 1.3 * sqrt( ( ((modelParameters.sy*pixelSize/2.0)^2 + pixelSize^2/12.0) / modelParameters.kValue) + ....
        ( 8 * pi * ( modelParameters.sy*pixelSize/2.0 )^4 * bg^2 / (pixelSize * modelParameters.kValue)^2 )  );
    
    nLoc = nLoc + 1;
    
    localizationAll(nLoc,:) = [xMPixels, yMPixels, xM, yM, struct2array(modelParameters), lpX, lpY];
    
    
end