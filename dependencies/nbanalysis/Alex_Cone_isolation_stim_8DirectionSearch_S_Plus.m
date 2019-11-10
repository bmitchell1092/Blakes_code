clear all
addpath('c:/users/mlab/documents')
%% CREATE CONE ISOLATION STIMULI
%=================================
% In sph2cart(Azimuth, Elevation, radius);

% S+ {Azimuth=90,     Elevation=45,     radius=0 to 0.9   }
% S- {Azimuth=270,    Elevation=-45,    radius=0 to 0.9   }

filename  = 'verboseLog.txt';
fileID    = fopen(filename,'w');
COM       = 'COM4';
% ******************************   S+   *************************************************************************************************
structConeIsolateStim.type = 'S+';

conespectraSampled = Stockman_Sharpe_cone_fundamentals(); % LMS spectra is between 392nm to 780nm in 4nm interval

%% Create centre surround stimuli
widthBackground  = 1200;
heightBackground = 800;
widthForeground  = 400;
heightForeground = 400;
backgroundGrey   = 128;
foregroundCoord_topLeftX = 400;
foregroundCoord_topLeftY = 200;
imageBack=backgroundGrey*ones(heightBackground,widthBackground,3);

%% gamma corect background stimuli

% adjust these values:
RgammaValue = 2.209671458402041;
GgammaValue = 2.50982067675356;
BgammaValue = 1.83058269411577;

imageBack(:,:,1) = ((imageBack(:,:,1))./255) .^ (1/RgammaValue) .*255;
imageBack(:,:,2) = ((imageBack(:,:,2))./255) .^ (1/GgammaValue) .*255;
imageBack(:,:,3) = ((imageBack(:,:,3))./255) .^ (1/BgammaValue) .*255;

structConeIsolateStim.backgroundRGBUncalibrated =[backgroundGrey backgroundGrey backgroundGrey];
structConeIsolateStim.backgroundRGBCalibrated   =[imageBack(1,1,1) imageBack(1,1,2) imageBack(1,1,3)];

figure; imshow(uint8(imageBack))
set(gcf, 'Position', get(1,'Screensize')); % Maximize figure.

%% GET SPECTROMETER READING OF BACKGROUND

obj = serial(COM);
fopen(obj);
spectrum=get_spectrometer_reading(obj);
fclose(obj);
specBackground_380_780_4nm= spectrum(4:end,2);

close all;

% Calculate background LMS response
backgroundL = sum(specBackground_380_780_4nm.*conespectraSampled(:,2));
backgroundM = sum(specBackground_380_780_4nm.*conespectraSampled(:,3));
backgroundS = sum(specBackground_380_780_4nm.*conespectraSampled(:,4));

structConeIsolateStim.backgroundLMSUncalibrated = [backgroundL backgroundM backgroundS];
nbytes = fprintf(fileID,'%1.6f %1.6f %1.6f\n',backgroundL, backgroundM, backgroundS);


%% change these values to obtain any number of cone isolation stimuli
%The following set of values will give you 9 cone isolation values
%If you want only one value increase the 'stimuliStep' to a large value and
%set the 'minRange' to a required value

stimuliStep = 0.1;
minRange    = 0.1;
maxRange    = 0.9;
arrayCoord  = minRange:stimuliStep:maxRange;

%% Steepest descent parameters
threshold       = 0.00002;   % Stopping criterion
maxIteration    = 2;         % max number of iterations
stepReduction   = 0.90;  %   % percentage reduction in step size between iterations (if no successful step is made at an iteration)

% azimuthCurrent=0*pi/180;
% elevationCurrent=45*pi/180;

totNumSearchDir = 9;  % 1=current, 2 to 5 for for directions
% coordLMS=zeros(totNumSearchDir,3);

azimStep        = 10*pi/180; % starting step =5 degrees
elevStep        = 10*pi/180; % starting step =5 degrees

arrayForegroundOptimised=zeros(size(arrayCoord,2),3);

%% steepest descent optimisation algorithm %%

for stimIndex = 1:1:size(arrayCoord,2)
    
    radiusCurrent = arrayCoord(stimIndex);
    
    % **********************   S+   *******************************
    
    % S+ cone isolation stimuli
    azimuthCurrent   = 90*pi/180;
    elevationCurrent = 45*pi/180;
    
    nbytes = fprintf(fileID,'      %1.2f \n', radiusCurrent);
    
    for iteration=1:1:maxIteration
        
        %2D search space (Azimuth and Elevation) But search in four directions
        coordLMS=zeros(totNumSearchDir,3);
        
        for dir=1:1:totNumSearchDir %radius is treated as independent variable(dim(1)=Azimuth and dim(2)=Elevation can vary at each step in)
            
            [stimIndex iteration dir]
            
            azimuthCandidate=azimuthCurrent;
            elevationCandidate=elevationCurrent;
            
            % search in 8 directions
            if dir==1
                % 1 is to calculate the metric at the current position
            end
            if dir==2
                azimuthCandidate=azimuthCurrent+azimStep;
            end
            if dir==3
                azimuthCandidate=azimuthCurrent-azimStep;
            end
            if dir==4
                elevationCandidate=elevationCurrent+elevStep;
            end
            if dir==5
                elevationCandidate=elevationCurrent+elevStep;
            end
            if dir==6
                azimuthCandidate=azimuthCurrent+azimStep;
                elevationCandidate=elevationCurrent+elevStep;
            end
            if dir==7
                azimuthCandidate=azimuthCurrent-azimStep;
                elevationCandidate=elevationCurrent+elevStep;
            end
            if dir==8
                azimuthCandidate=azimuthCurrent-azimStep;
                elevationCandidate=elevationCurrent-elevStep;
            end
            if dir==9
                azimuthCandidate=azimuthCurrent+azimStep;
                elevationCandidate=elevationCurrent-elevStep;
            end
            
            % convert to cartesian space
            [iXCoordinate, iYCoordinate,iZCoordinate ] = sph2cart(azimuthCandidate, elevationCandidate, radiusCurrent);
            
            ldrgyvImg(:,:,1)       = iZCoordinate;
            ldrgyvImg(:,:,2)       = iXCoordinate;
            ldrgyvImg(:,:,3)       = iYCoordinate;
            
            
            imgOutrgb              = ldrgyv2rgb_Monkey_Monitor(ldrgyvImg);
            imgOutrgb              = imgOutrgb.*255;
            
            imageForeground        = ones(heightForeground,widthForeground,3);
            imageBack              = backgroundGrey*ones(heightBackground,widthBackground,3);
            
            imageForeground(:,:,1) = imageForeground(:,:,1).*imgOutrgb(1);
            imageForeground(:,:,2) = imageForeground(:,:,2).*imgOutrgb(2);
            imageForeground(:,:,3) = imageForeground(:,:,3).*imgOutrgb(3);
            
            %Create centre surround display
            imageBack(foregroundCoord_topLeftY:foregroundCoord_topLeftY+heightForeground-1,foregroundCoord_topLeftX:foregroundCoord_topLeftX+widthForeground-1,:)=imageForeground;
            
            rgbForegroundCenterUncalibrated=imageBack(size(imageBack,1)/2,size(imageBack,2)/2,:);
            
            % Gamma correct foreground stimuli
            
            imageBack(:,:,1) = ((imageBack(:,:,1)./255) .^ (1/RgammaValue)).*255;
            imageBack(:,:,2) = ((imageBack(:,:,2)./255) .^ (1/GgammaValue)).*255;
            imageBack(:,:,3) = ((imageBack(:,:,3)./255) .^ (1/BgammaValue)).*255;
            
            figure;
            imshow(uint8(imageBack))
            set(gcf, 'Position', get(1,'Screensize')); % Maximize figure.
            
            % READ FOREGROUND SPECTRUM USING SPECTROMETER
            fopen(obj);
            spectrum=get_spectrometer_reading(obj);      % Get the spectrometer reading of the achromatic background and assign it to the variable 'data'   (B)
            fclose(obj);
            specForeground_380_780_4nm=spectrum(4:size(spectrum,1),2); %spectrometer gives 380nm to 780nm in 4nm interval
            close all;
            
            foregroundL = sum(specForeground_380_780_4nm.*conespectraSampled(:,2));
            foregroundM = sum(specForeground_380_780_4nm.*conespectraSampled(:,3));
            foregroundS = sum(specForeground_380_780_4nm.*conespectraSampled(:,4));
            
            coordLMS(dir,:) = [foregroundL foregroundM foregroundS];
            
        end
        
        % ******************************   S+   *****************************************************************************************************
        % for S+ cone isolation stimuli compare L and M (between foreground and background) *****************************   S+   ********************
        
        dist = sqrt((coordLMS(:,1)-backgroundL).^2+(coordLMS(:,2)-backgroundM).^2);
        
        [minDist, ind] = min(dist);
        
        nbytes = fprintf(fileID,'%1d, %1.8f      %1.4f  %1.6f %1.6f     ',ind, minDist, coordLMS(ind,:));
        
        % Break the loop if the distance between the foreground and the background is smaller than the threshold
        if(minDist<threshold)
            arrayForegroundOptimised(stimIndex,:)=coordLMS(ind,:);
            break;
        end
        
        if (ind==1) % if minimum is at current then reduce step size and repeat the search in all 8 directions in the search space
            azimStep=azimStep*stepReduction;
            elevStep=elevStep*stepReduction;
        else
            bestMatchingTargetLMS=[foregroundL foregroundM foregroundS];
            rgbForegroundCenterUncalibratedBestMatch=rgbForegroundCenterUncalibrated;
            rgbForegroundCenterCalibratedBestMatch=imageBack(size(imageBack,1)/2,size(imageBack,2)/2,:);
            % If there is a direction where the metric is small then move in that direction
            if ind==2
                azimuthCurrent=azimuthCurrent+azimStep;
            end
            if ind==3
                azimuthCurrent=azimuthCurrent-azimStep;
            end
            if ind==4
                elevationCurrent=elevationCurrent+elevStep;
            end
            if ind==5
                elevationCurrent=elevationCurrent-elevStep;
            end
            if ind==6
                azimuthCurrent=azimuthCurrent+azimStep;
                elevationCurrent=elevationCurrent+elevStep;
            end
            if ind==7
                azimuthCurrent=azimuthCurrent-azimStep;
                elevationCurrent=elevationCurrent+elevStep;
            end
            if ind==8
                azimuthCurrent=azimuthCurrent-azimStep;
                elevationCurrent=elevationCurrent-elevStep;
            end
            if ind==9
                azimuthCurrent=azimuthCurrent+azimStep;
                elevationCurrent=elevationCurrent-elevStep;
            end
            
        end
        
        nbytes = fprintf(fileID,'      %1.6f %1.6f       %1.6f %1.6f\n', azimuthCurrent, elevationCurrent, azimStep, elevStep);
        
        % At the end of max iteration
        if(iteration==maxIteration)
            arrayForegroundOptimised(stimIndex,:)=coordLMS(ind,:);
        end
        
    end
    
    structConeIsolateStim(stimIndex).foregroundDKLPolarCoordUncalibrated=[azimuthCurrent elevationCurrent radiusCurrent];
    structConeIsolateStim(stimIndex).foregroundLMSUncalibrated=bestMatchingTargetLMS;
    
    structConeIsolateStim(stimIndex).foregroundRGBUnCalibrated=rgbForegroundCenterUncalibratedBestMatch;
    structConeIsolateStim(stimIndex).foregroundRGBCalibrated=rgbForegroundCenterCalibratedBestMatch;
    
    nbytes = fprintf(fileID,'\n\n\n');
end

fclose(fileID);

% ******************************   S+   ****************************************************************************************************
%S+ cone isolation stimuli
save('structConeIsolateStim_S_Plus.mat','structConeIsolateStim');





