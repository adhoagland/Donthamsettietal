clear; close all

% Parameters to adjust

amtOfMovementToInitateBout = 1; % in distance units; the threshold amount of movement to initate bout
minDurationOfBout = 0.5; % in seconds; only count bouts that are longer than this amount
boutStopThresh = 50; % in seconds; the amount of time that has to pass before we call the bout finished 
smoothValue = 100; % 1 means no smoothing; usually 2-10 range works best, this represents the span or moving window range (lookup simple moving average formula)
errorValueThreshold = 50; % threshold to get rid of too high distance data that is incorrectly tracked
boutLengthToAverage = 1; % min length in seconds to compile the bout to average it - this doesn't do anything yet as of 3/31/20
seePlotsAsYouGoThroughData = 1; % 1 if you want to see the distance plots and how the program identifies bouts per mouse, 0 if you don't want to see graphs

redoMovementThresholds = 1; % 1 if you want to go through and select "amtOfMovementToInitateBout" for each mouse, 0 otherwise. If the variable doesn't already exist in the current directory then it will create it regardless

uvInds = 1; % indices of mice that are UV group 
blueInds = 2; % indices of mice that are blue group

% Don't need to understand anything below this point


disp('Select your excel data file')
[files, fileLocation]=uigetfile('.xlsx','Select a File','MultiSelect','on');
cd(fileLocation)

[behavData, text, raw]  = xlsread(files);

mouseNameIndices=zeros(size(text,2),1);
ledDelimiter=zeros(size(text,2),1);
mouseNames=cell(size(text,2),1);
ledWavelength=cell(size(text,2),1);

for ii = 1:size(text,2)
    mouseNameIndices(ii,1)=~isempty(text{2,ii});
    mouseNames{ii,1}=(text{2,ii});
    ledDelimiter(ii,1)=~isempty(text{1,ii});
    ledWavelength{ii,1}=(text{1,ii});
end

blueInd = find(ledDelimiter, 1, 'last' );
ledDelimiter(blueInd:end)=2;
ledDelimiter(1:blueInd-1)=1;

mouseNameLinInds = find(mouseNameIndices);

%% 

mouseThresholds=zeros(numel(mouseNameLinInds),1);

clear smoothedDistanceTraces boutTracesCompiled boutTracesIndsCompiled boutTracesOnDuration boutTracesOffDuration frameExposures percentMovingtime
for mouseNum = 1:numel(mouseNameLinInds)
    numel(mouseNameLinInds)
    mouseInd = mouseNameLinInds(mouseNum);

    disp(['You are currently on mouse ' text{2,mouseInd}])
    
    timepoints = behavData(:,mouseInd);
    distance = behavData(:,mouseInd+1);
    distance(isnan(distance))=[];
    
    trackErrorInds = find(distance>errorValueThreshold);
    if ~isempty(trackErrorInds)
        for errNum = 1:length(trackErrorInds)
            if trackErrorInds(errNum)==1
                distance(1)=distance(2); % if it's the first distance measurement that's wrong, take the second measurement
            else % otherwise take the previous distance measurement
                distance(trackErrorInds(errNum))=distance(trackErrorInds(errNum)-1);
            end
        end
    end
    
    distance = smooth(distance,7);

    if redoMovementThresholds==1 || ~exist('movementThresholds.mat')
    plot(distance)
    disp('Select movement threshold')
    [~,amtOfMovementToInitateBout,~] = ginput(1);
    close all
    mouseThresholds(mouseNum)=amtOfMovementToInitateBout;
    else
        if mouseNum==1
        load('movementThresholds.mat')
        end
    amtOfMovementToInitateBout=mouseThresholds(mouseNum);
    end

    frameExposureTime = mode(diff(timepoints));

    indicesAboveMovementThresh = distance>amtOfMovementToInitateBout; % change '>' to '>=' for greater than or equal to

    regionStats = regionprops(indicesAboveMovementThresh,'Area','PixelIdxList');

    regionLengths = [regionStats(:).Area];
    minDurationOfBoutInNumberOfFrames = round(minDurationOfBout/frameExposureTime);
    regionsTooShort = regionLengths<minDurationOfBoutInNumberOfFrames;
    regionStats(regionsTooShort)=[];

    revisedIndices = zeros(numel(distance),1);

    for regionNum = 1:numel(regionStats)
        revisedIndices(regionStats(regionNum).PixelIdxList)=1;
    end

    revisedIndicesZeroLocations = ~revisedIndices;
    zeroSpans = regionprops(revisedIndicesZeroLocations,'Area','PixelIdxList');
    zeroSpanLengths = [zeroSpans(:).Area];
    minOffDurationInNumberOfFrames = round(boutStopThresh/frameExposureTime);

    bouts2merge = find(zeroSpanLengths<minOffDurationInNumberOfFrames);

    finalBoutIndices = revisedIndices;
    for mergeRegNum = bouts2merge
        theseIndices = zeroSpans(mergeRegNum).PixelIdxList;
        revisedIndices(theseIndices)=1;
    end

    revisedProps = regionprops(logical(revisedIndices),'Area','PixelIdxList');
    finalBoutProps = regionprops(logical(finalBoutIndices),'Area','PixelIdxList');
    finalBoutPropsOffState = regionprops(~logical(finalBoutIndices),'Area','PixelIdxList');

    percentOnTime = mean(logical(finalBoutIndices));
    
    subplot(2,1,1)
    plot(distance,'m')
    xticks([round(15/frameExposureTime) round(30/frameExposureTime) round(45/frameExposureTime) round(60/frameExposureTime) round(75/frameExposureTime) round(90/frameExposureTime) round(105/frameExposureTime) round(120/frameExposureTime)])
    xticklabels({'15','30','45','60','75','90','105','120'})
    xlim([0 120/frameExposureTime])
%     ylim([0 15])
    ylabel('Distance')

    hold on

%     maxDist = max(distance)*1.02;
%     for lineSpanNum = 1:numel(revisedProps)
%         thisBoutLine = revisedProps(lineSpanNum).PixelIdxList;
%         yCoords = repmat(maxDist,numel(thisBoutLine),1);
%         plot(thisBoutLine,yCoords,'LineWidth',5)
%     end

    boutLengthToAverageInFrames = round(boutLengthToAverage/frameExposureTime);
    bouts2average = cell(numel(finalBoutProps),1);
    bouts2averageIndices = cell(numel(finalBoutProps),1);
    
    maxDist2 = max(distance)*1.07;
    for lineSpanNum = 1:numel(finalBoutProps)
        thisBoutLine = finalBoutProps(lineSpanNum).PixelIdxList;
        yCoords = repmat(maxDist2,numel(thisBoutLine),1);
        plot(thisBoutLine,yCoords,'LineWidth',5)

        thisBoutTrace = distance(thisBoutLine);
        bouts2average{lineSpanNum,1}=thisBoutTrace;
        bouts2averageIndices{lineSpanNum,1}=thisBoutLine;

    end
    minBoutLength =  round(boutLengthToAverage/frameExposureTime);
    
    boutCats = [];
    
    for cellNum = 1:numel(bouts2average)
        thisBout = bouts2average{cellNum};
        if bouts2averageIndices{lineSpanNum,1}(1)==1
            continue
        else
            if numel(thisBout) >= minBoutLength
                thisBoutApp = thisBout(1:minBoutLength)';
                boutCats=[boutCats; thisBoutApp];       
            end
        end
    end

    distanceThresholded = distance;
    distanceThresholded(~finalBoutIndices)=0;
    
    boutTracesCompiled{mouseNum,1}=boutCats;
    boutTracesOnDuration{mouseNum,1}=[finalBoutProps(:).Area]*frameExposureTime;
    boutTracesOffDuration{mouseNum,1}=[finalBoutPropsOffState(:).Area]*frameExposureTime;
    frameExposures{mouseNum,1}=frameExposureTime;
    percentMovingtime{mouseNum,1}=percentOnTime;
    smoothedDistanceTraces{mouseNum,1}=distanceThresholded;
    mouseIDs{mouseNum,1}=text{2,mouseInd};
    
    if seePlotsAsYouGoThroughData == 1
        disp('Press any key to advance to next mouse')
        waitforbuttonpress
    end
    close 
    
    
    
end

save('movementThresholds.mat','mouseThresholds');
save('compiledData.mat','mouseIDs','smoothedDistanceTraces','boutTracesCompiled','boutTracesOnDuration','boutTracesOffDuration','frameExposures','percentMovingtime');

%%%%
%% 

clear maxAmpsUV onDursUV offDursUV perActTimeUV nBoutsUV
for uvii = uvInds
    
    clear maxAmplitude
    theseBoutTraces = boutTracesCompiled{uvii};
    theseOnDurs = boutTracesOnDuration{uvii};
    theseOffDurs = boutTracesOffDuration{uvii};
    for subInd = 1:size(theseBoutTraces,1)
    maxAmplitude(subInd,1) = max(theseBoutTraces(subInd,:));
 
    end
    maxAmpsUV{uvii,1}=(maxAmplitude);
    onDursUV(uvii,1) = mean(theseOnDurs);
    offDursUV(uvii,1) = mean(theseOffDurs);
    perActTimeUV(uvii,1) = percentMovingtime{uvii};
    nBoutsUV(uvii,1)=size(theseBoutTraces,1);
end

clear maxAmpsBlue onDursBlue offDursBlue perAcperActTimeBluetTime nBoutsBlue
blueCnt= 1;
for blueii = blueInds
    clear maxAmplitude
    theseBoutTraces = boutTracesCompiled{blueii};
    theseOnDurs = boutTracesOnDuration{blueii};
    theseOffDurs = boutTracesOffDuration{blueii};
    for subInd = 1:size(theseBoutTraces,1)
    maxAmplitude(subInd,1) = max(theseBoutTraces(subInd,:));
    end
    maxAmpsBlue{blueCnt,1}=(maxAmplitude);
    onDursBlue(blueCnt,1) = mean(theseOnDurs);
    offDursBlue(blueCnt,1) = mean(theseOffDurs);
    perActTimeBlue(blueCnt,1) = percentMovingtime{blueii};
    nBoutsBlue(blueCnt,1)=size(theseBoutTraces,1);

    blueCnt = blueCnt+1;

end

meanAmpUV = mean(maxAmpsUV{:});
semAmpUV = std(maxAmpsUV{:}/sqrt(length(maxAmpsUV{:})));
meanAmpBlue = mean(maxAmpsBlue{:});
semAmpBlue = std(maxAmpsBlue{:}/sqrt(length(maxAmpsBlue{:})));

cc = linspecer(length(maxAmpsUV));
figure
hold on
superbar([meanAmpUV meanAmpBlue],'E',[semAmpUV semAmpBlue],'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',35,'BarWidth',0.4)
scatter(repmat(1.35,length(maxAmpsUV{:}),1),maxAmpsUV{:},250,'.','jitter','on', 'jitterAmount',0.1);
scatter(repmat(2.35,length(maxAmpsBlue{:}),1),maxAmpsBlue{:},250,'.','jitter','on', 'jitterAmount',.1);
set(gca,'fontsize',20)

set(gca,'XTick',[1 2])
xticklabels({'UV','Blue'})    
ylabel('Max Distance per Bout (cm)')

savefig('Max Dist per Bout.fig')
saveas(gcf,'Max Dist per Bout.tif')
saveas(gcf,'Max Dist per Bout.eps')

save('compiledData.mat','mouseIDs','maxAmpsUV','maxAmpsBlue','smoothedDistanceTraces','boutTracesCompiled','boutTracesOnDuration','boutTracesOffDuration','frameExposures','percentMovingtime');

%%%%% Bout duration

meanAmpUV = mean(onDursUV);
semAmpUV = std(onDursUV/sqrt(length(onDursUV)));
meanAmpBlue = mean(onDursBlue);
semAmpBlue = std(onDursBlue/sqrt(length(onDursBlue)));
pvalsAmp = ranksum(onDursUV,onDursBlue);

cc = linspecer(length(onDursUV));
figure
hold on
superbar([meanAmpUV meanAmpBlue],'E',[semAmpUV semAmpBlue],'P',[1 pvalsAmp],'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',35,'BarWidth',0.4)
scatter(repmat(1.35,length(onDursUV),1),onDursUV,250,'.','jitter','on', 'jitterAmount',0.1);
scatter(repmat(2.35,length(onDursBlue),1),onDursBlue,250,'.','jitter','on', 'jitterAmount',.1);
set(gca,'fontsize',20)

set(gca,'XTick',[1 2])
xticklabels({'UV','Blue'})    
ylabel('Bout duration (s)')

savefig('Bout duration.fig')
saveas(gcf,'Bout duration (s).tif')
saveas(gcf,'Bout duration (s).eps')

%%%%% interbout interval

meanAmpUV = mean(offDursUV);
semAmpUV = std(offDursUV/sqrt(length(offDursUV)));
meanAmpBlue = mean(offDursBlue);
semAmpBlue = std(offDursBlue/sqrt(length(offDursBlue)));
pvalsAmp = ranksum(offDursUV,offDursBlue);

cc = linspecer(length(offDursUV));
figure
hold on
superbar([meanAmpUV meanAmpBlue],'E',[semAmpUV semAmpBlue],'P',[1 pvalsAmp],'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',35,'BarWidth',0.4)
scatter(repmat(1.35,length(offDursUV),1),offDursUV,250,'.','jitter','on', 'jitterAmount',0.1);
scatter(repmat(2.35,length(offDursBlue),1),offDursBlue,250,'.','jitter','on', 'jitterAmount',.1);
set(gca,'fontsize',20)

set(gca,'XTick',[1 2])
xticklabels({'UV','Blue'})    
ylabel('Inter-bout interval (s)')

savefig('Inter-bout interval.fig')
saveas(gcf,'Inter-bout interval.tif')
saveas(gcf,'Inter-bout interval.eps')

%%%%% percent active time

meanAmpUV = mean(perActTimeUV);
semAmpUV = std(perActTimeUV/sqrt(length(perActTimeUV)));
meanAmpBlue = mean(perActTimeBlue);
semAmpBlue = std(perActTimeBlue/sqrt(length(perActTimeBlue)));
pvalsAmp = ranksum(perActTimeUV,perActTimeBlue);

cc = linspecer(length(perActTimeUV));
figure
hold on
superbar([meanAmpUV meanAmpBlue],'E',[semAmpUV semAmpBlue],'P',[1 pvalsAmp],'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',35,'BarWidth',0.4)
scatter(repmat(1.35,length(perActTimeUV),1),perActTimeUV,250,'.','jitter','on', 'jitterAmount',0.1);
scatter(repmat(2.35,length(perActTimeBlue),1),perActTimeBlue,250,'.','jitter','on', 'jitterAmount',.1);
set(gca,'fontsize',20)

set(gca,'XTick',[1 2])
xticklabels({'UV','Blue'})    
ylabel('Percent active time')

savefig('Percent active time.fig')
saveas(gcf,'Percent active time.tif')
saveas(gcf,'Percent active time.eps')

%%%%% number of bouts

meanAmpUV = mean(nBoutsUV);
semAmpUV = std(nBoutsUV/sqrt(length(nBoutsUV)));
meanAmpBlue = mean(nBoutsBlue);
semAmpBlue = std(nBoutsBlue/sqrt(length(nBoutsBlue)));
pvalsAmp = ranksum(nBoutsUV,nBoutsBlue);

cc = linspecer(length(nBoutsUV));
figure
hold on
superbar([meanAmpUV meanAmpBlue],'E',[semAmpUV semAmpBlue],'P',[1 pvalsAmp],'BarFaceColor', 'none', 'BarEdgeColor',cc,'PStarFontSize',35,'BarWidth',0.4)
scatter(repmat(1.35,length(nBoutsUV),1),nBoutsUV,250,'.','jitter','on', 'jitterAmount',0.1);
scatter(repmat(2.35,length(nBoutsBlue),1),nBoutsBlue,250,'.','jitter','on', 'jitterAmount',.1);
set(gca,'fontsize',20)

set(gca,'XTick',[1 2])
xticklabels({'UV','Blue'})    
ylabel('Mean number of Bouts')

savefig('Number of Bouts.fig')
saveas(gcf,'Number of Bouts.tif')
saveas(gcf,'Number of Bouts.eps')
