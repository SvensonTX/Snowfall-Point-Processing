%% Line scan statistics
% Generated from point data and snow mask created in LineScanProcessing

SnowData = cell(length(TLSlist),1);

for i = 1:length(TLSlist)

    SnowData{i} = PointData{i}(SnowMask{i,1},:);
    
end

% CleanIndex is robust line scans (i.e. disdrometer indicates snowfall)
CleanSnowData = SnowData(CleanIndex);
CleanPointData = PointData(CleanIndex);

%% Far point statistics
FarCountTotal = zeros(length(CleanSnowData),1);
Far1stCount = zeros(length(CleanSnowData),1);
Far2ndCount = zeros(length(CleanSnowData),1);

for k = 1:length(PointData)
    FarCountTotal(k) = sum(PointData{k}(:,2)>74);

    IndexFar1stReturns = PointData{k}(:,3)==1 & PointData{k}(:,2)>74;
    IndexFar2ndReturns = PointData{k}(:,3)==2 & PointData{k}(:,2)>74;
    Far1stCount(k) = sum(IndexFar1stReturns);
    Far2ndCount(k) = sum(IndexFar2ndReturns);
end

%% Near Snow Counts

AllNearCountTotal = zeros(length(SnowMask),1);

for k = 1:length(PointData)

    AllNearCountTotal(k) = length(PointData{k}(((PointData{k,1}(SnowMask{k},2)<=5))));

end



%% Snow Statistics

% Initialize variables
SnowCountTotal = zeros(length(CleanSnowData),1);
Snow1stCount = zeros(length(CleanSnowData),1);
Snow2ndCount = zeros(length(CleanSnowData),1);
Snow3rdCount = zeros(length(CleanSnowData),1);
Snow4thCount = zeros(length(CleanSnowData),1);

Returns1stCount = zeros(length(CleanSnowData),1);
Returns2ndCount = zeros(length(CleanSnowData),1);
Returns3rdCount = zeros(length(CleanSnowData),1);
Returns4thCount = zeros(length(CleanSnowData),1);

Return1stMax = zeros(length(CleanSnowData),1);  
Return1stMin = zeros(length(CleanSnowData),1);  
Return1stAvg = zeros(length(CleanSnowData),1);  
Return1stIMax =zeros(length(CleanSnowData),1);
Return1stVar = zeros(length(CleanSnowData),1);  
Return1stIVar =zeros(length(CleanSnowData),1); 
Return1stIMin =zeros(length(CleanSnowData),1);  
Return1stIAvg =zeros(length(CleanSnowData),1);  

Return2ndMax = zeros(length(CleanSnowData),1);  
Return2ndMin = zeros(length(CleanSnowData),1);  
Return2ndAvg = zeros(length(CleanSnowData),1);  
Return2ndIMax =zeros(length(CleanSnowData),1);
Return2ndVar = zeros(length(CleanSnowData),1);  
Return2ndIVar =zeros(length(CleanSnowData),1);
Return2ndIMin =zeros(length(CleanSnowData),1);  
Return2ndIAvg =zeros(length(CleanSnowData),1);  

Return3rdMax=  zeros(length(CleanSnowData),1);  
Return3rdMin = zeros(length(CleanSnowData),1);  
Return3rdAvg = zeros(length(CleanSnowData),1);  
Return3rdIMax =zeros(length(CleanSnowData),1); 
Return3rdVar = zeros(length(CleanSnowData),1);  
Return3rdIVar =zeros(length(CleanSnowData),1); 
Return3rdIMin =zeros(length(CleanSnowData),1);  
Return3rdIAvg =zeros(length(CleanSnowData),1);  

Return4thMax=  zeros(length(CleanSnowData),1);  
Return4thMin = zeros(length(CleanSnowData),1);  
Return4thAvg = zeros(length(CleanSnowData),1);  
Return4thIMax =zeros(length(CleanSnowData),1);
Return4thVar = zeros(length(CleanSnowData),1);  
Return4thIVar =zeros(length(CleanSnowData),1); 
Return4thIMin =zeros(length(CleanSnowData),1);  
Return4thIAvg =zeros(length(CleanSnowData),1);  

ReturnNum1stMax =  zeros(length(CleanSnowData),1);
ReturnNum1stMin =  zeros(length(CleanSnowData),1);
ReturnNum1stAvg =  zeros(length(CleanSnowData),1);
ReturnNum1stIMax = zeros(length(CleanSnowData),1);
ReturnNum1stVar =  zeros(length(CleanSnowData),1);
ReturnNum1stIVar = zeros(length(CleanSnowData),1);
ReturnNum1stIMin = zeros(length(CleanSnowData),1);
ReturnNum1stIAvg = zeros(length(CleanSnowData),1);

ReturnNum2ndMax= zeros(length(CleanSnowData),1);  
ReturnNum2ndMin= zeros(length(CleanSnowData),1);  
ReturnNum2ndAvg= zeros(length(CleanSnowData),1);  
ReturnNum2ndIMax =zeros(length(CleanSnowData),1); 
ReturnNum2ndVar= zeros(length(CleanSnowData),1);  
ReturnNum2ndIVar =zeros(length(CleanSnowData),1); 
ReturnNum2ndIMin =zeros(length(CleanSnowData),1); 
ReturnNum2ndIAvg =zeros(length(CleanSnowData),1); 

ReturnNum3rdMax= zeros(length(CleanSnowData),1);  
ReturnNum3rdMin= zeros(length(CleanSnowData),1);  
ReturnNum3rdAvg= zeros(length(CleanSnowData),1);  
ReturnNum3rdIMax =zeros(length(CleanSnowData),1);
ReturnNum3rdVar= zeros(length(CleanSnowData),1);  
ReturnNum3rdIVar =zeros(length(CleanSnowData),1); 
ReturnNum3rdIMin =zeros(length(CleanSnowData),1); 
ReturnNum3rdIAvg =zeros(length(CleanSnowData),1); 

ReturnNum4thMax= zeros(length(CleanSnowData),1);  
ReturnNum4thMin= zeros(length(CleanSnowData),1);  
ReturnNum4thAvg= zeros(length(CleanSnowData),1);  
ReturnNum4thIMax =zeros(length(CleanSnowData),1);
ReturnNum4thVar= zeros(length(CleanSnowData),1);  
ReturnNum4thIVar =zeros(length(CleanSnowData),1);
ReturnNum4thIMin =zeros(length(CleanSnowData),1); 
ReturnNum4thIAvg =zeros(length(CleanSnowData),1); 

for k = 1:length(CleanSnowData)
    
    SnowCountTotal(k) = length(CleanSnowData{k});

    Index1stReturns = CleanSnowData{k}(:,3)==1;
    Index2ndReturns = CleanSnowData{k}(:,3)==2;
    Index3rdReturns = CleanSnowData{k}(:,3)==3;
    Index4thReturns = CleanSnowData{k}(:,3)==4;
    Index1stReturnNumbers = CleanSnowData{k}(:,4)==1;
    Index2ndReturnNumbers = CleanSnowData{k}(:,4)==2;
    Index3rdReturnNumbers = CleanSnowData{k}(:,4)==3;
    Index4thReturnNumbers = CleanSnowData{k}(:,4)==4;

    % Return Counts
    Snow1stCount(k) = sum(Index1stReturns);
    Snow2ndCount(k) = sum(Index2ndReturns);
    Snow3rdCount(k) = sum(Index3rdReturns);
    Snow4thCount(k) = sum(Index4thReturns);
    
    %Return Total Counts
    Returns1stCount(k) = sum(Index1stReturnNumbers);
    Returns2ndCount(k) = sum(Index2ndReturnNumbers);
    Returns3rdCount(k) = sum(Index3rdReturnNumbers);
    Returns4thCount(k) = sum(Index4thReturnNumbers);


    %Snow Ranges
    SnowRanges = CleanSnowData{k}(:,2);
    
    %Snow Intensities
    SnowIntensities = CleanSnowData{k}(:,1);

    % 1st Returns
    Return1stMax(k) = max(SnowRanges(Index1stReturns));
    Return1stMin(k) = min(SnowRanges(Index1stReturns));
    Return1stAvg(k) = mean(SnowRanges(Index1stReturns));
    Return1stVar(k) = var(cast(SnowRanges(Index1stReturns),'single'));
    Return1stIMax(k) = max(SnowIntensities(Index1stReturns));
    Return1stIVar(k) = var(cast(SnowIntensities(Index1stReturns),'single'));
    Return1stIMin(k) = min(SnowIntensities(Index1stReturns));
    Return1stIAvg(k) = mean(SnowIntensities(Index1stReturns));
    
    % 2nd Returns
    Return2ndMax(k) = max(SnowRanges(Index2ndReturns));
    Return2ndMin(k) = min(SnowRanges(Index2ndReturns));
    Return2ndAvg(k) = mean(SnowRanges(Index2ndReturns));
    Return2ndIMax(k) = max(SnowIntensities(Index2ndReturns));
    Return2ndVar(k) = var(cast(SnowRanges(Index2ndReturns),'single'));
    Return2ndIVar(k) = var(cast(SnowIntensities(Index2ndReturns),'single'));
    Return2ndIMin(k) = min(SnowIntensities(Index2ndReturns));
    Return2ndIAvg(k) = mean(SnowIntensities(Index2ndReturns));
    

    % 3rd Returns
    if isempty(max(SnowRanges(Index3rdReturns))) 
        Return3rdMax(k) = 0;
        Return3rdAvg(k) = 0;
        Return3rdIMax(k) = 0;
        Return3rdVar(k) = 0;
        Return3rdIVar(k) = 0;
        Return3rdIMin(k) = 0;
        Return3rdIAvg(k) = 0;
    else
        Return3rdMax(k)= max(SnowRanges(Index3rdReturns));   
        Return3rdMin(k) = min(SnowRanges(Index3rdReturns));
        Return3rdAvg(k) = mean(SnowRanges(Index3rdReturns));
        Return3rdIMax(k) = max(SnowIntensities(Index3rdReturns));
        Return3rdVar(k) = var(cast(SnowRanges(Index3rdReturns),'single'));
        Return3rdIVar(k) = var(cast(SnowIntensities(Index3rdReturns),'single'));
        Return3rdIMin(k) = min(SnowIntensities(Index3rdReturns));
        Return3rdIAvg(k) = mean(SnowIntensities(Index3rdReturns));
    end
    
    % 4th Returns
    if isempty(max(SnowRanges(Index4thReturns))) 
        Return4thMax(k) = 0;
        Return4thAvg(k) = 0;
        Return4thIMax(k) = 0;
        Return4thVar(k) = 0;
        Return4thIVar(k) = 0;
        Return4thIMin(k) = 0;
        Return4thIAvg(k) = 0;
    else
        Return4thMax(k) = max(SnowRanges(Index4thReturns));   
        Return4thMin(k) = min(SnowRanges(Index4thReturns));
        Return4thAvg(k) = mean(SnowRanges(Index4thReturns));
        Return4thIMax(k) = max(SnowIntensities(Index4thReturns));
        Return4thVar(k) = var(cast(SnowRanges(Index4thReturns),'single'));
        Return4thIVar(k) = var(cast(SnowIntensities(Index4thReturns),'single'));
        Return4thIMin(k) = min(SnowIntensities(Index4thReturns));
        Return4thIAvg(k) = mean(SnowIntensities(Index4thReturns));
    end 
    
    
    % Total Return Numbers
    % 1 Return
    ReturnNum1stMax(k) = max(SnowRanges(Index1stReturnNumbers));
    ReturnNum1stMin(k) = min(SnowRanges(Index1stReturnNumbers));
    ReturnNum1stAvg(k) = mean(SnowRanges(Index1stReturnNumbers));
    ReturnNum1stIMax(k) = max(SnowIntensities(Index1stReturnNumbers));
    ReturnNum1stVar(k) = var(cast(SnowRanges(Index1stReturnNumbers),'single'));
    ReturnNum1stIVar(k) = var(cast(SnowIntensities(Index1stReturnNumbers),'single'));
    ReturnNum1stIMin(k) = min(SnowIntensities(Index1stReturnNumbers));
    ReturnNum1stIAvg(k) = mean(SnowIntensities(Index1stReturnNumbers));
    
    % 2 Returns
    ReturnNum2ndMax(k) = max(SnowRanges(Index2ndReturnNumbers));
    ReturnNum2ndMin(k) = min(SnowRanges(Index2ndReturnNumbers));
    ReturnNum2ndAvg(k) = mean(SnowRanges(Index2ndReturnNumbers));
    ReturnNum2ndIMax(k) = max(SnowIntensities(Index2ndReturnNumbers));
    ReturnNum2ndVar(k) = var(cast(SnowRanges(Index2ndReturnNumbers),'single'));
    ReturnNum2ndIVar(k) = var(cast(SnowIntensities(Index2ndReturnNumbers),'single'));
    ReturnNum2ndIMin(k) = min(SnowIntensities(Index2ndReturnNumbers));
    ReturnNum2ndIAvg(k) = mean(SnowIntensities(Index2ndReturnNumbers));
    

    % 3 Returns
    if isempty(max(SnowRanges(Index3rdReturnNumbers))) 
        ReturnNum3rdMax(k) = 0;
        ReturnNum3rdMin(k) = 0;
        ReturnNum3rdAvg(k) = 0;
        ReturnNum3rdVar(k) = 0;
        ReturnNum3rdIVar(k) = 0;
        ReturnNum3rdIMax(k) = 0;
        ReturnNum3rdIMin(k) = 0;
        ReturnNum3rdIAvg(k) = 0;
    else
        ReturnNum3rdMax(k) = max(SnowRanges(Index3rdReturnNumbers));
        ReturnNum3rdMin(k) = min(SnowRanges(Index3rdReturnNumbers));
        ReturnNum3rdAvg(k) = mean(SnowRanges(Index3rdReturnNumbers));
        ReturnNum3rdIMax(k) = max(SnowIntensities(Index3rdReturnNumbers));
        ReturnNum3rdVar(k) = var(cast(SnowRanges(Index3rdReturnNumbers),'single'));
        ReturnNum3rdIVar(k) = var(cast(SnowIntensities(Index3rdReturnNumbers),'single'));
        ReturnNum3rdIMin(k) = min(SnowIntensities(Index3rdReturnNumbers));
        ReturnNum3rdIAvg(k) = mean(SnowIntensities(Index3rdReturnNumbers));
    end
    
    % 4 Returns
    if isempty(max(SnowRanges(Index4thReturnNumbers))) 
        ReturnNum4thMax(k) = 0;
        ReturnNum4thMin(k) = 0;
        ReturnNum4thAvg(k) = 0;
        ReturnNum4thVar(k) = 0;
        ReturnNum4thIVar(k) = 0;
        ReturnNum4thIMax(k) = 0;
        ReturnNum4thIMin(k) = 0;
        ReturnNum4thIAvg(k) = 0;
    else
        ReturnNum4thMax(k) = max(SnowRanges(Index4thReturnNumbers));
        ReturnNum4thMin(k) = min(SnowRanges(Index4thReturnNumbers));
        ReturnNum4thAvg(k) = mean(SnowRanges(Index4thReturnNumbers));
        ReturnNum4thIMax(k) = max(SnowIntensities(Index4thReturnNumbers));
        ReturnNum4thVar(k) = var(cast(SnowRanges(Index4thReturnNumbers),'single'));
        ReturnNum4thIVar(k) = var(cast(SnowIntensities(Index4thReturnNumbers),'single'));
        ReturnNum4thIMin(k) = min(SnowIntensities(Index4thReturnNumbers));
        ReturnNum4thIAvg(k) = mean(SnowIntensities(Index4thReturnNumbers));
    end

end
%%

Intensity = CleanResults.DisIntensity + CleanResults.DisSnowIntensity

CleanDateTime = DateTime(CleanIndex)
CleanAVGsnowRange2 = AVGsnowRange2(CleanIndex)
SnowStats = table(CleanDateTime,Intensity,SnowCountTotal,Snow1stCount,Returns1stCount,Return1stMin,Return1stAvg,Return1stVar,Return1stMax,ReturnNum1stMin,ReturnNum1stAvg,ReturnNum1stVar,ReturnNum1stMax,Return1stIMin,Return1stIAvg,Return1stIVar,Return1stIMax,ReturnNum1stIMin,ReturnNum1stIAvg,ReturnNum1stIVar,ReturnNum1stIMax,...
                                                Snow2ndCount,Returns2ndCount,Return2ndMin,Return2ndAvg,Return2ndVar,Return2ndMax,ReturnNum2ndMin,ReturnNum2ndAvg,ReturnNum2ndVar,ReturnNum2ndMax,Return2ndIMin,Return2ndIAvg,Return2ndIVar,Return2ndIMax,ReturnNum2ndIMin,ReturnNum2ndIAvg,ReturnNum2ndIVar,ReturnNum2ndIMax,...
                                                Snow3rdCount,Returns3rdCount,Return3rdMin,Return3rdAvg,Return3rdVar,Return3rdMax,ReturnNum3rdMin,ReturnNum3rdAvg,ReturnNum3rdVar,ReturnNum3rdMax,Return3rdIMin,Return3rdIAvg,Return3rdIVar,Return3rdIMax,ReturnNum3rdIMin,ReturnNum3rdIAvg,ReturnNum3rdIVar,ReturnNum3rdIMax,...
                                                Snow4thCount,Returns4thCount,Return4thMin,Return4thAvg,Return4thVar,Return4thMax,ReturnNum4thMin,ReturnNum4thAvg,ReturnNum4thVar,ReturnNum4thMax,Return4thIMin,Return4thIAvg,Return4thIVar,Return4thIMax,ReturnNum4thIMin,ReturnNum4thIAvg,ReturnNum4thIVar,ReturnNum4thIMax...
                                                ,NearCountTotal,CleanAVGsnowRange2)%,FarCountTotal,Far1stCount,Far2ndCount,MidCountTotal)%FarCountTotal,Far1stCount,Far2ndCount)
SnowStats.Properties.VariableNames{1} = 'DateTime';
AllData = join(CleanResults,SnowStats)


%% Near Range Higher Spatial resolution
% Range Binning



RangeBin = 0:0.05:6;

formatOut = 'yyyymmdd HHMM-SS'
for k=1:335
    range1= SnowMask{k,1}((PointData{k,1}(:,3)==1));
    range2= SnowMask{k,1}((PointData{k,1}(:,3)==2));
    range3= SnowMask{k,1}((PointData{k,1}(:,3)==3));
    range4= SnowMask{k,1}((PointData{k,1}(:,3)==4));

    subplot(4,1,1)

    histogram(PointRanges{k,1}(range1,1),RangeBin)
%     HF = histfit(PointData{k,1}(range1,2),19,'exponential')

    title(strcat('Range Histogram for',' ',datestr(DateTime(k),formatOut)))
    subplot(4,1,2)
    histogram(PointRanges{k,1}(range2,1),RangeBin)
    subplot(4,1,3)
    histogram(PointRanges{k,1}(range3,1),RangeBin)
    subplot(4,1,4)
    histogram(PointRanges{k,1}(range4,1),RangeBin)


    saveas(gcf,strcat(datestr(DateTime(k),formatOut),' Line Scan Range Histogram - 5cm','.png'))

    close all
end


%% Spectralon Analysis
Results.MedSz = MedianSz;
Results.NearCountTotal = AllNearCountTotal;

Results.FarCountTotal = FarCountTotal;

Results.Far1stCount = Far1stCount;

Results.Far2ndCount = Far2ndCount;
Results.AVGsnowRange2 = AVGsnowRange2;

DiscrepantScans = Results(~CleanIndex,:);
UncertainResults = DiscrepantScans(DiscrepantScans.DisSnowIntensityAVG == 0,:);

% Time table of results
ResultsTimeTable=table2timetable(Results);

% First import InstrumentData.csv
InstrumentData = table2timetable(InstrumentData);

% SpectralonStats function imports statistics generated from
% 'ExtractTarget_withStatistics.py'
FrameSpec2ndReturn = table2timetable(SpectralonStats('SpectralonFrameStats_2ndReturns.csv'));
FrameSpec2ndReturn= rmmissing(FrameSpec2ndReturn);
FrameSpec1stReturn = table2timetable(SpectralonStats('SpectralonFrameStats_1stReturns.csv'));
FrameSpec1stReturn= rmmissing(FrameSpec1stReturn);
TreeFrame1stReturn = table2timetable(SpectralonStats('TreeTrunkFrameStats_1stReturns.csv'));
TreeFrame2ndReturn = table2timetable(SpectralonStats('TreeTrunkFrameStats_2ndReturns.csv'));

LineScanSpectralon = table2timetable(SpectralonStats('SpectralonStats.csv'));

FrameSnowSpectralon = table2timetable(SpectralonStats('SpectralonSnowFrameStats.csv'));
FrameSnowSpectralon(FrameSnowSpectralon.VarName3 == 0,:) = [];

% Sync different datasets together based on time

tol = minutes(0.4);

% Sync Frame Spectralon Scans with Instrument Data
tMatch_Results = FrameSpec1stReturn(withtol(InstrumentData.DateTime,tol),:).Date;
AllSpecMatch= InstrumentData(withtol(tMatch_Results, tol),:);
AllSpecMatch=retime(AllSpecMatch,tMatch_Results,'nearest');
AllSpecMatch=rmmissing(AllSpecMatch);
AllSpecMatch.Properties.VariableNames{1} = 'DateStr';
FrameSpecDataWithIns = synchronize(FrameSpec1stReturn,AllSpecMatch);

FrameSpecDataWithIns= rmmissing(FrameSpecDataWithIns);


% Sync Spectralon Scans with Line Scans
tMatch_Results = ResultsTimeTable(withtol(LineScanSpectralon.Date,tol),:).DateTime;
SpecMatch= LineScanSpectralon(withtol(tMatch_Results, tol),:);
SpecMatch=retime(SpecMatch,tMatch_Results,'nearest');
LineScanSpecData = synchronize(ResultsTimeTable,SpecMatch);
LineScanSpecData= rmmissing(LineScanSpecData);


% Sync Snow Frame Spectralon Scans with Line Scans
tMatch_Results = ResultsTimeTable(withtol(FrameSnowSpectralon.Date,tol),:).DateTime;
SpecMatch= FrameSnowSpectralon(withtol(tMatch_Results, tol),:);
SpecMatch=retime(SpecMatch,tMatch_Results,'nearest');
SpecMatch = removevars(SpecMatch,{'VarName4'});
SpecMatch=rmmissing(SpecMatch);
FrameScanSpecData = synchronize(ResultsTimeTable,SpecMatch);
FrameScanSpecData= rmmissing(FrameScanSpecData);

tol=minutes(5)
% Sync All Frame Spectralon Scans with Line Scans
tMatch_Results = ResultsTimeTable(withtol(FrameSnowSpectralon.Date,tol),:).DateTime;
AllSpecMatch= FrameSnowSpectralon(withtol(tMatch_Results, tol),:);
AllSpecMatch=retime(AllSpecMatch,tMatch_Results,'nearest');
AllSpecMatch=rmmissing(AllSpecMatch);
AllFrameSpecData = synchronize(ResultsTimeTable,AllSpecMatch);
AllFrameSpecData= rmmissing(AllFrameSpecData);


tol = minutes(0.4)
tMatch_Results = InstrumentData(withtol(ReFrameSpec.VarName1,tol),:).VarName1;
AllSpecMatch= ReFrameSpec(withtol(tMatch_Results, tol),:);
AllSpecMatch=retime(AllSpecMatch,tMatch_Results,'nearest');
AllSpecMatch=rmmissing(AllSpecMatch);
FrameSpecDataWithIns = synchronize(InstrumentData,AllSpecMatch);
FrameSpecDataWithIns= rmmissing(FrameSpecDataWithIns);



%% Range Binning

RangeBin = 1:1:50;

formatOut = 'yyyymmdd HHMM-SS'
for k=1:335
    
    range1= SnowMask{k,1}&((PointData{k,1}(:,3)==1));
    range2= SnowMask{k,1}&((PointData{k,1}(:,3)==2));
    range3= SnowMask{k,1}&((PointData{k,1}(:,3)==3));
    range4= SnowMask{k,1}&((PointData{k,1}(:,3)==4));
    subplot(4,1,1)

    histogram(PointData{k,1}(range1,2),RangeBin)
%     HF = histfit(PointData{k,1}(range1,2),19,'exponential')

    title(strcat('Range Histogram for',' ',datestr(DateTime(k),formatOut)))
    ylabel('1st Returns')

    subplot(4,1,2)
    histogram(PointData{k,1}(range2,2),RangeBin)
    ylabel('2nd Returns')

    subplot(4,1,3)
    histogram(PointData{k,1}(range3,2),RangeBin)
    ylabel('3rd Returns')

    subplot(4,1,4)
    histogram(PointData{k,1}(range4,2),RangeBin)
    ylabel('4th Returns')
    xlabel('Range (m)')

%     saveas(gcf,strcat(datestr(DateTime(k),formatOut),' Line Scan Range Histogram','.png'))
%     input('Enter')
    close all
end

%% e-folding length

RangeBin = 1:1:70;

eFoldFit = zeros(335,1);

formatOut = 'yyyymmdd HHMM-SS'
for k=1:335
    
    % Classified hydrometeors within 70 meter range
    range1= SnowMask{k,1} & ((PointData{k,1}(:,2)<70 ));

    figure()
    hold on 
    % Create exponential fit
    HF = histfit(PointData{k,1}(range1,2),70,'exponential');
    % Array of exponential fit values
    FitArray = HF(2,1).YData;
    % e-folding length
    eFold = HF(2,1).YData(1)/exp(1);
    % Index of fit array closest to e-folding length
    [~,eIdx] = min(abs(FitArray-eFold));
    
    plot(HF(2,1).XData(eIdx),HF(2,1).YData(eIdx),'g*','MarkerSize',12)
    eFoldFit(k) = HF(2,1).XData(eIdx);
    title(strcat('Range Histogram for','-',datestr(DateTime(k),formatOut)));
    ylabel('Classified Hydrometeor Counts');xlabel('Range (m)');legend('Range Counts','Exponential Fit','e-folding length');

    saveas(gcf,strcat(datestr(DateTime(k),formatOut),' Line Scan Exponential Fit','.png'))
    close all
end

% figure()
% scatter(eFoldFit,Results.DisSnowIntensity)
% title('E folding length');xlabel('E folding length (m)');ylabel('Disdrometer MOR (m)');
% figure()
% scatter(Results.AVGsnowRange2,Results.DisMOR)
% title('Average Snow Range');xlabel('Average Snow Range (m)');ylabel('Disdrometer MOR (m)');
% figure()
% scatter(eFoldFit,AVGsnowRange2)
% title('E folding length vs Average Snow Range');xlabel('E folding length (m)');ylabel('Average Snow Range (m)');

%% Streak Velocities

% Calculate streak velocities from point cloud (returns limited to 5 meters
% with a scalar field of linearity calculated from CloudCompare)
[StreakTable,StreakHist] = StreakVelocities('20191224-0600-01 - Cloud.csv')