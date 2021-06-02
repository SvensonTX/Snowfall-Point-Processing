% Batch line processing
tic % Start timer

% Creates lists from TLS and Disdromter CSV data within a directory
TLSlist = dir(strcat('D:\Research\OperationSnowflake\Data\LineDataCSV\IntensityReturns/*SOCS.csv'));
DisdroList = dir('D:\Research\OperationSnowflake\Data\LineDataCSV\IntensityReturns\SesameDisdroCSV/*00.csv');

% Initialize variables

% TLS Metrics
Filename = cell(length(TLSlist),1);
DateTime= datetime.empty(length(DisdroList),0);
CountAVG =zeros(length(TLSlist),1);
CountSTD = zeros(length(TLSlist),1);
AVGsnowRange = zeros(length(TLSlist),1);
AVGsnowRange2 = zeros(length(TLSlist),1);
VARsnowRange = zeros(length(TLSlist),1);
CycleDatas = cell(length(TLSlist),1);
CycleIntensities = cell(length(TLSlist),1);
SnowMask = cell(length(TLSlist),1);
PointData = cell(length(TLSlist),1);


% Disdrometer Data
DisIntensity= zeros(length(TLSlist),1);
DisSnowIntensity= zeros(length(TLSlist),1);
NumberOfParticles = zeros(length(TLSlist),1);
Radar = zeros(length(TLSlist),1);
DisIntensityAVG= zeros(length(DisdroList),1);
DisSnowIntensityAVG= zeros(length(DisdroList),1);
DisIntensityVAR= zeros(length(DisdroList),1);
DisSnowIntensityVAR= zeros(length(DisdroList),1);
DisMOR= zeros(length(DisdroList),1);
DisMORavg = zeros(length(DisdroList),1);
DisMORvar = zeros(length(DisdroList),1);
MedianSz = NaN(length(DisdroList),1);
MedianSp = NaN(length(DisdroList),1);

% Sesame Instrument Data
WSWindSpeed = zeros(length(DisdroList),1);
WSWindDir = zeros(length(DisdroList),1);
WSWindGust = zeros(length(DisdroList),1);
WSAirTemp = zeros(length(DisdroList),1);
WSRH = zeros(length(DisdroList),1);
CampAirTemp = zeros(length(DisdroList),1);
CampRH = zeros(length(DisdroList),1);



% Rotation for line scan plot
rotation = 30;
RM = [cosd(rotation) -sind(rotation); sind(rotation) cosd(rotation)];

for k = 57%1:length(TLSlist)
    
TLS = getfield(TLSlist,{k},'name');
Filename{k} = TLS;

% LoadSesameData is a matlab generated function to read CSV of instrument
% data generated from python script (InstrumentData.py)
DisdrometerData = LoadSesameData(getfield(DisdroList,{k},'name'));

% Extract Dates and Times from TLS filename
clear DateStr TimeStr Year Month Day Hour Minute Second
DateStr= string(extractBetween(TLS,1,8));
TimeStr = string(extractBetween(TLS,10,16));
Year = str2double(extractBetween(DateStr,1,4));
Month = str2double(extractBetween(DateStr,5,6));
Day = str2double(extractBetween(DateStr,7,8));
Hour = str2double(extractBetween(TimeStr,1,2));
Minute = str2double(extractBetween(TimeStr,3,4));
Second = str2double(extractBetween(TimeStr,6,7));
DateTime(k) = datetime(Year,Month,Day,Hour,Minute,Second);

% Read CVS point data 
% Assumes Data Already Read In Column Vectors
Points = csvread(strcat('D:\Research\OperationSnowflake\Data\LineDataCSV\IntensityReturns/',TLS));
PointInt = zeros(length(Points),1);
X = Points(:,1);
Y = Points(:,2);
Z = Points(:,3);
I = cast(Points(:,4),'uint16'); % Intensity
T = Points(:,5);                % Time
R = cast(Points(:,6),'uint8');  % Return number
N = cast(Points(:,7),'uint8');  % Total number of returns

% Calculate distance and angle for each point
dist = sqrt(X.^2 + Y.^2 + Z.^2);
theta = atan2(X,Z)*180/pi();

scancnt = 1; % Set line cycle number
clear delta
clear scannum

% Create a Vector scannum that will give line cycle number for each point
scannum = zeros(length(Points),1);
scannum(1) = 1;
delta = zeros(length(Points),1);
for i=2:length(theta)
    delta(i) = abs(theta(i)-theta(i-1));
    if  delta(i) > 45 %We have reached the end of the line scan cycle
        scancnt = scancnt+1; % Increase line scan cycle count
    end
    scannum(i) = scancnt; % Assign line scan count to point return
end

% Compute Maximum Observed Range in 0.05 degree increments
start_angle = min(theta);
stop_angle = max(theta);
incre = 0.05;

%max_range includes max range for each increment of angle.
count = 1;
max_range = zeros(length(start_angle+incre:incre:stop_angle),3);
for j=start_angle+incre:incre:stop_angle
    mask = theta < j & theta > j-incre;
    if any(mask) == 0
        clear mask
        continue
    else
        max_range(count,3) = max(dist(mask));
        max_range(count,1) = j-incre;
        max_range(count,2) = j;
        count = count + 1;
        clear mask;
    end
end

% Determine if it is snow.
clear snow
snow = zeros(length(Points),1);
for i=1:length(X)   
    for j=1:length(max_range)
        if theta(i) < max_range(j,2) && theta(i) > max_range(j,1)
            location = j;
        end
    end
    % Conditional statements to classify point return as hydrometeor.
    % If dist(i) is less than the max range minus a buffer AND within
    % 73 meteres AND theta(i) is above -47.5, then it is a snow return
    if dist(i) < max_range(location,3) - 0.5 && dist(i) < 73 || theta(i) > -47.5 
        snow(i) = 1;
    else        
        snow(i) = 0;
    end

end

% Now we have determined points that are snow.
clear check
CycleData = cell(scancnt,1);
% Save point data of each cycle for returns classified as hydrometeors
for i = 1:scancnt
    % Index snow in line cycle
    cycle = (scannum==i)&(snow==1);
    % Save point data for a specific line scan 
    CycleData{i} = [I(cycle),dist(cycle),R(cycle),N(cycle)];
end
CycleDatas{k} = CycleData;

clear idx
idx = zeros(length(Points),1);
% Find time of line cycle change
for i = 1:length(X)-1
    idx(i) = scannum(i+1) - scannum(i);
end
clear LineTimes
idx(end)=1;
LineScanIdx = find(idx);
LineTimes = T(LineScanIdx)';
% Relative time of cycle change from start of scan time
LineTimes = LineTimes-LineTimes(1);

clear mask_snow

% Extract disdrometer date and time (local time)
DisdroDate=datetime(table2array(DisdrometerData(6,3)),'InputFormat','yyyy-MM-dd');
DisdroTime = datetime(table2array(DisdrometerData(6,4)));

% Classified hydrometeor mask
mask_snow = snow == 1;

% Average and variance of classified hydrometeor range
AVGsnowRange(k) = mean(dist(mask_snow));
VARsnowRange(k) = var(dist(mask_snow));
% Save snow mask for given line scan
SnowMask{k} = mask_snow;

% Calculate average snow range with limiting extent of 70 meters to
% eliminate false positive hydrometeors classified at far ranges
SnowRanges = PointData{k}(SnowMask{k,1},2);
SnowRanges = SnowRanges(SnowRanges<70);
% Average snow range of classified hydrometeors within 70 meters
AVGsnowRange2(k) = mean(SnowRanges);

% Save point data for given line scan
PointData{k} = [I,dist,R,N];

% Create disdrometer size and speed histograms if available
try 
    % Input = date and time of disdrometer (local time)
    [FigSize,FigSpeed,MedSize,MedSpeed] = ParsivelSpectrums(DisdroDate.Year,DisdroDate.Month,DisdroDate.Day,DisdroTime.Hour,DisdroTime.Minute,0);
    set(gca,'FontSize',24)
    set(FigSize, 'WindowState', 'maximized');
    
    % % Save hydrometeor size histogram
    % saveas(FigSize,DateStr+' '+TimeStr+' Parsivel Size Spectrum'+'.png')

    set(FigSpeed, 'WindowState', 'maximized');
    
    % % Save hydrometeor speed histogram
	% saveas(FigSpeed,DateStr+' '+TimeStr+' Parsivel Speed Spectrum'+'.png')
    set(gca,'FontSize',24)
    
    % Median size
    MedianSz(k) = MedSize;
    % Median speed
    MedianSp(k) = MedSpeed;
    
catch
    % Warning if no parsivel data is available from disdrometer
    fprintf('No Parsivel Data for %s\n',DateTime(k))
end

 % Create line scan classification plot
try
    clear RotPoints TempPoints
    figure(k*2-1);
    % Rotate line scan points
    TempPoints = [X Z]';
    RotPoints = [RM*TempPoints]';
    
    hold on
    plot(RotPoints(mask_snow,1),RotPoints(mask_snow,2),'r.');
    plot(RotPoints(~mask_snow,1),RotPoints(~mask_snow,2),'b.')
    xlabel('X (m)')
    xlim([-140 0])
    ylabel('Z (m)')
    ylim([-20 80])
    title('Line Scan: ' + DateStr +' '+ TimeStr)
    % Add median disdrometer speed and size to graph
    text(-60,75,strcat('Median Speed = ',num2str(MedianSp(k)),' m/s'))
    text(-60,70,strcat('Median Size = ',num2str(MedianSz(k)),' mm'))
%     set(gca,'FontSize',24)
%     set(gcf, 'WindowState', 'maximized');
%     saveas(gcf,DateStr+' '+TimeStr+'Line Scan Corrected'+' .png')
    
catch
    % Warning if not able to plot line scan classification plot
    fprintf('Cannot plot %s Line Scan',DateTime(k))
end



% SnowCount is classified hydrometeors per line scan cycle
clear SnowCount
SnowCount = zeros(length(LineTimes),1);
for i=1:scancnt
    
    SnowCount(i) = sum((scannum == i) & (snow ==1));
end
% Average and standard deviation of classified hydrometeor counts per cycle
CountAVG(k) =mean(SnowCount);
CountSTD(k) = std(SnowCount);
%%

clear DDTime DDIntensity Times Start
Start=datetime(Year, Month, Day, Hour, Minute, Second);
Times = Start +seconds(LineTimes);
DDTime=datetime(table2array(DisdrometerData(4:8,2)));
% DisTime=datetime(Year,Month,Day,DDTime.Hour,DDTime.Minute,DDTime.Second)
DDIntensity=table2array(DisdrometerData(4:8,5));
DDSnowIntensity=table2array(DisdrometerData(4:8,11));

% Create line scan count plot with disdrometer intensities if available
try
    
    figure(k*2)
    hold on
    plot(Times,SnowCount)
    xlabel('Time')
    ylabel('Number of Snow Points Detected')
    title('Line Scan: ' + DateStr +' '+ TimeStr)
    yyaxis right
    stairs(DDTime(:,1),DDIntensity(:,1))
    stairs(DDTime(:,1),DDSnowIntensity(:,1))
    legend('Snow Count','Rain Intensity','Snow Intensity')
    ylabel('Intensity of Precipitation (mm/h)')
%     set(gca,'FontSize',24)
%     set(gcf, 'WindowState', 'maximized');
    saveas(gcf,DateStr+' '+TimeStr+' Snow Count'+' .png')
catch
    fprintf('Cannot plot %s Snow Count',DateTime(k))
end

% Close all figures
% close all




% Assign instrument data to variables for each line scan
DisMOR(k)=table2array(DisdrometerData(6,8));
DisMORavg(k)= mean(table2array(DisdrometerData(5:7,8)));
DisMORvar(k)= var(table2array(DisdrometerData(5:7,8)));
DisIntensityAVG(k)=mean(table2array(DisdrometerData(:,5)));
DisSnowIntensityAVG(k)=mean(table2array(DisdrometerData(:,11)));
DisIntensityVAR(k)=var(table2array(DisdrometerData(:,5)));
DisSnowIntensityVAR(k)=var(table2array(DisdrometerData(:,11)));
DisIntensity(k)= DDIntensity(3);
DisSnowIntensity(k)= DDSnowIntensity(3);
WSWindSpeed(k)=table2array(DisdrometerData(6,17));
WSWindDir(k)=table2array(DisdrometerData(6,18));
WSWindGust(k)=table2array(DisdrometerData(6,19));
WSAirTemp(k)=table2array(DisdrometerData(6,14));
WSRH(k)=table2array(DisdrometerData(6,15));
CampAirTemp(k)=table2array(DisdrometerData(6,43));
CampRH(k)=table2array(DisdrometerData(6,44));

NumberOfParticles(k) = table2array(DisdrometerData(6,10));
Radar(k) = table2array(DisdrometerData(6,7));
end
toc
%%

DateTime = DateTime';

% Compiled results
Results = table(DateTime,Filename,CountAVG,CountSTD,WSWindSpeed,WSWindDir,WSWindGust,WSAirTemp,WSRH,CampAirTemp,CampRH,DisMOR,DisMORavg,DisMORvar,DisIntensity,DisIntensityAVG,DisIntensityVAR,DisSnowIntensity,DisSnowIntensityAVG,DisSnowIntensityVAR,NumberOfParticles,Radar,AVGsnowRange,AVGsnowRange2,VARsnowRange);

% Index of disdrometer readings of no precipitation
CleanIndex = Results.Radar ~= -9.9990;

% Filter clean (robust) scans
CleanResults = Results(CleanIndex,:);

% Filter questionable (discrepant) scans
QuestionableResults = Results(~CleanIndex,:);
