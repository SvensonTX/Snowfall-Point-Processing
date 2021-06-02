function LoadSesameData = LoadSesameData(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   UNTITLED = IMPORTFILE(FILENAME) Reads data from text file FILENAME for
%   the default selection.
%
%   UNTITLED = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from rows
%   STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   Untitled = importfile('191127155900.csv', 2, 12);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2020/09/22 21:34:07

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Format for each line of text:
%   column1: double (%f)
%	column2: text (%s)
%   column3: datetimes (%{yyyy-MM-dd}D)
%	column4: datetimes (%{HH:mm:ss}D)
%   column5: double (%f)
%	column6: double (%f)
%   column7: double (%f)
%	column8: double (%f)
%   column9: double (%f)
%	column10: double (%f)
%   column11: double (%f)
%	column12: double (%f)
%   column13: double (%f)
%	column14: double (%f)
%   column15: double (%f)
%	column16: double (%f)
%   column17: double (%f)
%	column18: double (%f)
%   column19: double (%f)
%	column20: double (%f)
%   column21: double (%f)
%	column22: double (%f)
%   column23: double (%f)
%	column24: double (%f)
%   column25: double (%f)
%	column26: double (%f)
%   column27: double (%f)
%	column28: double (%f)
%   column29: double (%f)
%	column30: double (%f)
%   column31: double (%f)
%	column32: double (%f)
%   column33: double (%f)
%	column34: double (%f)
%   column35: double (%f)
%	column36: double (%f)
%   column37: double (%f)
%	column38: double (%f)
%   column39: double (%f)
%	column40: double (%f)
%   column41: double (%f)
%	column42: double (%f)
%   column43: double (%f)
%	column44: double (%f)
%   column45: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%q%q%{HH:mm:ss}D%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%*s%*s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
LoadSesameData = table(dataArray{1:end-1}, 'VariableNames', {'VarName1','DateTimemmddyyhhmmss','Date','Time','Intensityofprecipitationmmh','Precipitationsincestartmm','RadarreflectivitydBz','MORVisibilitym','SignalamplitudeofLaserband','Numberofdetectedparticles','Snowintensitymmh','RECORD','BattVolt','WS600AirTemp','WS600RH','WS600AirPressure','WS600WindSpeed','WS600WindDir','WS600WindGust','WS600WindQuality','WS600PrecipType','WS600Precip1min_hourlyrate','WS600Precip1min','WS600Precip15min','WS600PrecipWorkDay','WS600PrecipMidnight','WS600PrecipBoard','WS600PrecipWY','MetOnePrecip1min','MetOnePrecip15min','MetOnePrecipWorkDay','MetOnePrecipMidnight','MetOnePrecipBoard','MetOnePrecipWY','RMYoungPrecip1min','RMYoungPrecip15min','RMYoungPrecipWorkDay','RMYoungPrecipMidnight','RMYoungPrecipBoard','RMYoungPrecipWY','JuddSnowDepthBoard','JuddSnowDepthTotal','CampbellAirTemp','CampbellRH','SnowPillow'});

% For code requiring serial dates (datenum) instead of datetime, uncomment
% the following line(s) below to return the imported dates as datenum(s).

% Untitled.Date=datenum(Untitled.Date);
% Untitled.Time=datenum(Untitled.Time);
