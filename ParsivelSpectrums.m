function [FigSize, FigSpeed,MedSize,MedSpeed] = ParsivelSpectrums(year,month,day,hour,minute,second)

load('Plookup.mat') % Date of disdrometer record
load('Pdata.mat')   % Data from disdrometer, col 7 in cell

% Search rows for exact disdrometer date
rows = find( subset(:,1) == year & subset(:,2) == month & subset(:,3) == day & subset(:,4) == hour & subset(:,5) == minute & subset(:,6) == second);

FigSize = figure(1);
% Extract disdrometer data from cell
temp = ParsivelDataset{rows,7}{:,:};
% Array of size counts
SzCnts=sum(temp);
% Bar graph of size counts
bar(1:32,SzCnts)
xticks([1:32])
% Set x labels to disdrometer sizes
xticklabels({'0.062','0.187','0.312','0.437','0.562','0.687','0.812','0.937','1.062','1.187','1.375','1.625','1.875','2.125','2.375','2.75','3.25','3.75','4.25','4.75','5.5','6.5','7.5','8.5','9.5','11','13','15','17','18','21.5','24.5'})
xtickangle(90)
title(sprintf('Parsivel Size Spectrum: %d/%d/%d %d:%d:%d',year,month,day,hour,minute,second))
xlabel('Volume-equivalent Diameter (mm)')
ylabel('Count')
xlabelsz= xticklabels;


DisCnt = zeros(sum(SzCnts),1);
k=1;
for i = 1:length(SzCnts)
    for j=1:SzCnts(i)
        % Create array with size estimate for every particle
        DisCnt(k) = 1*str2double(cell2mat(xlabelsz(i)));
        k= k+1;
    end
end       
% Calculate median size from disdrometer
MedSize = median(DisCnt);

FigSpeed = figure(2);
% Array of speed counts
SpCnts=sum(temp,2);
% Bar graph of speed counts
bar(1:32,SpCnts)
xticks([1:32])
% Set x labels to disdrometer speeds
xticklabels({'0.05','0.15','0.25','0.35','0.45','0.55','0.65','0.75','0.85','0.95','1.1','1.3','1.5','1.7','1.9','2.2','2.6','3','3.4','3.8','4.4','5.2','6','6.8','7.6','8.8','10.4','12','13.6','15.2','17.6','20.8'})
xtickangle(90)
title(sprintf('Parsivel Speed Spectrum: %d/%d/%d %d:%d:%d',year,month,day,hour,minute,second))
xlabel('Speed (m/s)')
ylabel('Count')
xlabelspd=xticklabels;

DisCnt = zeros(sum(SpCnts),1);
k=1;
for i = 1:length(SpCnts)
    for j=1:SpCnts(i)
        % Create array with speed estimate for every particle
        DisCnt(k) = 1*str2double(cell2mat(xlabelspd(i)));
        k= k+1;
    end
end       
% Calculate median speed from disdrometer
MedSpeed = median(DisCnt);

end