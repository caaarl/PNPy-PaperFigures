colors =    [[0.,0.,0.]; ...
            [230., 159., 0.]; ...
            [86., 180., 233.]; ...
            [0., 158., 115.]...
            ]/255;

segmentLength = 15;
loadNames = {'sciaticCell.mat', 'vagusCell.mat'};

figure,
colorInds = [3,2];
for i = 1:2
    dataCell = importdata(loadNames{i}); % sciaticCell; %vagusCell; %  

    dataSegCell = section(dataCell, segmentLength);
    dirDiffNorms = dirDiffNorm(dataSegCell, segmentLength);

    
    subplot(2,1,i)
    histogram(dirDiffNorms, 'BinWidth', 0.05, 'EdgeColor', 'none', 'FaceColor', colors(colorInds(i),:))
    xlim([0,1])
    title(loadNames{i})
    
    xlabel('c (unitless)')
    ylabel('#')
end

% Print setup
width = 4.5; %cm
height = 9; %cm
fontSize = 8;
LineWidth = 1;
papersize = [21.6, 2.79]; % letter format get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print('fig12A','-depsc', '-r300', '-painters')