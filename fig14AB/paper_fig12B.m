colors =    [[0.,0.,0.]; ...
            [230., 159., 0.]; ...
            [86., 180., 233.]; ...
            [0., 158., 115.]...
            ]/255;

numPoints = 5000;
alphas = [0.2, 0.6, 1];

plotInd = 1;
axes = [];
figure,
for alpha = alphas

    a = 1;
    
    % normal distribution
    w = (randn(1, numPoints)/3);
    c = 2 * a * sin(0.5 * atan(w/a*alpha./(2.1-alpha)));
    
    ax = subplot(2,length(alphas),plotInd);
    axes = [axes, ax];
    histogram(abs(c), 'EdgeColor', 'none', 'FaceColor', colors(3,:), 'BinWidth', 0.05)
    title([num2str(alpha), ' normal'])
    
    % uniform distribution
    w = (rand(1, numPoints)*2-1)*alpha;
    c = 2 * a * sin(0.5 * atan(w/a*alpha./(2.1-alpha)));

    ax = subplot(2,length(alphas),plotInd + length(alphas));
    axes = [axes, ax];
    histogram(abs(c), 'EdgeColor', 'none', 'FaceColor', colors(1,:), 'BinWidth', 0.05)
    title([num2str(alpha), ' uniform'])
    
    plotInd = plotInd + 1;
    
    xlabel('c (unitless)')
    ylabel('#')
end

xlim([0,1])
linkaxes(axes, 'x')

% Print setup
width = 16; %cm
height = 9; %cm
fontSize = 8;
LineWidth = 1;
papersize = [21.6, 2.79]; % letter format get(gcf, 'PaperSize');
left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf,'PaperPosition', myfiguresize);
print('fig12B','-depsc', '-r300', '-painters')