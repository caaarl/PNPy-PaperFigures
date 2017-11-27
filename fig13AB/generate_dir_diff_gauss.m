numPoints = 1000;

plotInd = 1;
axes = [];
figure,
for alpha = [0.2, 0.4, 0.6, 0.8, 1]

    ax = subplot(2,3,plotInd);
    axes = [axes, ax];
    
    a = 1;
%     c = (rand(1, numPoints)*2-1)*alpha*a;
    c = (randn(1, numPoints)/3)*alpha*a;

    % b = 2 * a * sin(0.5 * atan(c/a));
    b = 2 * sin(0.5 * atan(c/a./(2-abs(c/a))));

    histogram(abs(b))
    title(num2str(alpha))
    hold on
    plotInd = plotInd + 1;
end

xlim([0,1])
linkaxes(axes, 'x')

% % Print setup
% width = 8; %cm
% height = 4.5; %cm
% fontSize = 8;
% LineWidth = 1;
% papersize = [21.6, 2.79]; % letter format get(gcf, 'PaperSize');
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(gcf,'PaperPosition', myfiguresize);
% print('normal_d_distr','-depsc')