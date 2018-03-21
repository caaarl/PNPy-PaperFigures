numPoints = 1000;

plotInd = 1;
axes = [];
figure,
for alpha = [0.2, 0.4, 0.6, 0.8, 1]

    ax = subplot(2,3,plotInd);
    axes = [axes, ax];
    
    a = 1;
    c = (rand(1, numPoints)*2-1)*alpha;
%     c = (randn(numPoints)/3)*alpha;

    b = 2 * a * sin(0.5 * atan(c/a));

    histogram(abs(b))
    title(num2str(alpha))
    hold on
    plotInd = plotInd + 1;
end

linkaxes(axes, 'x')