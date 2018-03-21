function dataSegCell = section(dataCell, segmentLen)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    interpStep = 1;

    dataSegCell = {};
    for dataInd = 1:length(dataCell)

        data = dataCell{dataInd};

        x = data(:,1);
        y = data(:,2);

        xInterp = min(x):interpStep:max(x);
        yInterp = interp1(x,y,xInterp);

        xSeg = [x(1)];
        ySeg = [y(1)];

        ind = 1;
        while true

            r = sqrt((xInterp(ind:end)-xSeg(end)).^2 + (yInterp(ind:end)-ySeg(end)).^2);
            r(r<segmentLen) = inf;
            [minVal, indLoc] = min(r);
            ind = indLoc + ind;

            if minVal == inf || ind > length(xInterp)
                break
            end

            xSeg = [xSeg, xInterp(ind)];
            ySeg = [ySeg, yInterp(ind)];

        end

        dataSeg = [xSeg; ySeg]';

        dataSegCell = horzcat(dataSegCell, {dataSeg});
    end

end

