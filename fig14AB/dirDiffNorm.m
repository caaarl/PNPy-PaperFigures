function dirDiffNorms = dirDiffNorm( dataSegCell, segLength)

    dirDiffNorms = [];
    for dataInd = 1:length(dataSegCell)

        dataSeg = dataSegCell{dataInd};

        dirs = diff(dataSeg,1);
        dirDiff = diff(dirs,1);

        dirDiffNorm = sqrt(sum(dirDiff.^2, 2));

        dirDiffNorms = [dirDiffNorms; dirDiffNorm];

    end
    
    dirDiffNorms = dirDiffNorms/segLength;

end

