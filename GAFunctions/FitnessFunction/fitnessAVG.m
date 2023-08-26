function fitvalue = fitnessAVG(genVec, scale,num)
genImage = uint8(makeMatrix3D_v2(genVec,scale,num));
propSize = abs(1200 - AverageGrainSize(genImage))/1200;
blackFrac = abs(0.5 - nnz(genVec == 0)/(num^3))/0.5;
whiteFrac = abs(0.3 - nnz(genVec == 255)/(num^3))/0.3;
single = singleSubVol(genVec, num);
fitvalue = propSize^2 + blackFrac^2 + whiteFrac^2 + single^2;
end
