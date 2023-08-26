function fitvalue = microQuality(genVec, scale, num, desFracNi, desFracYSZ, desAvgSize, desTauNi, desTauYSZ, voxelSize, power)

  microstructure = uint8(makeMatrix3D_v2(genVec,scale,num));

[avgSize, devSize] = AverageGrainSize(microstructure, voxelSize);
normAvgSize = abs(desAvgSize - avgSize)/desAvgSize;
normDevSize = devSize/desAvgSize;

if num ~= 8 % CHECK CONFIG - INTITIAL DETAIL
[tauNi, tauYSZ] = tortuosity(microstructure, voxelSize);

if isinf(tauNi) || isnan(tauNi)
    tauNi = desTauNi*(-1);
end

if isinf(tauYSZ) || isnan(tauYSZ)
    tauYSZ = desTauYSZ*(-1);
end

normTauNi = abs(desTauNi - tauNi)/desTauNi;
normTauYSZ = abs(desTauYSZ - tauYSZ)/desTauYSZ;

else
normTauNi = 0;
normTauYSZ = 0;
end
 

fracYSZ = abs(desFracYSZ - nnz(genVec == 127)/(num^3))/desFracYSZ;
fracNi = abs(desFracNi - nnz(genVec == 255)/(num^3))/desFracNi;

fitvalue = normAvgSize^power + normDevSize^power + normTauNi^power + normTauYSZ^power + 0.5*(fracYSZ + fracNi)^power;

end

