function [avgTauNi, avgTauYSZ] = tortuosity(structure3D, voxelSize)
PhaDir = [ [0,0,0]; [1,1,1]; [1,1,1]];
dimen = [voxelSize, voxelSize, voxelSize];
converted = convertForTau(structure3D);
tauResult = TauFactor('InLine', 1, 0, converted, PhaDir, dimen);

avgTauNi = (tauResult.Tau_W1.Tau + tauResult.Tau_W2.Tau + tauResult.Tau_W3.Tau)/3;
avgTauYSZ = (tauResult.Tau_G1.Tau + tauResult.Tau_G2.Tau + tauResult.Tau_G3.Tau)/3;

end

