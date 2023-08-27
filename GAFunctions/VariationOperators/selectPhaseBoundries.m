function boundries = selectPhaseBoundries(microstructureVector)
[~, voxelNumber] = size(microstructureVector);
belongToBoundry = zeros(1, voxelNumber);
    for voxel = 1:voxelNumber
        flag = 0;
        if voxel > 1
            if microstructureVector(1, voxel-1) ~= microstructureVector(1, voxel)
                flag = 1;
            end
        end
        if mod(voxel, num) ~= 0
            if microstructureVector(1, voxel+1) ~= microstructureVector(1, voxel)
                flag = 1;
            end
        end
        if  mod(voxel, num^2) > num
             if microstructureVector(1, voxel-num) ~= microstructureVector(1, voxel)
                 flag = 1;
             end
        end
        if  mod(voxel, num^2) < num^2-num
             if microstructureVector(1, voxel+num) ~= microstructureVector(1, voxel)
                 flag = 1;
             end
        end
        if voxel > num^2
            if microstructureVector(1, voxel-num^2) ~= microstructureVector(1,voxel)
                flag =1;
            end
        end
        if voxel < num^3-num^2
            if microstructureVector(1, voxel+num^2) ~= microstructureVector(1,voxel)
                flag = 1;
            end
        end
        belongToBoundry(voxelNumber) = flag;
    end
    boundries = microstructureVector(belongToBoundry == 1);
end

