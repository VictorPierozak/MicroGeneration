function converted = convertForTau(structure3D)
[structureSize, ~, ~] = size(structure3D);
converted = structure3D;
for x = 1:structureSize
    for y = 1:structureSize
        for z = 1:structureSize
            if converted(x,y,z) == 127
                converted(x,y,z) = 1;
            elseif converted(x,y,z) == 255
                converted(x,y,z) = 2;
            end
        end
    end
end
end

