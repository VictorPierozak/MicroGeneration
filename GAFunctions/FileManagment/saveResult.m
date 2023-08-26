function out = saveResult(structure, num, scale)
structure = makeMatrix3D(structure, scale, num );
structureSize = (num*scale);
structure = uint8(structure);
baseName = "Result/result_";
for n = 1:structureSize

slice = uint8(structure(:,:, n));
if n > 9
imwrite(slice, strcat(baseName, num2str(n), ".bmp"));
else
imwrite(slice, strcat(baseName,"0",num2str(n), ".bmp"));
end

end

out = 1;
end

