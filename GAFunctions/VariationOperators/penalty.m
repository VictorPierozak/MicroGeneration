function kill = penalty(genVec, scale, num, prevSize, prevBlack, prevWhite)
blackFrac = abs(0.4 - nnz(genVec == 0)/(num^3))/0.4;
whiteFrac = abs(0.3 - nnz(genVec == 255)/(num^3))/0.3;
greyFrac = abs(0.3 - nnz(genVec == 127)/(num^3))/0.3;

kill = 0;

if blackFrac > 1.1*prevBlack
    kill = 1;
end
if whiteFrac > 1.1*prevWhite
    kill = 1;
end
if greyFrac > 1.1*prevGrey
    kill = 1;
end

end

