function pop = initpop(popsize, varnum, collorArray)

[~, collorNum] = size(collorArray);
pop = fix( rand(popsize, varnum)*(collorNum)+1);
for i = 1:int64(collorNum)
pop(pop == i) = collorArray(i);
end

end


