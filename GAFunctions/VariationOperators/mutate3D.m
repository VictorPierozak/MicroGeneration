function mutationChildren = mutate3D(parents, options, nvars,... 
FitnessFcn, state, thisScore, thisPopulation, valueArray, scale, num, f)

killValue = 8;

toMutate = thisPopulation(parents,:);
[chNum, len] = size(toMutate);
m = initpop(chNum, len, valueArray);
mutNum = fix(rand(chNum,1)*(f) + 1 + num/2);
for i = 1:chNum
    mutated = toMutate(i, :);
    mutPos = rand(1, mutNum(i));
    mutPos = fix(mutPos*(nvars-1)+1);
    mutated(1, mutPos) = m(i, mutPos);
    %if FitnessFcn(mutated) < killValue
        toMutate(i,:) = mutated;
    %end
end

mutationChildren = toMutate;
end

