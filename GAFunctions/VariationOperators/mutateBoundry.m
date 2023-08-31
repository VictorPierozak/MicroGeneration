function mutationChildren = mutateBoundry(parents, options, nvars,... 
FitnessFcn, state, thisScore, thisPopulation)

mutationChildren = thisPopulation(parents,:);
[chNum, ~] = size(mutationChildren);
for i = 1:chNum
    boundries = selectPhaseBoundries(mutationChildren(i, :));
    mutationChildren(i,:) = applyMutationOnBoundries(boundries, mutationChildren(i,:));
end

end

