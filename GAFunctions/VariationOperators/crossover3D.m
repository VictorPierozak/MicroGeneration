function xoverKids = crossover3D(parents, options, nvars, FitnessFcn, ...
    unused,thisPopulation, scale, num, initialDetail)

toCross = thisPopulation(parents,:);
[~,parNum] = size(parents);
offsprings = zeros(parNum/2, nvars);

j = 1;
for i = 1:2:parNum
f1 = FitnessFcn(toCross(i,:));
f2 = FitnessFcn(toCross(i+1,:));
prob1 = 1 - f1/(f1+f2);
prob2 = 1 - f2/(f1+f2);
heritanceVec = rand(nvars,1);
if num > initialDetail
    heritanceVec = heritanceVec + (prob1-prob2)/2;
end
offsprings(j, heritanceVec > 0.5) = toCross(i, heritanceVec > 0.5);
offsprings(j, heritanceVec <= 0.5) = toCross(i+1, heritanceVec <= 0.5);
j=j+1;
end

xoverKids = offsprings;

end

