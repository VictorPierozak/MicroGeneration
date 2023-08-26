function mutationChildren = mutateBoundry(parents, options, nvars,... 
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
    flag = 0;
    for k = 1:mutNum(i)
        if mutPos > 1
            if mutated(1, mutPos-1) ~= mutated(1, mutPos)
                flag = 1;
            end
        end
        if mod(mutPos, num) ~= 0
            if mutated(1, mutPos+1) ~= mutated(1, mutPos)
                flag = 1;
            end
        end
        if  mod(mutPos, num^2) > num
             if mutated(1, mutPos-num) ~= mutated(1, mutPos)
                 flag = 1;
             end
        end
        if  mod(mutPos, num^2) < num^2-num
             if mutated(1, mutPos+num) ~= mutated(1, mutPos)
                 flag = 1;
             end
        end
        if mutPos > num^2
            if mutated(1, mutPos-num^2) ~= mutated(1,mutPos)
                flag =1;
            end
        end
        if mutPos < num^3-num^2
            if mutated(1, mutPos+num^2) ~= mutated(1,mutPos)
                flag = 1;
            end
        end
            if flag == 1
           mutated(1, mutPos) = m(i, mutPos);
            end
    end
    %if FitnessFcn(mutated) < killValue
        toMutate(i,:) = mutated;
    %end
end

mutationChildren = toMutate;
end

