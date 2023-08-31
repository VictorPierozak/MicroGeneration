function microstructureVec = applyMutationOnBoundries(boundryPosition, microstructureVec)

numOfMutation = size(boundryPosition,2);
for i = 1:numOfMutation
   if microstructureVec(boundryPosition(i)) == 0
      if rand() < 0.5
          microstructureVec(boundryPosition(i)) = 127;
      else
          microstructureVec(boundryPosition(i)) = 255;
      end
   elseif microstructureVec(boundryPosition(i)) == 127
      if rand() < 0.5
          microstructureVec(boundryPosition(i)) = 0;
      else
          microstructureVec(boundryPosition(i)) = 255;
      end
   else
       if rand() < 0.5
          microstructureVec(boundryPosition(i)) = 0;
       else
          microstructureVec(boundryPosition(i)) = 127;
       end
   end
end

end

