function extended = extend(vector, num, factor)
vector = repmat(vector, factor,1);
vector = reshape(vector, [], 1);
extended = zeros(num^3,1);
iter = 1;
for  z = 1:2*num^2:num^3
    for p = 1:2*num:num^2
        extended((z+p-1):(z+p+num-2),1) = vector(iter:(iter+num-1),1);
        extended((z+p-1+num):(z+p+2*num-2),1) = vector(iter:(iter+num-1),1);
        iter = iter + num;
    end
    extended((z+num^2):(z+2*num^2-1),1) = extended(z:(z+num^2-1));
end
extended = extended';

end

