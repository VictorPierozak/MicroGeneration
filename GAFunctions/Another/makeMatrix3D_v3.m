%num - how many times smaller is edge length of cube than whole structure
%scale - edge length of pixel of one cube (sub-volume)

function result = makeMatrix3D_v2(vector, scale, num)

if scale == 1
    %result = reshape(vector, num,num,num);
     v = 1;

    for p = 1:scale:num*scale

    for j = 1:scale:num*scale
        result(j:(j+scale-1),:,p:(p+scale-1)) = vector(:, v:(v+num*scale-1),:);
        v = v + num*scale;

    end
    end
else
    result = zeros(num);
    for r = 1:num
        for w = 1:num
            for k = 1:num
                result(k, r, w) = vector((num*num)*(r-1)+(num)*(w-1)+k);
            end
        end
    end
          

end

