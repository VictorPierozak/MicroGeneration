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
    A = ones(scale,scale,scale);
    result = zeros(scale*num);
    for r = 1:num
        for w = 1:num
            for k = 1:num
                result(((k-1)*scale+1):(k*scale), ((r-1)*scale+1):(r*scale), ((w-1)*scale+1):(w*scale)) = A.*vector((num*num)*(r-1)+(num)*(w-1)+k);
            end
        end
    end
          

end

