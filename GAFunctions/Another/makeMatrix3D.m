%num - how many times smaller is edge length of cube than whole structure
%scale - edge length of pixel of one cube (sub-volume)

function result = makeMatrix3D(vector, scale, num)

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
    vector = repmat(vector,scale,1);
    vector = reshape(vector,1,[]);
    vector = repmat(vector,scale,1);
    vector = repmat(vector,1,1,scale);
    result = zeros(num*scale, num*scale, num*scale);

    v = 1;

    for p = 1:scale:num*scale

    for j = 1:scale:num*scale
        result(j:(j+scale-1),:,p:(p+scale-1)) = vector(:, v:(v+num*scale-1),:);
        v = v + num*scale;

    end
    end

end

