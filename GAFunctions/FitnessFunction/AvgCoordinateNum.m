function coornum = AvgCoordinateNum(vec, num)
[~,sizevec] = size(vec);

coornum = 0;
    for i = 2:(sizevec-1)

    if vec(i-1) == vec(i)
        coornum = coornum + 1;
    end
    if vec(i+1) == vec(i)
        coornum =  coornum + 1;
    end

    if i > num && i < sizevec-num
        if vec(i-num) == vec(i) 
        coornum = coornum + 1;
        end
        if vec(i+num) == vec(i)
        coornum = coornum + 1;
        end
    end

    if i > num*num && i < sizevec-num*num
        if vec(i-num*num) == vec(i) 
        coornum = coornum + 1;
        end
        if vec(i+num*num) == vec(i)
        coornum = coornum + 1;
        end
    end

    end

coornum = coornum/num^3;

end

