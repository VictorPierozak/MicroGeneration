function msg = f_maxStallGeneration(prevmsg, max)

if prevmsg < max
    msg = prevmsg;
else
    msg = max;
end

end

