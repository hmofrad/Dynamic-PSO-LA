function tf = checkbnd(x,lbound,ubound)
    tf = true;
    for k=1:length(x)
        if(tf && x(k)>=lbound && x(k)<=ubound)
            tf = true;
        else
            tf = false;
        end
    end
