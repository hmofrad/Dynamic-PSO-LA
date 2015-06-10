function x = rndadd(x,f,nopdead,dead,gbest)
	ind = randi([1 size(x,1)],[1 nopdead]); %random selection of dead particles
    for i=1:length(ind)
        if (x(ind(i),end)    ~= gbest(end)) && dead(ind(i)) == 0
            dead(ind(i))      = 1;
            x(ind(i),1:end-1) = randi([-10 10],[1 size(x,2)-1]);
            x(ind(i),end)     = ff({f},x(ind(i),1:end-1));
        else
            dead(ind(i)) = 0;
        end
    end
end
