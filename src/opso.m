function [gbest_hist gbest] = opso(f,lbound,ubound,dim,nop,w,endgen)
    % f as fitness function
    % enter f as char like : #
    % lbound as lower bound of fitness function
    % lbound as upper bound of fitness function
    % dim as dimension
    % nop as number of particles
    % w as weigth interia
    % endgen as maximum generation number
    bnd = [lbound ubound];
    x = randi(bnd,[nop dim]); % initialize position
    spd = rand([nop dim]); % initial velocity
    x(:,end+1)=0;
    for i=1:nop
        x(i,end) = ff({f},x(i,1:end-1));
    end
    pbest = x; %initialize Best Particle Position
    [bst ind] = min(x(:,end));
    gbest = x(ind,:); % initialize global best position
    gbest_hist = gbest(end);
    c1 = 2; c2 = 2;
    for i=1:endgen
        for j=1:nop
            spd(j,:) = w.* spd(j,:)+c1.*rand(1,dim).*(pbest(j,1:end-1)-x(j,1:end-1))...
                                   +c2.*rand(1,dim).*(gbest(1:end-1) - x(j,1:end-1));
            ind = find(abs(spd)>0.3);
            spd(ind) = 0.3.*sign(spd(ind));
            x(j,1:end-1) = x(j,1:end-1)+spd(j,:);
            x(j,end) = ff({f},x(j,1:end-1));
            tf = checkbnd(x(j,1:end-1),lbound,ubound);
            if (tf && x(j,end)<pbest(j,end))
                pbest(j,:) = x(j,:);
            end
        end    
        [bst ind] = min(x(:,end));
        tf = checkbnd(x(ind,1:end-1),lbound,ubound);    
        if (tf && bst<gbest(end))
            gbest = x(ind,:);
        end 
        fprintf('iteration=%u,gbest=%e\n',i,gbest(end))
        gbest_hist = [gbest_hist gbest(end)];
    end
end
