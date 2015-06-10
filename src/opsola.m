function [O gbest] = opsola(f,lbound,ubound,dim,nop,w,endgen,alpha,beta)
    % f as fitness function
    % enter f as char like : #
    % lbound as lower bound of fitness function
    % lbound as upper bound of fitness function
    % dim as dimension
    % nop as number of particles
    % w as weigth interia
    % endgen as maximum generation number
    % alpha as reward impulse
    % beta  as penalty impulse
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
    O = gbest(end);
    c1 = 2; c2 = 2;
    % Learning Automata
    act1 = 1; % action 1  : global search
    act2 = 2; % action 2  : centralized global best search 
    action = [act1 act2];
    r = size(action,2);
    p = repmat(1/r,nop,r);
    for i=1:endgen
        for j=1:nop
            rw = randsrc(1000,1,[action; p(j,:)]); % roller wheel 
            act = rw(ceil(1000*rand()));  % selected action
            switch (act)
                case {1}     % global search
                    spd(j,:) = w.* spd(j,:)+c1.*rand(1,dim).*(pbest(j,1:end-1) - x(j,1:end-1))...
                                           +c2.*rand(1,dim).*(gbest(1:end-1)   - x(j,1:end-1));
                case {2}      % centralized global best search 
                    spd(j,:) =  c1.*rand(1,dim).*(pbest(j,1:end-1) - x(j,1:end-1))...
                               +c2.*rand(1,dim).*(gbest(1:end-1)   - x(j,1:end-1));
            end
            xtmp = x(j,end);
            ind = find(abs(spd)>0.3);
            spd(ind) = 0.3*sign(spd(ind));
            x(j,1:end-1) = x(j,1:end-1)+spd(j,:);
            x(j,end) = ff({f},x(j,1:end-1));
            if x(j,end) < xtmp
               p(j,act) = p(j,act) + alpha.*(1 - p(j,act)); % desired action
               p(j,action ~=act) = (1-alpha)*p(j,action ~=act);
            else
                p(j,act) = (1 - beta).*p(j,act); % non-desired action
                p(j,action ~=act) = (beta/(r-1))+(1-beta).*p(j,action ~=act);            
            end
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
        O = [O gbest(end)];
    end
end
