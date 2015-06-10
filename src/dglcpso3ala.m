function [O gbest] = dglcpso3ala(f,n,nn,lbound,ubound,dim,nop,w,a,b,endgen,alpha,beta)
    % f as fitness function
    % enter f as char like : #
    % n as neighborhood type
    % enter n as char like : 'fn'
    % nn as neighborhood particle number
    % lbound as lower bound of fitness function
    % lbound as upper bound of fitness function
    % dim as dimension
    % nop as number of particles
    % w as weigth interia
    % a as weight index
    % b as weight index
    % endgen as maximum generation number
    % alpha as reward impulse
    % beta as  penalty impulse
    % ACKLEY & ROSENBROCK
%     strgen = 2000;%f4
    strgen = 300;%f5
    deadgen = strgen:20:endgen; % dead particles generation
    nopdead = 5; % number of dead particle generation
    I = 1; % initial index of deadgen
    dead = zeros(nop,1); % dead matrix index
    c = 2;
    bnd = [lbound ubound];
    x = randi(bnd,[nop dim]); % initialize position
    spd = rand([nop dim]); % initial velocity
    x(:,end+1)=0;
    for i=1:nop
        x(i,end) = ff({f},x(i,1:end-1));
    end
%     maxval = ff({1},repmat(ubound,1,dim)); % MAX fittness value
%     maxval = 5*10^5;
    pbest = x; %initialize Best Particle Position
    [bst ind] = min(x(:,end));
    gbest = x(ind,:); % initialize global best position
    O = gbest(end);
    lbest = zeros(size(x)); % initialize local best position
    for i=1:nop
        lbest(i,:) = neighborhood({n},x,i,nn);
    end
    % Learning Automata
    act1 = 1; % action 1 : global and local combined search
    act2 = 2; % action 2 : centralized local best search 
    act3 = 3; % action 3 : centralized global best search
    action = [act1 act2 act3];
    r = size(action,2);
    p = repmat(1/r,nop,r);
    for i=1:endgen
        for j=1:nop
            rw = randsrc(1000,1,[action; p(j,:)]); % roller wheel 
            act = rw(ceil(1000*rand()));  % selected action
            switch (act)
                case {1}    % global and local combined search
                    spd(j,:) = w.* spd(j,:)+ rand(1,dim).*(a+1/(endgen+1-i)).*(pbest(j,1:end-1)-x(j,1:end-1))...
                                           +       (b-1/(endgen+1-i)).*(lbest(j,1:end-1)-x(j,1:end-1))...
                                           + rand(1,dim).*c                 .*(gbest(1:end-1)-x(j,1:end-1));
                case{2}     % centralized local best search 
                    spd(j,:) =             + rand(1,dim).*(a+1/(endgen+1-i)).*(pbest(j,1:end-1)-x(j,1:end-1))...
                                           +       (b-1/(endgen+1-i)).*(lbest(j,1:end-1)-x(j,1:end-1));
                case{3}     % centralized global best search 
                    spd(j,:) =               rand(1,dim).*(a+1/(endgen+1-i)).*(pbest(j,1:end-1)-x(j,1:end-1))...
                                           + rand(1,dim).*c                 .*(gbest(1:end-1)-x(j,1:end-1));            
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
        for j=1:nop
            tf = checkbnd(x(j,1:end-1),lbound,ubound);
            if tf
                lbesttmp = neighborhood({n},x,j,nn);
                if lbesttmp<lbest(j,end)
                    lbest(j,:) = lbesttmp;
                end
            end
        end
        [bst ind] = min(x(:,end));
        tf = checkbnd(x(ind,1:end-1),lbound,ubound);    
        if (tf && bst<gbest(end))
            gbest = x(ind,:);
        end
        fprintf('iteration=%u,gbest=%e\n',i,gbest(end))
        O = [O gbest(end)];
% ===============  ACKLEY & ROSENBROCK Specific =================
        if f == 4 || f == 5
            if i == deadgen(I)
                x = rndadd(x,f,nopdead,dead,gbest);
                I = I+1;
            end
        end
% ==============  GCA  ===============
    
%     hold on
%     axis([0 endgen 0 2.5*10^9])
% %     axis([3000 endgen 0 0.5])
%     title('convergence comparsion')
%     xlabel('iteration')
%     ylabel('gbest value')
%     plot(i,gbest(end),'--r','MarkerSize',3)    
%     if mod(i,50) == 0
%         plot(i,gbest(end),'--rs','MarkerSize',3)
%     end
%     drawnow
%     hold off

% ==============  LA GCA  ===============
% 	hold on
%     axis([0 endgen 0 1])
%     title('3-action probability selection variance')
%     xlabel('iteration')
%     ylabel('probability')
%     plot(i,p(1,1),'--k','MarkerSize',3)
%     plot(i,p(1,2),'--b','MarkerSize',3)
%     plot(i,p(1,3),'--r','MarkerSize',3)
%     if mod(i,50) == 0
%         plot(i,p(1,1),'--ks','MarkerSize',3)
%         plot(i,p(1,2),'--bh','MarkerSize',3)
%         plot(i,p(1,3),'--ro','MarkerSize',3)
%     end
% %     h = legend('combined search','local search','global search',1);
% %     set(h,'Interpreter','none')
%     drawnow
%     hold off
    
    end
end
