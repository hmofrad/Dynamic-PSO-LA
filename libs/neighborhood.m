function lbesttmp = neighborhood(c,x,i,nbr)
    nbr_type = c{1};
    nop = size(x,1);
    switch (nbr_type)
        case {'Fixed Neighborhood','fn'}
            nnbr = ceil(nop/nbr);   % number of neighborhood
            nind = ceil(i/nbr);     % current particle neighborhood index
            nlind = (nnbr-1)*nbr+1; % neighborhood last index
            if nnbr == nind
                [bst ind] = min(x(nlind:end,end));
                xtmp = x(nlind:end,:);
                lbesttmp =  xtmp(ind,:);
            else
                [bst ind] = min(x(nbr*(nind-1)+1:nbr*nind,end));
                xtmp = x(nbr*(nind-1)+1:nbr*nind,:);
                lbesttmp =  xtmp(ind,:);        
            end
        case {'Ring Neighborhood','rn'}
            if (i<=floor(nbr/2))
                [bst ind] = min([x(i-floor(nbr/2)+end:end,end);x(1:1+floor(nbr/2),end)]);
                xtmp = [x(i-floor(nbr/2)+end:end,:);x(1:1+floor(nbr/2),:)];
                lbesttmp =  xtmp(ind,:);    
            elseif (i+floor(nbr/2)>=nop)
                [bst ind] = min([x(i-floor(nbr/2):end,end);x(1:i+floor(nbr/2)-end,end)]);
                xtmp = [x(i-floor(nbr/2):end,end);x(1:i+floor(nbr/2)-end,end)];
                lbesttmp =  xtmp(ind,:);    
            else
                [bst ind] = min(x(i-floor(nbr/2):i+floor(nbr/2),end));
                xtmp = x(i-floor(nbr/2):i+floor(nbr/2),end);
                lbesttmp =  x(ind,:);
            end
        case {'Euclidean Neighborhood','en'}
            ecudist = @(x,y) sum((x-y).^2,2).^0.5; % compute euclidean distance
            xtmp = repmat(x(i,1:end-1),nop,1);
            d = ecudist(xtmp,x(:,1:end-1));
            [B ind] = sort(d);
            xtmp = x(ind(1:nbr),:);
            [C ind] = min(xtmp(:,end));
            lbesttmp =  x(ind,:);
    end
