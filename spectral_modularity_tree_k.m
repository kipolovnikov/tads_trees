function [numc, lbord, rbord, medsizes, mods] = spectral_modularity_tree_k(m, c, alpha)       

% k-modularity algorithm run on a sparse matrix
% with a cutoff for the minimal cluster size (l > 3)

    n = length(m);

    p = zeros(n, n);
    
    for i = 1:n
        for j = i+1:n
            p(i, j) = c/((abs(i-j))^alpha);
            p(j, i) = p(i, j);
        end
    end
    
    % calculating modularity matrix
    b = m;
    for i = 1:n
        b(i, i) = 0;
        for j = i+1:n
            b(i, j) = m(i, j) - p(i, j);
            b(j, i) = b(i, j);
        end
    end
    
    % hierarchical splitting till there is no more positive gain in modularity
    % score ($num$ equals -1)
    lbord = [1];
    rbord = [n];

    i = 1;
    lev = 0;
    medsizes = [];
    mod = [1];
    mods = [];
    meanmod = 0;
    while mod ~= [0]
        lev = lev + 1;
        medsizes(lev) = median(rbord-lbord)+1;
        lbord_new = [];
        rbord_new = [];
        mod = [0];
        k = 1;
        k1 = 1;
%         lev
%         lbord
%         rbord
        for i = 1:length(lbord)
           [num, md] = division_tree(b, lbord(i), rbord(i));
           %numv(i) = num;
           if num > 3 && rbord(i)-lbord(i)+1-num > 3
               mod(k1) = md;
               k1 = k1 + 1;
               lbord_new(k) = lbord(i);
               lbord_new(k+1) = lbord(i)+num;
               rbord_new(k) = lbord(i)+num-1;
               rbord_new(k+1) = rbord(i);

               k = k + 2;
           else
               lbord_new(k) = lbord(i);
               rbord_new(k) = rbord(i);

               k = k + 1;
           end
        end
        
        meanmod = mean(mod);
        mods(lev) = meanmod;
        lbord = [];
        rbord = [];
        lbord = lbord_new;
        rbord = rbord_new;

    end
    
        
    
    % CHECK
    numc = 0;
    for i = 1:length(lbord)
        len = rbord(i) - lbord(i) + 1;
        if len <= 3
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA!!!!!"
        end
    end

    
end