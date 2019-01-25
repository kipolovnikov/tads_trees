function [numc, lbord, rbord, medsizes, mods] = spectral_modularity_tree_1(m, g)       

% spectral modularity algorithm run on a sparse matrix
% with a cutoff for the minimal cluster size (l > 3)

    n = length(m);

    % strength of the weighted network
    w = sum(sum(m))/2;

    % vector with nodes degrees
    k = zeros(n, 1);
    for i = 1:n
        k(i) = sum(m(i, :));
    end
    
    % calculating modularity matrix
    b = m;
    for i = 1:n
        for j = i:n
            b(i, j) = m(i, j) - g*k(i)*k(j)/(2*w);
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
    mods = [];
    numv = [1];
    meanmod = 0;
    while allm(numv) == 1
        lev = lev + 1;
        medsizes(lev) = median(rbord-lbord);
        mods(lev) = meanmod;
        lbord_new = [];
        rbord_new = [];
        numv = [];
        mod = [];
        k = 1;
        for i = 1:length(lbord)
           [num, mod(i)] = division_tree(b, lbord(i), rbord(i));
           numv(i) = num;
           if num > 3 && num < rbord(i)-2
               lbord_new(k) = lbord(i);
               lbord_new(k+1) = lbord(i)+num;
               rbord_new(k) = lbord(i)+num-1;
               rbord_new(k+1) = rbord(i);
               k = k + 2;
           else
               lbord_new(k) = lbord(i);
               rbord_new(k) = rbord(i);
           end
        end
        meanmod = mean(mod);
        lbord = lbord_new;
        rbord = rbord_new;
    end
    
    
    % CHECK
    numc = 0;
    for i = 1:length(lbord)
        len = rbord(i) - lbord(i) + 1;
        if len <= 3
            printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA!!!!!")
        end
    end

    
end