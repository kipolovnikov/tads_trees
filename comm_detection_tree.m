function [num, lbord, rbord, medsizes, mods] = comm_detection_tree(m)

    % calculating optimal value of g
    g_opt = optimal_g(m);

    % corrected value of g_opt
    g_opt_corr = corrected_g_opt(m, g_opt);

    % obtaining number of clusters & positions of boundaries & median sizes
    % of clusters at different levels & modularity scores at diff levels
    [num, lbord, rbord, medsizes, mods] = spectral_modularity_tree(m, g_opt_corr);
    
end