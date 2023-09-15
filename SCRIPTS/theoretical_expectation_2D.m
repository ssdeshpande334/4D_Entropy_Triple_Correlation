function [scaled_expectation,expectation, p] = theoretical_expectation_2D(raster, neuron_window, time_window)
        motif_order = [
            1
            2 
            3 
            2
            3
            2
            3
            3
            3
            3
            3
            3
            3
            3];

    p = mean(raster,'all');
    scale = numel(raster);

    t_pm = (time_window);
    n_pm = (neuron_window);

    t_m = floor(t_pm / 2);
    t_p = ceil(t_pm / 2);

    triplet_count = [
        1  
        3*t_pm  
        t_pm*(t_pm-1)  
        3*n_pm  
        n_pm*(n_pm-1)  
        3*n_pm*t_pm  
        4*n_pm*t_p + 2*n_pm*t_m  
        4*n_pm*t_m + 2*n_pm*t_p  
        n_pm*t_p*(t_p-1) + 2*n_pm*t_m*t_p + n_pm*t_m*(t_m-1) 
        n_pm*(t_p)*(t_p-1) + n_pm*(t_m)*(t_m-1) + 2*n_pm*t_m*t_p 
        n_pm*t_m*(t_m-1) + 2*n_pm*t_p*t_m + n_pm*t_p*(t_p-1) 
        n_pm*(n_pm-1)*t_p + 2*n_pm*(n_pm-1)*t_m
        n_pm*(n_pm-1)*t_m + 2*n_pm*(n_pm-1)*t_p
        n_pm*(n_pm-1)*t_pm*(t_pm-1)
    ];


    scaled_expectation = numel(raster) .* (p .^ motif_order) .* triplet_count ./ scale;
    expectation = numel(raster) .* (p .^ motif_order) .* triplet_count;

end
