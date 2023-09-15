function [super_noise] = expectation_conditioned_on_constituent_parts_2D(actual, raster,neuron_window, time_window)
     
    [expected,unscaled_expected, p] = theoretical_expectation_2D(raster, neuron_window, time_window)
    
     super_noise = [
    
            expected(1),
            expected(2),  
            (expected(3) / (expected(2) * expected(1))) * (actual(2) * actual(1)),   
            expected(4),  
            (expected(5) / (expected(4) * expected(1))) * (actual(4) * actual(1)),  
            expected(6), 
            (expected(7) / sqrt(expected(2) * expected(4) * expected(6))) * sqrt(actual(2) * actual(4) * actual(6)),  
            (expected(8) / sqrt(expected(2) * expected(4) * expected(6))) * sqrt(actual(2) * actual(4) * actual(6)),  
            (expected(9) / sqrt(expected(2) * expected(6)^2)) * sqrt(actual(2) * actual(6)^2), 
            (expected(10) / sqrt(expected(2) * expected(6)^2)) * sqrt(actual(2) * actual(6)^2),  
            (expected(11) / sqrt(expected(2) * expected(6)^2)) * sqrt(actual(2) * actual(6)^2), 
            (expected(12) / sqrt(expected(4) * expected(6)^2)) * sqrt(actual(4) * actual(6)^2),  
            (expected(13) / sqrt(expected(4) * expected(6)^2)) * sqrt(actual(4) * actual(6)^2), 
            (expected(14) / (expected(6)^(3/2))) * (actual(6)^(3/2))  
        ];
     super_noise = super_noise';

end