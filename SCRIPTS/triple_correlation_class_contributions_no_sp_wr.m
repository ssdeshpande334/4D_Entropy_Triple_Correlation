function [tricorr_4D_distribution, class_contribution,class_count,contribution]= triple_correlation_class_contributions_no_sp_wr(arr, neuron_window, time_window)

    class_contribution = zeros(1,14);
    class_count = zeros(1,14);

    [N_neurons, N_times] = size(arr);

    tricorr_4D_distribution = zeros(neuron_window+1, time_window+1, neuron_window+1, time_window+1);

    max_time_lag = time_window / 2;

    for n1 = -floor((neuron_window / 2)):floor((neuron_window / 2)) 

        n1Lag = n1 + (neuron_window / 2) + 1;

    for n2 = -floor((neuron_window / 2)):floor((neuron_window / 2)) 

        n2Lag = n2 + (neuron_window / 2) + 1;

    for t1 = -floor((time_window / 2)):floor((time_window / 2))

        t1Lag = t1 + (time_window / 2) + 1;

    for t2 = -floor((time_window / 2)):floor((time_window / 2))

        t2Lag = t2 + (time_window / 2) + 1;

        contribution = 0; 
        
        class = network_motif_classification(n1, n2, t1, t2);
        class_count(1,class) = class_count(1,class) + 1;

                for r = (1-min(0,min(n1,n2))):(N_neurons-max(0,max(n1,n2)))

                for c = max_time_lag+1:N_times-max_time_lag           

                    contribution=contribution + arr(r,c)*arr(r+n1,c+t1)*arr(r+n2,c+t2);
                    
                    tricorr_contribution = arr(r,c)*arr(r+n1,c+t1)*arr(r+n2,c+t2);

                    tricorr_4D_distribution(n1Lag, t1Lag, n2Lag, t2Lag) = tricorr_4D_distribution(n1Lag, t1Lag, n2Lag, t2Lag) + tricorr_contribution; 
                    
                end
                end
                
        class_contribution(1,class) = contribution+class_contribution(1,class);
        
    end
    end
    end
    end

    
end
