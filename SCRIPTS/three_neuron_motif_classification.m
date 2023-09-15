function motif = three_neuron_motif_classification(n1, n2, t1, t2)
    % All neurons are distinct
    % Assume t1 <= t2
     vec = [0,t1,t2];
     n_distinct_times = length(unique(vec));

    if n_distinct_times == 1
        % All neurons are distinct, times are the same
        motif= 5;
    elseif n_distinct_times == 2
        if t1 == 0
            motif= 13;
        elseif t2 == 0
            motif= 12;
        elseif t1 == t2
            if t1 > 0
                motif= 12;
            else
                motif= 13;
            end
        else
            error("Shouldn't get here.")
        end 
    elseif n_distinct_times == 3
        motif= 14; 
    else
        error("Invalid number of same times: $n_distinct_times")
    end
end

