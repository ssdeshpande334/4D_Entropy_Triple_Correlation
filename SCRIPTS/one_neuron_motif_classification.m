function motif = one_neuron_motif_classification(n1, n2, t1, t2)
    % All neurons are the same
    % Assume t1 <= t2
    vec = [0,t1,t2];
    
    n_distinct_times = length(unique(vec));
    if n_distinct_times == 1
        % All neurons and times are the same
        motif = 1;
    elseif n_distinct_times == 2
        % All neurons are the same, two times
        motif = 2;
    elseif n_distinct_times == 3
        % All neurons are the same, all times distinct
        motif = 3;
    else
        disp('error: Invalid number of distinct times')
    end
end