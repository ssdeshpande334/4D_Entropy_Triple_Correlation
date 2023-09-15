function motif_class = network_motif_classification(n1, n2, t1, t2)
     vec = [0,n1,n2];
     n_distinct_neurons = length(unique(vec));
     
     list_vec = [n1, n2, t1, t2];
    if t1 < t2
        n1_new = n1;
        n2_new = n2;
        t1_new = t1;
        t2_new = t2;
    else
        n1_new = n2;
        n2_new = n1;
        t1_new = t2;
        t2_new = t1;
    end
    
    n1 = n1_new;
    n2 = n2_new;
    t1 = t1_new;
    t2 = t2_new;
    
    % Assume below that t1 <= t2

    if n_distinct_neurons == 1
        % All neurons are the same
        motif_class= one_neuron_motif_classification(n1, n2, t1, t2);
    elseif n_distinct_neurons == 2
        motif_class= two_neuron_motif_classification(n1, n2, t1, t2);
    elseif n_distinct_neurons == 3
        motif_class= three_neuron_motif_classification(n1, n2, t1, t2);
    else
        error('Invalid number of same neurons');
    end
end
