function motif = two_neuron_motif_classification(n1, n2, t1, t2)
    % Two neurons are the same (two distinct neurons involved)
    % Assume t1 <= t2
     vec = [0,t1,t2];
     n_distinct_times = length(unique(vec));
     
     if n_distinct_times == 1
        % Two neurons are the same, all times are the same
        motif= 4;
    elseif n_distinct_times == 2
        if (0 == t1 && 0 == n1) || (t1 == t2 && n1 == n2) || (0 == t2 && 0 == n2)
            % Base and first nodes are same; or first and second
            motif= 6;
        elseif (t1 == 0) || (t2 < 0)
            % Synchrony is first
            motif=7;
        elseif (t2 == 0) || (t1 > 0)
            % Synchrony is second
            motif= 8;
        end
    elseif n_distinct_times == 3
        if (n1 == 0)
            if (0 < t1)
                motif=11;
            elseif (t1 < 0)
                if t2 < 0
                    motif=10;
                elseif t2 > 0
                    motif=11;
                else
                    error("Shouldn't be here")
                end
            else
                error("Shouldn't be here")
            end
        elseif (n2 == 0)
            if (t2 < 0)
                motif=9;
            elseif (t2 > 0)
                if t1 > 0 
                    motif=10;
                elseif t1 < 0
                    motif=9;
                else
                    error("Shouldn't be here")
                end
            end
        elseif (n1 == n2)
            if (0 < t1)
                motif= 9;
            elseif (t1 < 0) &&(0 < t2)%(t1 < 0 || t2 > 0)
                motif=10;
            elseif (t2 < 0)
                motif=11;
            else
                error("Shouldn't be here")
            end
        else
            error("Shouldn't get here.")
        end
    else
        error("Invalid number of same times: $n_distinct_times")
    end
end

