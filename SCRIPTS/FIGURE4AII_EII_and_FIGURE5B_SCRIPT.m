clearvars
close all
clc

for iiiii = 1:100

    clearvars -except iiiii entropy_raster_4d_pure aovert_minus1_pure...
        neuron_window time_window
    disp(iiiii)
    % Set parameters
    num_rows = 50;
    num_cols = 50;
    num_spikes = 48;
    
    neuron_window = 8; 
    time_window = 6;
    
    row_min = neuron_window / 2 + 1;
    row_max = num_rows - neuron_window / 2;
    col_min = time_window / 2 + 1;
    col_max = num_cols - time_window / 2;
    

    % Initialize spike matrix
    spike_matrix = zeros(num_rows, num_cols);
    
    % Generate random spike indices within specified row and column ranges
    valid_indices = sub2ind([num_rows, num_cols], ...
        randi([row_min, row_max], num_spikes, 1), ...
        randi([col_min, col_max], num_spikes, 1));
    
    % Ensure unique spike locations
    while numel(unique(valid_indices)) < num_spikes
        valid_indices = sub2ind([num_rows, num_cols], ...
            randi([row_min, row_max], num_spikes, 1), ...
            randi([col_min, col_max], num_spikes, 1));
    end
    
    % Place spikes in the matrix
    spike_matrix(valid_indices) = 1;

    raster = spike_matrix;
    epoch_length =size(raster,2) - (time_window + 1);
    max_time_lag = ceil((time_window / 2));
    [N_neurons, N_times] = size(raster);
    
    post_end = 0;
    slice_start = post_end+1;
    
    pre_snippet_raster = raster(:,slice_start:slice_start+max_time_lag-1);
    pre_end = slice_start+max_time_lag-1;
    
    snippet_raster = raster(:, pre_end+1: pre_end+1+epoch_length);
    snip_end = pre_end+1+epoch_length;
    
    post_snippet_raster = raster(:, snip_end+1 : snip_end+1+max_time_lag-1);
    post_end = snip_end+1+max_time_lag-1;
    
    temp_snippet = cat(2,pre_snippet_raster, snippet_raster, post_snippet_raster);
    
    size(temp_snippet)
    [N_neurons, N_times] = size(temp_snippet);
    
    %compute tricorr
    [c3_4D_distribution, actual_contribution,class_count,contribution]= ...
        triple_correlation_class_contributions_no_sp_wr(temp_snippet, neuron_window, time_window);
    
    actual = actual_contribution./(numel(snippet_raster));
    [conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, snippet_raster,neuron_window, time_window);
    
    aovert_minus1 = (actual ./conditioned_expectation ) -1;
    
    c3_n1t1_distribution = c3_4D_distribution(:,:,neuron_window/2 + 1 ,time_window / 2 + 1);
    T_raster_n1t1=sum(sum(c3_n1t1_distribution));
    PDF_raster_n1t1=c3_n1t1_distribution./T_raster_n1t1;  
    
    c3_n2t2_distribution = c3_4D_distribution(neuron_window/2 + 1 ,time_window / 2 + 1,:,:);
    T_raster_n2t2=sum(sum(c3_n2t2_distribution));
    PDF_raster_n2t2=c3_n2t2_distribution./T_raster_n2t2;
    
    T_raster_4d=sum(sum(sum(sum(c3_4D_distribution))));
    PDF_raster_4d=c3_4D_distribution./T_raster_4d;
    
    
    %Compute the 4D Entropy
    bins_4d = [neuron_window+1 time_window+1 neuron_window+1 time_window+1];
    entropy_raster_4d = 0;
    
    for i=1:bins_4d(1)
        for j=1:bins_4d(2)
            for k = 1:bins_4d(3)
                for m = 1:bins_4d(4)
                    temp_4d = PDF_raster_4d;
                    if temp_4d(i,j,k,m)~=0        
                        entropy_raster_4d=entropy_raster_4d-temp_4d(i,j,k,m)*log2(temp_4d(i,j,k,m));       
                    end
                end
            end
    
        end
    end
    
    entropy_raster_4d_pure(iiiii,1) = entropy_raster_4d;
    aovert_minus1_pure(iiiii,:) = aovert_minus1;
end

raster_pure = raster;

%
for iiiii = 1:100

    clearvars -except iiiii entropy_raster_4d_pure aovert_minus1_pure...
        entropy_raster_4d_3A aovert_minus1_3A raster_pure...
        neuron_window time_window

    disp(iiiii)
    % Set parameters
    num_rows = 50;
    num_cols = 50;
    num_spikes = 48;
   
    
    row_min = neuron_window / 2 + 1;
    row_max = num_rows - neuron_window / 2;
    col_min = time_window / 2 + 1;
    col_max = num_cols - time_window / 2;
    
    % Initialize spike matrix
    spike_matrix = zeros(num_rows, num_cols);
    
    % Generate random spike indices within specified row and column ranges
    valid_indices = sub2ind([num_rows, num_cols], ...
        randi([row_min, row_max], num_spikes, 1), ...
        randi([col_min, col_max], num_spikes, 1));
    
    % Ensure unique spike locations
    while numel(unique(valid_indices)) < num_spikes
        valid_indices = sub2ind([num_rows, num_cols], ...
            randi([row_min, row_max], num_spikes, 1), ...
            randi([col_min, col_max], num_spikes, 1));
    end
    
    % Place spikes in the matrix
    spike_matrix(valid_indices) = 1;
    
    % Find empty columns for the synchronous spikes
    empty_cols = find(sum(spike_matrix(row_min:row_max, col_min:col_max)) == 0);
    if isempty(empty_cols)
        error('No empty columns available for placing synchronous spikes.');
    end
    
    % Randomly select a column from the empty columns
    sync_col = col_min + empty_cols(randi(length(empty_cols))) - 1;
    
    % Find empty rows for placing the synchronous spikes
    found = false;
    for i = row_min:(row_max - 2)
        if all(spike_matrix(i:i+2, sync_col) == 0)
            sync_rows = i:i+2;
            found = true;
            break;
        end
    end
    
    if ~found
        error('No empty rows available for placing synchronous spikes.');
    end
    
    % Place the synchronous spikes in the matrix
    spike_matrix(sync_rows, sync_col) = 1;

    raster = spike_matrix;
    epoch_length =size(raster,2) - (time_window + 1);
    max_time_lag = ceil((time_window / 2));
    [N_neurons, N_times] = size(raster);
    
    post_end = 0;
    slice_start = post_end+1;
    
    pre_snippet_raster = raster(:,slice_start:slice_start+max_time_lag-1);
    pre_end = slice_start+max_time_lag-1;
    
    snippet_raster = raster(:, pre_end+1: pre_end+1+epoch_length);
    snip_end = pre_end+1+epoch_length;
    
    post_snippet_raster = raster(:, snip_end+1 : snip_end+1+max_time_lag-1);
    post_end = snip_end+1+max_time_lag-1;
    
    temp_snippet = cat(2,pre_snippet_raster, snippet_raster, post_snippet_raster);
    
    size(temp_snippet)
    [N_neurons, N_times] = size(temp_snippet);
    
    %compute tricorr
    [c3_4D_distribution, actual_contribution,class_count,contribution]= ...
        triple_correlation_class_contributions_no_sp_wr(temp_snippet, neuron_window, time_window);
    
    actual = actual_contribution./(numel(snippet_raster));
    [conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, snippet_raster,neuron_window, time_window);
    
    aovert_minus1 = (actual ./conditioned_expectation ) -1;
    
    c3_n1t1_distribution = c3_4D_distribution(:,:,neuron_window/2 + 1 ,time_window / 2 + 1);
    T_raster_n1t1=sum(sum(c3_n1t1_distribution));
    PDF_raster_n1t1=c3_n1t1_distribution./T_raster_n1t1;  
    
    c3_n2t2_distribution = c3_4D_distribution(neuron_window/2 + 1 ,time_window / 2 + 1,:,:);
    T_raster_n2t2=sum(sum(c3_n2t2_distribution));
    PDF_raster_n2t2=c3_n2t2_distribution./T_raster_n2t2;
    
    T_raster_4d=sum(sum(sum(sum(c3_4D_distribution))));
    PDF_raster_4d=c3_4D_distribution./T_raster_4d;
    
    
    %Compute the 4D Entropy
    bins_4d = [neuron_window+1 time_window+1 neuron_window+1 time_window+1];
    entropy_raster_4d = 0;
    
    for i=1:bins_4d(1)
        for j=1:bins_4d(2)
            for k = 1:bins_4d(3)
                for m = 1:bins_4d(4)
                    temp_4d = PDF_raster_4d;
                    if temp_4d(i,j,k,m)~=0        
                        entropy_raster_4d=entropy_raster_4d-temp_4d(i,j,k,m)*log2(temp_4d(i,j,k,m));       
                    end
                end
            end
    
        end
    end
    
    entropy_raster_4d_3A(iiiii,1) = entropy_raster_4d;
    aovert_minus1_3A(iiiii,:) = aovert_minus1;
end

raster_3A = raster;

%%
for iiiii = 1:100
    clearvars -except iiiii entropy_raster_4d_pure aovert_minus1_pure...
        entropy_raster_4d_3A aovert_minus1_3A...
        entropy_raster_4d_3B aovert_minus1_3B...
        raster_pure raster_3A...
        neuron_window time_window

    disp(iiiii)
    
    % Set parameters
    num_rows = 50;
    num_cols = 50;
    total_spikes = 48; % Total number of spikes
    
    
    % Define boundaries
    row_min = neuron_window / 2 + 1;
    row_max = num_rows - neuron_window / 2;
    col_min = time_window / 2 + 1;
    col_max = num_cols - time_window / 2;
    
    % Define minimum separations
    min_row_sep = (neuron_window / 2) + 1;
    min_col_sep = (time_window / 2) + 1;
    
    % Create spike raster
    spike_raster = zeros(num_rows, num_cols);
    
    % Define section dimensions
    section_rows = floor(num_rows / 3);
    section_cols = floor(num_cols / 4);
    
    % Place 3 spikes in synchrony in each section within the boundaries
    for i = 1:3
        for j = 1:4
            % Define section boundaries
            row_start = (i - 1) * section_rows + 1;
            row_end = i * section_rows;
            col_start = (j - 1) * section_cols + 1;
            col_end = j * section_cols;
            
            % Randomly select a column within the section boundaries
            col = randi([max(col_start + min_col_sep, col_min), min(col_end - min_col_sep, col_max)]);
    
            % Randomly select a row within the section boundaries
            row = randi([max(row_start + min_row_sep, row_min), min(row_end - min_row_sep - 2, row_max - 2)]); % Ensure there's enough space for 3 consecutive rows within the boundaries
            
            % Place 3-spike synchrony
            spike_raster(row:row+2, col) = 1;
        end
    end
    
    % Poisson spike process parameters
    remaining_spikes = total_spikes - 36; % Remaining spikes needed
    
    % Place Poisson spikes within boundaries
    while sum(spike_raster(:)) < total_spikes
        % Generate random row and column within boundaries
        row = randi([row_min, row_max]);
        col = randi([col_min, col_max]);
        
        % Check if the bin is empty
        if spike_raster(row, col) == 0
            % Add a spike
            spike_raster(row, col) = 1;
        end
    end
        
    raster = spike_raster;
    epoch_length = size(raster, 2) - (time_window + 1);
    max_time_lag = ceil(time_window / 2);
    [N_neurons, N_times] = size(raster);
    
    post_end = 0;
    slice_start = post_end + 1;
    
    pre_snippet_raster = raster(:, slice_start:slice_start + max_time_lag - 1);
    pre_end = slice_start + max_time_lag - 1;
    
    snippet_raster = raster(:, pre_end + 1: pre_end + 1 + epoch_length);
    snip_end = pre_end + 1 + epoch_length;
    
    post_snippet_raster = raster(:, snip_end + 1: snip_end + 1 + max_time_lag - 1);
    post_end = snip_end + 1 + max_time_lag - 1;
    
    temp_snippet = cat(2, pre_snippet_raster, snippet_raster, post_snippet_raster);
    
    size(temp_snippet)
    [N_neurons, N_times] = size(temp_snippet);
    
    % Compute tricorr
    [c3_4D_distribution, actual_contribution, class_count, contribution] = ...
        triple_correlation_class_contributions_no_sp_wr(temp_snippet, neuron_window, time_window);
    
    actual = actual_contribution ./ numel(snippet_raster);
    [conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, snippet_raster, neuron_window, time_window);
    
    aovert_minus1 = (actual ./ conditioned_expectation) - 1;
    
    c3_n1t1_distribution = c3_4D_distribution(:,:,neuron_window / 2 + 1, time_window / 2 + 1);
    T_raster_n1t1 = sum(sum(c3_n1t1_distribution));
    PDF_raster_n1t1 = c3_n1t1_distribution ./ T_raster_n1t1;  
    
    c3_n2t2_distribution = c3_4D_distribution(neuron_window / 2 + 1, time_window / 2 + 1, :, :);
    T_raster_n2t2 = sum(sum(c3_n2t2_distribution));
    PDF_raster_n2t2 = c3_n2t2_distribution ./ T_raster_n2t2;
    
    T_raster_4d = sum(sum(sum(sum(c3_4D_distribution))));
    PDF_raster_4d = c3_4D_distribution ./ T_raster_4d;
    
    % Compute the 4D Entropy
    bins_4d = [neuron_window + 1, time_window + 1, neuron_window + 1, time_window + 1];
    entropy_raster_4d = 0;
    
    for i = 1:bins_4d(1)
        for j = 1:bins_4d(2)
            for k = 1:bins_4d(3)
                for m = 1:bins_4d(4)
                    temp_4d = PDF_raster_4d;
                    if temp_4d(i, j, k, m) ~= 0        
                        entropy_raster_4d = entropy_raster_4d - temp_4d(i, j, k, m) * log2(temp_4d(i, j, k, m));       
                    end
                end
            end
        end
    end
    
    entropy_raster_4d_3B(iiiii, 1) = entropy_raster_4d;
    aovert_minus1_3B(iiiii, :) = aovert_minus1;
end

raster_3B = raster;

%%
for iiiii = 1:100
    clearvars -except iiiii entropy_raster_4d_pure aovert_minus1_pure...
        entropy_raster_4d_3A aovert_minus1_3A...
        entropy_raster_4d_3B aovert_minus1_3B...
        entropy_raster_4d_3C aovert_minus1_3C...
        raster_pure raster_3A raster_3B...
        neuron_window time_window

    disp(iiiii)
    
    % Set parameters
    num_rows = 50;
    num_cols = 50;
    total_spikes = 48; % Total number of spikes
   
    
    % Define boundaries
    row_min = neuron_window / 2 + 1;
    row_max = num_rows - neuron_window / 2;
    col_min = time_window / 2 + 1;
    col_max = num_cols - time_window / 2;
    
    % Define minimum separations
    min_row_sep = (neuron_window / 2) + 1;
    min_col_sep = (time_window / 2) + 1;
    
    % Create spike raster
    spike_raster = zeros(num_rows, num_cols);
    
    % Define section dimensions
    section_rows = floor(num_rows / 3);
    section_cols = floor(num_cols / 2);
    
    % Place 3 spikes in synchrony in each section within the boundaries
    for i = 1:3
        for j = 1
            % Define section boundaries
            row_start = (i - 1) * section_rows + 1;
            row_end = i * section_rows;
            col_start = (j - 1) * section_cols + 1;
            col_end = j * section_cols;
            
            % Randomly select a column within the section boundaries
            col = randi([col_start + min_col_sep, col_end - min_col_sep]);
    
            % Randomly select a row within the section boundaries
            row = randi([row_start + min_row_sep, row_end - min_row_sep - 2]); % Ensure there's enough space for 3 consecutive rows within the boundaries
            
            % Place 3-spike synchrony
            spike_raster(row:row+2, col) = 1;
        end
    end
    
    % Poisson spike process parameters
    lambda = total_spikes - 9; % Remaining spikes needed
    
    % Place Poisson spikes within boundaries
    while sum(spike_raster(:)) < total_spikes
        % Generate random row and column within boundaries
        row = randi([row_min, row_max]);
        col = randi([col_min, col_max]);
        
        % Check if the bin is empty
        if spike_raster(row, col) == 0
            % Add a spike with Poisson probability
            if rand < lambda / (num_rows * num_cols)
                spike_raster(row, col) = 1;
            end
        end
    end
    
    raster = spike_raster;
    epoch_length = size(raster, 2) - (time_window + 1);
    max_time_lag = ceil(time_window / 2);
    [N_neurons, N_times] = size(raster);
    
    post_end = 0;
    slice_start = post_end + 1;
    
    pre_snippet_raster = raster(:, slice_start:slice_start + max_time_lag - 1);
    pre_end = slice_start + max_time_lag - 1;
    
    snippet_raster = raster(:, pre_end + 1: pre_end + 1 + epoch_length);
    snip_end = pre_end + 1 + epoch_length;
    
    post_snippet_raster = raster(:, snip_end + 1: snip_end + 1 + max_time_lag - 1);
    post_end = snip_end + 1 + max_time_lag - 1;
    
    temp_snippet = cat(2, pre_snippet_raster, snippet_raster, post_snippet_raster);
    
    size(temp_snippet)
    [N_neurons, N_times] = size(temp_snippet);
    
    % Compute tricorr
    [c3_4D_distribution, actual_contribution, class_count, contribution] = ...
        triple_correlation_class_contributions_no_sp_wr(temp_snippet, neuron_window, time_window);
    
    actual = actual_contribution ./ numel(snippet_raster);
    [conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, snippet_raster, neuron_window, time_window);
    
    aovert_minus1 = (actual ./ conditioned_expectation) - 1;
    
    c3_n1t1_distribution = c3_4D_distribution(:,:,neuron_window / 2 + 1, time_window / 2 + 1);
    T_raster_n1t1 = sum(sum(c3_n1t1_distribution));
    PDF_raster_n1t1 = c3_n1t1_distribution ./ T_raster_n1t1;  
    
    c3_n2t2_distribution = c3_4D_distribution(neuron_window / 2 + 1, time_window / 2 + 1, :, :);
    T_raster_n2t2 = sum(sum(c3_n2t2_distribution));
    PDF_raster_n2t2 = c3_n2t2_distribution ./ T_raster_n2t2;
    
    T_raster_4d = sum(sum(sum(sum(c3_4D_distribution))));
    PDF_raster_4d = c3_4D_distribution ./ T_raster_4d;
    
    % Compute the 4D Entropy
    bins_4d = [neuron_window + 1, time_window + 1, neuron_window + 1, time_window + 1];
    entropy_raster_4d = 0;
    
    for i = 1:bins_4d(1)
        for j = 1:bins_4d(2)
            for k = 1:bins_4d(3)
                for m = 1:bins_4d(4)
                    temp_4d = PDF_raster_4d;
                    if temp_4d(i, j, k, m) ~= 0        
                        entropy_raster_4d = entropy_raster_4d - temp_4d(i, j, k, m) * log2(temp_4d(i, j, k, m));       
                    end
                end
            end
        end
    end
    
    entropy_raster_4d_3C(iiiii, 1) = entropy_raster_4d;
    aovert_minus1_3C(iiiii, :) = aovert_minus1;
end

raster_3C = raster;

%%
for iiiii = 1:100
    clearvars -except iiiii entropy_raster_4d_pure aovert_minus1_pure...
        entropy_raster_4d_3A aovert_minus1_3A...
        entropy_raster_4d_3B aovert_minus1_3B...
        entropy_raster_4d_3C aovert_minus1_3C...
        entropy_raster_4d_3D aovert_minus1_3D...
        raster_pure raster_3A raster_3B raster_3C...
        neuron_window time_window


    disp(iiiii)
    
        % Set parameters
    num_rows = 50;
    num_cols = 50;
    total_spikes = 48; % Total number of spikes
    
    
    % Define boundaries
    row_min = neuron_window / 2 + 1;
    row_max = num_rows - neuron_window / 2;
    col_min = time_window / 2 + 1;
    col_max = num_cols - time_window / 2;
    
    % Define minimum separations
    min_row_sep = (neuron_window / 2) + 1;
    min_col_sep = (time_window / 2) + 1;
    
    % Create spike raster
    spike_raster = zeros(num_rows, num_cols);
    
    % Define section dimensions
    section_rows = floor(num_rows / 3);
    section_cols = floor(num_cols / 3);
    
    % Place 3 spikes in synchrony in each section within the boundaries
    for i = 1:3
        for j = 1:3
            % Define section boundaries
            row_start = (i - 1) * section_rows + 1;
            row_end = i * section_rows;
            col_start = (j - 1) * section_cols + 1;
            col_end = j * section_cols;
            
            % Randomly select a column within the section boundaries
            col = randi([col_start + min_col_sep, col_end - min_col_sep]);
    
            % Randomly select a row within the section boundaries
            row = randi([row_start + min_row_sep, row_end - min_row_sep - 2]); % Ensure there's enough space for 3 consecutive rows within the boundaries
            
            % Place 3-spike synchrony
            spike_raster(row:row+2, col) = 1;
        end
    end
    
    % Poisson spike process parameters
    lambda = total_spikes - 27; % Remaining spikes needed
    
    % Place Poisson spikes within boundaries
    while sum(spike_raster(:)) < total_spikes
        % Generate random row and column within boundaries
        row = randi([row_min, row_max]);
        col = randi([col_min, col_max]);
        
        % Check if the bin is empty
        if spike_raster(row, col) == 0
            % Add a spike with Poisson probability
            if rand < lambda / (num_rows * num_cols)
                spike_raster(row, col) = 1;
            end
        end
    end

    raster = spike_raster;
    epoch_length = size(raster, 2) - (time_window + 1);
    max_time_lag = ceil(time_window / 2);
    [N_neurons, N_times] = size(raster);
    
    post_end = 0;
    slice_start = post_end + 1;
    
    pre_snippet_raster = raster(:, slice_start:slice_start + max_time_lag - 1);
    pre_end = slice_start + max_time_lag - 1;
    
    snippet_raster = raster(:, pre_end + 1: pre_end + 1 + epoch_length);
    snip_end = pre_end + 1 + epoch_length;
    
    post_snippet_raster = raster(:, snip_end + 1: snip_end + 1 + max_time_lag - 1);
    post_end = snip_end + 1 + max_time_lag - 1;
    
    temp_snippet = cat(2, pre_snippet_raster, snippet_raster, post_snippet_raster);
    
    size(temp_snippet)
    [N_neurons, N_times] = size(temp_snippet);
    
    % Compute tricorr
    [c3_4D_distribution, actual_contribution, class_count, contribution] = ...
        triple_correlation_class_contributions_no_sp_wr(temp_snippet, neuron_window, time_window);
    
    actual = actual_contribution ./ numel(snippet_raster);
    [conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, snippet_raster, neuron_window, time_window);
    
    aovert_minus1 = (actual ./ conditioned_expectation) - 1;
    
    c3_n1t1_distribution = c3_4D_distribution(:,:,neuron_window / 2 + 1, time_window / 2 + 1);
    T_raster_n1t1 = sum(sum(c3_n1t1_distribution));
    PDF_raster_n1t1 = c3_n1t1_distribution ./ T_raster_n1t1;  
    
    c3_n2t2_distribution = c3_4D_distribution(neuron_window / 2 + 1, time_window / 2 + 1, :, :);
    T_raster_n2t2 = sum(sum(c3_n2t2_distribution));
    PDF_raster_n2t2 = c3_n2t2_distribution ./ T_raster_n2t2;
    
    T_raster_4d = sum(sum(sum(sum(c3_4D_distribution))));
    PDF_raster_4d = c3_4D_distribution ./ T_raster_4d;
    
    % Compute the 4D Entropy
    bins_4d = [neuron_window + 1, time_window + 1, neuron_window + 1, time_window + 1];
    entropy_raster_4d = 0;
    
    for i = 1:bins_4d(1)
        for j = 1:bins_4d(2)
            for k = 1:bins_4d(3)
                for m = 1:bins_4d(4)
                    temp_4d = PDF_raster_4d;
                    if temp_4d(i, j, k, m) ~= 0        
                        entropy_raster_4d = entropy_raster_4d - temp_4d(i, j, k, m) * log2(temp_4d(i, j, k, m));       
                    end
                end
            end
        end
    end
    
    entropy_raster_4d_3D(iiiii, 1) = entropy_raster_4d;
    aovert_minus1_3D(iiiii, :) = aovert_minus1;
end

raster_3D = raster;

%%
for iiiii = 1:100

    clearvars -except iiiii entropy_raster_4d_pure aovert_minus1_pure...
    entropy_raster_4d_3A aovert_minus1_3A...
    entropy_raster_4d_3B aovert_minus1_3B...
    entropy_raster_4d_3C aovert_minus1_3C...
    entropy_raster_4d_3D aovert_minus1_3D...
    entropy_raster_4d_3E aovert_minus1_3E...
    raster_pure raster_3A raster_3B raster_3C raster_3D raster_3E...
        neuron_window time_window

% Set parameters
num_rows = 50;
num_cols = 50;
num_sections = 16; % Number of sections
num_consecutive_spikes_per_section = 3; % 3-spike in consecutive rows in each section

% Define boundaries
row_min = neuron_window / 2 + 1;
row_max = num_rows - neuron_window / 2;
col_min = time_window / 2 + 1;
col_max = num_cols - time_window / 2;

% Create spike raster
spike_raster = zeros(num_rows, num_cols);

% Define the column ranges with at least 5 columns of separation
col_ranges = [
    col_min, floor(num_cols/4) - 5;
    floor(num_cols/4) + 5, 2*floor(num_cols/4) - 5;
    2*floor(num_cols/4) + 5, 3*floor(num_cols/4) - 5;
    3*floor(num_cols/4) + 5, col_max;
];

% Define the row ranges
row_ranges = [
    row_min, floor(num_rows/4);
    floor(num_rows/4)+1, 2*floor(num_rows/4);
    2*floor(num_rows/4)+1, 3*floor(num_rows/4);
    3*floor(num_rows/4)+1, row_max;
];

% To track columns with consecutive spikes
cols_with_consecutive_spikes = [];

% Place 3-spike consecutive rows set in each section, grouped by sets of 4 sections sharing the same column
for g = 1:4
    placed = false;
    while ~placed
        % Randomly select a column within the allowed range for the current group
        col = randi([col_ranges(g, 1), col_ranges(g, 2)]);
        
        % Check if the column is already used
        if ismember(col, cols_with_consecutive_spikes)
            continue;
        end
        
        rows_valid = true;
        row_positions = zeros(4, 1);

        for s = 1:4
            section_idx = (g - 1) * 4 + s;
            % Randomly select a position for the first spike of the consecutive set
            row = randi([row_ranges(s, 1), row_ranges(s, 2)]);
            row_positions(s) = row;

            % Check if the consecutive spikes fit within the boundaries
            if spike_raster(row, col) == 1 || spike_raster(row + 1, col) == 1 || spike_raster(row + 2, col) == 1
                rows_valid = false;
                break;
            end
        end

        if rows_valid
            for s = 1:4
                row = row_positions(s);
                spike_raster(row, col) = 1; % First spike
                spike_raster(row + 1, col) = 1; % Second spike
                spike_raster(row + 2, col) = 1; % Third spike
            end
            cols_with_consecutive_spikes = [cols_with_consecutive_spikes; col]; % Track this column
            placed = true;
        end
    end
end

    raster = spike_raster;
    epoch_length = size(raster, 2) - (time_window + 1);
    max_time_lag = ceil(time_window / 2);
    [N_neurons, N_times] = size(raster);
    
    post_end = 0;
    slice_start = post_end + 1;
    
    pre_snippet_raster = raster(:, slice_start:slice_start + max_time_lag - 1);
    pre_end = slice_start + max_time_lag - 1;
    
    snippet_raster = raster(:, pre_end + 1: pre_end + 1 + epoch_length);
    snip_end = pre_end + 1 + epoch_length;
    
    post_snippet_raster = raster(:, snip_end + 1: snip_end + 1 + max_time_lag - 1);
    post_end = snip_end + 1 + max_time_lag - 1;
    
    temp_snippet = cat(2, pre_snippet_raster, snippet_raster, post_snippet_raster);
    
    size(temp_snippet)
    [N_neurons, N_times] = size(temp_snippet);
    
    % Compute tricorr
    [c3_4D_distribution, actual_contribution, class_count, contribution] = ...
        triple_correlation_class_contributions_no_sp_wr(temp_snippet, neuron_window, time_window);
    
    actual = actual_contribution ./ numel(snippet_raster);
    [conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, snippet_raster, neuron_window, time_window);
    
    aovert_minus1 = (actual ./ conditioned_expectation) - 1;
    
    c3_n1t1_distribution = c3_4D_distribution(:,:,neuron_window / 2 + 1, time_window / 2 + 1);
    T_raster_n1t1 = sum(sum(c3_n1t1_distribution));
    PDF_raster_n1t1 = c3_n1t1_distribution ./ T_raster_n1t1;  
    
    c3_n2t2_distribution = c3_4D_distribution(neuron_window / 2 + 1, time_window / 2 + 1, :, :);
    T_raster_n2t2 = sum(sum(c3_n2t2_distribution));
    PDF_raster_n2t2 = c3_n2t2_distribution ./ T_raster_n2t2;
    
    T_raster_4d = sum(sum(sum(sum(c3_4D_distribution))));
    PDF_raster_4d = c3_4D_distribution ./ T_raster_4d;
    
    % Compute the 4D Entropy
    bins_4d = [neuron_window + 1, time_window + 1, neuron_window + 1, time_window + 1];
    entropy_raster_4d = 0;
    
    for i = 1:bins_4d(1)
        for j = 1:bins_4d(2)
            for k = 1:bins_4d(3)
                for m = 1:bins_4d(4)
                    temp_4d = PDF_raster_4d;
                    if temp_4d(i, j, k, m) ~= 0        
                        entropy_raster_4d = entropy_raster_4d - temp_4d(i, j, k, m) * log2(temp_4d(i, j, k, m));       
                    end
                end
            end
        end
    end

    entropy_raster_4d_3E(iiiii, 1) = entropy_raster_4d;
    aovert_minus1_3E(iiiii, :) = aovert_minus1;
end

raster_3E = raster;
    %%
    close all
entropy_raster_4d_3E_rep = entropy_raster_4d_3E; 

motifV_pure = aovert_minus1_pure(:,6);
motifV_3A = aovert_minus1_3A(:,6);
motifV_3C = aovert_minus1_3C(:,6);
motifV_3D = aovert_minus1_3D(:,6);
motifV_3E = aovert_minus1_3E(:,6);

motif4_pure = aovert_minus1_pure(:,5);
motif4_3A = aovert_minus1_3A(:,5);
motif4_3C = aovert_minus1_3C(:,5);
motif4_3D = aovert_minus1_3D(:,5);
motif4_3E = aovert_minus1_3E(:,5);

motif1_pure = aovert_minus1_pure(:,2);
motif1_3A = aovert_minus1_3A(:,2);
motif1_3C = aovert_minus1_3C(:,2);
motif1_3D = aovert_minus1_3D(:,2);
motif1_3E = aovert_minus1_3E(:,2);

motif3_pure = aovert_minus1_pure(:,4);
motif3_3A = aovert_minus1_3A(:,4);
motif3_3C = aovert_minus1_3C(:,4);
motif3_3D = aovert_minus1_3D(:,4);
motif3_3E = aovert_minus1_3E(:,4);

%% separate figures for constituent and third-order motif classes
dummy_val = repelem(-20,100);
dummy_val = dummy_val';
% Combine data for each group
group1 = {entropy_raster_4d_pure, motif4_pure,dummy_val,dummy_val};
group2 = {entropy_raster_4d_3A,  motif4_3A,dummy_val,dummy_val};
group3 = {entropy_raster_4d_3C, motif4_3C,dummy_val,dummy_val};
group4 = {entropy_raster_4d_3D, motif4_3D,dummy_val,dummy_val};
group5 = {entropy_raster_4d_3E, motif4_3E,dummy_val,dummy_val};

combined_data = [group1, group2, group3, group4, group5];

figure;
h = boxplotGroup(combined_data, ...
    'GroupLines', true, 'Colors', 'gm');
ylim([-1.5 18])
yline(0,'k')
title('Figure 5Bi')

group1 = {motif1_pure, motif3_pure, motifV_pure,dummy_val,dummy_val};
group2 = {motif1_3A, motif3_3A, motifV_3A, dummy_val,dummy_val};
group3 = {motif1_3C, motif3_3C, motifV_3C, dummy_val,dummy_val};
group4 = {motif1_3D, motif3_3D, motifV_3D, dummy_val,dummy_val};
group5 = {motif1_3E, motif3_3E, motifV_3E, dummy_val,dummy_val};

% Combine all groups
clear combined_data
combined_data = [group1, group2, group3, group4, group5];

figure;
h = boxplotGroup(combined_data, ...
    'GroupLines', true, 'Colors', 'bkrgg');
ylim([-1.5 18])
hold on;
yline(0,'k')
title('Figure 5Bii')

%%
figure; hold on
markersize = 50;

hold on; ylabel('')
hold on;
raster  = raster_pure;
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);

        end;
    end;
end;
axis off
box off
title('Figure 4Aii')

figure;
markersize = 50;
hold on; ylabel('')
hold on;
raster = raster_3A; 
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);

        end;
    end;
end;
axis off
box off
title('Figure 4Bii')

figure;
markersize = 50;
hold on; ylabel('')
hold on;
raster = raster_3C; 
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);

        end;
    end;
end;
axis off
box off
title('Figure 4Cii')

figure;
markersize = 50;
hold on; ylabel('')
hold on;
raster = raster_3D; 
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);

        end;
    end;
end;
axis off
box off
title('Figure 4Dii')

figure;
markersize = 50;
hold on; ylabel('')
hold on;
raster = raster_3E; 
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);

        end;
    end;
end;
axis off
box off
title('Figure 4Eii')