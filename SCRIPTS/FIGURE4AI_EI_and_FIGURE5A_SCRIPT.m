clearvars
close all
clc

% SCRIPT FOR FIGURE 4AI-EI & FIGURE 5
% Motif class II

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

%%
for iiiii = 1:100

    clearvars -except iiiii entropy_raster_4d_pure aovert_minus1_pure...
        entropy_raster_4d_9A aovert_minus1_9A...
        entropy_raster_4d_9C aovert_minus1_9C...
        entropy_raster_4d_9D aovert_minus1_9D...
        raster_pure raster_9A raster_9B raster_9C raster_9D raster_9E...
        neuron_window time_window

    disp(iiiii)
% Set parameters
num_rows = 50;
num_cols = 50;
num_sections = 3; % Number of sections
num_consecutive_spikes_per_section = 3; % 3-spike in consecutive columns in each section
num_poisson_spikes = 48 - num_sections * num_consecutive_spikes_per_section; % Remaining spikes to be placed randomly

% Define boundaries
row_min = neuron_window / 2 + 1;
row_max = num_rows - neuron_window / 2;
col_min = time_window / 2 + 1;
col_max = num_cols - time_window / 2;

% Create spike raster
spike_raster = zeros(num_rows, num_cols);

% Divide the raster into 3 sections
sections = [
    1, 1, floor(num_rows/3), floor(num_cols/2);
    floor(num_rows/3)+1, 1, 2*floor(num_rows/3), floor(num_cols/2);
    2*floor(num_rows/3)+1, 1, num_rows, floor(num_cols/2);
];

% Place one 3-spike consecutive columns set in each section
for s = 1:num_sections
    placed = false;
    while ~placed
        % Get section boundaries
        row_start = sections(s, 1);
        col_start = sections(s, 2);
        row_end = sections(s, 3);
        col_end = sections(s, 4);

        % Randomly select a position for the first spike of the consecutive set
        row = randi([max(row_min, row_start), min(row_max, row_end)]);
        col = randi([max(col_min, col_start), min(col_max-2, col_end)]); % Ensure there's space for the consecutive spikes

        % Check if the consecutive spikes fit within the boundaries
        if col <= col_end-2
            % Place the consecutive spikes
            if spike_raster(row, col) == 0 && spike_raster(row, col+1) == 0 && spike_raster(row, col+2) == 0
                spike_raster(row, col) = 1;     % First spike
                spike_raster(row, col+1) = 1;   % Second spike
                spike_raster(row, col+2) = 1;   % Third spike
                placed = true;
            end
        end
    end
end

% Add Poisson-distributed spikes
% Ensure we only place Poisson spikes in empty spots
remaining_spikes = num_poisson_spikes;
while remaining_spikes > 0
    % Randomly select a position for a Poisson spike
    row = randi([row_min, row_max]);
    col = randi([col_min, col_max]);

    % Check if the position is empty
    if spike_raster(row, col) == 0
        spike_raster(row, col) = 1;
        remaining_spikes = remaining_spikes - 1;
    end
end

% Use the created raster for further analysis
raster = spike_raster;
epoch_length = size(raster, 2) - (time_window + 1);
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

post_end = 0;
slice_start = post_end + 1;

pre_snippet_raster = raster(:, slice_start:slice_start + max_time_lag - 1);
pre_end = slice_start + max_time_lag - 1;

snippet_raster = raster(:, pre_end + 1: pre_end + 1 + epoch_length);
snip_end = pre_end + 1 + epoch_length;

post_snippet_raster = raster(:, snip_end + 1 : snip_end + 1 + max_time_lag - 1);
post_end = snip_end + 1 + max_time_lag - 1;

temp_snippet = cat(2, pre_snippet_raster, snippet_raster, post_snippet_raster);

size(temp_snippet)
[N_neurons, N_times] = size(temp_snippet);

% Compute tricorr
[c3_4D_distribution, actual_contribution, class_count, contribution] = ...
    triple_correlation_class_contributions_no_sp_wr(temp_snippet, neuron_window, time_window);

actual = actual_contribution./(numel(snippet_raster));
[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, snippet_raster, neuron_window, time_window);

aovert_minus1 = (actual ./ conditioned_expectation) - 1;

c3_n1t1_distribution = c3_4D_distribution(:,:,neuron_window/2 + 1 ,time_window / 2 + 1);
T_raster_n1t1 = sum(sum(c3_n1t1_distribution));
PDF_raster_n1t1 = c3_n1t1_distribution ./ T_raster_n1t1;  

c3_n2t2_distribution = c3_4D_distribution(neuron_window/2 + 1 ,time_window / 2 + 1,:,:);
T_raster_n2t2 = sum(sum(c3_n2t2_distribution));
PDF_raster_n2t2 = c3_n2t2_distribution ./ T_raster_n2t2;

T_raster_4d = sum(sum(sum(sum(c3_4D_distribution))));
PDF_raster_4d = c3_4D_distribution ./ T_raster_4d;

% Compute the 4D Entropy
bins_4d = [neuron_window+1, time_window+1, neuron_window+1, time_window+1];
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

entropy_raster_4d_9A(iiiii, 1) = entropy_raster_4d;
aovert_minus1_9A(iiiii, :) = aovert_minus1;
end

raster_9A = raster;

%%
for iiiii = 1:100

    clearvars -except iiiii entropy_raster_4d_pure aovert_minus1_pure...
        entropy_raster_4d_9A aovert_minus1_9A...
        entropy_raster_4d_9C aovert_minus1_9C...
        entropy_raster_4d_9D aovert_minus1_9D...
        raster_pure raster_9A raster_9B raster_9C raster_9D raster_9E...
        neuron_window time_window

    disp(iiiii)
% Set parameters
num_rows = 50;
num_cols = 50;
num_sections = 9; % Number of sections
num_consecutive_spikes_per_section = 3; % 3-spike in consecutive columns in each section
num_poisson_spikes = 48 - num_sections * num_consecutive_spikes_per_section; % Remaining spikes to be placed randomly

% Define boundaries
row_min = neuron_window / 2 + 1;
row_max = num_rows - neuron_window / 2;
col_min = time_window / 2 + 1;
col_max = num_cols - time_window / 2;

% Create spike raster
spike_raster = zeros(num_rows, num_cols);

% Divide the raster into 9 sections (3x3 grid)
sections = [
    1, 1, floor(num_rows/3), floor(num_cols/3);
    1, floor(num_cols/3)+1, floor(num_rows/3), 2*floor(num_cols/3);
    1, 2*floor(num_cols/3)+1, floor(num_rows/3), num_cols;
    floor(num_rows/3)+1, 1, 2*floor(num_rows/3), floor(num_cols/3);
    floor(num_rows/3)+1, floor(num_cols/3)+1, 2*floor(num_rows/3), 2*floor(num_cols/3);
    floor(num_rows/3)+1, 2*floor(num_cols/3)+1, 2*floor(num_rows/3), num_cols;
    2*floor(num_rows/3)+1, 1, num_rows, floor(num_cols/3);
    2*floor(num_rows/3)+1, floor(num_cols/3)+1, num_rows, 2*floor(num_cols/3);
    2*floor(num_rows/3)+1, 2*floor(num_cols/3)+1, num_rows, num_cols;
];

% Place one 3-spike consecutive columns set in each section
for s = 1:num_sections
    placed = false;
    while ~placed
        % Get section boundaries
        row_start = sections(s, 1);
        col_start = sections(s, 2);
        row_end = sections(s, 3);
        col_end = sections(s, 4);

        % Randomly select a position for the first spike of the consecutive set
        row = randi([max(row_min, row_start), min(row_max, row_end)]);
        col = randi([max(col_min, col_start), min(col_max-2, col_end)]); % Ensure there's space for the consecutive spikes

        % Check if the consecutive spikes fit within the boundaries
        if col <= col_end-2
            % Place the consecutive spikes
            if spike_raster(row, col) == 0 && spike_raster(row, col+1) == 0 && spike_raster(row, col+2) == 0
                spike_raster(row, col) = 1;     % First spike
                spike_raster(row, col+1) = 1;   % Second spike
                spike_raster(row, col+2) = 1;   % Third spike
                placed = true;
            end
        end
    end
end

% Add Poisson-distributed spikes
% Ensure we only place Poisson spikes in empty spots
remaining_spikes = num_poisson_spikes;
while remaining_spikes > 0
    % Randomly select a position for a Poisson spike
    row = randi([row_min, row_max]);
    col = randi([col_min, col_max]);

    % Check if the position is empty
    if spike_raster(row, col) == 0
        spike_raster(row, col) = 1;
        remaining_spikes = remaining_spikes - 1;
    end
end

% Use the created raster for further analysis
raster = spike_raster;
epoch_length = size(raster, 2) - (time_window + 1);
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

post_end = 0;
slice_start = post_end + 1;

pre_snippet_raster = raster(:, slice_start:slice_start + max_time_lag - 1);
pre_end = slice_start + max_time_lag - 1;

snippet_raster = raster(:, pre_end + 1: pre_end + 1 + epoch_length);
snip_end = pre_end + 1 + epoch_length;

post_snippet_raster = raster(:, snip_end + 1 : snip_end + 1 + max_time_lag - 1);
post_end = snip_end + 1 + max_time_lag - 1;

temp_snippet = cat(2, pre_snippet_raster, snippet_raster, post_snippet_raster);

size(temp_snippet)
[N_neurons, N_times] = size(temp_snippet);

% Compute tricorr
[c3_4D_distribution, actual_contribution, class_count, contribution] = ...
    triple_correlation_class_contributions_no_sp_wr(temp_snippet, neuron_window, time_window);

actual = actual_contribution./(numel(snippet_raster));
[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, snippet_raster, neuron_window, time_window);

aovert_minus1 = (actual ./ conditioned_expectation) - 1;

c3_n1t1_distribution = c3_4D_distribution(:,:,neuron_window/2 + 1 ,time_window / 2 + 1);
T_raster_n1t1 = sum(sum(c3_n1t1_distribution));
PDF_raster_n1t1 = c3_n1t1_distribution ./ T_raster_n1t1;  

c3_n2t2_distribution = c3_4D_distribution(neuron_window/2 + 1 ,time_window / 2 + 1,:,:);
T_raster_n2t2 = sum(sum(c3_n2t2_distribution));
PDF_raster_n2t2 = c3_n2t2_distribution ./ T_raster_n2t2;

T_raster_4d = sum(sum(sum(sum(c3_4D_distribution))));
PDF_raster_4d = c3_4D_distribution ./ T_raster_4d;

% Compute the 4D Entropy
bins_4d = [neuron_window+1, time_window+1, neuron_window+1, time_window+1];
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

entropy_raster_4d_9C(iiiii, 1) = entropy_raster_4d;
aovert_minus1_9C(iiiii, :) = aovert_minus1;
end

raster_9C = raster;

%%
for iiiii = 1:100

    clearvars -except iiiii entropy_raster_4d_pure aovert_minus1_pure...
        entropy_raster_4d_9A aovert_minus1_9A...
        entropy_raster_4d_9C aovert_minus1_9C...
        entropy_raster_4d_9D aovert_minus1_9D...
        raster_pure raster_9A raster_9B raster_9C raster_9D raster_9E...
        neuron_window time_window

    disp(iiiii)
% Set parameters
num_rows = 50;
num_cols = 50;
num_sections = 12; % Number of sections
num_consecutive_spikes_per_section = 3; % 3-spike in consecutive columns in each section
num_poisson_spikes = 48 - num_sections * num_consecutive_spikes_per_section; % Remaining spikes to be placed randomly

% Define boundaries
row_min = neuron_window / 2 + 1;
row_max = num_rows - neuron_window / 2;
col_min = time_window / 2 + 1;
col_max = num_cols - time_window / 2;

% Create spike raster
spike_raster = zeros(num_rows, num_cols);

% Divide the raster into 12 sections
sections = [
    1, 1, floor(num_rows/4), floor(num_cols/3);
    1, floor(num_cols/3)+1, floor(num_rows/4), 2*floor(num_cols/3);
    1, 2*floor(num_cols/3)+1, floor(num_rows/4), num_cols;
    floor(num_rows/4)+1, 1, 2*floor(num_rows/4), floor(num_cols/3);
    floor(num_rows/4)+1, floor(num_cols/3)+1, 2*floor(num_rows/4), 2*floor(num_cols/3);
    floor(num_rows/4)+1, 2*floor(num_cols/3)+1, 2*floor(num_rows/4), num_cols;
    2*floor(num_rows/4)+1, 1, 3*floor(num_rows/4), floor(num_cols/3);
    2*floor(num_rows/4)+1, floor(num_cols/3)+1, 3*floor(num_rows/4), 2*floor(num_cols/3);
    2*floor(num_rows/4)+1, 2*floor(num_cols/3)+1, 3*floor(num_rows/4), num_cols;
    3*floor(num_rows/4)+1, 1, num_rows, floor(num_cols/3);
    3*floor(num_rows/4)+1, floor(num_cols/3)+1, num_rows, 2*floor(num_cols/3);
    3*floor(num_rows/4)+1, 2*floor(num_cols/3)+1, num_rows, num_cols;
];

% Place one 3-spike consecutive columns set in each section
for s = 1:num_sections
    placed = false;
    while ~placed
        % Get section boundaries
        row_start = sections(s, 1);
        col_start = sections(s, 2);
        row_end = sections(s, 3);
        col_end = sections(s, 4);

        % Randomly select a position for the first spike of the consecutive set
        row = randi([max(row_min, row_start), min(row_max, row_end)]);
        col = randi([max(col_min, col_start), min(col_max-2, col_end)]); % Ensure there's space for the consecutive spikes

        % Check if the consecutive spikes fit within the boundaries
        if col <= col_end-2
            % Place the consecutive spikes
            if spike_raster(row, col) == 0 && spike_raster(row, col+1) == 0 && spike_raster(row, col+2) == 0
                spike_raster(row, col) = 1;     % First spike
                spike_raster(row, col+1) = 1;   % Second spike
                spike_raster(row, col+2) = 1;   % Third spike
                placed = true;
            end
        end
    end
end

% Add Poisson-distributed spikes
% Ensure we only place Poisson spikes in empty spots
remaining_spikes = num_poisson_spikes;
while remaining_spikes > 0
    % Randomly select a position for a Poisson spike
    row = randi([row_min, row_max]);
    col = randi([col_min, col_max]);

    % Check if the position is empty
    if spike_raster(row, col) == 0
        spike_raster(row, col) = 1;
        remaining_spikes = remaining_spikes - 1;
    end
end

% Use the created raster for further analysis
raster = spike_raster;
epoch_length = size(raster, 2) - (time_window + 1);
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

post_end = 0;
slice_start = post_end + 1;

pre_snippet_raster = raster(:, slice_start:slice_start + max_time_lag - 1);
pre_end = slice_start + max_time_lag - 1;

snippet_raster = raster(:, pre_end + 1: pre_end + 1 + epoch_length);
snip_end = pre_end + 1 + epoch_length;

post_snippet_raster = raster(:, snip_end + 1 : snip_end + 1 + max_time_lag - 1);
post_end = snip_end + 1 + max_time_lag - 1;

temp_snippet = cat(2, pre_snippet_raster, snippet_raster, post_snippet_raster);

size(temp_snippet)
[N_neurons, N_times] = size(temp_snippet);

% Compute tricorr
[c3_4D_distribution, actual_contribution, class_count, contribution] = ...
    triple_correlation_class_contributions_no_sp_wr(temp_snippet, neuron_window, time_window);

actual = actual_contribution./(numel(snippet_raster));
[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, snippet_raster, neuron_window, time_window);

aovert_minus1 = (actual ./ conditioned_expectation) - 1;

c3_n1t1_distribution = c3_4D_distribution(:,:,neuron_window/2 + 1 ,time_window / 2 + 1);
T_raster_n1t1 = sum(sum(c3_n1t1_distribution));
PDF_raster_n1t1 = c3_n1t1_distribution ./ T_raster_n1t1;  

c3_n2t2_distribution = c3_4D_distribution(neuron_window/2 + 1 ,time_window / 2 + 1,:,:);
T_raster_n2t2 = sum(sum(c3_n2t2_distribution));
PDF_raster_n2t2 = c3_n2t2_distribution ./ T_raster_n2t2;

T_raster_4d = sum(sum(sum(sum(c3_4D_distribution))));
PDF_raster_4d = c3_4D_distribution ./ T_raster_4d;

% Compute the 4D Entropy
bins_4d = [neuron_window+1, time_window+1, neuron_window+1, time_window+1];
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

entropy_raster_4d_9D(iiiii, 1) = entropy_raster_4d;
aovert_minus1_9D(iiiii, :) = aovert_minus1;
end

raster_9D = raster;

%%
for iiiii = 1:100

    clearvars -except iiiii entropy_raster_4d_pure aovert_minus1_pure...
        entropy_raster_4d_9A aovert_minus1_9A...
        entropy_raster_4d_9C aovert_minus1_9C...
        entropy_raster_4d_9D aovert_minus1_9D...
        entropy_raster_4d_9E aovert_minus1_9E...
        raster_pure raster_9A raster_9B raster_9C raster_9D raster_9E...
        neuron_window time_window

% Set parameters
num_rows = 50;
num_cols = 50;
num_sections = 16; % Number of sections
num_consecutive_spikes_per_section = 3; % 3-spike in consecutive columns in each section
num_poisson_spikes = 48 - num_sections * num_consecutive_spikes_per_section; % Remaining spikes to be placed randomly

% Set parameters
num_rows = 50;
num_cols = 50;
num_sections = 16; % Number of sections
num_consecutive_spikes_per_section = 3; % 3-spike in consecutive columns in each section
num_poisson_spikes = 48 - num_sections * num_consecutive_spikes_per_section; % Remaining spikes to be placed randomly

% Define boundaries
row_min = neuron_window / 2 + 1;
row_max = num_rows - neuron_window / 2;
col_min = time_window / 2 + 1;
col_max = num_cols - time_window / 2;

% Create spike raster
spike_raster = zeros(num_rows, num_cols);

% Divide the raster into 16 sections
sections = [
    1, 1, floor(num_rows/4), floor(num_cols/4);
    1, floor(num_cols/4)+1, floor(num_rows/4), 2*floor(num_cols/4);
    1, 2*floor(num_cols/4)+1, floor(num_rows/4), 3*floor(num_cols/4);
    1, 3*floor(num_cols/4)+1, floor(num_rows/4), num_cols;
    floor(num_rows/4)+1, 1, 2*floor(num_rows/4), floor(num_cols/4);
    floor(num_rows/4)+1, floor(num_cols/4)+1, 2*floor(num_rows/4), 2*floor(num_cols/4);
    floor(num_rows/4)+1, 2*floor(num_cols/4)+1, 2*floor(num_rows/4), 3*floor(num_cols/4);
    floor(num_rows/4)+1, 3*floor(num_cols/4)+1, 2*floor(num_rows/4), num_cols;
    2*floor(num_rows/4)+1, 1, 3*floor(num_rows/4), floor(num_cols/4);
    2*floor(num_rows/4)+1, floor(num_cols/4)+1, 3*floor(num_rows/4), 2*floor(num_cols/4);
    2*floor(num_rows/4)+1, 2*floor(num_cols/4)+1, 3*floor(num_rows/4), 3*floor(num_cols/4);
    2*floor(num_rows/4)+1, 3*floor(num_cols/4)+1, 3*floor(num_rows/4), num_cols;
    3*floor(num_rows/4)+1, 1, num_rows, floor(num_cols/4);
    3*floor(num_rows/4)+1, floor(num_cols/4)+1, num_rows, 2*floor(num_cols/4);
    3*floor(num_rows/4)+1, 2*floor(num_cols/4)+1, num_rows, 3*floor(num_cols/4);
    3*floor(num_rows/4)+1, 3*floor(num_cols/4)+1, num_rows, num_cols;
];

% Define the row ranges with at least 5 rows of separation
row_ranges = [
    row_min, floor(num_rows/4) - 5;
    floor(num_rows/4) + 5, 2*floor(num_rows/4) - 5;
    2*floor(num_rows/4) + 5, 3*floor(num_rows/4) - 5;
    3*floor(num_rows/4) + 5, row_max;
];

% To track rows with consecutive spikes
rows_with_consecutive_spikes = [];

% Place 3-spike consecutive columns set in each section, grouped by sets of 4 sections sharing the same row
for g = 1:4
    placed = false;
    while ~placed
        % Randomly select a row within the allowed range for the current group
        row = randi([row_ranges(g, 1), row_ranges(g, 2)]);
        col_positions = zeros(4, 1);
        cols_valid = true;

        for s = 1:4
            section_idx = (g - 1) * 4 + s;
            % Get section boundaries
            row_start = sections(section_idx, 1);
            col_start = sections(section_idx, 2);
            row_end = sections(section_idx, 3);
            col_end = sections(section_idx, 4);

            % Randomly select a position for the first spike of the consecutive set
            col = randi([max(col_min, col_start), min(col_max - 2, col_end)]); % Ensure there's space for the consecutive spikes
            col_positions(s) = col;

            % Check if the consecutive spikes fit within the boundaries
            if col > col_end - 2 || spike_raster(row, col) == 1 || spike_raster(row, col + 1) == 1 || spike_raster(row, col + 2) == 1
                cols_valid = false;
                break;
            end
        end

        if cols_valid
            for s = 1:4
                col = col_positions(s);
                spike_raster(row, col) = 1; % First spike
                spike_raster(row, col + 1) = 1; % Second spike
                spike_raster(row, col + 2) = 1; % Third spike
            end
            rows_with_consecutive_spikes = [rows_with_consecutive_spikes; row]; % Track this row
            placed = true;
        end
    end
end

% Add Poisson-distributed spikes
% Ensure we only place Poisson spikes in empty spots and not in the same rows as the consecutive spikes
remaining_spikes = num_poisson_spikes;
while remaining_spikes > 0
    % Randomly select a position for a Poisson spike
    row = randi([row_min, row_max]);
    col = randi([col_min, col_max]);

    % Check if the position is empty and not in the rows with consecutive spikes
    if spike_raster(row, col) == 0 && ~ismember(row, rows_with_consecutive_spikes)
        spike_raster(row, col) = 1;
        remaining_spikes = remaining_spikes - 1;
    end
end

% Use the created raster for further analysis
raster = spike_raster;
epoch_length = size(raster, 2) - (time_window + 1);
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

post_end = 0;
slice_start = post_end + 1;

pre_snippet_raster = raster(:, slice_start:slice_start + max_time_lag - 1);
pre_end = slice_start + max_time_lag - 1;

snippet_raster = raster(:, pre_end + 1: pre_end + 1 + epoch_length);
snip_end = pre_end + 1 + epoch_length;

post_snippet_raster = raster(:, snip_end + 1 : snip_end + 1 + max_time_lag - 1);
post_end = snip_end + 1 + max_time_lag - 1;

temp_snippet = cat(2, pre_snippet_raster, snippet_raster, post_snippet_raster);

size(temp_snippet)
[N_neurons, N_times] = size(temp_snippet);

% Compute tricorr
[c3_4D_distribution, actual_contribution, class_count, contribution] = ...
    triple_correlation_class_contributions_no_sp_wr(temp_snippet, neuron_window, time_window);

actual = actual_contribution./(numel(snippet_raster));
[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, snippet_raster, neuron_window, time_window);

aovert_minus1 = (actual ./ conditioned_expectation) - 1;

c3_n1t1_distribution = c3_4D_distribution(:,:,neuron_window/2 + 1 ,time_window / 2 + 1);
T_raster_n1t1 = sum(sum(c3_n1t1_distribution));
PDF_raster_n1t1 = c3_n1t1_distribution ./ T_raster_n1t1;  

c3_n2t2_distribution = c3_4D_distribution(neuron_window/2 + 1 ,time_window / 2 + 1,:,:);
T_raster_n2t2 = sum(sum(c3_n2t2_distribution));
PDF_raster_n2t2 = c3_n2t2_distribution ./ T_raster_n2t2;

T_raster_4d = sum(sum(sum(sum(c3_4D_distribution))));
PDF_raster_4d = c3_4D_distribution ./ T_raster_4d;

% Compute the 4D Entropy
bins_4d = [neuron_window+1, time_window+1, neuron_window+1, time_window+1];
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

entropy_raster_4d_9E(iiiii, 1) = entropy_raster_4d;
aovert_minus1_9E(iiiii, :) = aovert_minus1;
end

raster_9E = raster;

%%
close all

motifV_pure = aovert_minus1_pure(:,6);
motifV_9A = aovert_minus1_9A(:,6);
motifV_9C = aovert_minus1_9C(:,6);
motifV_9D = aovert_minus1_9D(:,6);
motifV_9E = aovert_minus1_9E(:,6);

%
motif2_pure = aovert_minus1_pure(:,3);
motif2_9A = aovert_minus1_9A(:,3);
motif2_9C = aovert_minus1_9C(:,3);
motif2_9D = aovert_minus1_9D(:,3);
motif2_9E = aovert_minus1_9E(:,3);

motif1_pure = aovert_minus1_pure(:,2);
motif1_9A = aovert_minus1_9A(:,2);
motif1_9C = aovert_minus1_9C(:,2);
motif1_9D = aovert_minus1_9D(:,2);
motif1_9E = aovert_minus1_9E(:,2);

motif3_pure = aovert_minus1_pure(:,4);
motif3_9A = aovert_minus1_9A(:,4);
motif3_9C = aovert_minus1_9C(:,4);
motif3_9D = aovert_minus1_9D(:,4);
motif3_9E = aovert_minus1_9E(:,4);

%% separate figures for constituent and third-order motif classes
dummy_val = repelem(-20,100);
dummy_val = dummy_val';

% Combine data for each group
group1 = {entropy_raster_4d_pure, motif2_pure,dummy_val,dummy_val};
group2 = {entropy_raster_4d_9A,  motif2_9A,dummy_val,dummy_val};
group3 = {entropy_raster_4d_9C, motif2_9C,dummy_val,dummy_val};
group4 = {entropy_raster_4d_9D, motif2_9D,dummy_val,dummy_val};
group5 = {entropy_raster_4d_9E, motif2_9E,dummy_val,dummy_val};

combined_data = [group1, group2, group3, group4, group5];
% Plot grouped boxplots with custom colors
figure;
h = boxplotGroup(combined_data, ...
    'GroupLines', true, 'Colors', 'gm');
ylim([-1.5 18])
yline(0,'k')
title('Figure 5Ai')

group1 = {motif1_pure, motif3_pure, motifV_pure,dummy_val,dummy_val};
group2 = {motif1_9A, motif3_9A, motifV_9A, dummy_val,dummy_val};
group3 = {motif1_9C, motif3_9C, motifV_9C, dummy_val,dummy_val};
group4 = {motif1_9D, motif3_9D, motifV_9D, dummy_val,dummy_val};
group5 = {motif1_9E, motif3_9E, motifV_9E, dummy_val,dummy_val};

% Combine all groups
clear combined_data
combined_data = [group1, group2, group3, group4, group5];

% Labels for primary and secondary groups
% Plot grouped boxplots with custom colors
figure;
h = boxplotGroup(combined_data, ...
    'GroupLines', true, 'Colors', 'rkbgg');
ylim([-1.5 18])
hold on;
yline(0,'k')
title('Figure 5Aii')

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
title('Figure 4Ai')
axis off
box off

figure;
markersize = 50;
hold on; ylabel('')
hold on;
raster = raster_9A; 
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);

        end;
    end;
end;
title('Figure 4Bi')
axis off
box off

figure;
markersize = 50;
hold on; ylabel('')
hold on;
raster = raster_9C; 
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);

        end;
    end;
end;
title('Figure 4Ci')
axis off
box off

figure;
markersize = 50;
hold on; ylabel('')
hold on;
raster = raster_9D; 
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);

        end;
    end;
end;
title('Figure 4Di')
axis off
box off

figure;
markersize = 50;
hold on; ylabel('')
hold on;
raster = raster_9E; 
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);

        end;
    end;
end;
title('Figure 4Ei')
axis off
box off
