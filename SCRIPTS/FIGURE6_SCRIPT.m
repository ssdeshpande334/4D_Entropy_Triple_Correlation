% Script to generate Figure 6
% Triple correlation analysis with raster + box plots, regenerating raster each iteration
clearvars
close all
clc

% Parameters
num_neurons = 10;       % Number of neurons
duration = 500;         % Duration of spike train (samples)
sample_rate = 500;  % Samples per second
spike_rate_hz = 15; % Spikes per second
rate = spike_rate_hz / sample_rate; % Probability per bin
motif_prob = 0.20;      % Probability of embedding a 3-spike motif
noise_ratio = 0.3;      % Noise level for structured cases
neuron_window = 10;     % Number of neurons in window for triple correlation
time_windows = [50, 100, 200]; % Short, middle, long time lags
num_iterations = 10;   % Number of times to regenerate raster & compute triple correlation

% Compute triple correlation and collect statistics
cases = {'poisson', 'motif_poisson', 'noisy_motif', 'random_motif'};

for c = 1:length(cases)

    % Storage for runs of motif contributions
    motif_contributions = cell(length(time_windows), 1);

    for iter = 1:num_iterations
        disp(iter)
        switch cases{c}
            case 'poisson'
                temp_snippet = generate_poisson_raster(num_neurons, duration, round(num_neurons * duration * rate / num_neurons));
            case 'motif_poisson'
                temp_snippet = generate_motif_poisson_raster(num_neurons, duration, round(num_neurons * duration * rate / num_neurons), motif_prob);
            case 'noisy_motif'
                temp_snippet = generate_noisy_motif_raster(num_neurons, duration, round(num_neurons * duration * rate / num_neurons), motif_prob, noise_ratio);
            case 'random_motif'
                temp_snippet = generate_pure_motif_raster(num_neurons, duration, round(num_neurons * duration * rate / num_neurons), motif_prob);
        end
        
        % Loop over time windows
        for t_idx = 1:length(time_windows)
            time_window = time_windows(t_idx);
            flank_size = round(time_window / 2); % Flank size is half of the time window
           
            % Flank the raster 
            zero_flank = zeros(num_neurons, flank_size);
            flanked_snippet = [zero_flank, temp_snippet, zero_flank];

            % Compute triple correlation
            tic
            [~, actual_contribution, ~, ~] = ...
                triple_correlation_class_contributions_no_sp_wr(flanked_snippet, neuron_window, time_window);
            toc
            % Compute conditioned expectation
            actual = actual_contribution ./ numel(temp_snippet);
            conditioned_expectation = expectation_conditioned_on_constituent_parts_2D(actual, temp_snippet, neuron_window, time_window);
            aovert_minus1 = (actual ./ conditioned_expectation) - 1;

            % Extract motif class contributions (columns 1, 2, and 3) and store them
            motif_contributions{t_idx}(iter, :) = (aovert_minus1(:, 1:3));
        end
    end

    % Raster Plot (Only from last iteration)
    figure;
    hold on;
    
    % Extract raster from the last iteration
    raster = temp_snippet;
    
    % Find all spike locations
    [neuron_idx, time_idx] = find(raster);
    
    % Identify motif spikes (red) and noisy spikes (black)
    motif_times = []; motif_neurons = [];
    noise_times = []; noise_neurons = [];
    
    for i = 1:num_neurons
        for t = 1:duration - 2
            % Check for 3-spike motif
            if raster(i, t) == 1 && raster(i, t+1) == 1 && raster(i, t+2) == 1
                motif_times = [motif_times, t, t+1, t+2];
                motif_neurons = [motif_neurons, i, i, i];
            end
        end
    end
    
    % Identify noisy spikes (spikes that are not part of 3-spike motifs)
    for i = 1:length(time_idx)
        if ~ismember(time_idx(i), motif_times)
            noise_times = [noise_times, time_idx(i)];
            noise_neurons = [noise_neurons, neuron_idx(i)];
        end
    end
    
    % Plot noisy spikes in black
    scatter(noise_times, noise_neurons, 20, 'k', 'filled'); 
    
    % Plot motif spikes in red
    scatter(motif_times, motif_neurons, 20, 'r', 'filled');
    
    % Formatting
    title(['Raster: ' strrep(cases{c}, '_', ' ')]);
    xlabel('Time (ms)');
    ylabel('Neuron');
    set(gca, 'Color', 'w'); % White background
    xlim([0 duration]); % Show only first 100 ms
    ylim([0 num_neurons + 1]);
    hold off;

    % Box Plots for Time Windows
    figure; hold on; 
    for t_idx = 1:length(time_windows)
        subplot(1, 3, t_idx);
        boxplot(motif_contributions{t_idx}, 'Labels', {'Class 0', 'Class 1', 'Class 2'});
        ylabel('Contribution');
        title(['Lag = ' num2str(time_windows(t_idx)) ' ms']);
        yline(0,'k-')
    end
end

%%
% Compute SNR only for Case 2 (motif_poisson) and Case 3 (noisy_motif)
snr_values = struct();
snr_avg = struct(); % Store average SNR per case

% Only compute SNR for these cases
snr_cases = {'motif_poisson', 'noisy_motif'};

for c = 1:length(snr_cases)
    case_name = snr_cases{c};
    snr_list = zeros(num_iterations, 1); % Store SNR per iteration

    for iter = 1:num_iterations
        % Generate a new raster for this iteration
        switch case_name
            case 'motif_poisson'
                raster = generate_motif_poisson_raster(num_neurons, duration, round(num_neurons * duration * rate / num_neurons), motif_prob);
            case 'noisy_motif'
                raster = generate_noisy_motif_raster(num_neurons, duration, round(num_neurons * duration * rate / num_neurons), motif_prob, noise_ratio);
        end

        % Count total spikes
        total_spikes = sum(raster(:));

        % Count motif spikes (spikes part of a 3-spike motif)
        motif_spikes = 0;
        for i = 1:num_neurons
            for t = 1:duration - 2
                if raster(i, t) == 1 && raster(i, t+1) == 1 && raster(i, t+2) == 1
                    motif_spikes = motif_spikes + 3; % Count all 3 spikes in the motif
                end
            end
        end

        % Compute Poisson-only spikes (assumed to be total spikes minus structured motifs)
        poisson_spikes = total_spikes - motif_spikes;

        % Compute SNR (motif spikes / poisson spikes)
        if poisson_spikes > 0
            snr_list(iter) = 10 * log10(motif_spikes / poisson_spikes);        
        else
            snr_list(iter) = NaN; % Avoid division by zero
        end
    end

    % Store the SNR values for this case
    snr_values.(case_name) = snr_list;

    % Compute and store the average SNR for this case
    snr_avg.(case_name) = nanmean(snr_list); % Use nanmean to ignore NaN values
end

% Display SNR values and average SNR
disp('SNR Values (Motif Spikes / Poisson Spikes) per iteration:');
disp(snr_values);

disp('Average SNR for motif_poisson and noisy_motif:');
disp(snr_avg);


%% Function Definitions

% Poisson spike raster with exact spike count
function raster = generate_poisson_raster(num_neurons, duration, spike_target)
    raster = zeros(num_neurons, duration);
    for i = 1:num_neurons
        spike_indices = randperm(duration, spike_target);
        raster(i, spike_indices) = 1;
    end
end

% Poisson +  motif class II while keeping spike count fixed
function raster = generate_motif_poisson_raster(num_neurons, duration, spike_target, motif_prob)
    raster = zeros(num_neurons, duration);
    for i = 1:num_neurons
        spikes = zeros(1, duration);
        count = 0;
        while count < spike_target
            if rand < motif_prob && count + 3 <= spike_target
                t = randi([1, duration - 2]);
                spikes([t, t+1, t+2]) = 1; % Adjacent  motif class II
                count = count + 3;
            else
                t = randi([1, duration]);
                if spikes(t) == 0
                    spikes(t) = 1;
                    count = count + 1;
                end
            end
        end
        raster(i, :) = spikes;
    end
end

% Structured motifs embedded in noise, ensuring spike count is fixed
function raster = generate_noisy_motif_raster(num_neurons, duration, spike_target, motif_prob, noise_ratio)
    raster = generate_motif_poisson_raster(num_neurons, duration, spike_target, motif_prob);
    
    % Controlled noise addition: flip a fraction of spikes randomly
    total_spikes = sum(raster(:));
    noise_spikes = round(noise_ratio * total_spikes); % Define noise level
    noise_indices = randperm(numel(raster), noise_spikes);
    
    % Flip only selected spike positions
    raster(noise_indices) = ~raster(noise_indices);
end

% Pure motif class II raster with fixed spike count
function raster = generate_pure_motif_raster(num_neurons, duration, spike_target, motif_prob)
    raster = zeros(num_neurons, duration);
    for i = 1:num_neurons
        spikes = zeros(1, duration);
        count = 0;
        while count < spike_target
            t = randi([1, duration - 2]);
            if all(spikes([t, t+1, t+2]) == 0) % Ensure motifs do not overlap
                spikes([t, t+1, t+2]) = 1; % Adjacent  motif class II
                count = count + 3;
            end
        end
        raster(i, :) = spikes;
    end
end