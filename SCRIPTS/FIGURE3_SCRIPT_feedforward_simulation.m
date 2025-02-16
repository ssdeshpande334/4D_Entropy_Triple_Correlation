%% Script for Figure 3 of 4D Entropy Paper

clearvars
close all
clc

% STEP 1: GENERATE A SPIKE RASTER

set(0,'defaultAxesFontSize',20)

C=50;  % columns
R=50;  % rows

p = zeros(R,C);

str = 'Raster: Feedforward (Motif Class XIII)';
for n=5:10:35
    for t=5:10:35
         p(n,t)=1;
         p(n+2,t+1) = 1;
         p(n+1,t+3) = 1;
    end
end

raster = p;
neuron_window = 10; %size(raster,1)*2 -2;
time_window = 8;
epoch_length =size(raster,2) - (time_window + 1);
max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

figure(1); hold on
markersize = 50;
hold on; ylabel('')
hold on;
for k=1:size(raster,1);
    for tt=1:size(raster,2); 
        if raster(k,tt)~=0; 
            plot(tt,raster(k,tt)-2*k,'k.','MarkerSize',markersize);

        end;
    end;
end;
xlabel('Time')
ylabel('Neuron #')
title('Figure 3A Raster: Feedforward (Motif Class XIII)','fontsize',40)
axis off
box off

% figure;
% heatmap(p,'CellLabelColor','none');
% colormap(flipud(hot))
% colorbar off
% xlabel('Time bin #')
% ylabel('Neuron #')


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
temp_snippet = raster;
[N_neurons, N_times] = size(temp_snippet);

%compute tricorr
[c3_4D_distribution, actual_contribution,class_count,contribution]= ...
    triple_correlation_class_contributions_no_sp_wr(temp_snippet, neuron_window, time_window);

actual = actual_contribution./(numel(temp_snippet));
[conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual, temp_snippet,neuron_window, time_window);

aovert_minus1 = (actual ./conditioned_expectation ) -1;

c3_n1t1_distribution = c3_4D_distribution(:,:,neuron_window/2 + 1 ,time_window / 2 + 1);
T_raster_n1t1=sum(sum(c3_n1t1_distribution));
PDF_raster_n1t1=c3_n1t1_distribution./T_raster_n1t1;  

c3_n2t2_distribution = c3_4D_distribution(neuron_window/2 + 1 ,time_window / 2 + 1,:,:);
T_raster_n2t2=sum(sum(c3_n2t2_distribution));
PDF_raster_n2t2=c3_n2t2_distribution./T_raster_n2t2;

T_raster_4d=sum(sum(sum(sum(c3_4D_distribution))));
PDF_raster_4d=c3_4D_distribution./T_raster_4d;

figure;
h = heatmap(PDF_raster_n1t1); %,'CellLabelColor','none');
h.CellLabelFormat = '%.2f';
h.FontSize = 40;
%h.GridVisible = 'off';
colormap("parula")
xlabel('temporal lag, t_1')
ylabel('spatial lag, n_1')
title({'Figure 3C'})
for i = 1:numel(h.XDisplayLabels)
    h.XDisplayLabels{i} = str2num(h.XDisplayLabels{i}) - ((time_window/2) + 1);
end
for i = 1:numel(h.YDisplayLabels)
    h.YDisplayLabels{i} = str2num(h.YDisplayLabels{i}) - ((neuron_window/2) + 1);
end
clim([0 0.35])

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


%% now compute for 100 iterations of a surrogate raster
n_iterations = 100;
for i = 1:n_iterations
    disp(i)
    
    [pre_noise_raster{i,1}] = generate_noise_raster(raster(:,1:max_time_lag));
    [curr_noise_raster{i,1}] = generate_noise_raster(raster(:,max_time_lag+1:N_times-max_time_lag));
    [post_noise_raster{i,1}] = generate_noise_raster(raster(:,N_times-max_time_lag+1:N_times));
    [test_raster{i,1}] = generate_noise_raster(raster);
    noise_raster{i,1} = cat(2,pre_noise_raster{i,1},curr_noise_raster{i,1},post_noise_raster{i,1});
    tic

    %compute tricorr - no spatial wrapping
    [c3_noise_4D_distribution, noise_contribution,class_count,contribution]= ...
        triple_correlation_class_contributions_no_sp_wr(noise_raster{i,1}, neuron_window, time_window);
    toc

    c3_n1t1_noise_distribution(:,:,i) = c3_noise_4D_distribution(:,:,neuron_window/2 + 1 ,time_window / 2 + 1);
    T_noise_n1t1=sum(sum(c3_n1t1_noise_distribution));
    PDF_noise_n1t1(:,:,i)=c3_n1t1_noise_distribution(:,:,i)./T_raster_n1t1;  
    
    all_noise_iterations_c3_n1t1_noise_distribution{i,1} = c3_n1t1_noise_distribution(:,:,i);
    all_noise_iterations_PDF_n1t1_noise_distribution{i,1} = PDF_noise_n1t1(:,:,i);

    all_noise_iterations_c3_noise_4D_distribution{i,1} = c3_noise_4D_distribution;

    T_noise_4d=sum(sum(sum(sum(c3_noise_4D_distribution))));
    PDF_noise_4d(:,:,:,:,i)=c3_noise_4D_distribution./T_noise_4d;

    all_noise_iterations_PDF_noise_4D_distribution{i,1} = PDF_noise_4d(:,:,:,:,i);
    PDF_noise_4D_distribution = PDF_noise_4d(:,:,:,:,i);
    PDF_noise_n1t1_distribution = PDF_noise_n1t1(:,:,i);

    %Compute the 4D Entropy
    bins_4d = [neuron_window+1 time_window+1 neuron_window+1 time_window+1];
    entropy_noise_4d(i,1) = 0;
    
    for ii=1:bins_4d(1)
        for j=1:bins_4d(2)
            for k = 1:bins_4d(3)
                for m = 1:bins_4d(4)
                    temp_4d = PDF_noise_4d(:,:,:,:,i);
                    if temp_4d(ii,j,k,m)~=0        
                        entropy_noise_4d(i,1)=entropy_noise_4d(i,1)-temp_4d(ii,j,k,m)*log2(temp_4d(ii,j,k,m));       
                    end
                end
            end
    
        end
    end

end

figure; hold on
markersize = 50;
hold on; ylabel('')
hold on;
for k=1:size(noise_raster{1,1},1);
    for tt=1:size(noise_raster{1,1},2); 
        if noise_raster{1,1}(k,tt)~=0; 
            plot(tt,noise_raster{1,1}(k,tt)-2*k,'k.','MarkerSize',markersize);

        end;
    end;
end;
xlabel('Time')
ylabel('Neuron #')
title('Figure 3B Surrogate Raster','fontsize',40)
axis off
box off
%%
figure;
h = heatmap(mean(PDF_noise_n1t1,3)); %,'CellLabelColor','none');
h.FontSize = 40;
h.CellLabelFormat = '%.3f';
%h.GridVisible = 'off';
colormap("parula")
xlabel('temporal lag, t_1')
ylabel('spatial lag, n_1')
title({'Figure 3D Surrogate PDF'})
for i = 1:numel(h.XDisplayLabels)
    h.XDisplayLabels{i} = str2num(h.XDisplayLabels{i}) - ((time_window/2) + 1);
end
for i = 1:numel(h.YDisplayLabels)
    h.YDisplayLabels{i} = str2num(h.YDisplayLabels{i}) - ((neuron_window/2) + 1);
end
clim([0 0.35])

%% compute the maximum entropy (based on a uniform distribution)
% Define the ranges for each dimension
n_range = -(neuron_window/2):(neuron_window/2); 
t_range = -(time_window/2):(time_window/2);

% Determine the sizes of the dimensions
n1 = length(n_range); 
n2 = length(n_range); 
t1 = length(t_range); 
t2 = length(t_range); 

% Calculate the total number of elements in the 4D array
total_elements = n1 * n2 * t1 * t2;

% Calculate the uniform probability for each element
uniform_probability = 1 / total_elements;

% Initialize the 4D array with uniform probabilities
PDF_uniform_4d = ones(n1, t1, n2, t2) * uniform_probability;

% Verify that the sum of all probabilities is 1
sum_probabilities = sum(PDF_uniform_4d(:)); 

%Compute the 4D maximum Entropy
bins_4d = [neuron_window+1 time_window+1 neuron_window+1 time_window+1];
entropy_maximum_4d = 0;

for i=1:bins_4d(1)
    for j=1:bins_4d(2)
        for k = 1:bins_4d(3)
            for m = 1:bins_4d(4)
                temp_4d = PDF_uniform_4d;
                if temp_4d(i,j,k,m)~=0        
                    entropy_maximum_4d=entropy_maximum_4d-temp_4d(i,j,k,m)*log2(temp_4d(i,j,k,m));       
                end
            end
        end

    end
end


%% plot the entropies
figure;
hold on
set(0,'defaultAxesFontSize',30)

hold on;
boxplot(entropy_raster_4d,'Colors','r')
boxplot(entropy_noise_4d,'Colors','b')
boxplot(entropy_maximum_4d,'Colors','g')
ylim([0 15])
title('Figure 3E (Network: red, surrogate: blue, uniform/maximum: green)')
ylabel('4D entropy')
