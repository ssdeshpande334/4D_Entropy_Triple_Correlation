%% Run triple correlation and compute the 4D spatiotemporal distribution, 
% PDF, and associated 4D Shannon's entropy from the rat cortical cultures
% spiking rasters

%datafilename = ''
raster_temp = load(datafilename);
raster_temp = raster_temp.raster_3D;

raster = reshape(raster_temp,size(raster_temp,1) * size(raster_temp,2),size(raster_temp,3));

% get file details, directory path, filename and extension - here it
% shoud be .mat
[fpath, fname, fext] = fileparts(datafilename);
%let us create log file with same name, it will be text file
logfilename = append(fname, '_log.txt');
fidlog = fopen(logfilename, 'at+');

neuron_window = size(raster,1)*2 -2;
time_window = 50;
n_iterations = 1;
epoch_length = 500;
n_epochs = 30;
start_time = 1000;

max_time_lag = ceil((time_window / 2));
[N_neurons, N_times] = size(raster);

size_snippet = 2*max_time_lag + epoch_length + 1;
total_length = size(raster,2);

init = 1;
post_end = start_time;

for ii = init:n_epochs
    tic
    disp(ii)

    if (post_end + 1 + size_snippet) < size(raster,2)

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
        
        tic
        %compute tricorr - no spatial wrapping
        [c3_4D_distribution, actual_contribution,class_count,contribution]= ...
            triple_correlation_class_contributions_no_sp_wr(temp_snippet, neuron_window, time_window);
        toc
        
        actual{ii,1} = actual_contribution./(numel(snippet_raster));
        [conditioned_expectation] = expectation_conditioned_on_constituent_parts_2D(actual{ii,1}, snippet_raster,neuron_window, time_window);
        
        aovert_minus1{ii,1} = (actual{ii,1} ./conditioned_expectation ) -1;
        
        T_raster_4d=sum(sum(sum(sum(c3_4D_distribution))));
        PDF_raster_4d=c3_4D_distribution./T_raster_4d;
        
        % Compute the entropy of the 4D PDF

        temp_4d = PDF_raster_4d + eps;
        entropy_4d_raster(ii,1) = -sum(temp_4d(:) .* log2(temp_4d(:)));
        
    end

end
clear raster_temp c3_4D_distribution PDF_raster_4d P_marginals

resultfile = append('RESULTS_', fname);
save(resultfile);
