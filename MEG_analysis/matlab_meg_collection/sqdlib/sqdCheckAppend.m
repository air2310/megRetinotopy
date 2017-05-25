%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: sqdCheckAppend.m
%
% Purpose: Check to see whether a file is the result
%   of appending files together. Used automatically
%   by sqdAppend.m
%
% Inputs: 
%   append_file: The appended together file.
%   files: The files that did the appending.
%
% Outputs: None.
%
% Usage: sqdCheckAppend('New.sqd',{'File1.sqd','File2.sqd'})
%
% Author: Doug Bemis
% Date: 6/3/10
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sqdCheckAppend(append_file,files)

disp('Getting info...');
step_size = 50000;
append_info = sqdhandle(append_file);
checked = 0;
for i = 1:length(files)
    
    % See if we've got enough
    disp(['Checking ' files{i} '...']);
    info = sqdhandle(files{i});
    if get(append_info,'SamplesAvailable') < get(info,'SamplesAvailable')+checked
        error('Not enough samples to check. Too many files?');
    end
    
    % If so, check it out
    sample_start = 1;
    total_samples = get(info,'SamplesAvailable');
    while sample_start <= total_samples

        % Set the end
        if (sample_start + step_size - 1) < total_samples
             sample_end = (sample_start + step_size - 1);
        else
            sample_end = total_samples;
        end

        % Read and check
        old_data = sqdread(files{i},'SAMPLES',[sample_start sample_end]);
        new_data = sqdread(append_file,'SAMPLES',[checked+sample_start checked+sample_end]);

        % Check
        % NOTE: Seems to sometimes have a very tiny error.
        % TODO: See if this is important...
        if ~isempty(find(abs(old_data - new_data) > 0.00000000001,1))
            [x y] = find(abs(old_data - new_data) > 0.00000000001,1);
            error(['Found a discrepancy in channel ' num2str(y) ' at Old: ' num2str(checked+sample_start+x) '; New: ' num2str(sample_start+x) '. Exiting...']); 
        end

        % Update
        disp([num2str(sample_end / total_samples * 100) '%...']);
        sample_start = sample_end + 1;
    end
    checked = total_samples + checked;
end
if get(append_info,'SamplesAvailable') ~= checked
    error('Not all samples checked. Maybe some more files?');
end
disp('All done.');
