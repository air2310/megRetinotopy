%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: sqdAppend.m
%
% Purpose: Append .sqd files together. 
%
% Inputs: 
%   * new_file: The new file to create.
%   * files: The files to append (in the order
%       to be appended.
% 
% Outputs: None.
%
% Usage: sqdAppend('New.sqd',{'File1.sqd','File2.sqd'})
%
% Author:Doug Bemis
% Date: 3/1/09
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sqdAppend(new_file,files)

% This is slightly more complicated because
% MATLAB runs out of memory.

% Write each one
expected_samples = 0;
for i = 1:length(files)
    disp(['Appending ' files{i} '...']);
    disp('Getting info...');
    info = sqdhandle(files{i});
    
    % These things have to be the same for now
    if i == 1
       
        % Get it started
        disp('Starting the file...');
        start_num = 10;
        data = sqdread(files{i},'SAMPLES',[1 start_num]);
        putdata(info,new_file,'Data',data);
        new_info = sqdhandle(new_file);
        sample_start = 1+start_num;

    else
        if (get(info,'ChannelCount') ~= get(new_info,'ChannelCount')) || ...
            (get(info,'SampleRate') ~= get(new_info,'SampleRate')) || ...
            (get(info,'AcquisitionType') ~= get(new_info,'AcquisitionType')) || ...
            (~strcmp(get(info,'Datatype'),get(new_info,'Datatype')))

            error('Infos are not matching. Exiting...');
        end
    
        % Looks like this one's only a warning
        if (get(info,'RawOffset') ~= get(new_info,'RawOffset'))
            disp('WARNING: RawOffsets do not match...');
        end
        
        sample_start = 1;
    end
    
    % Simply go through and apppend
    % Don't want to run out of memory
    step_size = 50000;
    total_samples = get(info,'SamplesAvailable');
    disp('Writing data...');
    while sample_start <= total_samples

        % Set the end
        if (sample_start + step_size - 1) < total_samples
             sample_end = (sample_start + step_size - 1);
        else
            sample_end = total_samples;
        end

        % Read and append
        data = sqdread(files{i},'SAMPLES',[sample_start sample_end]);
        err = putdata(new_info,new_file,'ACTION','Append','DATA',data);

        % Check
        if err ~= (get(info,'ChannelCount') * (sample_end - sample_start + 1))
            error('Something failed writing. Exiting...'); 
        end

        % Update
        disp([num2str(sample_end / total_samples * 100) '%...']);
        sample_start = sample_end + 1;
    end
    expected_samples = expected_samples + get(info,'SamplesAvailable');
end

% Check, thoroughly
sqdCheckAppend(new_file,files);
