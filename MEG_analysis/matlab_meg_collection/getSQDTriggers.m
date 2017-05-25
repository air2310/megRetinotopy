%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: getSQDTriggers.m
%
% Finds all the triggers from
% a raw MEG160 data file.
%
% Arguments: 
%   * filename: The name of the .sqd file
%   to read from.
%   * channels: A numerical list of channels
%   to read the triggers from. (Usually 160:167)
%   * out_file: The file to print the triggers
%   in a format for BESA. (Leave as '' to omit).
%
% Author: Alec Marantz
% Date: 7/28/08
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function triggers = getSQDTriggers(filename, channels, out_file)

% Read in the file
[data,info] = sqdread(filename,'Channels',channels);

% Find large changes in amplitude
% between successive time samples
[R,C] = find( diff(data) > 1500);

% Transform the column index into
% a binary digit
for i = 1:length(C)
    C(i) = 2^(C(i)-1);
end;

% Sort by time
A = sortrows ([R C], 1);

% Now, there can be some imprecision
% in the trigger markers, so, we'll
% say anything within trigger_width
% samples is part of the same trigger,
% and that each trigger bit can only 
% fire once per trigger.
trigger_width = 5;
next_trig_start = 1;
triggers = [];
while next_trig_start <= size(A,1)

    % Get the window for the current trigger
    t_start = next_trig_start;
    while next_trig_start <= size(A,1) && ...
            A(next_trig_start,1) <= A(t_start,1) + trigger_width
        next_trig_start = next_trig_start + 1; 
    end
    t_end = next_trig_start - 1;
    
    % Add up the unique triggers
    t_values = [];
    for i = t_start:t_end
        if isempty(find(t_values == A(i,2),1))
            t_values(length(t_values)+1) = A(i,2); 
        end
    end
    
    % And set the trigger
    triggers(size(triggers,1)+1,:) = [A(t_start,1) sum(t_values)];
end

% Transform time into msecs
triggers(:,1) = triggers(:,1) * (1000/info.SampleRate);

% If there's an out file, convert them to BESA format
if ~isempty(out_file)
    
    % Open it up
    evt_file = fopen(out_file,'w');
    fprintf(evt_file,'Tmu        \t\tCode\tTriNo\tComnt\n');
    
            
    % Write the triggers to the file
    for i = 1:size(triggers,1)
        fprintf(evt_file, '%d\t1\t%d\tComment - TODO\n', [triggers(i,1)*1000 triggers(i,2)]);
    end
    
    % And done
    fclose(evt_file);
end
