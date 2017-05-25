function triggers = all_triggerd(filename, channels)
[data,info] = sqdread(filename,'Channels',channels);
[R,C] = find( -diff(data) > (15000/(info.InputGain * info.OutputGain)));
for i = 1:length(C)
    C(i) = 2^(C(i)-1);
end;
A = sortrows ([R C], 1);
G = find(diff(A(:,1))>0);
counter = 1;
trigcount = 1;
triggers(trigcount,1) = A(G(1), 1);
triggers(trigcount,2) = sum(A(counter:G(1), 2));
for i = 2:(length(G))
    if A(G(i),1) - triggers(trigcount,1) > 1
        trigcount = trigcount + 1;
        triggers(trigcount,1) = A(G(i),1);
        triggers(trigcount,2) = sum(A(counter:G(i), 2));
    end;
        counter = G(i) +1;
end;
if A(length(A), 1) - triggers(trigcount,1) > 1
    trigcount = trigcount + 1;
    triggers(trigcount,1) = A(length(A),1);
    triggers(trigcount,2) = sum(A(counter:length(A),2));
end;
triggers(:,1) = triggers(:,1) * (1000/info.SampleRate);