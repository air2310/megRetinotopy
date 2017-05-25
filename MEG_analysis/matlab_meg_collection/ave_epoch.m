function avepoc = ave_epoch(datastack)
for i=1:length(datastack(:,1,1))
    for j=1:length(datastack(1,:,1))
        avepoc(i,j)=mean(datastack(i,j,:));
    end;
end;
plot (avepoc, 'DisplayName','avepoc', 'YDataSource', 'avepoc'); figure(gcf)
end