function avepoc = ave_epochA(data)
for i=1:length(data(:,1))
    avepoc(i)=mean(data(i,:));
end
plot (avepoc, 'DisplayName','avepoc', 'YDataSource', 'avepoc'); figure(gcf)
end