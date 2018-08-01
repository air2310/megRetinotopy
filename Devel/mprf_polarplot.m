function mprf_polarplot(Amp,ph,hold)

 if notDefined('hold')
     hold=0;
 end

 NaNs = find(isnan(Amp));
 if ~isempty(NaNs)
  %myWarnDlg('contains epochs with nans.  These epochs are being ignored.');
  notNaNs = find(~isnan(Amp));
  Amp = Amp(notNaNs);
  ph = ph(notNaNs);
end
 
% Window header
headerStr = ['Phase vs. Amplitude'];
set(gcf,'Name',headerStr); 

% Plot it
fontSize = 6;
symbolSize = 3;

 % polar plot
  x = Amp.*cos(ph);
  y = Amp.*sin(ph);

  % polar plot params
  params.grid = 'on';
  params.line = 'off';
  params.gridColor = [0.7,0.7,0.7];
  params.fontSize = fontSize;
  params.symbol = 'o';
  params.size = symbolSize;
  params.color = 'w';
  params.fillColor = 'w';
  params.maxAmp = max(Amp(:));
  params.ringTicks = [0:0.2:1].*max(Amp);
  params.gridLineWidth   = 0.01;
  params.lineWidth       = 0.2;
  % Use 'polarPlot' to set up grid'
  %clf
  if hold == 1
  else
      polarPlot(0,params);
  end
  % finish plotting it
  h=plot(x,y,'bo','MarkerSize',symbolSize);
  set(h,'MarkerFaceColor','b')



end