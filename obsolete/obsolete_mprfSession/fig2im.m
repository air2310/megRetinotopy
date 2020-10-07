function [] = fig2im(output_type,quality,which_figs,figNames,renderer)
% Prints matlab figures to file format.
% output_type:  output format, default is .pdf
% quality:      quality of the figure. Higher numbers is higher quality.
%               default is 600
% which_figs:   figure handel (usually number of figure), can be vector. 
%               If not defined, will attempt to convert all open figures. 
%               Assumes that figure handles are 1 : n_figures in that case.        
% figNames:     cell with size which_figs for names of the figure after
%               saving. Usage is recommended (speeds up saving)
% renderer:     Type of rendered used by print. Default 'painters'
% All inputs are optional

% 11-2011, BK:  wrote it




if notDefined('figNames')
    useFigName =1;
else
    if iscell(figNames)
        useFigName = 0;
    else
        error('figNames must be cell')
    end
end

if notDefined('output_type')
    output_type = 'pdf';
end

if notDefined('quality')
    quality = '600';
end

if notDefined('which_figs')
    which_figs = sort(get(0,'children'))';
end

if notDefined('renderer')
    renderer = 'painters';
end

renderer = ['-' lower(renderer)];

if useFigName == 0
    if size(figNames) == size(which_figs)
    else
        error('which_figs and figNames must be the same size')
    end
end

output_type = ['-d' output_type];
quality = ['-r' quality];


for b = which_figs;
    figure(b)
    if useFigName == 1
        tmp = get(gcf);
        imName = tmp.FileName;
        ind = imName == '.';
        imName = imName(~ind);
        
    else
        imName = figNames{find(which_figs==b)};
    end
    
    if strcmpi(renderer,'-opengl')
        ps = get(gcf,'PaperSize');
        
    elseif strcmpi(renderer,'-painters')
        ps = get(gcf,'PaperSize');
        
    end
    
    %fillPage(gcf,'margins',[0 0 0 0],'paperunits','centimeters','papersize',[2 2]);
    fillPage(gcf,'margins',[0 0 0 0],'papersize',[100 28.125]);
    print(renderer, output_type,quality, imName);
end


end