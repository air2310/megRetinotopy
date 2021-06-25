function [] = makeSupplementFigureS1()
% Function to create the Supplement Figure S1: MEG head plots with
% SSVEF 10 Hz coherence and split half reliability for individual subjects.

plotAverage         = false;
plotSupplementalFig = true;
    
makeFigure3AB_SSVEFCoherence(1:10,plotAverage,plotSupplementalFig)
makeFigure3C_SSVEFReliability(1:10,plotAverage,plotSupplementalFig)

return