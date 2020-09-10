function [] = makeSupplementFigureS1()

plotAverage         = false;
plotSupplementalFig = true;
    
makeFigure3AB_SSVEFCoherence(1:10,plotAverage,plotSupplementalFig)
makeFigure3C_SSVEFReliability(1:10,plotAverage,plotSupplementalFig)

return