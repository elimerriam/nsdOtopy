% fig2.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: fig2()
%   by: zvi roth
%   date: 7/29/2022
%   purpose: Plot scatter plots for Figure 2
%   uses files created by: getVoxPref.m

close all
clear all
tic

figRoi=1;

toSavePdf = 0;
imgFormat = 'jpg';
subjects = [1:8];

ifig=0;
nrois = 4;

imgScaling = 0.5;
global interpSz; interpSz= 714*imgScaling;
global backgroundSz; backgroundSz= 1024*imgScaling;
global degPerPix; degPerPix = 8.4/(interpSz*imgScaling);

nhistbins = 30;
histAlpha = 0.5;

eccMin = 0.1;
eccMax = 10;
nbins = 20;
binBorders = logspace(log10(eccMin),log10(eccMax),nbins+1);
nbins = length(binBorders)-1;
for i=2:length(binBorders)
    binCenters(i-1) = (binBorders(i)+binBorders(i-1))/2;
end
binBorders;

%bins for pRF angle
angMin = 0;
angMax = 2*pi;
angBinBorders = linspace(angMin,angMax,nbins+1);
nbins = length(angBinBorders)-1;
for i=2:length(angBinBorders)
    angBinCenters(i-1) = (angBinBorders(i)+angBinBorders(i-1))/2;
end
angBinBorders;

%bins for R^2
rMin = 0;
rMax = 1;
rBinBorders = linspace(rMin,rMax,nbins+1);
nbins = length(rBinBorders)-1;
for i=2:length(rBinBorders)
    rBinCenters(i-1) = (rBinBorders(i)+rBinBorders(i-1))/2;
end
rBinBorders;

%scatter parameters
markersize = 1;
edgeAlpha = 0.3;%0.07
markerColor = [0 0 0];
prfThresh = 0;

prffolder = ['~/NSD/prfsample/'];
figFolder = ['/Users/rothzn/Documents/MATLAB/NSD/figures/'];

allOri = cell(1,nrois);
allLevVig = cell(1,nrois);
allLevFull = cell(1,nrois);
allPrfR2 = cell(1,nrois);
allNsdCorr = cell(1,nrois);
allNsdOriCorr = cell(1,nrois);
allSynthCorr = cell(1,nrois);
allSynthOriCorr = cell(1,nrois);
allPrfX = cell(1,nrois);
allPrfY = cell(1,nrois);
allPrfEcc = cell(1,nrois);
allPrfAng = cell(1,nrois);
allSynthOri = cell(1,nrois);
allSynthLevVig = cell(1,nrois);
allSynthLevFull = cell(1,nrois);
allImprovCorr = cell(1,nrois);
allNsdOriR2 = cell(1,nrois);
allNsdR2 = cell(1,nrois);

for isub=1:length(subjects)
    subnum = subjects(isub);
    load([prffolder 'voxModelPref_sub' num2str(isub) '.mat'],'allRoiPrf','roiLevVig','roiLevFull',...
        'roiOri','roiNsdCorr','roiNsdOriCorr','roiNsdOriR2','roiNsdR2',...
        'roiSynthCorr','roiSynthOriCorr','roiSynthOri',...
        'roiSynthLevVig','roiSynthLevFull', ...
        'nsdSynthImprov_pval', 'nsdSynthImprov_corr','roiOriDeviation',...
        'visRoiData','roiNames','combinedRoiNames','roiInd','prefAnalysis');
    subAnalysis(isub) = prefAnalysis;
    subNsdSynthImprov_corr(isub,:,:) = nsdSynthImprov_corr;
    subNsdSynthImprov_pval(isub,:,:) = nsdSynthImprov_pval;
    for iroi=1:nrois
        allPrfX{iroi} = [allPrfX{iroi}; allRoiPrf{iroi}.x];
        allPrfY{iroi} = [allPrfY{iroi}; allRoiPrf{iroi}.y];
        allPrfEcc{iroi} = [allPrfEcc{iroi}; allRoiPrf{iroi}.ecc];
        allPrfAng{iroi} = [allPrfAng{iroi}; allRoiPrf{iroi}.ang];
        allPrfR2{iroi} = [allPrfR2{iroi}; allRoiPrf{iroi}.r2];
        allOri{iroi} =  [allOri{iroi} roiOri{iroi}];
        allLevVig{iroi} = [allLevVig{iroi} roiLevVig{iroi}];
        allLevFull{iroi} = [allLevFull{iroi} roiLevFull{iroi}];
        allNsdCorr{iroi} = [allNsdCorr{iroi} roiNsdCorr{iroi}];
        allNsdOriCorr{iroi} = [allNsdOriCorr{iroi} roiNsdOriCorr{iroi}];
        allNsdOriR2{iroi} = [allNsdOriR2{iroi} roiNsdOriR2{iroi}];
        allNsdR2{iroi} = [allNsdR2{iroi} roiNsdR2{iroi}];
        allSynthCorr{iroi} = [allSynthCorr{iroi} roiSynthCorr{iroi}];
        allSynthOriCorr{iroi} = [allSynthOriCorr{iroi} roiSynthOriCorr{iroi}];
        allSynthOri{iroi} = [allSynthOri{iroi}; roiSynthOri{iroi}'];
        allSynthLevVig{iroi} = [allSynthLevVig{iroi} roiSynthLevVig{iroi}'];
        allSynthLevFull{iroi} = [allSynthLevFull{iroi} roiSynthLevFull{iroi}'];
        allImprovCorr{iroi} = [allImprovCorr{iroi} roiOriDeviation{iroi}];
    end
end

ifig=ifig+1; h=figure(ifig); clf;
rows=2;
cols = 4;
isubplot=0;
iroi=figRoi;
isplit=3;

constColor = [133,149,225]/255;
fullColor = [224,123,145]/255;%pinkish red
r2Color = [156,222,214]/255;
linewidth=2;
%%
%histogram of R2 for full and constrained
isubplot=isubplot+1; subplot(rows,cols,isubplot);
histogram(allNsdR2{iroi}(isplit,:),nhistbins,'faceColor',constColor,'faceAlpha',histAlpha,'Normalization','probability'); hold on;
histogram(allNsdOriR2{iroi}(isplit,:),nhistbins,'faceColor',fullColor,'faceAlpha',histAlpha,'Normalization','probability'); hold on;
xlabel('\itR^2');
ylabel('\itproportion of voxels');
legend('\itconstrained','\itfull');
axis square
mean(allNsdR2{iroi}(isplit,:));
mean(allNsdOriR2{iroi}(isplit,:));
pR2improvMedian = ranksum(allNsdOriR2{iroi}(isplit,:),allNsdR2{iroi}(isplit,:))
[corrR2, pCorrR2] = corr(allNsdCorr{iroi}(isplit,:)',allNsdOriCorr{iroi}(isplit,:)');

% scatter plot of full R2 vs. constrained R2
isubplot=isubplot+1; subplot(rows,cols,isubplot);
scatter(allNsdR2{iroi}(isplit,:),allNsdOriR2{iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 0.2]); ylim([0 0.2]);
axis square
xlabel('\itconstrained R^2');
ylabel('\itfull R^2');
box on

% NSD constrained R2  vs. pRF R2
isubplot=isubplot+1; subplot(rows,cols,isubplot);
scatter(allPrfR2{iroi}./100,allNsdR2{iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 1]);  ylim([0 0.2]); %ylim([0 1]);
axis square
xlabel('\itpRF R^2');
ylabel('\itconstrained R^2');
box on
hold on
temp = allNsdR2{iroi}(isplit,:);
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{iroi}./100>rBinBorders(ibin) & allPrfR2{iroi}./100<=rBinBorders(ibin+1)));
end
p=plot(rBinCenters,binData,'color',constColor,'linestyle','-','linewidth',linewidth);
[corrConstR2PrfR2, pCorrConstR2PrfR2] = corr(allNsdR2{iroi}(isplit,:)',allPrfR2{iroi}./100,'type','Pearson');

% NSD full R2  vs. pRF R2
isubplot=isubplot+1; subplot(rows,cols,isubplot);
scatter(allPrfR2{iroi}./100,allNsdOriR2{iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 1]);  ylim([0 0.2]);
axis square
xlabel('\itpRF R^2');
ylabel('\itfull R^2');
box on
hold on
temp = allNsdOriR2{iroi}(isplit,:);
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{iroi}./100>rBinBorders(ibin) & allPrfR2{iroi}./100<=rBinBorders(ibin+1)));
end
p=plot(rBinCenters,binData,'color',fullColor,'linestyle','-','linewidth',linewidth);
[corrFullR2PrfR2, pCorrFullR2PrfR2] = corr(allNsdOriR2{iroi}(isplit,:)',allPrfR2{iroi}./100,'type','Pearson');


% NSD R2 improvement vs. pRF R2
isubplot=isubplot+1; subplot(rows,cols,isubplot);
improv = allNsdOriR2{iroi}(isplit,:) - allNsdR2{iroi}(isplit,:);
scatter(allPrfR2{iroi}./100,allNsdOriR2{iroi}(isplit,:) - allNsdR2{iroi}(isplit,:),markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 1]);
axis square
xlabel('\itpRF R^2');
ylabel('\itfull R^2 - const. R^2');
box on
hold on
temp = improv;
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{iroi}./100>rBinBorders(ibin) & allPrfR2{iroi}./100<=rBinBorders(ibin+1)));
end
p=plot(rBinCenters,binData,'color',r2Color,'linestyle','-','linewidth',linewidth);
ylim([-0.01 0.03]);
[corrImprovPrfR2, pCorrImprovPrfR2] = corr(temp',allPrfR2{iroi}./100,'type','Pearson');


% scatter plot: full model improvement vs. pRF eccentricity
isubplot=isubplot+1; subplot(rows,cols,isubplot);
improv = allNsdOriR2{iroi}(isplit,:) - allNsdR2{iroi}(isplit,:);
scatter(allPrfEcc{iroi},temp,markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
axis square
xlabel('\itpRF eccentricity (deg)');
ylabel('\itfull R^2 - const. R^2');
box on
hold on
for ibin=1:nbins
    binData(ibin) = mean(temp(allPrfR2{iroi}>prfThresh &  allPrfEcc{iroi}>=binBorders(ibin) & allPrfEcc{iroi}<binBorders(ibin+1)));
end
p=plot(binCenters,binData,'color',r2Color,'linestyle','-','linewidth',linewidth);
xlim([0 10]);
ylim([-0.01 0.03]);

% scatter plot: pRF R^2 vs. pRF eccentricity
isubplot=isubplot+1; subplot(rows,cols,isubplot);
scatter(allPrfEcc{iroi},allPrfR2{iroi}/100,markersize,markerColor,'MarkerFaceAlpha',0,'MarkerEdgeAlpha',edgeAlpha);
xlim([0 10]);
axis square
xlabel('\itpRF eccentricity (deg)');
ylabel('\itpRF R^2');
box on
hold on
for ibin=1:nbins
    binData(ibin) = mean(allPrfR2{iroi}(allPrfR2{iroi}>prfThresh & allPrfEcc{iroi}>=binBorders(ibin) & allPrfEcc{iroi}<binBorders(ibin+1))/100);
end
plot(binCenters,binData,'color',r2Color,'linestyle','-','linewidth',linewidth);
ylim([0 1]);
xlim([0 10]);

%%
set(gcf,'position',[150 180 940 420]);
if toSavePdf
    savepdf(h, [figFolder 'fig2.pdf']);
    saveas(h, [figFolder 'fig2.' imgFormat]);
end

toc
