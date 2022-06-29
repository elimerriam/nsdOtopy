% fig3.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: fig3()
%   by: zvi roth
%   date: 7/29/2022
%   purpose: Plot preferred orientations for Figure 3
%   uses files created by: getVoxPref.m

close all
clear all
tic

toSavePdf = 0;

imgFormat = 'jpg';
subjects = [1:8];
% subjects = [7];
ifig=0;
nrois = 4;

imgScaling = 0.5;
global interpSz; interpSz= 714*imgScaling;
global backgroundSz; backgroundSz= 1024*imgScaling;
global degPerPix; degPerPix = 8.4/interpSz;


%scatter parameters
markersize = 1;
edgeAlpha = 0.3;%0.07
markerColor = [0 0 0];
prfThresh = 0;

prffolder = ['~/NSD/prfsample/'];
figFolder = ['/Users/rothzn/Documents/MATLAB/NSD/figures/'];

allOri = cell(1,nrois);
allResidOri = cell(1,nrois);
allResidOriOri = cell(1,nrois);
allPredOri = cell(1,nrois);
allPredOriOri = cell(1,nrois);

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
allOriDeviation = cell(1,nrois);
allVertDeviation = cell(1,nrois);
allCardDeviation = cell(1,nrois);
allNsdOriR2 = cell(1,nrois);
allNsdR2 = cell(1,nrois);
allSubInd = cell(1,nrois);

for isub=1:length(subjects)
    subnum = subjects(isub);
    load([prffolder 'voxModelPref_sub' num2str(subnum) '.mat'],'allRoiPrf','roiLevVig','roiLevFull',...
        'roiOri','roiNsdCorr','roiNsdOriCorr','roiNsdOriR2','roiNsdR2',...
        'residOri','residOriOri','predOri','predOriOri',...
        'roiSynthCorr','roiSynthOriCorr','roiSynthOri',...
        'roiSynthLevVig','roiSynthLevFull', ...
        'nsdSynthImprov_pval', 'nsdSynthImprov_corr','roiOriDeviation','roiVertDeviation','roiCardDeviation',...
        'visRoiData','roiNames','combinedRoiNames','roiInd','prefAnalysis','nsplits');
    
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
        allResidOri{iroi} =  [allResidOri{iroi} residOri{iroi}];
        allResidOriOri{iroi} =  [allResidOriOri{iroi} roiOri{iroi}];
        allPredOri{iroi} =  [allPredOri{iroi} residOriOri{iroi}];
        allPredOriOri{iroi} =  [allPredOriOri{iroi} predOriOri{iroi}];
        
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
        allOriDeviation{iroi} = [allOriDeviation{iroi} roiOriDeviation{iroi}];
        allVertDeviation{iroi} = [allVertDeviation{iroi} roiVertDeviation{iroi}];
        allCardDeviation{iroi} = [allCardDeviation{iroi} roiCardDeviation{iroi}];
        allSubInd{iroi} = [allSubInd{iroi}; subnum*ones(size(roiOriDeviation{iroi},2),1)];
    end
end

%%
ifig=ifig+1; h=figure(ifig); clf;
rows=2;
cols=3;
isplit = nsplits;
isubplot=0;

iroi=1;

%% preferred ORIENTATION
isubplot=isubplot+1;
subplot(rows,cols, isubplot);
plotOriLines(allOri{iroi}(isplit,:), allPrfX{iroi}, allPrfY{iroi}, allPrfEcc{iroi},(3*allNsdR2{iroi}(isplit,:)));

xlabel('\itx position (deg)');
ylabel('\ity position (deg)');

%% preferred ORIENTATION - single splits
cols=2*cols;
isubplot=cols;
for isplit=1:2
    isubplot = isubplot+1;
    subplot(rows,cols, isubplot);
    plotOriLines(allOri{iroi}(isplit,:), allPrfX{iroi}, allPrfY{iroi}, allPrfEcc{iroi},(3*allNsdR2{iroi}(3,:)));
    
    xlabel('\itx position (deg)');
    if isplit==1
        ylabel('\ity position (deg)');
    else
        set(gca,'yTick',[]);
        set(gca,'yTicklabels',{});
    end
end
isplit=3;

%% orientation color legend
cols = cols/2;
isubplot=2;
subplot(rows,cols, isubplot);
noris = 16;
oris = linspace(0,2*pi,noris+1);
oris = oris(1:end-1);
R = 0.3;
[oriX, oriY] = pol2cart(pi/2-oris,R);
linelength = 0.1;
linewidth = 4;
cMap = turbo(256);

angleColor = zeros(100,100);
[X,Y] = meshgrid(-100:100, -100:100);
[TH,R] = cart2pol(X,Y);
TH(R>100) = NaN;
imagesc(mod(pi/2+TH,pi));
colormap turbo
axis square

axis square
axis off


%% scale legend
isubplot=isubplot+1;
subplot(rows,cols, isubplot);
nlines = 6;
orientations = (pi/4)*ones(1,nlines);
x = 100*linspace(-1, 1, nlines);
y = ones(1,nlines);
snr = linspace(0.01, 0.2, nlines);
snr = 3*snr;
snr = snr*200;
lineWidth = 0.01*snr;
lineLength = 0.4*snr;
cMap = turbo(256);
for iline=1:nlines
    drawOriLine(x(iline), y(iline), pi/2-orientations(iline), lineLength(iline), lineWidth(iline), cMap(1+floor((orientations(iline))*255/(pi)),:));
    hold on
end
xlim([-interpSz interpSz]); ylim([-interpSz interpSz]);

axis square
axis off

%%
set(gcf,'position',[150 180 3*250 rows*210]);
h.Units = 'centimeters';
h.PaperSize=[10 3.5];
if toSavePdf
    print('-painters','-dpdf',[figFolder 'fig3_left']);
end

toc
%%
function prefOri = plotOriLines(prefOri, prfX, prfY, prfEcc,r2)
r2 = (r2)*200;
minWidth = 0.001;
r2(isnan(r2)) = minWidth;
r2(r2<minWidth) = minWidth;
numvox = length(prefOri);
lineWidth = 0.01*r2;
lineLength = 0.4*r2;
cMap = turbo(256);

for ivox=1:numvox
    %if the coefficients are NaN, don't plot
    if ~isnan(prefOri(ivox))
        h=drawOriLine(prfX(ivox), prfY(ivox), pi/2-prefOri(ivox), lineLength(ivox), lineWidth(ivox), cMap(1+floor((prefOri(ivox))*255/(pi)),:));
        hold on
    end
end

global interpSz;% = 714;
global backgroundSz;% = 1024;
global degPerPix;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
linecolor = [0.5 0.5 0.5];
line([0 0], [-backgroundSz backgroundSz],'color',linecolor);
line([-backgroundSz backgroundSz],[0 0], 'color',linecolor);
line([-interpSz/2 -interpSz/2], [-interpSz/2 interpSz/2], 'color',linecolor);
line([interpSz/2 interpSz/2], [-interpSz/2 interpSz/2], 'color',linecolor);
line([-interpSz/2 interpSz/2],[interpSz/2 interpSz/2], 'color',linecolor);
line([-interpSz/2 interpSz/2],[-interpSz/2 -interpSz/2], 'color',linecolor);
xlim([-interpSz interpSz]); ylim([-interpSz interpSz]);
set(gca,'xTick',[-interpSz/2 0 interpSz/2]);
set(gca,'xTicklabels',{-degPerPix*interpSz/2, 0, degPerPix*interpSz/2});
set(gca,'yTick',[-interpSz/2 0 interpSz/2]);
set(gca,'yTicklabels',{-degPerPix*interpSz/2, 0, degPerPix*interpSz/2});
box on
axis square

end
