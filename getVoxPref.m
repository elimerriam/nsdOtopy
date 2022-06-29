% getVoxPref.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: getVoxPref(1,1)
%   by: zvi roth
%   date: 7/29/2022
%   purpose: extract preferred orientation from regression weights, and sum
%   across partitions
%   uses files created by: regressPrfSplit.m
%   creates files used by: fig$$.m


function getVoxPref(isub,numregions)

%uses data from regressPrfSplit.m

%chooses preferred orientation/levels by averaging across filters
mrQuit
close all
% clear all
global interpSz;
global backgroundSz;
global degPerPix;
global prefAnalysis;

prefAnalysis = 3;

toSavePdf = 0;

numOrientations = 8;
figFolder = ['/Users/rothzn/Documents/MATLAB/NSD/figures/'];
nperms=1000;
prffolder = ['~/NSD/prfsample/'];

interpSz= 714;
backgroundSz= 1024;
bandpass = 1; bandMin = 1; bandMax = 7;
bandpassStr = '';
if bandpass
    bandpassStr = ['_bandpass' num2str(bandMin) 'to' num2str(bandMax)];
end

imgScaling = 0.5;
interpSz= 714*imgScaling;
backgroundSz= 1024*imgScaling;
degPerPix = 8.4/(714*imgScaling);

gratings = load(['~/NSD/gratings/gratings.mat'],'cpds','angles','freqs','numOrientations','numLevels',...
    'sumOriEnergy','modelOriEnergy','normResp','backgroundSize','imgScaling');
for iregion=1:numregions
    visualRegion = iregion;%V1,V2,V3,V4
    load(fullfile(prffolder,['regressPrfSplit' bandpassStr '_v' num2str(visualRegion) '_sub' num2str(isub)  '.mat']), ...
        'nsd', 'synth',...
        'numLevels', 'numOrientations','rois','nvox','roiPrf','nsplits');
    
    if length(rois)>1 %combine across ventral and dorsal ROIs
        oldNsd = nsd;
        nsd.voxResidual{1} = [];
        nsd.voxOriResidual{1} = [];
        nsd.voxResidualSplit{1} = [];
        nsd.voxOriResidualSplit{1} = [];
        nsd.r2{1} = [];
        nsd.r2ori{1} = [];
        nsd.r2split{1} = [];
        nsd.r2oriSplit{1} = [];
        nsd.voxCoef{1} = [];
        nsd.voxOriCoef{1} = [];
        nsd.voxPredOriCoef{1} = [];
        nsd.voxOriPredOriCoef{1} = [];
        nsd.voxResidOriCoef{1} = [];
        nsd.voxOriResidOriCoef{1} = [];
        
        nsd.voxPredOriR2{1} = [];
        nsd.voxOriPredOriR2{1} = [];
        nsd.voxResidOriR2{1} = [];
        nsd.voxOriResidOriR2{1} = [];
        
        nsd.pearsonRori{1} = [];
        nsd.pearsonR{1} = [];
        
        nsd.roiInd{1} = [];
        
        oldSynth = synth;
        synth.pearsonR{1} = [];
        synth.pearsonRori{1} = [];
        synth.r2ori{1} = [];
        synth.voxResidual{1} = [];
        synth.voxOriResidual{1} = [];
        synth.voxCoef{1} = [];
        synth.voxOriCoef{1} = [];
        
        for iroi=1:length(rois)
            nsd.voxResidual{1} = cat(2,nsd.voxResidual{1},oldNsd.voxResidual{iroi});
            nsd.voxOriResidual{1} = cat(2,nsd.voxOriResidual{1},oldNsd.voxOriResidual{iroi});
            nsd.voxResidualSplit{1} = cat(2,nsd.voxResidualSplit{1},oldNsd.voxOriResidual{iroi});
            nsd.voxOriResidualSplit{1} = cat(2,nsd.voxOriResidualSplit{1},oldNsd.voxOriResidual{iroi});
            
            nsd.pearsonRori{1} = cat(2,nsd.pearsonRori{1},oldNsd.pearsonRori{iroi});
            nsd.pearsonR{1} = cat(2,nsd.pearsonR{1},oldNsd.pearsonR{iroi});
            nsd.r2{1} = cat(2,nsd.r2{1},oldNsd.r2{iroi});
            nsd.r2ori{1} = cat(2,nsd.r2ori{1},oldNsd.r2ori{iroi});
            nsd.r2split{1} = cat(2,nsd.r2split{1},oldNsd.r2split{iroi});
            nsd.r2oriSplit{1} = cat(2,nsd.r2oriSplit{1},oldNsd.r2oriSplit{iroi});
            nsd.voxCoef{1} = cat(2,nsd.voxCoef{1},oldNsd.voxCoef{iroi});
            nsd.voxOriCoef{1} = cat(2,nsd.voxOriCoef{1},oldNsd.voxOriCoef{iroi});
            nsd.voxPredOriCoef{1} = cat(2,nsd.voxPredOriCoef{1},oldNsd.voxPredOriCoef{iroi});
            nsd.voxOriPredOriCoef{1} = cat(2,nsd.voxOriPredOriCoef{1},oldNsd.voxOriPredOriCoef{iroi});
            nsd.voxResidOriCoef{1} = cat(2,nsd.voxResidOriCoef{1},oldNsd.voxResidOriCoef{iroi});
            nsd.voxOriResidOriCoef{1} = cat(2,nsd.voxOriResidOriCoef{1},oldNsd.voxOriResidOriCoef{iroi});
            
            nsd.voxPredOriR2{1} = cat(2,nsd.voxPredOriR2{1},oldNsd.voxPredOriR2{iroi});
            nsd.voxOriPredOriR2{1} = cat(2,nsd.voxOriPredOriR2{1},oldNsd.voxOriPredOriR2{iroi});
            nsd.voxResidOriR2{1} = cat(2,nsd.voxResidOriR2{1},oldNsd.voxResidOriR2{iroi});
            nsd.voxOriResidOriR2{1} = cat(2,nsd.voxOriResidOriR2{1},oldNsd.voxOriResidOriR2{iroi});
            
            nsd.roiInd{1} = cat(1,nsd.roiInd{1}, oldNsd.roiInd{iroi});
            
            synth.voxResidual{1} = cat(2,synth.voxResidual{1},oldSynth.voxResidual{iroi});
            synth.voxOriResidual{1} = cat(2,synth.voxOriResidual{1},oldSynth.voxOriResidual{iroi});
            synth.pearsonRori{1} = cat(2,synth.pearsonRori{1},oldSynth.pearsonRori{iroi});
            synth.pearsonR{1} = cat(2,synth.pearsonR{1},oldSynth.pearsonR{iroi});
            synth.voxCoef{1} = cat(1,synth.voxCoef{1},oldSynth.voxCoef{iroi});
            synth.voxOriCoef{1} = cat(1,synth.voxOriCoef{1},oldSynth.voxOriCoef{iroi});
            
        end
        oldPrf = roiPrf; clear roiPrf;
        roiPrf{1}.ecc=[];
        roiPrf{1}.ang=[];
        roiPrf{1}.sz=[];
        %         roiPrf{1}.exponent=[];
        %         roiPrf{1}.gain=[];
        roiPrf{1}.r2=[];
        roiPrf{1}.x=[];
        roiPrf{1}.y=[];
        for iroi=1:length(rois)
            roiPrf{1}.ecc = cat(1,roiPrf{1}.ecc,oldPrf{iroi}.ecc);
            roiPrf{1}.ang = cat(1,roiPrf{1}.ang,oldPrf{iroi}.ang);
            roiPrf{1}.sz = cat(1,roiPrf{1}.sz,oldPrf{iroi}.sz);
            %             roiPrf{1}.exponent = cat(1,roiPrf{1}.exponent,oldPrf{iroi}.exponent);
            %             roiPrf{1}.gain = cat(1,roiPrf{1}.gain,oldPrf{iroi}.gain);
            roiPrf{1}.r2 = cat(1,roiPrf{1}.r2,oldPrf{iroi}.r2);
            roiPrf{1}.x = cat(1,roiPrf{1}.x,oldPrf{iroi}.x);
            roiPrf{1}.y = cat(1,roiPrf{1}.y,oldPrf{iroi}.y);
        end
        rois = 1;
    end
    
    %% AVERAGE SPLITS
    nsd.voxCoef{1}(nsplits+1,:,:) = mean(nsd.voxCoef{1},1);
    nsd.voxOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxOriCoef{1},1);
    nsd.voxPredOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxPredOriCoef{1},1);
    nsd.voxOriPredOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxOriPredOriCoef{1},1);
    nsd.voxResidOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxResidOriCoef{1},1);
    nsd.voxOriResidOriCoef{1}(nsplits+1,:,:) = mean(nsd.voxOriResidOriCoef{1},1);
    
    nsd.pearsonRori{1}(nsplits+1,:) = mean(nsd.pearsonRori{1},1);
    nsd.pearsonR{1}(nsplits+1,:) = mean(nsd.pearsonR{1},1);
    nsd.r2{1}(nsplits+1,:) = mean(nsd.r2{1},1);
    nsd.r2ori{1}(nsplits+1,:) = mean(nsd.r2ori{1},1);
    nsd.r2split{1}(nsplits+1,:) = mean(nsd.r2split{1},1);
    nsd.r2oriSplit{1}(nsplits+1,:) = mean(nsd.r2oriSplit{1},1);
    
    synth.voxResidual{1}(nsplits+1,:,:) = mean(synth.voxResidual{1},1);
    synth.voxOriResidual{1}(nsplits+1,:,:) = mean(synth.voxOriResidual{1},1);
    synth.pearsonRori{1}(nsplits+1,:) = mean(synth.pearsonRori{1},1);
    synth.pearsonR{1}(nsplits+1,:) = mean(synth.pearsonR{1},1);
    
    nsplits = nsplits+1;
    
    %% COMPUTE PREFERRED LEVEL - CENTER OF MASS
    clear vigPrefLevel fullPrefLevel residPrefLevel residOriPrefLevel predCoefLevel predOriCoefLevel
    for  iroi=1:length(rois)
        for isplit=1:nsplits
            %constrained model
            vigCoef = squeeze(nsd.voxCoef{iroi}(isplit,:,1:end-1));
            vigPrefLevel{iroi}(isplit,:) = gratingPrefFreq(vigCoef,gratings.sumOriEnergy);
            
            %full orientation model
            fullCoef = squeeze(nsd.voxOriCoef{iroi}(isplit,:,1:end-1));
            fullPrefLevel{iroi}(isplit,:) = gratingPrefFreq(fullCoef,gratings.modelOriEnergy);
            
            %full model on residuals of constrained model
            fullCoef = squeeze(nsd.voxResidOriCoef{iroi}(isplit,:,1:end-1));
            residPrefLevel{iroi}(isplit,:) = gratingPrefFreq(fullCoef,gratings.modelOriEnergy);
            
            %full model on residual of full model
            fullCoef = squeeze(nsd.voxOriResidOriCoef{iroi}(isplit,:,1:end-1));
            residOriPrefLevel{iroi}(isplit,:) = gratingPrefFreq(fullCoef,gratings.modelOriEnergy);
            
            %full model on constrained prediction
            fullCoef = squeeze(nsd.voxPredOriCoef{iroi}(isplit,:,1:end-1));
            predCoefLevel{iroi}(isplit,:) = gratingPrefFreq(fullCoef,gratings.modelOriEnergy);
            
            %full model on full model prediction
            fullCoef = squeeze(nsd.voxOriPredOriCoef{iroi}(isplit,:,1:end-1));
            predOriCoefLevel{iroi}(isplit,:) = gratingPrefFreq(fullCoef,gratings.modelOriEnergy);
            
        end
        %%%%%%%%%%%%%% SYNTH
        %constrained model
        synthVigCoef = squeeze(synth.voxCoef{iroi}(:,1:end-1));
        synthVigPrefLevel{iroi} = gratingPrefFreq(synthVigCoef,gratings.sumOriEnergy);
        
        %full orientation model
        synthFullCoef = squeeze(synth.voxOriCoef{iroi}(:,1:end-1));
        synthFullPrefLevel{iroi} = gratingPrefFreq(synthFullCoef,gratings.modelOriEnergy);
    end
    
    
    %% COMPUTE PREFERRED ORIENTATION - CIRCULAR CENTER OF MASS
    
    clear fullPrefOri residPrefOri residOriPrefOri predCoefOri predOriCoefOri oriDeviation vertDeviation cardDeviation
    clear fullOriModul fullPrefAmp fullAntiAmp
    for  iroi=1:length(rois)%rois=1
        for isplit=1:nsplits
            %full orientation model
            fullCoef = squeeze(nsd.voxOriCoef{iroi}(isplit,:,1:end-1));
            fullPrefOri{iroi}(isplit,:) = gratingPrefOri(fullCoef,gratings.modelOriEnergy);
            
            %full model on residuals of constrained model
            fullCoef = squeeze(nsd.voxResidOriCoef{iroi}(isplit,:,1:end-1));
            residPrefOri{iroi}(isplit,:) = gratingPrefOri(fullCoef,gratings.modelOriEnergy);
            
            %full model on residual of full model
            fullCoef = squeeze(nsd.voxOriResidOriCoef{iroi}(isplit,:,1:end-1));
            residOriPrefOri{iroi}(isplit,:) = gratingPrefOri(fullCoef,gratings.modelOriEnergy);
            
            %full model on constrained prediction
            fullCoef = squeeze(nsd.voxPredOriCoef{iroi}(isplit,:,1:end-1));
            predCoefOri{iroi}(isplit,:) = gratingPrefOri(fullCoef,gratings.modelOriEnergy);
            
            %full model on full model prediction
            fullCoef = squeeze(nsd.voxOriPredOriCoef{iroi}(isplit,:,1:end-1));
            predOriCoefOri{iroi}(isplit,:) = gratingPrefOri(fullCoef,gratings.modelOriEnergy);
            
        end
        %cross-validated measure of orientation modulation
        for isplit=1:2
            fullCoef = squeeze(nsd.voxOriCoef{iroi}(isplit,:,1:end-1));
            [fullOriModul{iroi}(isplit,:), fullPrefAmp{iroi}(isplit,:), fullAntiAmp{iroi}(isplit,:)] = gratingOriModulation(fullPrefOri{iroi}(3-isplit,:),fullCoef,gratings.modelOriEnergy);
        end
        fullOriModul{iroi}(nsplits,:) = mean(fullOriModul{iroi}(1:2,:),1);
        fullPrefAmp{iroi}(nsplits,:) = mean(fullPrefAmp{iroi}(1:2,:),1);
        fullAntiAmp{iroi}(nsplits,:) = mean(fullAntiAmp{iroi}(1:2,:),1);
        
        %%%%%%%%%%%%% SYNTH
        %full orientation model
        synthFullCoef = squeeze(synth.voxOriCoef{iroi}(:,1:end-1));
        synthFullPrefOri{iroi} = gratingPrefOri(synthFullCoef,gratings.modelOriEnergy);
    end
    
    %% SCATTER PLOT ORIENTATION DEVIATION FROM RADIAL VS ECCENTRICITY
    for  iroi=1:length(rois)%rois=1
        for isplit=1:nsplits
            %multiply by 2 for circular dist, then divide by 2
            oriDeviation{iroi}(isplit,:) = 0.5*circ_dist(2*(pi/2-fullPrefOri{iroi}(isplit,:)),2*(mod(roiPrf{iroi}.ang*pi/180,pi)'));
        end
    end
    %% SCATTER PLOT ORIENTATION DEVIATION FROM VERTICAL VS ECCENTRICITY
    for  iroi=1:length(rois)%rois=1
        for isplit=1:nsplits
            vertDeviation{iroi}(isplit,:) = 0.5*circ_dist(2*fullPrefOri{iroi}(isplit,:),0);
        end
    end
    
    %% SCATTER PLOT ORIENTATION DEVIATION FROM CARDINAL VS ECCENTRICITY
    for  iroi=1:length(rois)%rois=1
        for isplit=1:nsplits
            vertAngle1 = roiPrf{iroi}.ang*pi/180>=pi/4 & roiPrf{iroi}.ang*pi/180<3*pi/4;
            vertAngle2 = roiPrf{iroi}.ang*pi/180>=5*pi/4 & roiPrf{iroi}.ang*pi/180<7*pi/4;
            vertAngle = vertAngle1 | vertAngle2;
            horizAngle1 = roiPrf{iroi}.ang*pi/180<pi/4;
            horizAngle2 = roiPrf{iroi}.ang*pi/180>=3*pi/4 & roiPrf{iroi}.ang*pi/180<5*pi/4;
            horizAngle3 = roiPrf{iroi}.ang*pi/180>=7*pi/4;
            horizAngle = horizAngle1 | horizAngle2 | horizAngle3;
            %multiply by 2 for circular dist to work on pi cycle instead of 2pi.
            cardDeviation{iroi}(isplit,vertAngle) = 0.5*circ_dist(2*fullPrefOri{iroi}(isplit,vertAngle),0);
            cardDeviation{iroi}(isplit,horizAngle) = 0.5*circ_dist(2*fullPrefOri{iroi}(isplit,horizAngle),pi);%using pi instead of pi/2
        end
    end
    
    %% Distribution of cross-validated SYNTHETIC  R^2 across voxels
    for  iroi=1:length(rois)%rois=1
        for isplit=1:nsplits
            [nsdSynthImprov_corr(iregion,isplit), nsdSynthImprov_pval(iregion,isplit)] = corr((nsd.pearsonRori{iroi}(isplit,:)-nsd.pearsonR{iroi}(isplit,:))', (synth.pearsonRori{iroi}(isplit,:) - synth.pearsonR{iroi}(isplit,:))', 'rows','complete');
        end
    end
    
    %%
    %save preferred orientation and level for this ROI
    allRoiPrf{iregion} = roiPrf{iroi};%iroi=1
    roiLevVig{iregion} = vigPrefLevel{iroi};
    roiLevFull{iregion} = fullPrefLevel{iroi};
    roiOri{iregion} = fullPrefOri{iroi};
    residOri{iregion} = residPrefOri{iroi};
    residOriOri{iregion} = residOriPrefOri{iroi};
    predOri{iregion} = predCoefOri{iroi};
    predOriOri{iregion} = predOriCoefOri{iroi};
    oriModulation{iregion} = fullOriModul{iroi};
    oriPrefAmp{iregion} = fullPrefAmp{iroi};
    oriAntiAmp{iregion} = fullAntiAmp{iroi};
    
    
    roiInd{iregion} = nsd.roiInd{iroi};
    roiNsdCorr{iregion} = nsd.pearsonR{iroi};
    roiNsdOriCorr{iregion} = nsd.pearsonRori{iroi};
    roiNsdOriR2{iregion} = nsd.r2oriSplit{iroi};
    roiNsdR2{iregion} = nsd.r2split{iroi};
    
    roiNsdOriPredOriR2{iregion} = nsd.voxOriPredOriR2{iroi};
    roiNsdOriResidOriR2{iregion} = nsd.voxOriResidOriR2{iroi};
    roiNsdPredOriR2{iregion} = nsd.voxPredOriR2{iroi};
    roiNsdResidOriR2{iregion} = nsd.voxResidOriR2{iroi};
    
    roiSynthCorr{iregion} = synth.pearsonR{iroi};
    roiSynthOriCorr{iregion} = synth.pearsonRori{iroi};
    roiSynthOri{iregion} = synthFullPrefOri{iroi};
    roiSynthLevVig{iregion} = synthVigPrefLevel{iroi};
    roiSynthLevFull{iregion} = synthFullPrefLevel{iroi};
    roiOriDeviation{iregion} = oriDeviation{iroi};
    roiVertDeviation{iregion} = vertDeviation{iroi};
    roiCardDeviation{iregion} = cardDeviation{iroi};
end


%save all ROIs to create overlay
roifolder = ['~/NSD/sub' num2str(isub) '_betas_func1pt8mm/'];
visualRoisFile = fullfile(roifolder,'prf-visualrois.nii');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
combinedRoiNames = {'V1','V2','V3','hV4'};

save([prffolder 'voxModelPref_sub' num2str(isub) '.mat'],'allRoiPrf','roiLevVig','roiLevFull',...
    'roiOri','roiNsdCorr','roiNsdOriCorr','roiNsdOriR2','roiNsdR2',...
    'roiNsdResidOriR2','roiNsdOriResidOriR2','roiNsdPredOriR2','roiNsdOriPredOriR2',...
    'residOri','residOriOri','predOri','predOriOri',...
    'oriModulation','oriPrefAmp','oriAntiAmp',...
    'roiSynthCorr','roiSynthOriCorr','roiSynthOri',...
    'roiSynthLevVig','roiSynthLevFull', ...
    'nsdSynthImprov_pval', 'nsdSynthImprov_corr','roiOriDeviation','roiVertDeviation','roiCardDeviation',...
    'visRoiData','roiNames','combinedRoiNames','roiInd','prefAnalysis','nsplits');

%%
    function prefAngle = gratingPrefOri(fullCoef,modelOriEnergy);
        [numFreqs, numAngles, numLevels, numOrientations] = size(modelOriEnergy);
        numvox = size(fullCoef,1);
        coefMat = reshape(fullCoef,numvox,numLevels,numOrientations);%vox x levels x orientations
        coefVec = reshape(fullCoef,numvox,numLevels*numOrientations);%vox x coefficients
        modelEnergy = reshape(modelOriEnergy,numFreqs*numAngles,numLevels*numOrientations);
        modelEnergy = modelEnergy';%(ilev*orientation, isf*iangle)
        
        voxGratingResp = coefVec*modelEnergy;%vox,isf*iangle
        angles = linspace(0,180,numAngles+1);
        angles = angles(1:numAngles);
        switch prefAnalysis
            case 1
                % 1 - simple max
                [m maxInd] = max(voxGratingResp,[],2);%vox
                [prefFreqNum, prefAngleNum] = ind2sub([numFreqs, numAngles],maxInd);%vox
                prefAngle = angles(prefAngleNum)*pi/180;
            case 2
                % 2 - average across gratings spatial frequency, and then max
                voxAngleResp = squeeze(mean(reshape(voxGratingResp,numvox, numFreqs, numAngles),2));
                [m prefAngleNum] = max(voxAngleResp,[],2);%ivox
                prefAngle = angles(prefAngleNum)*pi/180;
            case 3
                % 3 - average across gratings spatial frequency, and then circular weighted mean
                voxAngleResp = squeeze(mean(reshape(voxGratingResp,numvox, numFreqs, numAngles),2));
                %subtract minimum response (might be negative)
                voxAngleResp = voxAngleResp - min(voxAngleResp,[],2);
                theta = linspace(0,2*pi,numAngles+1);%for circular calculation
                theta = theta(1:end-1);
                for ivox=1:numvox
                    prefAngle(ivox) = circ_mean(theta',voxAngleResp(ivox,:)');
                end
                prefAngle = mod(prefAngle,2*pi);%from [-pi, pi] to [0 2pi]
                prefAngle = prefAngle./2;%range 0 to pi.
        end
    end

%%
    function [oriModul, prefOriAmp, antiPrefOriAmp, prefOriInd] = gratingOriModulation(prefOri,fullCoef,modelOriEnergy);
        [numFreqs, numAngles, numLevels, numOrientations] = size(modelOriEnergy);
        numvox = size(fullCoef,1);
        coefMat = reshape(fullCoef,numvox,numLevels,numOrientations);%vox x levels x orientations
        coefVec = reshape(fullCoef,numvox,numLevels*numOrientations);%vox x coefficients
        modelEnergy = reshape(modelOriEnergy,numFreqs*numAngles,numLevels*numOrientations);
        modelEnergy = modelEnergy';%(ilev*orientation, isf*iangle)
        
        voxGratingResp = coefVec*modelEnergy;%vox,isf*iangle
        angles = linspace(0,180,numAngles+1);
        angles = angles(1:numAngles);
        prefOriInd = zeros(numvox,1);
        minDist = zeros(numvox,1);
        prefOriAmp = zeros(numvox,1);
        antiPrefOriAmp = zeros(numvox,1);
        
        randorder = randperm(numvox);
        for ivox=1:numvox
            [minDist(ivox), prefOriInd(ivox)] = min(abs(circ_dist(angles*pi/180,prefOri(ivox))));%multiplied by pi/180 to be circular around 2pi
        end
        antiPrefOriInd = 1+mod(prefOriInd - 1 - length(angles)/2,length(angles));
        
        %how to deal with frequency: 1 - average across frequencies. 2 -
        %use preferred frequency
        
        voxGratingResp = reshape(voxGratingResp, numvox,numFreqs,numAngles);
        for ivox=1:numvox
            prefOriAmp(ivox) = mean(voxGratingResp(ivox,:,prefOriInd(ivox)),2);%mean across frequencies
            antiPrefOriAmp(ivox) = mean(voxGratingResp(ivox,:,antiPrefOriInd(ivox)),2);
        end
        oriModul = (prefOriAmp - antiPrefOriAmp)./(prefOriAmp + antiPrefOriAmp);
        
    end

%%
    function prefFreqNum = gratingPrefFreq(fullCoef,modelOriEnergy);
        %         global prefAnalysis
        [numFreqs, numAngles, numLevels, numOrientations] = size(modelOriEnergy);
        numvox = size(fullCoef,1);
        coefMat = reshape(fullCoef,numvox,numLevels,numOrientations);%vox x levels x orientations
        coefVec = reshape(fullCoef,numvox,numLevels*numOrientations);%vox x coefficients
        modelEnergy = reshape(modelOriEnergy,numFreqs*numAngles,numLevels*numOrientations);
        modelEnergy = modelEnergy';%(ilev*orientation, isf*iangle)
        
        voxGratingResp = coefVec*modelEnergy;%vox,isf*iangle
        switch prefAnalysis
            case 1
                % 1 - simple max
                [m maxInd] = max(voxGratingResp,[],2);%vox
                [prefFreqNum, prefAngleNum] = ind2sub([numFreqs, numAngles],maxInd);%vox
            case 2
                % 2 - average across gratings orientation, and then max
                voxFreqResp = squeeze(mean(reshape(voxGratingResp,numvox, numFreqs, numAngles),3));
                [m prefFreqNum] = max(voxFreqResp,[],2);%ivox
            case 3
                % 3 - average across gratings orientation, and then weighted mean
                voxFreqResp = squeeze(mean(reshape(voxGratingResp,numvox, numFreqs, numAngles),3));
                %subtract minimum response
                voxFreqResp = voxFreqResp - min(voxFreqResp,[],2);
                prefFreqNum = (voxFreqResp * [1:numFreqs]')./sum(voxFreqResp,2);
        end
    end
end