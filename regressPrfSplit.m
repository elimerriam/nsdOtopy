% regressPrfSplit.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: regressPrfSplit(1,1)
%   by: zvi roth
%   date: 7/29/2022
%   purpose: Perform linear regression on the response amplitudes for each voxel with filter output values as predictors
%   uses files created by: prfSampleModel.m, prfSampleModel_synth.m
%   creates files used by: getVoxPref.m

function regressPrfSplit(isub,visualRegions)

tic

%%%%%%%%%%%%%%%%%%%%%%%%
nsessionsSub = [40 40 32 30 40 32 40 30];
nsessions=nsessionsSub(isub);
nsplits=2;
bandpass = 1; bandMin = 1; bandMax = 7;
boxfolder = ['~/Box/NSD/prfsample' noBackgroundStr '/'];
prffolder = ['~/NSD/prfsample' noBackgroundStr '/'];
betasfolder = ['~/NSD/sub' num2str(isub) '_betas_func1pt8mm/'];
% stimfilename = fullfile(folder,'nsdsynthetic_colorstimuli_subj01.hdf5');
nsdfolder = '~/NSD/';
roifolder = ['~/NSD/sub' num2str(isub) '_betas_func1pt8mm/'];
visualRoisFile = fullfile(roifolder,'prf-visualrois.nii');%V1v, V1d, V2v, V2d, V3v, V3d, and hV4
visRoiData = niftiread(visualRoisFile);
roiNames = {'V1v','V1d','V2v','V2d','V3v','V3d','hV4'};
visRoiData = visRoiData(:);

for visualRegion=visualRegions
    
    load(fullfile(boxfolder,['prfSampleStim_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
        'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg',...
        'roiPrf');
    %if prf sampling was done with the nonlinear CSS prf, then we want to
    %define the weights for the constrained model as a sum of the
    %orientation model across orientations:
    for roinum=1:length(rois)
        prfSampleLev{roinum} = squeeze(sum(prfSampleLevOri{roinum},4));
    end
    
    if bandpass
        for roinum=1:length(rois)
            prfSampleLevOri{roinum} = prfSampleLevOri{roinum}(:,:,bandMin:bandMax,:);
            prfSampleLev{roinum} = prfSampleLev{roinum}(:,:,bandMin:bandMax);
        end
        numLevels = bandMax-bandMin+1;
    end
    
    for roinum=1:length(rois); iroi = rois(roinum); roiBetas{roinum}=[]; end
    for isession=1:nsessions
        betasfilename = fullfile(betasfolder,['betas_session' num2str(isession,'%02.f') '.nii']);
        betas = niftiread(betasfilename);
        betas = cast(betas,'double');
        betas = betas/300;
        betas=reshape(betas,[],size(betas,4));
        for roinum=1:length(rois)
            iroi = rois(roinum);
            roiBetas{roinum} = [roiBetas{roinum} betas(visRoiData==iroi,:)];
            roiInd{roinum} = find(visRoiData==iroi);
        end
    end
    
    nsdDesignFilename = fullfile(nsdfolder, 'nsd_expdesign.mat');
    nsdDesign = load(nsdDesignFilename);
    nsdDesign.masterordering;%for each of 30000 trials, what is corresponding image (out of 10000 images)
    subDesign = nsdDesign.subjectim(isub,nsdDesign.masterordering);%for each of 30000 trials, what is corresponding image (out of 73000 images)
    [imgTrials, imgNum] = ismember(subDesign, allImgs);%logical array
    
    %if less than 40 sessions, only use image trials that were actually presented
    imgTrials = imgTrials(1:size(roiBetas{roinum},2));
    imgNum = imgNum(1:size(roiBetas{roinum},2));
    
    imgTrialSum = cumsum(imgTrials);
    splitImgTrials = repmat(imgTrials,2,1);
    
    midImg = ceil(median(imgNum));
    splitImgTrials(1,imgNum<midImg) = zeros;
    splitImgTrials(2,imgNum>=midImg) = zeros;
    maxNumTrials = max(sum(splitImgTrials,2));
    
    r2 = cell(length(rois),1);
    r2ori = cell(length(rois),1);
    r2oriSplit = cell(length(rois),1);
    r2split = cell(length(rois),1);
    voxOriResidOriR2 = cell(length(rois),1);
    voxOriPredOriR2 = cell(length(rois),1);
    voxResidOriR2 = cell(length(rois),1);
    voxPredOriR2 = cell(length(rois),1);
    
    pearsonRori = cell(length(rois),1);
    pearsonR = cell(length(rois),1);
    
    for roinum=1:length(rois)
        nvox(roinum) = size(roiBetas{roinum},1);
        voxOriResidual{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        voxResidual{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        voxOriResidualSplit{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        voxResidualSplit{roinum} = NaN(nsplits, nvox(roinum),maxNumTrials);
        voxOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        voxCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels+1);
        voxPredOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        voxPredCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels+1);
        voxOriPredOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        voxResidOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        voxOriResidOriCoef{roinum} = zeros(nsplits, nvox(roinum),numLevels*numOrientations+1);
        
        %get model coefficients for each voxel, within each split
        for isplit=1:nsplits
            imgTrials = splitImgTrials(isplit,:);
            numTrials = sum(imgTrials);
            for ivox=1:nvox(roinum)
                voxBetas = roiBetas{roinum}(ivox,imgTrials>0)';
                voxPrfSample = squeeze(prfSampleLev{roinum}(imgNum(imgTrials>0),ivox,:));
                %add constant predictor
                voxPrfSample(:,end+1) = ones;
                voxCoef{roinum}(isplit,ivox,:) = voxPrfSample\voxBetas;
                
                voxPrfOriSample = squeeze(prfSampleLevOri{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
                
                %add constant predictor
                voxPrfOriSample(:,end+1) = ones;
                voxOriCoef{roinum}(isplit,ivox,:) = voxPrfOriSample\voxBetas;%check vox 144 in first ROI
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                %regress vignetting predicted timecourse with orientation model
                voxPred = squeeze(voxCoef{roinum}(isplit,ivox,:))'*voxPrfSample';
                voxPredOriCoef{roinum}(isplit,ivox,:) = regress(voxPred',voxPrfOriSample);
                
                voxOriPred = squeeze(voxOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxOriPredOriCoef{roinum}(isplit,ivox,:) = regress(voxOriPred',voxPrfOriSample);
                
                %compute within-split residuals
                voxOriResidual{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxResidual{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxCoef{roinum}(isplit,ivox,:))'*voxPrfSample';
                
                %regress residuals of vignetting model with full model
                voxResidOriCoef{roinum}(isplit,ivox,:) = regress(squeeze(voxResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:)))),voxPrfOriSample);
                %regress residuals of orientation model with orientation model
                voxOriResidOriCoef{roinum}(isplit,ivox,:) = regress(squeeze(voxOriResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:)))),voxPrfOriSample);
                
                %r2 within split
                r2{roinum}(isplit,ivox) = rsquared(voxResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                r2ori{roinum}(isplit,ivox) = rsquared(voxOriResidual{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                %resid r2 within split
                %residuals of full orientation model:
                oriResidOriResidual = voxBetas' - squeeze(voxOriResidOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxOriResidOriR2{roinum}(isplit,ivox) = rsquared(oriResidOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                %prediction of full orientation model:
                oriPredOriResidual = voxBetas' - squeeze(voxOriPredOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxOriPredOriR2{roinum}(isplit,ivox) = rsquared(oriPredOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                
                %residuals of constrained model:
                residOriResidual = voxBetas' - squeeze(voxResidOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxResidOriR2{roinum}(isplit,ivox) = rsquared(residOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                
                %prediction of constrained model:
                predOriResidual = voxBetas' - squeeze(voxPredOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxPredOriR2{roinum}(isplit,ivox) = rsquared(predOriResidual(1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
            end
        end
        
        for isplit=1:nsplits
            imgTrials = splitImgTrials(isplit,:);
            numTrials = sum(imgTrials);
            for ivox=1:nvox(roinum)
                
                voxBetas = roiBetas{roinum}(ivox,imgTrials>0)';
                
                voxPrfSample = squeeze(prfSampleLev{roinum}(imgNum(imgTrials>0),ivox,:));
                %add constant predictor
                voxPrfSample(:,end+1) = ones;
                voxPrfOriSample = squeeze(prfSampleLevOri{roinum}(imgNum(imgTrials>0),ivox,:,:));
                voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
                
                %add constant predictor
                voxPrfOriSample(:,end+1) = ones;
                
                voxOriResidualSplit{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxOriCoef{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfOriSample';
                voxResidualSplit{roinum}(isplit,ivox,1:numTrials) = voxBetas' - squeeze(voxCoef{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfSample';
                
                %r2 between splits
                r2split{roinum}(isplit,ivox) = rsquared(voxResidualSplit{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                r2oriSplit{roinum}(isplit,ivox) = rsquared(voxOriResidualSplit{roinum}(isplit,ivox,1:sum(splitImgTrials(isplit,:))), roiBetas{roinum}(ivox,imgTrials>0));
                
                %corr between splits
                pearsonRori{roinum}(isplit,ivox) = corr(voxBetas,(squeeze(voxOriCoef{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfOriSample')');
                pearsonR{roinum}(isplit,ivox) = corr(voxBetas,(squeeze(voxCoef{roinum}(nsplits-isplit+1,ivox,:))'*voxPrfSample')');
                
            end
        end
        
        
    end
    
    nsd.voxResidual = voxResidual;
    nsd.voxOriResidual = voxOriResidual;
    nsd.voxResidualSplit = voxResidualSplit;
    nsd.voxOriResidualSplit = voxOriResidualSplit;
    nsd.r2 = r2;
    nsd.r2ori = r2ori;
    nsd.r2split = r2split;
    nsd.r2oriSplit = r2oriSplit;
    nsd.pearsonRori = pearsonRori;
    nsd.pearsonR = pearsonR;
    nsd.imgNum = imgNum;
    nsd.splitImgTrials = splitImgTrials;
    nsd.voxCoef = voxCoef;
    nsd.voxOriCoef = voxOriCoef;
    nsd.voxPredOriCoef = voxPredOriCoef;
    nsd.voxOriPredOriCoef = voxOriPredOriCoef;
    nsd.voxResidOriCoef = voxResidOriCoef;
    nsd.voxOriResidOriCoef = voxOriResidOriCoef;
    nsd.voxPredOriR2 = voxPredOriR2;
    nsd.voxOriPredOriR2 = voxOriPredOriR2;
    nsd.voxResidOriR2 = voxResidOriR2;
    nsd.voxOriResidOriR2 = voxOriResidOriR2;
    nsd.roiInd = roiInd;
    
    
    %% SYNTHETIC STIMULI
    load(fullfile(prffolder,['prfSampleSynth_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']),'prfSampleLevOri','prfSampleLev',...
        'rois','allImgs','numLevels','numOrientations','interpImgSize','backgroundSize','pixPerDeg','imgScaling',...
        'roiPrf');
    if bandpass
        for roinum=1:length(rois)
            prfSampleLevOri{roinum} = prfSampleLevOri{roinum}(:,:,bandMin:bandMax,:);
            prfSampleLev{roinum} = prfSampleLev{roinum}(:,:,bandMin:bandMax);
        end
        numLevels = bandMax-bandMin+1;
    end
    
    allImgs = 105:216;%1:220; ONLY GRATINGS AND SPIRALS
    synthfolder = ['~/NSD/sub' num2str(isub) '_synth_func1pt8mm/'];
    betasfilename = fullfile(synthfolder,'betas_nsdsynthetic.nii');
    betas = niftiread(betasfilename);
    betas = cast(betas,'double');
    betas = betas/300;
    
    synthDesignFilename = fullfile(nsdfolder, 'nsdsynthetic_expdesign.mat');
    synthDesign = load(synthDesignFilename);
    synthDesign.masterordering;
    
    [imgTrials, imgNum] = ismember(synthDesign.masterordering, allImgs);%logical array
    betas=reshape(betas,[],size(betas,4));
    nimgs=length(allImgs);
    imgTrialSum = cumsum(imgTrials);
    totalTrials = max(imgTrialSum);%246
    
    clear roiBetas pearsonRori pearsonR voxOriResidual voxResidual
    for roinum=1:length(rois)
        iroi = rois(roinum);
        roiBetas{roinum} = betas(visRoiData==iroi,:);
    end
    r2 = cell(length(rois),1);
    r2ori = cell(length(rois),1);
    voxOriResidual=cell(length(rois),1);
    voxResidual=cell(length(rois),1);
    voxOriCoef = cell(length(rois),1);
    voxCoef = cell(length(rois),1);
    for roinum=1:length(rois)
        nvox(roinum) = size(roiBetas{roinum},1);
        voxOriResidual{roinum} = zeros(nsplits, nvox(roinum),nimgs);
        voxResidual{roinum} = zeros(nsplits, nvox(roinum),nimgs);
        %compute R^2 for both models, within splits, and between splits
        for isplit=1:nsplits
            
            %average across trials of each image
            clear imgBetas imgPrfSample imgPrfOriSample
            for iimg=1:nimgs
                imgBetas(:,iimg) = mean(roiBetas{roinum}(:,imgNum==iimg),2);
            end
            
            for ivox=1:nvox(roinum)
                voxBetas = imgBetas(ivox,:)';
                voxPrfSample = squeeze(prfSampleLev{roinum}(:,ivox,:));
                voxPrfSample(:,end+1) = ones;
                voxCoef{roinum}(ivox,:) = voxPrfSample\voxBetas;
                
                voxPrfOriSample = squeeze(prfSampleLevOri{roinum}(:,ivox,:,:));
                voxPrfOriSample = reshape(voxPrfOriSample,[],numLevels*numOrientations);
                voxPrfOriSample(:,end+1) = ones;
                voxOriCoef{roinum}(ivox,:) = voxPrfOriSample\voxBetas;%check vox 144 in first ROI
                voxOriResidual{roinum}(isplit,ivox,:) = voxBetas' - squeeze(nsd.voxOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample';
                voxResidual{roinum}(isplit,ivox,:) = voxBetas' - squeeze(nsd.voxCoef{roinum}(isplit,ivox,:))'*voxPrfSample';
                pearsonRori{roinum}(isplit,ivox) = corr(voxBetas,(squeeze(nsd.voxOriCoef{roinum}(isplit,ivox,:))'*voxPrfOriSample')');
                pearsonR{roinum}(isplit,ivox) = corr(voxBetas,(squeeze(nsd.voxCoef{roinum}(isplit,ivox,:))'*voxPrfSample')');
            end
        end
    end
    
    synth.roiBetas = roiBetas;
    synth.prfSampleLevOri = prfSampleLevOri;
    synth.prfSampleLev = prfSampleLev;
    synth.imgTrials = imgTrials;
    synth.imgNum = imgNum;
    synth.pearsonR = pearsonR;
    synth.pearsonRori = pearsonRori;
    synth.voxOriCoef = voxOriCoef;
    synth.voxCoef = voxCoef;
    synth.voxResidual = voxResidual;
    synth.voxOriResidual = voxOriResidual;
    
    %% SAVE RESULTS
    bandpassStr = ['_bandpass' num2str(bandMin) 'to' num2str(bandMax)];
    save(fullfile(prffolder,['regressPrfSplit' bandpassStr '_v' num2str(visualRegion) '_sub' num2str(isub) '.mat']), ...
        'nsd', 'synth',...
        'numLevels', 'numOrientations','rois','nvox','roiPrf','nsplits');
    toc
end

end


function r2 = rsquared(Xresid, Xorig)
r2 = 1 - sum((Xresid).^2)/sum((Xorig - mean(Xorig)).^2);
end