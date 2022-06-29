%synthStim.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: synthStim()
%   by: zvi roth
%   date: 7/29/2022
%   purpose: run synthetic stimuli through the models, get filter energy responses,
%   and save output
%   creates files used by: prfSampleModel_synth.m
% uses the steerable pyramid: https://github.com/elimerriam/stimulusVignetting

% In brief, the nsdsynthetic experiment involved 8 runs, alternating between fixation runs and memory runs.
% During the fixation runs, the subject performed a fixation task on a central dot;
% during the memory runs, the subject performed a 1-back task on the presented images.
% Each image was presented for 2 s and was followed by a 2-s gap before the next trial.
% There are a total of 284 distinct images.
% The stimuli and trial ordering in the 8 runs was exactly the same for all subjects (including the 1-back trials).
% In each run, there were a total of 93 stimulus trials (as well as blank trials).
% The 93 stimulus trials consisted of 83 regular stimulus trials plus 10 special 1-back trials (in which the presented stimulus was identical to the previously presented stimulus).
% There were a total of 93 x 8 = 744 stimulus trials conducted in the scan session.
% There is a separate .hdf5 file created for each of the NSD subjects (based on color calibration that was tailored to each subject).
% After concatenating the 220 images with the 64 images, the result is 284 images.
% These images are shown on a gray background with RGB value (126,110,108).
% There are 714 rows and 1360 columns.
% The reason for the non-square shape of the image is that the word images extend far to the left and far to the right.

% <masterordering> is 1 x 744 with the sequence of trials (indices relative to the 284 images)
% <stimpattern> is 1 session x 8 runs x 107 trials.
% elements are 0/1 indicating when stimulus trials actually occur.
close all
clear all

interpImgSize = 714;
backgroundSize = 1024;
imgScaling = 0.5;
%%
% construct quad frequency filters
numOrientations = 8;
bandwidth = 1;
dims = [backgroundSize backgroundSize];
dims = dims*imgScaling;
numLevels = maxLevel(dims,bandwidth);
[freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numOrientations, bandwidth);

%%
backgroundColor(1,1,:) = uint8([126,110,108]);
nsdfolder = '~/NSD/';
synthDesignFilename = fullfile(nsdfolder, 'nsdsynthetic_expdesign.mat');
synthDesign = load(synthDesignFilename);
synthDesign.masterordering;

betasfolder = '~/NSD/sub1_synth_func1pt8mm/';
betasfilename = fullfile(betasfolder,'betas_nsdsynthetic.nii');

stimfolder = '~/NSD/syntheticStim/';
stimfilename = fullfile(stimfolder,'nsdsynthetic_stimuli.hdf5');%[3 1360 714 220]
stiminfo = h5info(stimfilename);

imgSizeX = stiminfo.Datasets.Dataspace.Size(2);%1360;
imgSizeY = stiminfo.Datasets.Dataspace.Size(3);%714;
numImgs = stiminfo.Datasets.Dataspace.Size(4);%220



nImgs = 1;

%105-128: radial gratings (spokes)
%129-152: spirals
%153-176: angular gratings (rings)
%177-216: spirals
tic
fixPoint(1,1,:) = uint8([255 0 0]);
for imgNum=1:220
    origImg = h5read(stimfilename,'/imgBrick/',[1 1 1 imgNum],[3 imgSizeX imgSizeY nImgs]);
    origImg = permute(origImg,[3 2 1]);%[425,425,3]
    
    
    %add red semi-transparent fixation point
    origImg(imgSizeY/2-8:imgSizeY/2+8,imgSizeX/2-8:imgSizeX/2+8,:) = ...
        (origImg(imgSizeY/2-8:imgSizeY/2+8,imgSizeX/2-8:imgSizeX/2+8,:) + ...
        repmat(fixPoint,17,17,1))/2;
    
    %%add background
    bigImg = repmat(backgroundColor,backgroundSize,backgroundSize,1);
    
    %     bigImg(1+backgroundSize/2-interpImgSize/2 : backgroundSize/2+interpImgSize/2, 1+backgroundSize/2-interpImgSize/2 : backgroundSize/2+interpImgSize/2,:) = origImg(:,imgSizeX/2-imgSizeY/2+1:imgSizeX/2+imgSizeY/2,:);
    bigImg(1+backgroundSize/2-interpImgSize/2 : backgroundSize/2+interpImgSize/2, :,:) = origImg(:,imgSizeX/2-backgroundSize/2+1:imgSizeX/2+backgroundSize/2,:);
    
    
    % imagesc(mean(bigImg,3));
    imshow(bigImg);
    % imagesc(cast(bigImg,'double'));
    axis image
    
    %change to grayscale
    % for now, simply by averaging across RGB channels
    bigImg = mean(bigImg,3);
    bigImg = single(bigImg);
    
    %DOWNSAMPLE
    bigImg = imresize(bigImg,imgScaling);
    
    %% pass image through steerable pyramid
    pyramidfolder = '~/NSD/syntheticStim/pyramid/';
    pyramidfilename = ['pyrImg' num2str(imgNum) '.mat'];
    [pyr, pind] = buildQuadBands(bigImg, freqRespsImag, freqRespsReal);
    sumOri = cell(numLevels,1);
    modelOri = cell(numLevels,1);
    for ilev = 1:numLevels
        % loop over levels and orientations of the pyramid
        % initialize output
        sumOri{ilev}(:,:) = zeros(dims(1), dims(2),'single');
        modelOri{ilev} = zeros(numOrientations, dims(1), dims(2),'single');
        for orientation = 1:numOrientations
            if normResp
                nEnergies = normEnergies(pyr,pind,numOrientations,0.1);
                thisBand = abs(accessSteerBand(nEnergies,pind,numOrientations,ilev,orientation));
            else
                thisBand = abs(accessSteerBand(pyr, pind, numOrientations,ilev, orientation)).^2;
            end
            sumOri{ilev}(:,:) = sumOri{ilev}(:,:) + thisBand;
            modelOri{ilev}(orientation,:,:) = thisBand;
        end
    end
    save(fullfile(pyramidfolder, pyramidfilename), 'interpImgSize','backgroundSize','numOrientations','imgScaling',...
        'bandwidth','dims','bigImg','sumOri','modelOri','numLevels','normResp');
end
toc

