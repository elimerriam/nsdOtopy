% makeGratings.m
%
% associated with the following publication: Roth, ZN, Kay, K, and Merriam, EP (2022).
% Massive natural scene sampling reveals reliable coarse-scale orientation tuning in human V1
% DOI:
%
%   usage: makeGratings()
%   by: zvi roth
%   date: 10/25/2023
%   purpose: Create gratings used to measure orientation preference of
%   voxel model fits.
%   gratings saved as gratings.mat are used by getVoxPref.m


interpImgSize = 714;
backgroundSize = 1024;
interpDeg = 8.4;
imgScaling = 0.5;
width = backgroundSize*imgScaling;
pixPerDeg = interpImgSize*imgScaling/interpDeg;
totalDeg = width/pixPerDeg; %including gray background

numFreqs = 30;
numAngles = 32;
minCycles = 1;
maxCycles = backgroundSize*imgScaling/2;
minCpd = minCycles/totalDeg;
maxCpd = maxCycles/totalDeg;
cpds = logspace(log10(minCpd),log10(maxCpd),numFreqs);

angles = linspace(0,180,numAngles+1);
angles = angles(1:numAngles);

% from mglMakeGrating.m
phase = 0;
phase = pi*phase/180;

% get a grid of x and y coordinates that has
% the correct number of pixels
x = -(width-1)/2:(width-1)/2;
y = -(width-1)/2:(width-1)/2;
[xMesh,yMesh] = meshgrid(x,y);

freqs = cpds/pixPerDeg;
allGratings = zeros(length(freqs),length(angles),width,width);

for isf=1:length(freqs)
    sf = freqs(isf);
    for iangle=1:length(angles)
        angle = angles(iangle);
        
        % calculate orientation
        angle = pi*angle/180;
        
        a=cos(angle)*sf*2*pi;
        b=sin(angle)*sf*2*pi;
        % compute grating
        im = cos(a*xMesh+b*yMesh+phase);
        allGratings(isf,iangle,:,:) = im;

        %% pass image through steerable pyramid
        normResp=0;
        % construct quad frequency filters
        numOrientations = 8;
        bandwidth = 1;
        dims = [backgroundSize backgroundSize];
        dims = dims*imgScaling;
        numLevels = maxLevel(dims,bandwidth);
        [freqRespsImag, freqRespsReal, pind] = makeQuadFRs(dims, numLevels, numOrientations, bandwidth);
        
        %%
        [pyr, pind] = buildQuadBands(im, freqRespsImag, freqRespsReal);
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
                
                %sum energy
                sumOriEnergy(isf,iangle,ilev) = sum(sumOri{ilev}(:));
                modelOriEnergy(isf,iangle,ilev,orientation) = sum(thisBand(:));
            end
        end
    end
end

save(['~/NSD/gratings/gratings.mat'],'cpds','angles','freqs','numOrientations','numLevels',...
    'sumOriEnergy','modelOriEnergy','normResp','backgroundSize','imgScaling');
save(['~/NSD/gratings/allGratings.mat'],'cpds','angles','freqs','numOrientations','numLevels',...
    'sumOriEnergy','modelOriEnergy','normResp','backgroundSize','imgScaling','allGratings','-v7.3');