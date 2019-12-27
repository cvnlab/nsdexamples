%% Example 6: Basic loading and inspection of NSD betas

%% Introduction

% In this script, we go through an example of loading and inspecting the NSD betas.
%
% Skills/concepts:
% - Thinking about data units
% - Data visualization
% - Simple ROI definition
% - Non-trivial indexing into the NSD data



%% General setup

% define
stimfile = '~/nsd/nsddata_stimuli/stimuli/nsd/nsd_stimuli.hdf5';
expfile = '~/nsd/nsddata/experiments/nsd/nsd_expdesign.mat';
subjix = 3;                % which NSD subject we are analyzing
nsess = 32;                % how many NSD sessions are available
betaver = 'betas_fithrf';  % which beta version to load



%% Define ROI

% We are going to examine NSD betas for one region of interest (ROI).
% Specifically, right hemisphere (RH) fusiform face area (FFA).
% Here, we go through one approach for defining that region.

% load in t-values for faces vs. nonfaces from the floc experiment
a1 = load_untouch_nii(sprintf('~/nsd/nsddata/ppdata/subj%02d/func1pt8mm/floc_facestval.nii.gz',subjix));

% in ITK-SNAP, look at the NIFTI file and determine a small box in 
% the appropriate anatomical location and containing high t-values (e.g. > 5)
lrix = 58:66;  % from left to right
paix = 31:54;  % from posterior to anterior
isix = 25:32;  % from inferior to superior

% note that this a rough-and-ready approach to defining the ROI. it is
% sufficient for the sake of this example, but a more accurate approach 
% would be to examine the data on the cortical surface.

% create a binary volume with the box
boxvol = zeros(size(a1.img));
boxvol(lrix,paix,isix) = 1;

% define the ROI as within the box AND t > 5
mask = boxvol & a1.img > 5;

% to aid in loading the data, compute a tight fitting box and associated indices
[d1,d2,d3,ii] = computebrickandindices(mask);



%% Load in the betas

% load data
data = [];  % voxels x 750 trials x sessions
for p=1:nsess
  fprintf('sess %d...',p);
  file0 = sprintf('~/nsd/nsddata_betas/ppdata/subj%02d/func1pt8mm/%s/betas_session%02d.mat',subjix,betaver,p);
  a1 = matfile(file0);
  temp = double(a1.betas(d1,d2,d3,:))/300;  % convert to double and then convert to percent signal change
  temp = squish(temp,3);                    % flatten voxels
  temp = temp(ii,:);                        % extract the voxels we want
  data(:,:,p) = temp;                       % record
end



%% Inspect the data

% visualize data in original units (percent signal change)
figure; hold on;
imagesc(reshape(data,size(data,1),[]),[-10 10]);
axis ij tight;
colormap(cmapsign4); colorbar;
xlabel('Trial');
ylabel('Voxel');
title('Response (% BOLD)');
%%

% Notice the high heterogeneity across voxels.
%%

% visualize each voxel's mean PSC
figure; hold on;
bar(mean(mean(data,2),3));
xlabel('Voxel');
ylabel('Average response (% BOLD)');
%%

% This reinforces the point: voxels have highly different overall BOLD responses.
%%

% on a per-voxel basis, z-score the betas obtained in each session.
% this is a fairly drastic measure, but ensures stationarity across time
% and comparable units across voxels.
dataZ = calczscore(data,2);

% visualize data again in z-score units
figure; hold on;
imagesc(reshape(dataZ,size(dataZ,1),[]),[-3 3]);
axis ij tight;
colormap(cmapsign4); colorbar;
xlabel('Trial');
ylabel('Voxel');
title('Response (z-score units)');
%%



%% Load in experiment information

% load
exp1 = load(expfile);
theorder = exp1.masterordering(1:750*nsess);  % the trials that we have data for
uniqueix = union(theorder,[]);                % unique indices into the 10k
length(theorder)  % total number of trials
length(uniqueix)  % total number of unique images
%%



%% Visualize ROI-averaged responses

% define
numtodo = 100;  % number of distinct images to plot responses for

% massage dimensionality
dataZ = reshape(dataZ,size(dataZ,1),[]);  % voxels x trials*sessions

% make the plot
todo = picksubset(uniqueix,numtodo);  % pick a small subset to plot
versions = {'Regular' 'Shuffled'};
for ver=1:2
  figureprep([100 100 1000 300],1); hold on;
  avgresp = [];
  for p=1:length(todo)
    switch ver
    case 1
      ix = find(theorder==todo(p));              % which trials correspond to the image
    case 2
      ix = find(permutedim(theorder==todo(p)));  % SHUFFLE!
    end
    yy = mean(dataZ(:,ix),1);      % compute ROI average for each trial
    scatter(repmat(p,[1 length(yy)]),yy,'ro');
    avgresp(p) = mean(yy);
  end
  h = plot(1:length(todo),avgresp,'ko-');
  set(h,'MarkerFaceColor','k');
  ax = axis;
  axis([ax(1:2) -2 2]);
  xlabel('Image');
  ylabel('Response (z-score units)');
  title(versions{ver});
end
%%

% In the 'Regular' plot, we are looking for small within-image variability
% compared to the across-image variability. The 'Shuffled' plot looks slightly
% different compared to the 'Regular' plot, indicating that the image-evoked 
% signal in the data is fairly weak. This is expected given the noise in fMRI
% data, the fact that an aggressive rapid event-related design was used in NSD 
% (where responses to successive trials overlap substantially), and the small
% number of trials per distinct image that was used in the design.



%% Look at best and worst images

% compute the mean across trials only for those images with all 3 trials
newdata = [];       % voxels x images with the trial-averaged response
newdatastim = [];   % 1 x images with indices into the 10k
for p=1:length(uniqueix)
  ix = find(theorder==uniqueix(p));
  if length(ix)==3
    newdata(:,end+1) = mean(dataZ(:,ix),2);
    newdatastim(end+1) = uniqueix(p);
  end
end

% sort ROI-averaged response in descending order
[ss,ssix] = sort(mean(newdata,1),'descend');

% plot best and worst images
strs = {'Best' 'Worst'};
for flag=1:2
  figureprep([100 100 1000 300],1);
  for p=1:5
    
    % figure out 73k ID
    switch flag
    case 1
      id73k = exp1.subjectim(subjix,newdatastim(ssix(p)));
    case 2
      id73k = exp1.subjectim(subjix,newdatastim(ssix(end-p+1)));
    end

    % get image (425 x 425 x 3, uint8)
    im = permute(h5read(stimfile,'/imgBrick',[1 1 1 id73k],[3 425 425 1]),[3 2 1]);
    
    % plot it
    subplot(1,5,p);
    imshow(im);
    if p==1
      title(strs{flag});
    end

  end
end
%%

% As expected, the images that most strongly drove the response tend to have faces.


