%% Example 5: Inspect data at scale

%% Introduction

% In this script, we demonstrate two different ways to quickly and effectively
% inspect large amounts of data. The first way is to use a high framerate
% movie, thus exploiting the rapidity of our visual system. The second way is 
% to use an image to visualize many time series at once.
%
% Skills/concepts:
% - Manipulation of images, volumes, and movies
% - Image normalization and color ranges
% - aseg (from FreeSurfer)



%% Make a movie of brain volumes over time

% In this example, we will inspect all of the pre-processed fMRI time-series data
% from the first NSD session from subj01. Knowing how to use low-level routines
% to create the movie is useful as it allows you to customize the movie to your needs.

% define
files = '~/nsd/nsddata_timeseries/ppdata/subj01/func1pt8mm/timeseries/timeseries_session01_*.nii.gz';
outfile = 'testmovie';  % we will add file extensions later
fps = 72;               % frames per second for the movie

% find the files
files0 = matchfiles(files);

% load the data from each run and construct orthographic image inspections
im = uint8([]);
cnt = 1;
for p=1:length(files0)
  
  % load data and convert to double
  data0 = load_untouch_nii(files0{p});
  data0 = double(data0.img);

  % for each volume, make an image
  for r=1:size(data0,4)
  
    % pull out middle slice in all three dimensions
    imA = [];  % 120 x 120 x 3
    imA = cat(3,imA,placematrix(zeros(120,120),             squeeze(data0(:,:,round(end/2),r))));
    imA = cat(3,imA,placematrix(zeros(120,120),rotatematrix(squeeze(data0(:,round(end/2),:,r)),1,2,1)));
    imA = cat(3,imA,placematrix(zeros(120,120),rotatematrix(squeeze(data0(round(end/2),:,:,r)),1,2,1)));
    
    % normalize, make a mosaic, and convert to uint8
    im(:,:,cnt) = uint8(255*makeimagestack(imA,[0 2000],1,[1 3]));
    cnt = cnt + 1;
    
    % the first ten volumes of each run get marked by a little white square
    if r <= 10
      im(1:4,1:4,end) = 255;
    end
    
  end

end

% visualize one image
figure; imshow(im(:,:,1));
%%

% write all of the images to a single .mov file
imagesequencetomovie(im,sprintf('%s.mov',outfile),fps);

% make a compressed version of the movie (mpeg-4 format)
unix_wrapper(sprintf('HandBrakeCLI -i %s.mov -q 2 --strict-anamorphic -o %s.m4v',outfile,outfile));

% note that instead of QTWriter (used by imagesequencetomovie.m), we could perhaps use 
% other utilities such as built-in MATLAB toolbox functions and/or ffmpeg.



%% Use a carpet plot (Power NeuroImage 2017) to quickly inspect time-series data

% In this example, we will use a carpet plot to visualize the same data as above.

% define
files = '~/nsd/nsddata_timeseries/ppdata/subj01/func1pt8mm/timeseries/timeseries_session01_*.nii.gz';
asegfile = '~/nsd/nsddata/ppdata/subj01/func1pt8mm/aseg.nii.gz';

% define more
grayix = [42 3];  % gray matter  (rh is 42, lh is 3)
wmix =   [41 2];  % white matter (rh is 41, lh is 2)
csfix =  24;      % CSF
maxplot = 500;    % maximum number of voxels to plot from each compartment

% use aseg to decide which voxels to extract
aseg = load_untouch_nii(asegfile);
grayvx = picksubset(find(ismember(aseg.img,grayix)),maxplot);
wmvx   = picksubset(find(ismember(aseg.img,wmix)),maxplot);
csfvx  = picksubset(find(ismember(aseg.img,csfix)),maxplot);

% calculate row bounds
rowbounds = [length(grayvx) length(wmvx) length(csfvx)];
rowbounds = cumsum(rowbounds);

% find the time-series files
files0 = matchfiles(files);

% for each file, collect the data
colbounds = [];
alldata = [];
for p=1:length(files0)
  
  % load data, convert to double, flatten voxels
  data0 = load_untouch_nii(files0{p});
  data0 = squish(double(data0.img),3);
  
  % extract data and z-score each voxel
  temp = cat(1,data0(grayvx,:),data0(wmvx,:),data0(csfvx,:));
  temp = calczscore(temp,2);  % there are other ways we could normalize...
  
  % record
  alldata = cat(2,alldata,temp);

  % calculate column bounds
  colbounds = [colbounds size(data0,2)];

end
colbounds = cumsum(colbounds);

% make the visualization
figureprep([100 100 800 600],1); hold on;
imagesc(alldata,[-3 3]);
colormap(gray);
axis ij;
axis([.5 size(alldata,2)+.5 .5 size(alldata,1)+.5]);
straightline(rowbounds+.5,'h','r-');
straightline(colbounds+.5,'v','r-');
xlabel('Volume number');
ylabel('Voxel');
%%

% In the plot, the three row compartments are gray matter, white matter, and CSF,
% and the column compartments correspond to different fMRI runs.


