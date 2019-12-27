%% Example 9: Some example analyses of the behavioral data

%% Introduction

% In this script, we perform a few simple analyses of the behavioral data
% that were collected during the NSD experiment. These data are rich and 
% extensive, as they reflect up to 30,000 trials from a given subject.
%
% Skills/concepts:
% - Various plotting techniques
% - Bootstrapping



%% Load data

% Load behavioral .tsv file. A few hints:
%  2 = session
%  7 = time
%  9 = iscorrect
% 10 = RT
% 19 = missingdata
thedata = {};
for subjix=1:8
  file0 = sprintf('~/nsd/nsddata/ppdata/subj%02d/behav/responses.tsv',subjix);
  a1 = importdata(file0);
  thedata{subjix} = a1.data;
end
thedata
%%



%% Histogram of RTs

% plot histogram of RTs
bins =  0:50:4200;  % use the same bins for every subject
figureprep([100 100 1000 400],1);
for subjix=1:8
  subplot(2,4,subjix); hold on;
  hist(thedata{subjix}(:,10),bins);
  straightline(0:500:4200,'v','r:');
  xlabel('Reaction time (ms)');
  ylabel('Frequency');
  title(sprintf('Subject %d',subjix));
end
%%



%% RTs as a function of time

% define
subjix = 1;

% Plot RT as a function of time
figure; hold on;
scatter(thedata{subjix}(:,7),thedata{subjix}(:,10),'ro');  % note that NaNs just disappear
xlabel('Number of days');
ylabel('Reaction time (ms)');
%%

% The figure is hard to interpret given the very large
% number of dots. Thus, we will also plot a summary metric.
%%

% Plot the median RT in each session
allsess = unique(thedata{subjix}(:,2));  % all sessions
time0 = [];  % 1 x N with the mean time
md = [];     % 1 x N with the median RT
for p=1:length(allsess)
  ix = find(thedata{subjix}(:,2) == allsess(p));  % trials to consider
  ix = ix(isfinite(thedata{subjix}(ix,10)));      % only consider those with valid data
  time0(p) = mean(thedata{subjix}(ix,7));         % mean time
  md(p) = median(thedata{subjix}(ix,10));         % median RT
end
plot(time0,md,'ko-','LineWidth',2);
%%

% Note that the second black data point is a little funny because it 
% reflects data pooled from two split sessions (frankenstein session).
% Reaction times were fairly stable throughout the experiment.
%%



%% Calculate percent correct in each session, bootstrapping to get reliability

% Here, we will use bootstrapping to estimate the reliability of 
% the percent correct obtained in each scan session. The use of
% bootstrapping is a bit overkill, as we could alternatively just
% use parametric error estimates from the binomial distribution.
% But, it is nonetheless useful to demonstrate how bootstrapping 
% can be implemented.

% do it
numboot = 100;
pctcorrect = NaN*zeros(8,40,numboot);  % initialize with NaNs
for subjix=1:8
  for p=1:40
    ix = thedata{subjix}(:,2)==p;  % find trials
    if sum(ix ~= 0)                % some subjects did not complete all 40
      for boot=1:numboot
      
        % extract data
        subjdata = thedata{subjix}(ix,:);  % trials x columns

        % perform bootstrap sampling (see also bootstrp.m)
        n = size(subjdata,1);              % how many trials are there?
        bootix = ceil(n*rand(1,n));        % generate bootstrap indices
        bootdata = subjdata(bootix,:);     % create the sample

        % calculate percent correct
        isok = bootdata(:,19)==0;  % which rows have valid data?
        pctcorrect(subjix,p,boot) = mean(bootdata(isok,9)==1) * 100;  % NaNs are treated as false here

      end
    end
  end
end

% visualize
figure; hold on;
cmap0 = jet(8);
h1 = []; h2 = [];
for subjix=1:8
  pp0 = prctile(pctcorrect(subjix,:,:),[16 84],3);  % 1 x 40 x 2 (68% confidence interval)
  md0 = median(pctcorrect(subjix,:,:),3);           % 1 x 40
  h1(subjix) = errorbar3(1:40,md0,permute(pp0,[3 2 1]),'v',(cmap0(subjix,:)+[1 1 1])/2);
  h2(subjix) = plot(md0,'-','LineWidth',2,'Color',cmap0(subjix,:));
end
uistack(h2,'top');  % ensure the median lines are on top
xlabel('Session number');
ylabel('Percent correct');
%%


