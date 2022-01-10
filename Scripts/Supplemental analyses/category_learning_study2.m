clear all;
clc;

%% import data
filename = ['Category learning_Study 2_cleaned.xlsx'];
rawData = importdata(filename);
subjectIDlist = unique(rawData.data(:,1));

% check number of trials per participant and drop participants < 144
for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index = find(rawData.data(:,1)==subjectID);
    subjectData = rawData.data(index,:);
    numTrials(i_subject) = size(subjectData,1);
end
subjectsToDrop = subjectIDlist(numTrials<144);
subjectIDlist(numTrials<144) = [];
subjectsToDrop_index = ismember(rawData.data(:,1),subjectsToDrop);
rawData.textdata(subjectsToDrop_index,:) = [];
rawData.data(subjectsToDrop_index,:) = [];

%% derive factors from task output
textData = rawData.textdata(2:end,:);
for i_row = 1:size(rawData.data,1)
    procedure = char(textData(i_row,4));
    if strcmp(procedure(1:4),'Face') == 1
        modality(i_row) = 1;
    elseif strcmp(procedure(1:4),'Body') == 1
        modality(i_row) = 2;
    else
        modality(i_row) = 3; % voice
    end
    if strcmp(procedure(end-7:end),'Mismatch') == 1
        match(i_row) = 0;
    else
        match(i_row) = 1; % match
    end
    scenario = char(textData(i_row,5));
    if strcmp(scenario(1:5),'gluck') == 1
        emotion(i_row) = 2;
    elseif strcmp(scenario(1:5),'itosh') == 1
        emotion(i_row) = 3;
    elseif strcmp(scenario(1:5),'greng') == 1
        emotion(i_row) = 4;
    else
        emotion(i_row) = 1; % fear
    end
    trialNumber = rawData.data(i_row,3);
    if trialNumber < 49
        block(i_row) = 1;
    elseif trialNumber > 48 && trialNumber < 97
        block(i_row) = 2;
    else
        block(i_row) = 3;
    end
end

%% sort trials based on emotion category
sortData = [rawData.data transpose(emotion)];
sortData = sortrows(sortData,[1 11 3]);
for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index = find(sortData(:,1)==subjectID);
    subjectData = sortData(index,:);
    trialCount_subject = [1:1:sum(subjectData(:,11)==1) 1:1:sum(subjectData(:,11)==2) 1:1:sum(subjectData(:,11)==3) 1:1:sum(subjectData(:,11)==4)];
    if i_subject == 1
        trialCount = trialCount_subject;
    else
        trialCount = [trialCount trialCount_subject];
    end
    clear trialCount_subject
end
sortData = [sortData transpose(trialCount)];
sortData = sortrows(sortData,[1 3]);

%% compile variables of interest
data =  [rawData.data(:,1:2) transpose(block) sortData(:,end) transpose(match) transpose(emotion) transpose(modality) rawData.data(:,9)];
variables = {'PPID','Condition','Block','TrialCount','Match','Emotion','Modality','Accuracy'};
data_Table = array2table(data,'VariableNames',variables);
%writetable(data_Table,'Study2.xlsx');

%% aggregate accuracy per subject on desired factors/variables
for i_subject = 1:length(subjectIDlist)
    subjectData = [];
    subjectID = subjectIDlist(i_subject);
    index = find(data(:,1)==subjectID);
    subjectData = data(index,:);
    % condition
    condition(i_subject) = nanmean(subjectData(:,2));
    % overall
    overall(i_subject) = nanmean(subjectData(:,8));
    % block
    block1(i_subject) = nanmean(subjectData(subjectData(:,3)==1,8));
    block2(i_subject) = nanmean(subjectData(subjectData(:,3)==2,8));
    block3(i_subject) = nanmean(subjectData(subjectData(:,3)==3,8));
    % trialCount
    trialCount1(i_subject) = nanmean(subjectData(subjectData(:,4)==1,8));
    trialCount36(i_subject) = nanmean(subjectData(subjectData(:,4)==max(subjectData(:,4)),8));
    trialCount1to4(i_subject) = nanmean(subjectData(subjectData(:,4)<=4,8));
    trialCount32to36(i_subject) = nanmean(subjectData(subjectData(:,4)>=32,8));
    % emotion
    fear(i_subject) = nanmean(subjectData(subjectData(:,6)==1,8));
    gluck(i_subject) = nanmean(subjectData(subjectData(:,6)==2,8));
    itosh(i_subject) = nanmean(subjectData(subjectData(:,6)==3,8));
    kreng(i_subject) = nanmean(subjectData(subjectData(:,6)==4,8));
    % modality
    face(i_subject) = nanmean(subjectData(subjectData(:,7)==1,8));
    body(i_subject) = nanmean(subjectData(subjectData(:,7)==2,8));
    voice(i_subject) = nanmean(subjectData(subjectData(:,7)==3,8));
end

%% compile aggregated data
aggData = [transpose(condition) transpose(block1) transpose(block2) transpose(block3)];
variables = {'condition','block1','block2','block3'};
aggData_Table = array2table(aggData,'VariableNames',variables);
aggData_Table.condition = categorical(aggData_Table.condition);

%% run mixed-design ANOVA
blocks = table([1 2 3]','VariableNames',{'blocks'});
rm = fitrm(aggData_Table,'block1-block3~condition','WithinDesign',blocks);
ranovatbl_within = ranova(rm); % within-subjects effect of trial bin
ranovatbl_between = ranova(rm,'WithinModel','blocks'); % between-subjects effect of condition
partialEtaSq_within = table2array(ranovatbl_within(1,1))/(table2array(ranovatbl_within(1,1))+table2array(ranovatbl_within(3,1)));
partialEtaSq_between = table2array(ranovatbl_between(2,1))/(table2array(ranovatbl_between(2,1))+table2array(ranovatbl_between(3,1)));
% pairwise comparisons for block
[~,p_12,~,stats_12] = ttest(aggData(:,2),aggData(:,3));
[~,p_13,~,stats_13] = ttest(aggData(:,2),aggData(:,4));
[~,p_23,~,stats_23] = ttest(aggData(:,3),aggData(:,4));
p_12 = p_12*3; % Bonferroni correction for number of comparisons (3)
p_13 = p_13*3;
p_23 = p_23*3;
% descriptive statistics for block
mBlock1 = mean(aggData(:,2));
sdBlock1 = std(aggData(:,2));
mBlock2 = mean(aggData(:,3));
sdBlock2 = std(aggData(:,3));
mBlock3 = mean(aggData(:,4));
sdBlock3 = std(aggData(:,4));
% descriptive statistics for condition
aggData(:,6) = (aggData(:,2)+aggData(:,3)+aggData(:,4))/3;
aggData_label = aggData(aggData(:,1)==1,:);
aggData_noLabel = aggData(aggData(:,1)==2,:);
mLabel = mean(aggData_label(:,6));
sdLabel = std(aggData_label(:,6));
mNoLabel = mean(aggData_noLabel(:,6));
sdNoLabel = std(aggData_noLabel(:,6));
N = size(aggData,1);
N_label = size(aggData_label,1);
N_noLabel = size(aggData_noLabel,1);

%% generate plot
label_toPlot = mean(aggData_label(:,2:4));
noLabel_toPlot = mean(aggData_noLabel(:,2:4));
label_errorBars = std(aggData_label(:,2:4));
noLabel_errorBars = std(aggData_noLabel(:,2:4));
figure; 
plot_label = errorbar(label_toPlot,label_errorBars,'-b');
hold on;
plot_noLabel = errorbar(noLabel_toPlot,noLabel_errorBars,'-k');
xlim([0.5 3.5]);
xticks([1 2 3]);
xlabel({'trial block'});
ylim([0.5 1]);
ylabel({'proportion correct'});
hold off;
legend('label','no label');

%% generate effect size estimates for meta-analysis
% between-subjects effect of condition, following Goh et al. (2016)
d_between_num = mLabel-mNoLabel;
d_between_denom = sqrt(((N_label-1)*sdLabel^2+(N_noLabel-1)*sdNoLabel^2)/(N-2));
d_between = d_between_num/d_between_denom;
se_between = sqrt(N/(N_label*N_noLabel)+(d_between^2/(2*N)));

% within-subjects (repeated measure) effect of block (1 vs 3)
%r_13 = corr(aggData(:,2),aggData(:,4));
d_within_num = mBlock3-mBlock1;
d_within_denom = stats_13.sd; % following Gibbons et al. (1993); Morris & DeShon (2002)
d_within = d_within_num/d_within_denom;
se_within = sqrt(1/N+(d_within^2/(2*N))); % adapted from Goh et al (2016); Morris & DeShon (2002)

%% descriptive statistics for overall performance
mOverall = mean(overall);
sdOverall = std(overall);

%% descriptive statistics for trial
% first trial
mTrial1 = mean(rawData.data(rawData.data(:,3)==1,9));
sdTrial1 = std(rawData.data(rawData.data(:,3)==1,9));
% last trial
mTrial144 = mean(rawData.data(rawData.data(:,3)==144,9));
sdTrial144 = std(rawData.data(rawData.data(:,3)==144,9));

%% descriptive statistics for trial count
% first trial per emotion
mTrialCount1 = mean(trialCount1);
sdTrialCount1 = std(trialCount1);
% last trial per emotion
mTrialCount36 = mean(trialCount36);
sdTrialCount36 = std(trialCount36);

% first 4 trials per emoion
mTrialCount1to4 = mean(trialCount1to4);
sdTrialCount1to4 = std(trialCount1to4);
% last 4 trials per emotion
mTrialCount32to36 = mean(trialCount32to36);
sdTrialCount32to36 = std(trialCount32to36);

%% descriptive statistics for emotion and modality
mFear = mean(fear);
sdFear = std(fear);
mGluck = mean(gluck);
sdGluck = std(gluck);
mItosh = mean(itosh);
sdItosh = std(itosh);
mKreng = mean(kreng);
sdKreng = std(kreng);

mFace = mean(face);
sdFace = std(face);
mBody = mean(body);
sdBody = std(body);
mVoice = mean(voice);
sdVoice = std(voice);