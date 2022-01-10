clear all;
clc;

%% import data
filename = ['Category learning_Study 3_cleaned.xlsx'];
rawData = importdata(filename);
subjectIDlist = unique(rawData.data(:,1));

%% derive factors from task output
textData = rawData.textdata(2:end,:);
for i_row = 1:size(textData,1)
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
    elseif strcmp(procedure(end-3:end),'Test') == 1
        match(i_row) = 2; % test block (no feedback)
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
end

%% sort trials into bins based on emotion category
binAssignments = [ones(10,1); 2*ones(10,1); 3*ones(10,1); 4*ones(10,1)];
trialBin = repmat(binAssignments,2*length(subjectIDlist),1);
trialCount = repmat(transpose(1:1:40),2*length(subjectIDlist),1);
sortData = [rawData.data transpose(emotion)];
sortData = sortrows(sortData,[1 11 3]);
sortData = [sortData trialBin trialCount];
sortData = sortrows(sortData,[1 3]);

%% compile variables of interest
data =  [rawData.data(:,1:2) sortData(:,end-1:end) transpose(match) transpose(emotion) transpose(modality) rawData.data(:,9)];
variables = {'PPID','Condition','TrialBin','TrialCount','Match','Emotion','Modality','Accuracy'};
data_Table = array2table(data,'VariableNames',variables);
%writetable(data_Table,'Study3.xlsx');

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
    % trialBin
    bin1(i_subject) = nanmean(subjectData(subjectData(:,3)==1,8));
    bin2(i_subject) = nanmean(subjectData(subjectData(:,3)==2,8));
    bin3(i_subject) = nanmean(subjectData(subjectData(:,3)==3,8));
    bin4(i_subject) = nanmean(subjectData(subjectData(:,3)==4,8));
    % trialCount
    trialCount1(i_subject) = nanmean(subjectData(subjectData(:,4)==1,8));
    trialCount40(i_subject) = nanmean(subjectData(subjectData(:,4)==40,8));
    trialCount1to4(i_subject) = nanmean(subjectData(subjectData(:,4)<=4,8));
    trialCount36to40(i_subject) = nanmean(subjectData(subjectData(:,4)>=36,8));
    % modality
    face(i_subject) = nanmean(subjectData(subjectData(:,7)==1,8));
    body(i_subject) = nanmean(subjectData(subjectData(:,7)==2,8));
    voice(i_subject) = nanmean(subjectData(subjectData(:,7)==3,8));
end

%% compile aggregated data
aggData = [transpose(condition) transpose(bin1) transpose(bin2) transpose(bin3) transpose(bin4)];
variables = {'condition','bin1','bin2','bin3','bin4'};
aggData_Table = array2table(aggData,'VariableNames',variables);
aggData_Table.condition = categorical(aggData_Table.condition);

%% run mixed-design ANOVA
bins = table([1 2 3 4]','VariableNames',{'bins'});
rm = fitrm(aggData_Table,'bin1-bin4~condition','WithinDesign',bins);
ranovatbl_within = ranova(rm); % within-subjects effect of trial bin
ranovatbl_between = ranova(rm,'WithinModel','bins'); % between-subjects effect of condition
partialEtaSq_within = table2array(ranovatbl_within(1,1))/(table2array(ranovatbl_within(1,1))+table2array(ranovatbl_within(3,1)));
partialEtaSq_between = table2array(ranovatbl_between(2,1))/(table2array(ranovatbl_between(2,1))+table2array(ranovatbl_between(3,1)));
% pairwise comparisons for trial bin
[~,p_12,~,stats_12] = ttest(aggData(:,2),aggData(:,3));
[~,p_13,~,stats_13] = ttest(aggData(:,2),aggData(:,4));
[~,p_14,~,stats_14] = ttest(aggData(:,2),aggData(:,5));
[~,p_23,~,stats_23] = ttest(aggData(:,3),aggData(:,4));
[~,p_24,~,stats_24] = ttest(aggData(:,3),aggData(:,5));
[~,p_34,~,stats_34] = ttest(aggData(:,4),aggData(:,5));
p_12 = p_12*6; % Bonferroni correction for number of comparisons (6)
p_13 = p_13*6;
p_14 = p_14*6;
p_23 = p_23*6;
p_24 = p_24*6;
p_34 = p_34*6;
% descriptive statistics for trial bin
mBin1 = mean(aggData(:,2));
sdBin1 = std(aggData(:,2));
mBin2 = mean(aggData(:,3));
sdBin2 = std(aggData(:,3));
mBin3 = mean(aggData(:,4));
sdBin3 = std(aggData(:,4));
mBin4 = mean(aggData(:,5));
sdBin4 = std(aggData(:,5));
% descriptive statistics for condition
aggData(:,6) = (aggData(:,2)+aggData(:,3)+aggData(:,4)+aggData(:,5))/4;
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
label_toPlot = mean(aggData_label(:,2:5));
noLabel_toPlot = mean(aggData_noLabel(:,2:5));
label_errorBars = std(aggData_label(:,2:5));
noLabel_errorBars = std(aggData_noLabel(:,2:5));
figure; 
plot_label = errorbar(label_toPlot,label_errorBars,'-b');
hold on;
plot_noLabel = errorbar(noLabel_toPlot,noLabel_errorBars,'-k');
xlim([0.5 4.5]);
xticks([1 2 3 4]);
xlabel({'trial bin'});
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

% within-subjects (repeated measure) effect of trial bin (1 vs 4)
%r_14 = corr(aggData(:,2),aggData(:,5));
d_within_num = mBin4-mBin1;
d_within_denom = stats_14.sd; % following Gibbons et al. (1993); Morris & DeShon (2002)
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
mTrial80 = mean(rawData.data(rawData.data(:,3)==80,9));
sdTrial80 = std(rawData.data(rawData.data(:,3)==80,9));

%% descriptive statistics for trial count
% first trial per emotion
mTrialCount1 = mean(trialCount1);
sdTrialCount1 = std(trialCount1);
% last trial per emotion
mTrialCount40 = mean(trialCount40);
sdTrialCount40 = std(trialCount40);

% first 4 trials per emoion
mTrialCount1to4 = mean(trialCount1to4);
sdTrialCount1to4 = std(trialCount1to4);
% last 4 trials per emotion
mTrialCount36to40 = mean(trialCount36to40);
sdTrialCount36to40 = std(trialCount36to40);

%% descriptive statistics for emotion and modality
mFace = mean(face);
sdFace = std(face);
mBody = mean(body);
sdBody = std(body);
mVoice = mean(voice);
sdVoice = std(voice);