clear all;
clc;

%% import data
filename = ['Category learning_Study 5_cleaned.xlsx'];
rawData = importdata(filename);
subjectIDlist = unique(rawData.data(:,1));

%% derive factors from task output
textData = rawData.textdata(2:end,:);
for i_row = 1:size(textData,1)
    correctKey = char(textData(i_row,6));
    correctResponse = char(textData(i_row,10));
    if strcmp(correctKey,'A') == 1 && strcmp(correctResponse,'a') == 1
        match(i_row) = 1;
    elseif strcmp(correctKey,'K') == 1 && strcmp(correctResponse,'k') == 1
        match(i_row) = 1;
    else
        match(i_row) = 0; % mismatch
    end
    modality(i_row) = 1; % face
end

%% sort trials into bins based on emotion category
binAssignments = [ones(16,1); 2*ones(16,1); 3*ones(16,1); 4*ones(16,1); 5*ones(16,1); 6*ones(16,1)];
trialBin = repmat(binAssignments,2*length(subjectIDlist),1);
trialCount = repmat(transpose(1:1:96),2*length(subjectIDlist),1);
sortData = rawData.data;
sortData = sortrows(sortData,[1 9 3]);
sortData = [sortData trialBin trialCount];
sortData = sortrows(sortData,[1 3]);

%% compile variables of interest
data =  [rawData.data(:,1:2) sortData(:,end-1:end) transpose(match) rawData.data(:,9) transpose(modality) rawData.data(:,12)];
variables = {'PPID','Condition','TrialBin','TrialCount','Match','Emotion','Modality','Accuracy'};
data_Table = array2table(data,'VariableNames',variables);
writetable(data_Table,'Study5.xlsx');

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
    bin5(i_subject) = nanmean(subjectData(subjectData(:,3)==5,8));
    bin6(i_subject) = nanmean(subjectData(subjectData(:,3)==6,8));
    % trialCount
    trialCount1(i_subject) = nanmean(subjectData(subjectData(:,4)==1,8));
    trialCount96(i_subject) = nanmean(subjectData(subjectData(:,4)==96,8));
    trialCount1to4(i_subject) = nanmean(subjectData(subjectData(:,4)<=4,8));
    trialCount92to96(i_subject) = nanmean(subjectData(subjectData(:,4)>=92,8));
end

%% compile aggregated data
aggData = [transpose(condition) transpose(bin1) transpose(bin2) transpose(bin3) transpose(bin4) transpose(bin5) transpose(bin6)];
variables = {'condition','bin1','bin2','bin3','bin4','bin5','bin6'};
aggData_Table = array2table(aggData,'VariableNames',variables);
aggData_Table.condition = categorical(aggData_Table.condition);

%% run mixed-design ANOVA
bins = table([1 2 3 4 5 6]','VariableNames',{'bins'});
rm = fitrm(aggData_Table,'bin1-bin6~condition','WithinDesign',bins);
ranovatbl_within = ranova(rm); % within-subjects effect of trial bin
ranovatbl_between = ranova(rm,'WithinModel','bins'); % between-subjects effect of condition
partialEtaSq_within = table2array(ranovatbl_within(1,1))/(table2array(ranovatbl_within(1,1))+table2array(ranovatbl_within(3,1)));
partialEtaSq_between = table2array(ranovatbl_between(2,1))/(table2array(ranovatbl_between(2,1))+table2array(ranovatbl_between(3,1)));
% pairwise comparisons for trial bin
[~,p_12,~,stats_12] = ttest(aggData(:,2),aggData(:,3));
[~,p_13,~,stats_13] = ttest(aggData(:,2),aggData(:,4));
[~,p_14,~,stats_14] = ttest(aggData(:,2),aggData(:,5));
[~,p_15,~,stats_15] = ttest(aggData(:,2),aggData(:,6));
[~,p_16,~,stats_16] = ttest(aggData(:,2),aggData(:,7));
[~,p_23,~,stats_23] = ttest(aggData(:,3),aggData(:,4));
[~,p_24,~,stats_24] = ttest(aggData(:,3),aggData(:,5));
[~,p_25,~,stats_25] = ttest(aggData(:,3),aggData(:,6));
[~,p_26,~,stats_26] = ttest(aggData(:,3),aggData(:,7));
[~,p_34,~,stats_34] = ttest(aggData(:,4),aggData(:,5));
[~,p_35,~,stats_35] = ttest(aggData(:,4),aggData(:,6));
[~,p_36,~,stats_36] = ttest(aggData(:,4),aggData(:,7));
[~,p_45,~,stats_45] = ttest(aggData(:,5),aggData(:,6));
[~,p_46,~,stats_46] = ttest(aggData(:,5),aggData(:,7));
[~,p_56,~,stats_56] = ttest(aggData(:,6),aggData(:,7));
p_12 = p_12*15; % Bonferroni correction for number of comparisons (15)
p_13 = p_13*15;
p_14 = p_14*15;
p_15 = p_15*15;
p_16 = p_16*15;
p_23 = p_23*15;
p_24 = p_24*15;
p_25 = p_25*15;
p_26 = p_26*15;
p_34 = p_34*15;
p_35 = p_35*15;
p_36 = p_36*15;
p_45 = p_45*15;
p_46 = p_46*15;
p_56 = p_56*15;
% descriptive statistics for trial bin
mBin1 = mean(aggData(:,2));
sdBin1 = std(aggData(:,2));
mBin2 = mean(aggData(:,3));
sdBin2 = std(aggData(:,3));
mBin3 = mean(aggData(:,4));
sdBin3 = std(aggData(:,4));
mBin4 = mean(aggData(:,5));
sdBin4 = std(aggData(:,5));
mBin5 = mean(aggData(:,6));
sdBin5 = std(aggData(:,6));
mBin6 = mean(aggData(:,7));
sdBin6 = std(aggData(:,7));
% descriptive statistics for condition
aggData(:,8) = (aggData(:,2)+aggData(:,3)+aggData(:,4)+aggData(:,5)+aggData(:,6)+aggData(:,7))/6;
aggData_label = aggData(aggData(:,1)==1,:);
aggData_noLabel = aggData(aggData(:,1)==2,:);
mLabel = mean(aggData_label(:,8));
sdLabel = std(aggData_label(:,8));
mNoLabel = mean(aggData_noLabel(:,8));
sdNoLabel = std(aggData_noLabel(:,8));
N = size(aggData,1);
N_label = size(aggData_label,1);
N_noLabel = size(aggData_noLabel,1);

%% generate plot
label_toPlot = mean(aggData_label(:,2:7));
noLabel_toPlot = mean(aggData_noLabel(:,2:7));
label_errorBars = std(aggData_label(:,2:7));
noLabel_errorBars = std(aggData_noLabel(:,2:7));
figure; 
plot_label = errorbar(label_toPlot,label_errorBars,'-b');
hold on;
plot_noLabel = errorbar(noLabel_toPlot,noLabel_errorBars,'-k');
xlim([0.5 6.5]);
xticks([1 2 3 4 5 6]);
xlabel({'trial bin'});
ylim([0.5 1]);
ylabel({'proportion correct'});
hold off;
legend('label','no label','Location','southeast');

%% generate effect size estimates for meta-analysis
% between-subjects effect of condition, following Goh et al. (2016)
d_between_num = mLabel-mNoLabel;
d_between_denom = sqrt(((N_label-1)*sdLabel^2+(N_noLabel-1)*sdNoLabel^2)/(N-2));
d_between = d_between_num/d_between_denom;
se_between = sqrt(N/(N_label*N_noLabel)+(d_between^2/(2*N)));

% within-subjects (repeated measure) effect of trial bin (1 vs 6)
d_within_num = mBin6-mBin1;
d_within_denom = stats_16.sd; % following Gibbons et al. (1993); Morris & DeShon (2002)
%r_16 = corr(aggData(:,2),aggData(:,7));
%d_within_denom = sqrt((sdBin1^2+sdBin6^2)-(2*r_16*sdBin1*sdBin6));
d_within = d_within_num/d_within_denom;
se_within = sqrt(1/N+(d_within^2/(2*N))); % adapted from Goh et al (2016); Morris & DeShon (2002)

% standard errors based on Morris & DeShon's (2002) formula for variance
% n_between = (N_label*N_noLabel)/N;
% df_between = N-2;
% cdf_between = 1-3/(4*df_between-1);
% var_between = (1/n_between)*(df_between/(df_between-2))*(1+n_between*d_between^2)-(d_between^2/cdf_between^2);
% se_between = sqrt(var_between/N)

% df_within = N-1;
% cdf_within = 1-3/(4*df_within-1);
% var_within = (1/N)*(df_within/(df_within-2))*(1+N*d_within^2)-(d_within^2/cdf_within^2);
% se_within = sqrt(var_within/N)

% this website is helpful for general checking: https://www.psychometrica.de/effect_size.html
% 95%CI = X +/- 1.96*SE

%% descriptive statistics for overall performance
mOverall = mean(overall);
sdOverall = std(overall);

%% descriptive statistics for trial
% first trial
mTrial1 = mean(rawData.data(rawData.data(:,3)==1,12));
sdTrial1 = std(rawData.data(rawData.data(:,3)==1,12));
% last trial
mTrial192 = mean(rawData.data(rawData.data(:,3)==192,12));
sdTrial192 = std(rawData.data(rawData.data(:,3)==192,12));

%% descriptive statistics for trial count
% first trial per emotion
mTrialCount1 = mean(trialCount1);
sdTrialCount1 = std(trialCount1);
% last trial per emotion
mTrialCount96 = mean(trialCount96);
sdTrialCount96 = std(trialCount96);

% first 4 trials per emoion
mTrialCount1to4 = mean(trialCount1to4);
sdTrialCount1to4 = std(trialCount1to4);
% last 4 trials per emotion
mTrialCount92to96 = mean(trialCount92to96);
sdTrialCount92to96 = std(trialCount92to96);