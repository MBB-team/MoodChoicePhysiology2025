%% Analysis of behaviour: emotion ratings and decision-making

% General settings
    data_directory = 'C:\Users\rheerema\OneDrive\Experiment data\MoodChoicePhysiology2024'; %Fill in the directory where the data is stored here
    load('participants.mat') %Load the participants table
    allRatings = struct; %rating results
    choice_mdlfree = struct; %choice behaviour results
    
% Loop through participants
    for ppt = 1:size(participants,1)
        %Load the data
            disp(['PPT #' num2str(ppt)])
            load([data_directory filesep participants.dataset{ppt} filesep 'AllData.mat'])
            
        %% Ratings analysis
            %Settings
                which_emotions = 1:5;
                which_ratings = {'RateHappy','RateSad','RateAngry','RateFear','Mood'};
            %Compute mood score
                AllData.affect.Mood = AllData.affect.RateHappy-AllData.affect.RateSad;
                AllData.affect.Mood(ismember(AllData.affect.Condition,{'anger','fear'})) = NaN;
                AllData.affect.Mood = nanzscore(AllData.affect.Mood);
            %Gater results in "allRatings" structure
                for emo = 1:length(which_emotions)
                    for R = 1:length(which_ratings)
                        select = AllData.affect.condition==which_emotions(emo);
                        try
                            allRatings.(which_ratings{R})(ppt,emo) = nanmean(AllData.affect.(which_ratings{R})(select));
                        catch
                            allRatings.(which_ratings{R})(ppt,emo) = NaN;
                        end
                    end
                end
            %Mean and SEM per study (for visualization)
                if ppt == size(participants,1)
                    include_studies = [0,1]; %[exploratory, confirmatory]
                    confirmatory = double(participants.study>2); %binary: confirmatory
                    for R = 1:length(which_ratings)
                        rating_data = allRatings.(which_ratings{R});
                        study_data = cell(2,3); %Three conditions (top row: means, bottom row: SEM)
                        for cond = which_emotions
                            cond_means = NaN(1,2); cond_SEM = cond_means;
                            for study = 1:2
                                cond_means(study) = nanmean(rating_data(confirmatory==include_studies(study),cond));
                                cond_SEM(study) = nanstd(rating_data(confirmatory==include_studies(study),cond))/...
                                    sqrt(sum(confirmatory==include_studies(study)));
                            end
                            study_data{1,cond} = cond_means;
                            study_data{2,cond} = cond_SEM;
                        end
                        allRatings.perStudy.(which_ratings{R}) = study_data;
                    end
                end
            
        %% Choice data
            %Exclude choice types from datasets that are dissimilar to the same type in other studies
                include_choicetypes = true(1,3);   
                for type = 1:3
                    if participants.study(ppt) == 1 && type == 2 %Excude probability discounting from study 1 
                        include_choicetypes(type) = 0;
                    elseif participants.study(ppt) == 2 && type == 1 %Exclude delay discounting from study 2
                        include_choicetypes(type) = 0;
                    end
                end
            %Choice rate per emotion
                for emo = 1:5
                    select = ismember(AllData.trialinfo.choicetype,find(include_choicetypes)) & AllData.trialinfo.condition == emo;
                    choice_mdlfree.choiceRate(ppt,emo) = nanmean(AllData.trialinfo.choiceLL(select));
                end
            %Logistic regression against mood                
                for type = 1:3 %delay, risk, effort
                    %inclusion criteria
                        i_select = AllData.trialinfo.choicetype==type;
                        type_choicerate = nanmean(AllData.trialinfo.choiceLL(i_select)); %can't regress if insufficient variability in choice behaviour
                        if ~ismember(type,find(include_choicetypes)) || type_choicerate < 0.05 || type_choicerate > 0.95
                            choice_mdlfree.beta_mood(ppt,type) = NaN;
                            continue
                        end
                    %regress
                        fit_model = fitglm(AllData.trialinfo.mood(i_select),AllData.trialinfo.choiceLL(i_select),'Distribution','binomial');
                        choice_mdlfree.beta_mood(ppt,type) = fit_model.Coefficients.Estimate(2);
                end
            %Response time
                %Get data (only happy/sad/neutral)
                    select = ismember(AllData.trialinfo.choicetype,find(include_choicetypes)) & ismember(AllData.trialinfo.condition,[1,2,5]);
                    RT = AllData.trialinfo.RT; %Raw response time in seconds
                    RT(~select) = NaN;
                %Preprocessing
                    RT_upperlimit = 10; %For trimming: upper limit RT
                    RT_lowerlimit = 0.5; %For trimming: lower limit RT (already enforced experimentally in studies 3 and 4)
                    n_std_exclude = 3; %Remove outliers: +/- n standard deviations from the mean
                    RT(RT > RT_upperlimit) = NaN; %Trim out choices more than 10s
                    RT(RT < RT_lowerlimit) = NaN; %Trim out choices less than 0.5s
                    RT(RT > nanmean(RT)+n_std_exclude*nanstd(RT) | RT < nanmean(RT)-n_std_exclude*nanstd(RT)) = NaN; %Choices more or less than 3 standard deviations above or below the mean
                    fit_model = fitglm(AllData.trialinfo.trial,RT); %Regress out trial number
                    RT = fit_model.Residuals.Raw;
                    RT = nanzscore(RT); %Standardize
                %Get average per emotion X chosen option
                    choice_mdlfree.RT_perEmo_SS(ppt,:) = [nanmean(RT(AllData.trialinfo.condition==1 & AllData.trialinfo.choiceLL==0)),...
                        nanmean(RT(AllData.trialinfo.condition==2 & AllData.trialinfo.choiceLL==0)),...
                        nanmean(RT(AllData.trialinfo.condition==5 & AllData.trialinfo.choiceLL==0))];
                    choice_mdlfree.RT_perEmo_LL(ppt,:) = [nanmean(RT(AllData.trialinfo.condition==1 & AllData.trialinfo.choiceLL==1)),...
                        nanmean(RT(AllData.trialinfo.condition==2 & AllData.trialinfo.choiceLL==1)),...
                        nanmean(RT(AllData.trialinfo.condition==5 & AllData.trialinfo.choiceLL==1))];
    end %for ppt
    
% Store
    save([cd filesep 'Results' filesep 'InductionRatings'],'allRatings')
    save([cd filesep 'Results' filesep 'choiceModelFree'],'choice_mdlfree')
