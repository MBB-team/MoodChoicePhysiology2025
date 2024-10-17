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
                %All ratings per condition
                    for emo = 1:length(which_emotions)
                        select = AllData.affect.condition==which_emotions(emo);
                        for R = 1:length(which_ratings)
                            try
                                allRatings.(which_ratings{R})(ppt,emo) = nanmean(AllData.affect.(which_ratings{R})(select));
                            catch
                                allRatings.(which_ratings{R})(ppt,emo) = NaN;
                            end
                        end
                    end
                %Target rating 
                    if strcmp(participants.experiment{ppt},'exploratory')
                        ratings = [AllData.affect.RateHappy,AllData.affect.RateSad,AllData.affect.RateAngry,AllData.affect.RateFear];
                        i_nontarget = setdiff(1:4,emo);
                    elseif strcmp(participants.experiment{ppt},'confirmatory')
                        ratings = [AllData.affect.RateHappy,AllData.affect.RateSad];
                        i_nontarget = setdiff(1:2,emo);
                    end
                    %Difference w.r.t. mean of other ratings
                        for emo = 1:4 %no target rating for neutral
                            try
                                allRatings.deltaTarget(ppt,emo) = nanmean(ratings(AllData.affect.condition==emo,emo) - nanmean(ratings(AllData.affect.condition==emo,i_nontarget),2));
                            catch %no anger/fear ratings in confirmatory studies
                                allRatings.deltaTarget(ppt,emo) = NaN;
                            end
                        end
                    %Target rating decline
                        targetRating = NaN(size(ratings,1),1);
                        for emo = 1:size(ratings,2) 
                            targetRating(AllData.affect.condition==emo) = ratings(AllData.affect.condition==emo,emo);
                        end
                        fit_model = fitglm(AllData.affect.induction,targetRating);
                        allRatings.beta_target_trl(ppt,1) = fit_model.Coefficients.Estimate(2);
                    %Regression against same rating on subsequent induction
                        ratings = nanzscore(ratings); %must standardize, otherwise correlations will always be positive
                        targetRating = NaN(size(ratings,1),1); nextRating = targetRating;
                        for trl = 1:length(targetRating)-1
                            if AllData.affect.condition(trl) < 5
                                targetRating(trl) = ratings(trl,AllData.affect.condition(trl));
                                nextRating(trl) = ratings(trl+1,AllData.affect.condition(trl));
                            end
                        end
                        fit_model = fitglm(targetRating,nextRating);
                        allRatings.beta_target_next(ppt,1) = fit_model.Coefficients.Estimate(2);
                %Correlation between rated happiness and sadness
                    allRatings.R_happySad(ppt,1) = RH_Corr(AllData.affect.RateHappy,AllData.affect.RateSad);
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
            %Choice rate
                %Per emotion
                    for emo = 1:5
                        select = ismember(AllData.trialinfo.choicetype,find(include_choicetypes)) & AllData.trialinfo.condition == emo;
                        choice_mdlfree.choiceRate(ppt,emo) = nanmean(AllData.trialinfo.choiceLL(select));
                    end
                %Regressed against trial number
                    fit_model = fitglm(AllData.trialinfo.trial,AllData.trialinfo.choiceLL,'Distribution','binomial');
                    choice_mdlfree.beta_trialno(ppt,1) = fit_model.Coefficients.Estimate(2);
            %Logistic regression of choices against ratings               
                %Per cost type
                    for type = 1:3 %delay, risk, effort
                        %inclusion criteria
                            i_select = AllData.trialinfo.choicetype==type;
                            type_choicerate = nanmean(AllData.trialinfo.choiceLL(i_select)); %can't regress if insufficient variability in choice behaviour
                            if ~ismember(type,find(include_choicetypes)) || type_choicerate < 0.05 || type_choicerate > 0.95
                                choice_mdlfree.beta_mood(ppt,type) = NaN;
                                choice_mdlfree.beta_RateAngry(ppt,type) = NaN;
                                choice_mdlfree.beta_RateFear(ppt,type) = NaN;
                                continue
                            end
                        %regress against mood (happy/sad/neutral inductions only)
                            fit_model = fitglm(AllData.trialinfo.mood(i_select),AllData.trialinfo.choiceLL(i_select),'Distribution','binomial');
                            choice_mdlfree.beta_mood(ppt,type) = fit_model.Coefficients.Estimate(2);
                        %regress against anger and fear ratings
                            if strcmp(participants.experiment(ppt),'exploratory')
                                ratingData_choice = NaN(size(AllData.trialinfo,1),2);
                                for emo = 1:2
                                    if emo == 1; rating = 'RateAngry'; else; rating = 'RateFear'; end
                                    ratingData_induction = nanzscore(AllData.affect.(rating));
                                    for trl = 1:length(ratingData_choice)
                                        ratingData_choice(trl,emo) = ratingData_induction(AllData.trialinfo.induction(trl));
                                    end
                                    fit_model = fitglm(ratingData_choice(i_select,emo),AllData.trialinfo.choiceLL(i_select),'Distribution','binomial');
                                    choice_mdlfree.(['beta_' rating])(ppt,type) = fit_model.Coefficients.Estimate(2);
                                end
                            else
                                choice_mdlfree.beta_RateAngry(ppt,type) = NaN;
                                choice_mdlfree.beta_RateFear(ppt,type) = NaN;
                            end
                    end %for cost type
                %Across cost types
                    i_select = ismember(AllData.trialinfo.choicetype,find(include_choicetypes));
                    if nanmean(AllData.trialinfo.choiceLL(i_select)) >= 0.05 && nanmean(AllData.trialinfo.choiceLL(i_select)) <= 0.95
                        %Regress against mood
                            fit_model = fitglm(AllData.trialinfo.mood(i_select),AllData.trialinfo.choiceLL(i_select),'Distribution','binomial');
                            choice_mdlfree.beta_mood(ppt,4) = fit_model.Coefficients.Estimate(2);
                        %Regress against anger and fear
                            if strcmp(participants.experiment(ppt),'exploratory')
                                fit_model = fitglm(ratingData_choice(i_select,1),AllData.trialinfo.choiceLL(i_select),'Distribution','binomial');
                                choice_mdlfree.beta_RateAngry(ppt,4) = fit_model.Coefficients.Estimate(2);
                                fit_model = fitglm(ratingData_choice(i_select,2),AllData.trialinfo.choiceLL(i_select),'Distribution','binomial');
                                choice_mdlfree.beta_RateFear(ppt,4) = fit_model.Coefficients.Estimate(2);
                            else
                                choice_mdlfree.beta_RateAngry(ppt,4) = NaN;
                                choice_mdlfree.beta_RateFear(ppt,4) = NaN;
                            end
                    else
                        choice_mdlfree.beta_mood(ppt,4) = NaN;
                        choice_mdlfree.beta_RateAngry(ppt,4) = NaN;
                        choice_mdlfree.beta_RateFear(ppt,4) = NaN;
                    end
            %Response time
                %Get data (only happy/sad/neutral)
                    select = ismember(AllData.trialinfo.choicetype,find(include_choicetypes)) & ismember(AllData.trialinfo.condition,[1,2,5]);
                    RT = AllData.trialinfo.RT; %Raw response time in seconds
                    RT(~select) = NaN;
                %Preprocessing
                    RT_upperlimit = 10; %For trimming: upper limit RT
                    RT_lowerlimit = 0.75; %For trimming: lower limit RT (already enforced experimentally in studies 3 and 4)
                    n_std_exclude = 3; %Remove outliers: +/- n standard deviations from the mean
                    RT(RT > RT_upperlimit) = NaN; %Trim out choices more than 10s
                    RT(RT < RT_lowerlimit) = NaN; %Trim out choices less than 0.75s
                    RT(RT > nanmean(RT)+n_std_exclude*nanstd(RT) | RT < nanmean(RT)-n_std_exclude*nanstd(RT)) = NaN; %#ok<*NANSTD> %Choices more or less than 3 standard deviations above or below the mean
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
                %Regress against chosen option and mood
                    fit_model = fitglm([AllData.trialinfo.mood,1-AllData.trialinfo.choiceLL],RT,'interactions');
                    choice_mdlfree.beta_RT(ppt,:) = fit_model.Coefficients.Estimate(2:end)';
    end %for ppt
    
% Store
    save([cd filesep 'Results' filesep 'InductionRatings'],'allRatings')
    save([cd filesep 'Results' filesep 'choiceModelFree'],'choice_mdlfree')
