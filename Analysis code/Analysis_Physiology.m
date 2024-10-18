%% Analysis of physiology

% General settings
    data_directory = 'C:\Users\rheerema\OneDrive\Experiment data\MoodChoicePhysiology2024'; %Fill in the directory where the data is stored here
    load('participants.mat') %Load the participants table
    load('choiceModelBased_noMood.mat') %Load the winning choice model that does not include mood
    phys = {'EDA','pupil','zygomaticus','corrugator'}; %physiological signals
    samples_EDA = 10:110; %Induction epoch for skin conductance data
    samples_EMG = 150:650; %Induction epoch for facial musculature data
    samples_pupil = 1:601; %Induction epoch for facial musculature data
    physiology_averages_HSN = struct; %results for happy/sad/neutral
    physiology_averages_AFN = struct; %results for anger/fear/neutral
    study_sessions = {kron(1:4,ones(1,15)),kron(1:3,ones(1,25)),kron(1:3,ones(1,20)),kron(1:2,ones(1,30))}; %division of experiment into sessions, per study
    include_choicetypes = {[1,3],[2,3],1:3,1:3}; %cost types to include, per study
    physiology_correlations = NaN(7,7,size(participants,1)); %predefine correlation matrix
    R_mood_rated_proxy = NaN(size(participants,1),1); %predefine the correlation coefficient array
    beta_mood_residuals = cell(size(participants,1),3); %predefine regression coefficient matrix
    binned_physiology = struct; %predefine the structure that contains binned data
    n_bins = 9; %number of bins

% Loop through participants
    for ppt = 1:size(participants,1)
        %Load the data
            if participants.biopac(ppt)
                disp(['PPT #' num2str(ppt)])
                load([data_directory filesep participants.dataset{ppt} filesep 'AllData.mat'])
                sessionNumber = study_sessions{participants.study(ppt)};
                sessionNumber = sessionNumber(AllData.affect.induction)'; %Session number within experiment
            else
                continue
            end
        %Get time-resolved averages per induction condition
            for emo = 1:5
                if ismember(emo,[1,2,5]) %happy/sad/neutral
                    for j = 1:length(phys)
                        %Select data
                            i_cond = AllData.affect.condition == emo;
                            switch phys{j}
                                case 'pupil'
                                    phys_data = cell2mat(cellfun(@(x)(x(:,samples_pupil)),...
                                        AllData.pupil.induction(i_cond & ~cellfun(@isempty,AllData.pupil.induction)),'UniformOutput',false));
                                case 'EDA'
                                    phys_data = cell2mat(cellfun(@(x)(x(:,samples_EDA)),...
                                        AllData.EDA.induction(i_cond & ~cellfun(@isempty,AllData.EDA.induction)),'UniformOutput',false));
                                case 'zygomaticus'
                                    phys_data = cell2mat(cellfun(@(x)(x(:,samples_EMG)),...
                                        AllData.EMG.zygomaticus(i_cond & ~cellfun(@isempty,AllData.EMG.zygomaticus)),'UniformOutput',false));
                                case 'corrugator'
                                    phys_data = cell2mat(cellfun(@(x)(x(:,samples_EMG)),...
                                        AllData.EMG.corrugator(i_cond & ~cellfun(@isempty,AllData.EMG.corrugator)),'UniformOutput',false));
                            end
                        %Get average
                            if isempty(phys_data)
                                physiology_averages_HSN.(phys{j}){ppt,emo} = [];
                            else
                                physiology_averages_HSN.(phys{j}){ppt,emo} = nanmean(phys_data);
                            end
                    end
                elseif ismember(emo,[3,4]) && strcmp(participants.experiment(ppt),'exploratory') %anger/fear/neutral in exploratory studies
                    for j = 1:length(phys)    
                        for ii = 1:2
                            %Select data
                                i_cond = AllData.affect.condition == emo; i_emo = emo;
                                if ii==2; i_cond = AllData.affect.condition == 5; i_emo = 5;
                                end
                                switch phys{j}
                                    case 'pupil'
                                        phys_data = cell2mat(cellfun(@(x)(x(:,samples_pupil)),...
                                            AllData.pupil.AFN(i_cond & ~cellfun(@isempty,AllData.pupil.AFN)),'UniformOutput',false));
                                    case 'EDA'
                                        phys_data = cell2mat(cellfun(@(x)(x(:,samples_EDA)),...
                                            AllData.EDA.AFN.induction(i_cond & ~cellfun(@isempty,AllData.EDA.AFN.induction)),'UniformOutput',false));
                                    case 'zygomaticus'
                                        phys_data = cell2mat(cellfun(@(x)(x(:,samples_EMG)),...
                                            AllData.EMG.AFN.zygomaticus(i_cond & ~cellfun(@isempty,AllData.EMG.AFN.zygomaticus)),'UniformOutput',false));
                                    case 'corrugator'
                                        phys_data = cell2mat(cellfun(@(x)(x(:,samples_EMG)),...
                                            AllData.EMG.AFN.corrugator(i_cond & ~cellfun(@isempty,AllData.EMG.AFN.corrugator)),'UniformOutput',false));
                                end %switch   
                            %Get average
                                if isempty(phys_data)
                                    physiology_averages_AFN.(phys{j}){ppt,i_emo} = [];
                                else
                                    physiology_averages_AFN.(phys{j}){ppt,i_emo} = nanmean(phys_data);
                                end
                        end %for ii
                    end %for j
                end %if HSN or AFN
            end %for emo
        %Get single-number physiology scores per induction
            physiology = table;
            for j = 1:length(phys)
                %Select data
                    switch phys{j}
                        case 'pupil'
                            i_trials = ismember(AllData.affect.condition, [1,2,5]) & ~cellfun(@isempty,AllData.pupil.induction);
                            phys_data = cell2mat(cellfun(@(x)(x(:,samples_pupil)), AllData.pupil.induction(i_trials),'UniformOutput',false));
                            phys_data = nanmedian(phys_data,2); %#ok<*NANMEDIAN>
                        case 'EDA'
                            i_trials = ismember(AllData.affect.condition, [1,2,5]);
                            phys_data = AllData.EDA.pspm_amplitudes(i_trials);
                        case 'zygomaticus'
                            i_trials = ismember(AllData.affect.condition, [1,2,5]) & ~cellfun(@isempty,AllData.EMG.zygomaticus);
                            phys_data = cell2mat(cellfun(@(x)(x(:,samples_EMG)),AllData.EMG.zygomaticus(i_trials),'UniformOutput',false));
                            phys_data = nanmedian(phys_data,2); %#ok<*NANMEDIAN>
                        case 'corrugator'
                            i_trials = ismember(AllData.affect.condition, [1,2,5]) & ~cellfun(@isempty,AllData.EMG.corrugator);
                            phys_data = cell2mat(cellfun(@(x)(x(:,samples_EMG)),AllData.EMG.corrugator(i_trials),'UniformOutput',false));
                            phys_data = nanmedian(phys_data,2); %#ok<*NANMEDIAN>
                    end
                %Data correction
                    if any(i_trials)
                        fit_model = fitglm(sessionNumber(i_trials),phys_data,'CategoricalVars',true);
                        phys_data = fit_model.Residuals.Raw; %Correct for trial number and session number
                        phys_data = nanzscore(phys_data); %Standardize
                    end
                %Store
                    physiology.(phys{j}) =  NaN(length(AllData.affect.condition),1);
                    physiology.(phys{j})(i_trials) = phys_data;
            end %for j
        %Get behavioural measures for correlating
            rated_mood = nanzscore(AllData.affect.RateHappy-AllData.affect.RateSad);
            rated_mood(~ismember(AllData.affect.condition,[1,2,5])) = NaN;
            if strcmp(participants.experiment(ppt),'exploratory')
                ratings = [AllData.affect.RateHappy,AllData.affect.RateSad,AllData.affect.RateAngry,AllData.affect.RateFear,NaN(size(AllData.affect.RateHappy))];
            elseif strcmp(participants.experiment(ppt),'confirmatory')
                ratings = [AllData.affect.RateHappy,AllData.affect.RateSad,NaN(length(AllData.affect.RateHappy),3)];
            end
            target_rating = NaN(size(AllData.affect.induction));
            choice_rate = NaN(size(AllData.affect.induction));
            for ind = 1:length(choice_rate)
                if ismember(AllData.affect.condition(ind),[1,2,5])
                    ii_trials = ismember(AllData.trialinfo.choicetype,include_choicetypes{participants.study(ppt)}) & ...
                        AllData.trialinfo.induction==AllData.affect.induction(ind);
                    choice_rate(ind) = nanmean(AllData.trialinfo.choiceLL(ii_trials));
                end
                target_rating(ind) = ratings(ind,AllData.affect.condition(ind));
            end
            correlation_data = [physiology.pupil,physiology.EDA,physiology.zygomaticus,physiology.corrugator,...
                rated_mood,choice_rate,target_rating];
            physiology_correlations(:,:,ppt) = RH_Corr(correlation_data);
        %Get choice trial info
            ii_trials = ismember(AllData.trialinfo.choicetype,include_choicetypes{participants.study(ppt)}) & ...
                ismember(AllData.trialinfo.condition,[1,2,5]);
            trialinfo = AllData.trialinfo(ii_trials,:);
        %Get physiological arousal proxy measure
            proxy_arousal = mean([physiology.pupil,physiology.EDA],2);
            proxy_arousal = RH_Normalize(proxy_arousal);
            trialinfo.proxy_arousal = NaN(size(trialinfo.induction));
            for trl = 1:size(trialinfo,1)
                trialinfo.proxy_arousal(trl) = proxy_arousal(trialinfo.induction(trl));
            end
        %Get physiological valence proxy measure
            proxy_valence = physiology.zygomaticus-physiology.corrugator;
            fit_model = fitglm(proxy_arousal,proxy_valence);
            proxy_valence = fit_model.Residuals.Raw; %orthogonalise valence w.r.t. arousal
            proxy_valence = tanh(proxy_valence); %rescale
            trialinfo.proxy_valence = NaN(size(trialinfo.induction));
            for trl = 1:size(trialinfo,1)
                trialinfo.proxy_valence(trl) = proxy_valence(trialinfo.induction(trl));
            end
        %Bin the valence and arousal proxies based on rated mood
            [sorted_mood,I] = sort(rated_mood);
            sorted_proxy_arousal = proxy_arousal(I);
            sorted_proxy_valence = proxy_valence(I);
            n_per_bin = ceil(sum(~isnan(sorted_proxy_valence))/n_bins);
            if n_per_bin * n_bins > length(sorted_proxy_valence)
                n_per_bin = floor(sum(~isnan(sorted_proxy_valence))/n_bins);
            end
            for bin = 1:n_bins
                i_bin = (bin-1)*n_per_bin + (1:n_per_bin);
                binned_physiology.rated_mood(ppt,bin) = nanmean(sorted_mood(i_bin));
                binned_physiology.proxy_arousal(ppt,bin) = nanmean(sorted_proxy_arousal(i_bin));
                binned_physiology.proxy_valence(ppt,bin) = nanmean(sorted_proxy_valence(i_bin));                
            end
        %Get physiological mood proxy measure and its correlation with the rated mood
            trialinfo.proxy_mood = trialinfo.proxy_valence.*trialinfo.proxy_arousal./trialinfo.ind_trialno;
            if ~all(isnan(trialinfo.proxy_mood))
                R_mood_rated_proxy(ppt) = RH_Corr(trialinfo.proxy_mood,trialinfo.mood);
            end
        %Get choice model residuals and regress against physiological mood proxy
            trialinfo.residuals = winning_model_data.residuals{ppt}';
            choicerate = nanmean(trialinfo.choiceLL);
            if choicerate > 0.05 && choicerate < 0.95 %otherwise regression won't work
                fit_model_rated = fitglm(trialinfo.mood,trialinfo.residuals);
                beta_mood_residuals{ppt,1} = fit_model_rated.Coefficients.Estimate';
                if ~all(isnan(trialinfo.proxy_mood))
                    %Predict residuals with mood proxy
                        fit_model_proxy = fitglm(trialinfo.proxy_mood,trialinfo.residuals);
                        beta_mood_residuals{ppt,2} = fit_model_proxy.Coefficients.Estimate';
                    %Predict residuals with mood proxy and rated mood score
                        fit_model_both = fitglm([trialinfo.mood,trialinfo.proxy_mood],trialinfo.residuals);
                        beta_mood_residuals{ppt,3} = fit_model_both.Coefficients.Estimate';
                else
                    beta_mood_residuals{ppt,2} = NaN(1,2);
                    beta_mood_residuals{ppt,3} = NaN(1,3);
                end
            end
        %Bin the choice model residuals based on the physiological mood proxy
            [sorted_proxy,I] = sort(trialinfo.proxy_mood);
            sorted_residuals = trialinfo.residuals(I);
            n_per_bin = ceil(sum(~isnan(sorted_proxy))/n_bins);
            if n_per_bin * n_bins > length(sorted_proxy)
                n_per_bin = floor(sum(~isnan(sorted_proxy))/n_bins);
            end
            for bin = 1:n_bins
                i_bin = (bin-1)*n_per_bin + (1:n_per_bin);
                binned_physiology.proxy_mood(ppt,bin) = nanmean(sorted_proxy(i_bin));
                binned_physiology.residuals(ppt,bin) = nanmean(sorted_residuals(i_bin));
            end
    end %for ppt

% Correct missing data
    fields = fieldnames(binned_physiology);
    for f = 1:length(fields)
        binned_physiology.(fields{f})(~participants.biopac,:) = NaN;
    end

% Save
    save([cd filesep 'Results\physiology_averages'],'physiology_averages_AFN','physiology_averages_HSN')
    save([cd filesep 'Results\physiology_analysis'],'physiology_correlations','beta_mood_residuals','binned_physiology','R_mood_rated_proxy')
