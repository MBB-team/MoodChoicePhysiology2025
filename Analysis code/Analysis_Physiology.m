%% Analysis of physiology

% General settings
    data_directory = 'C:\Users\Roeland\OneDrive\Experiment data\MoodChoicePhysiology2024'; %Fill in the directory where the data is stored here
    load('participants.mat') %Load the participants table
    phys = {'EDA','pupil','zygomaticus','corrugator'}; %physiological signals
    samples_EDA = 10:110; %Induction epoch for skin conductance data
    samples_EMG = 150:650; %Induction epoch for facial musculature data
    samples_pupil = 1:601; %Induction epoch for facial musculature data
    physiology_averages_HSN = struct; %results for happy/sad/neutral
    physiology_averages_AFN = struct; %results for anger/fear/neutral

% Loop through participants
    for ppt = 1:size(participants,1)
        %Load the data
            if participants.biopac(ppt)
                disp(['PPT #' num2str(ppt)])
                load([data_directory filesep participants.dataset{ppt} filesep 'AllData.mat'])
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
                    switch participants.study(ppt)
                        case 1; s = kron(1:4,ones(1,15));
                        case 2; s = kron(1:3,ones(1,25));
                        case 3; s = kron(1:3,ones(1,20));
                        case 4; s = kron(1:2,ones(1,30));
                    end
                    sessionNumber = s(AllData.affect.induction)'; %Session number within experiment
                    fit_model = fitglm([AllData.affect.induction(i_trials),sessionNumber(i_trials)],phys_data,'CategoricalVars',[false,true]);
                    phys_data = fit_model.Residuals.Raw;
                %Standardize
                    phys_data = nanzscore(phys_data);
                %Store
                    physiology.(phys{j}) =  NaN(length(AllData.affect.condition),1);
                    physiology.(phys{j})(i_trials) = phys_data;
            end %for j
        %Get behavioural measures for correlating
            rated_mood = nanzscore(AllData.affect.RateHappy-AllData.affect.RateSad);
            ratings = nanzscore([AllData.affect.RateHappy,AllData.affect.RateSad,AllData.affect.RateAngry,AllData.affect.RateFear,NaN(size(AllData.affect.RateHappy))]);
            target_rating = NaN(size(AllData.affect.induction));
            choice_rate = NaN(size(AllData.affect.induction));
            for ind = 1:length(choice_rate)
                if ismember(AllData.affect.condition(ind),[1,2,5])
                    choice_rate(ind) = nanmean(AllData.trialinfo.choiceLL(AllData.trialinfo.induction==AllData.affect.induction(ind)));
                end
                target_rating(ind) = ratings(ind,AllData.affect.condition(ind));
            end
        %Get physiological arousal summary measure
%             results(ppt).physArousal = mean([physiology.pupil,physiology.EDA],2);
%             arousal = NaN(size(trialinfo,1),1);
%             for i = 1:length(score)
%                 arousal(trialinfo.induction==Affectdata.induction(i)) = score(i);
%             end
        %Get mood proxy measure from EMG
%             valence = NaN(size(trialinfo,1),1);
%             EMGdelta = EMG_MakeDeltaSignal(Affectdata,[],participant_list(ppt).dataset,participant_list(ppt).study); 
%             for i = 1:length(EMGdelta)
%                 valence(trialinfo.induction==Affectdata.induction(i)) = EMGdelta(i);
%             end
        %Orthogonalize
%             [~,fit] = RH_GLM(arousal,valence);
%             valence = fit.Residuals.Raw;
%         %Only take the sign
%             valence = tanh(valence);
%             results(ppt).arousal = arousal;
%             results(ppt).valence = valence;
        %Core affect signal
%             CoreAffect = valence.*arousal./trialinfo.ind_trialno;
        %Get residuals
%             res = ...
%             betas = RH_GLM(CoreAffect,res);
%             analysis.beta_signal(ppt,:) = betas';
%             results(ppt).CoreAffect = CoreAffect;
        %Bin the signals according to valence
%         n_bins = 9;
%         [sortVal,I] = sort(nanzscore(trialinfo.valence));
%         sortPhysAr = arousal(I);
%         sortEMG = valence(I);
%         n_per_bin = ceil(sum(~isnan(sortVal))/n_bins);
%         if n_per_bin * n_bins > length(sortVal)
%             n_per_bin = floor(sum(~isnan(sortVal))/n_bins);
%         end
%         for bin = 1:n_bins
%             i_bin = (bin-1)*n_per_bin + (1:n_per_bin);
%             analysis.binned_PhysArousal(ppt,bin) = nanmean(sortPhysAr(i_bin));
%             analysis.binned_valence(ppt,bin) = nanmean(sortVal(i_bin));
%             analysis.binned_EMGvalence(ppt,bin) = nanmean(sortEMG(i_bin));
%         end
            
    end %for ppt

% Correlate with beta_mood
    %...

% Save
%     save([cd filesep 'Results\physiology_averages'],'physiology_averages_AFN','physiology_averages_HSN')