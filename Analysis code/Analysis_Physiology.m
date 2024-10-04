%% Analysis of physiology

% General settings
    data_directory = 'C:\Users\rheerema\OneDrive\Experiment data\MoodChoicePhysiology2024'; %Fill in the directory where the data is stored here
    load('participants.mat') %Load the participants table
    phys = {'EDA','pupil','zygomaticus','corrugator'}; %physiological signals
    samples_EDA = 1:101; %Induction epoch for skin conductance data
    samples_EMG = 151:651; %Induction epoch for facial musculature data
    samples_pupil = 1:601; %Induction epoch for facial musculature data
    EMG_kernel = 20; %Extra smoothing kernel for EMG (400ms)
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
                                    phys_data = cell2mat(AllData.EMG.zygomaticus(i_cond & ~cellfun(@isempty,AllData.EMG.zygomaticus)));
                                    if ~isempty(phys_data)
                                        for k = 1:size(phys_data,1)
                                            if ~all(isnan(phys_data(k,:)))
                                                phys_data(k,:) = ft_preproc_smooth(phys_data(k,:),EMG_kernel); %extra smoothing
                                            end
                                        end
                                        phys_data = phys_data(:,samples_EMG);
                                    end
                                case 'corrugator'
                                    phys_data = cell2mat(AllData.EMG.corrugator(i_cond & ~cellfun(@isempty,AllData.EMG.corrugator)));
                                    if ~isempty(phys_data)
                                        for k = 1:size(phys_data,1)
                                            if ~all(isnan(phys_data(k,:)))
                                                phys_data(k,:) = ft_preproc_smooth(phys_data(k,:),EMG_kernel); %extra smoothing
                                            end
                                        end
                                        phys_data = phys_data(:,samples_EMG);
                                    end
                            end
                        %Get average
                            if isempty(phys_data)
                                physiology_averages_HSN.(phys{j}){ppt,emo} = [];
                            else
                                phys_data = phys_data-phys_data(:,1); %baseline-correction
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
                                        phys_data = cell2mat(AllData.EMG.AFN.zygomaticus(i_cond & ~cellfun(@isempty,AllData.EMG.AFN.zygomaticus)));
                                        if ~isempty(phys_data)
                                            for k = 1:size(phys_data,1)
                                                if ~all(isnan(phys_data(k,:)))
                                                    phys_data(k,:) = ft_preproc_smooth(phys_data(k,:),EMG_kernel); %extra smoothing
                                                end
                                            end
                                            phys_data = phys_data(:,samples_EMG);
                                        end
                                    case 'corrugator'
                                        phys_data = cell2mat(AllData.EMG.AFN.corrugator(i_cond & ~cellfun(@isempty,AllData.EMG.AFN.corrugator)));
                                        if ~isempty(phys_data)
                                            for k = 1:size(phys_data,1)
                                                if ~all(isnan(phys_data(k,:)))
                                                    phys_data(k,:) = ft_preproc_smooth(phys_data(k,:),EMG_kernel); %extra smoothing
                                                end
                                            end
                                            phys_data = phys_data(:,samples_EMG);
                                        end
                                end %switch   
                            %Get average
                                if isempty(phys_data)
                                    physiology_averages_AFN.(phys{j}){ppt,i_emo} = [];
                                else
                                    phys_data = phys_data-phys_data(:,1); %baseline-correction
                                    physiology_averages_AFN.(phys{j}){ppt,i_emo} = nanmean(phys_data);
                                end
                        end %for ii
                    end %for j
                end %if HSN or AFN
            end %for emo
        %Get single-number physiology scores per induction
            %...
    end %for ppt


% Save
    save([cd filesep 'Results\physiology_averages'],'physiology_averages_AFN','physiology_averages_HSN')