
% Data
    data = squeeze(physiology_correlations(3:4,5,:))';
    is_percentage = false;
    effect_name = 'R';
    split_studies = false;
% Run test
    if split_studies
        ii_exp = strcmp(participants.experiment,'exploratory');
    else
        ii_exp = true(length(data),1);
    end
    [~,p1,CI1,stats1] = ttest(data(ii_exp,:));
    [~,p2,CI2,stats2] = ttest(data(~ii_exp,:));
    M_expl = nanmean(data(ii_exp,:))';
    M_conf = nanmean(data(~ii_exp,:))';
% Reformat percentages
    if is_percentage
        M_expl = round(100*M_expl,1);
        M_conf = round(100*M_conf,1);
        CI1 = round(100*CI1,1);
        CI2 = round(100*CI2,1);
    end
% Results table
    T = table;
    T.t_expl = round(stats1.tstat',1);
    T.df_expl = stats1.df'; 
    T.([effect_name '_expl']) = M_expl;
    T.p_expl = p1';
    T.CI_expl = CI1';
    T.t_conf = round(stats2.tstat',1);
    T.df_conf = stats2.df';
    T.p_conf = p2';
    T.([effect_name '_conf']) = M_conf;
    T.CI_conf = CI2';
% Display
    clc
    disp(T)
    writetable(T,'results')