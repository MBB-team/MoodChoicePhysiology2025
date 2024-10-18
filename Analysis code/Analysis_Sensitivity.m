% Sensitivity analysis

% Get beta_mood parameter
    load('choiceModelBased_withMood.mat')
    i_winning_model = 1;
    beta_mood = [Inversionresults(i_winning_model).muPhi.betaMood]'; %The effect of interest

% Sample size vs. minimum detectable effect
    mu = mean(beta_mood); %mean of the effect of interest
    sigma = std(beta_mood); %standard deviation of effect of interest
    Cohen_d = mu/sigma; %Cohen's D
    N = length(beta_mood); %current sample size
    nn = 2:150; %Sample as an independent variable
    MDE_80 = sampsizepwr('t',[mu, sigma],[],0.8,nn) - mu; %Minimal detectable effect at 80% power
    MDE_90 = sampsizepwr('t',[mu, sigma],[],0.9,nn) - mu; %Minimal detectable effect at 90% power
    N_min_80 = nn(find(MDE_80<mu,1,'first')); %Minimal sample to find current effect size at 80% power
    N_min_90 = nn(find(MDE_90<mu,1,'first')); %Minimal sample to find current effect size at 90% power
    
% Visualize
    figure; 
    subplot(1,3,1); hold on; box on
    colors = [0,0.447,0.741;0.850,0.325,0.0980;0.929,0.694,0.125;0.494,0.184,0.556];
    plot(nn,MDE_90,'linewidth',2,'color',[0.1 0.1 0.1])
    plot(nn,MDE_80,'linewidth',2,'color',[0.6 0.6 0.6])
    scatter(N,MDE_90(nn==N),60,colors(3,:),'filled')
    scatter(N,MDE_80(nn==N),60,colors(4,:),'filled')
    plot(0:150,ones(1,151)*mu,'--','color',colors(2,:))
    plot(N_min_90*ones(1,2),[0,mu],':','color',colors(3,:))
    plot(N_min_80*ones(1,2),[0,mu],':','color',colors(4,:))
    xlabel('Sample Size'); xticks(0:15:150); ylabel('Minimal detectable effect'); ylim([0 0.3])
    legend({'','',['MDE = ' num2str(round(MDE_90(nn==N),3)) ' (90% power)'],...
        ['MDE = ' num2str(round(MDE_80(nn==N),3)) ' (80% power)'],...
        'Mean \beta_{Mood}',...
        ['N = ' num2str(N_min_90) ' (90% power)'],...
        ['N = ' num2str(N_min_80) ' (80% power)']},'Location','NorthEast')
        
% Recovery analysis
    load('recoveryAnalysis.mat')
    %Visualize betaMood
        subplot(1,3,2); hold on; box on
        scatter(results.bM_sim,results.bM_inv,10,'filled','MarkerEdgeColor','k');
        axis([-2,2,-2,2]); xticks(-2:1:2); yticks(-2:1:2); %axis square; 
        xlabel('Simulated \beta_{M}'); ylabel('Fitted \beta_{M}'); set(gca,'fontsize',12)
    %Recovery of betaMood vs. other parameters
        subplot(1,3,3); hold on; box on
        bar(all_beta(end,:));
        labels = {'intc.','kD','kR','kE','\beta_{0,D}','\beta_{0,R}','\beta_{0,E}','\beta_{1,D}','\beta_{1,R}','\beta_{1,E}',...
            '\gamma_D','\gamma_R','\gamma_E','\beta_{M}'};
        xticks(1:length(labels)); set(gca,'XTickLabel',labels,'fontsize',12); xtickangle(90)
        ylabel('Regression weight'); ylim([-0.05,0.9]);