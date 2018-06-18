clear;clc;
close all;

load('..\..\entrainment\data\eeg_label.mat')
load('..\..\fieldtrip_eeg_clean\mat\acticap-64ch-standard2_ferrara.mat');
load('..\..\fieldtrip_eeg_clean\mat\acticap-64ch-standard2_ferrara_neighb.mat');
load('..\data\partialCoh\Delayed_partialCohopp-removedFirst-0.5.mat')
subject_name = {'Alice','Lucrezia','Elena','Jonluca','Manu','Sara','Marco','Elisa','Pasquale','Linda','Leonardo','Gianluca1','Federica','Silvia'...
    ,'Andrea','Giorgia','Laura','Daniel','Giada','Pagani','Silvia2',...
    'Elenora','Martina','Tommaso','Francesca'};
feature = {'envelop';'jawaopening';'lipaparature';'TTCD';'TBCD';'TMCD';'lipProtrusion'};
% 
% delay = 1:11;
% P = [];
% for d =1:length(delay)
%     A = Data(:,1,d);
%     B = Data(:,2,d);
%     %     figure;
%     %     ksdensity(A)
%     %     hold on
%     %     ksdensity(B)
%     [pval, t_orig, crit_t, est_alpha, seed_state]=mult_comp_perm_t2(A,B,1000,-1);
%     P = [P;pval];
%     
%     
%     
% end
% [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(P)

%%
load('..\data\partialCoh\partial-removedFirst-0.5all.mat')
subj = 25;
delay = -0.5:0.1:0.5;
%dd=6
A = squeeze(Data(:,1,:,:));
B = squeeze(Data(:,2,:,:));

ori=[];null=[];
for i= 1:subj
    ori{i}.cohspctrm = (squeeze(A(i,:,:)));% - squeeze(B(i,:,:));
    null{i}.cohspctrm = (squeeze(B(i,:,:)));% ori{i}.cohspctrm;
    %null{i}.cohspctrm(:) = 0;

    ori{i}.label = label;
    null{i}.label = label;
    
    ori{i}.freq = delay;
    null{i}.freq = ori{i}.freq;
end

cfg = [];
cfg.method    = 'montecarlo';
cfg.statistic = 'ft_statfun_depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg.tail             = 0;
cfg.clustertail      = cfg.tail;
cfg.alpha            = 0.025;
cfg.correcttail      = 'prob';  % 'prob''alpha'
cfg.numrandomization = 1000;
cfg.neighbours       = neighbours;
cfg.parameter   = 'cohspctrm';
cfg.channel          = label;
%     cfg.avgoverfreq = 'yes';

design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i) = i;
end
design(2,1:subj)        = 1;
design(2,subj+1:2*subj) = 2;

cfg.design   = design;
cfg.uvar     = 1;
cfg.ivar     = 2;

[statA] = ft_freqstatistics(cfg,ori{:},null{:});

cfg=[];
cfg.layout  = lay;
cfg.parameter = 'stat';
cfg.subplotsize = [1 1];
% cfg.alpha                     = 0.08;
ft_clusterplot(cfg,statA);

%%
D = [];
for d=1:11
    S = [];
    for s=1:25
        [a,b] = max(A(s,:,d));
        S = [S;b];
    end
    D = [D S];
end
d = D;
a = unique(d);
out = [a,histc(d(:),a)];
out = sortrows(out,2);
a = label(out(:,1));

D=[];
aa = squeeze(mean(A,1));
for d=1:11
    [a,b] = max(aa(:,d));
    D = [D;b];
end
a = label(D);
















