function stat=coherence_stat(PLV,PLV_S,freq_band,a,label,neighbours)


subj=size(PLV,1);
ori=[];null=[];
for i= 1:subj
    ori{i}.cohspctrm =squeeze(PLV(i,:,a)) - squeeze(PLV_S(i,:,a));
    null{i}.cohspctrm = squeeze(PLV_S(i,:,a));% ori{i}.cohspctrm;
    null{i}.cohspctrm(:) = 0;
       
    ori{i}.label = label;
    null{i}.label = label;
    
    ori{i}.freq = freq_band(a);
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
cfg.correcttail      = 'alpha';  % 'prob''alpha'
cfg.numrandomization = 500;
cfg.neighbours       = neighbours;
cfg.parameter   = 'cohspctrm';
cfg.channel          = label;
% cfg.avgoverfreq = 'yes';

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

[stat] = ft_freqstatistics(cfg,ori{:},null{:});


end