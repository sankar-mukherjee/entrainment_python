%% coherence between speech and all egg channels using diiferent metric value
clc; clear; close all;


data_path = 'C:\Users\SMukherjee\Desktop\projects\listen_italian_motor_entrainment\analysis\python\data\partialCoh\';

trailLen = 2;
removedFirst = 0;
feature = {'envelop';'jawaopening';'lipaparature';'TTCD';'TBCD';'TMCD';'lipProtrusion'};
condition = {'hyper','normal','hypo','all'};
delay = 0;

subject_name = {'Alice','Lucrezia','Elena','Jonluca','Manu','Sara','Marco','Elisa','Pasquale','Linda','Leonardo','Gianluca1','Federica','Silvia','Andrea','Giorgia','Laura','Daniel','Giada','Pagani','Silvia2',...
    'Elenora','Martina','Tommaso','Francesca'};

remove_feature_contribution = {{'envelop'};{'jawaopening'};{'envelop'};{'lipaparature'};};
freq_band = [1:40];

data = [];
data.Feature = [];
data.RemovedFeature = [];
data.Condition = [];
data.Delay = [];
data.Data = [];
data.noTrials = [];
data.Subject = [];

k=1;
for s = 1:length(subject_name)    
    for c=1:length(condition)
        for d=1:length(delay)
            fiff_file = [data_path 'partialCoh-trailLen-' num2str(trailLen) '-removedFirst-' num2str(removedFirst) 's-condition-' condition{c} '-delay-' num2str(delay(d)) 's-' subject_name{s} '.fif'];
            
            
            cfg = [];
            cfg.dataset = fiff_file;
            A = ft_preprocessing(cfg);
            
            cfg            = [];
            cfg.output     = 'powandcsd';
            cfg.method     = 'mtmfft';
            cfg.taper      =  'hanning';
            cfg.foi        = freq_band;
            cfg.pad        = 'nextpow2';
            freqA          = ft_freqanalysis(cfg, A);
            
            cfgg            = [];
            cfgg.method     = 'coh';
            cfgg.complex     = 'abs';
            cfg.channelcmb = { 'All'};
            
            for ff=1:length(feature)
                cfgg.partchannel = feature{ff};
                a = ft_connectivityanalysis(cfgg, freqA);
                
                for chComb=60:65
                    x = strsplit(a.label{chComb},'\');
                    
                    partialCoh=squeeze(a.cohspctrm(chComb,1:59,:));
                    
                    
                    data.Feature{k} = x{1};
                    data.RemovedFeature{k} = x{2};
                    data.Condition{k} = condition{c};
                    data.Delay{k} = num2str(delay(d));
                    data.Data{k} = partialCoh;
                    data.noTrials{k} = length(A.trial);
                    data.Subject{k} = subject_name{s};
                    k=k+1;
                end
            end
        end
    end
end
freq=freqA.freq;
save([data_path 'PartialCoherence.mat'],'data','freq');












