clear;clc;
close all;

load('..\..\entrainment\data\eeg_label.mat')
load('..\..\fieldtrip_eeg_clean\mat\acticap-64ch-standard2_ferrara.mat');
load('..\..\fieldtrip_eeg_clean\mat\acticap-64ch-standard2_ferrara_neighb.mat');
load('..\data\SurrogateCoherence\SurrogateCoherence-0.5.mat')
subject_name = {'Alice','Lucrezia','Elena','Jonluca','Manu','Sara','Marco','Elisa','Pasquale','Linda','Leonardo','Gianluca1','Federica','Silvia','Andrea','Giorgia','Laura','Daniel','Giada','Pagani','Silvia2',...
    'Elenora','Martina','Tommaso','Francesca'};
feature = {'envelop';'jawaopening';'lipaparature';'TTCD';'TBCD';'TMCD';'lipProtrusion'};
feature = 'envelop';

delay = 0:0.1:1;
freq_band = 1:40;
target_freq = 1:3;

D = [];
for d =1:length(delay)
    
    dd = num2str(delay(d));
    
    PLVA = [];PLV_SA=[];
    for s = 1:length(subject_name)
        a = find(contains(data.Subject,subject_name{s}));
        b = find(contains(data.Delay,dd));
        c = find(contains(data.Feature,feature));
        a = intersect(a,b);
        a = intersect(a,c);
        
        PLVA{s} = data.Data{a};
        PLV_SA{s} = data.Surrogate{a};
    end
    PLVA = cat(3,PLVA{:});
    PLVA = permute(PLVA,[3 1 2]);
    PLV_SA = cat(3,PLV_SA{:});
    PLV_SA = permute(PLV_SA,[3 1 2]);
    
    statA=coherence_stat(PLVA,PLV_SA,freq_band,target_freq,label,neighbours);
    
    if(isfield(statA,'posclusters'))
        if not(isempty(statA.posclusters))
            if(statA.posclusters(1).prob<=0.05)
                D{d,1} = statA.posclusters(1).prob;
                D{d,2} = dd;
            end
        end
    end
end
cfg=[];
cfg.layout  = lay;
cfg.parameter = 'stat';
cfg.subplotsize = [1 2];
ft_clusterplot(cfg,statA);
