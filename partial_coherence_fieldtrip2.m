%% coherence between speech and all egg channels using diiferent metric value
clc; clear; close all;

data_path = '..\data\partialCoh\';

trailLen = 2;
removedFirst = 0.5;
feature = {'envelop';'jawaopening';'lipaparature';'TTCD';'TBCD';'TMCD';'lipProtrusion'};
condition = {'Hyper','Normal','Hypo','All'};
delay = 0:0.1:1;

subject_name = {'Alice','Lucrezia','Elena','Jonluca','Manu','Sara','Marco','Elisa','Pasquale','Linda','Leonardo','Gianluca1','Federica','Silvia','Andrea','Giorgia','Laura','Daniel','Giada','Pagani','Silvia2',...
    'Elenora','Martina','Tommaso','Francesca'};

freq_band = [1,3;4,6];

for s = 1:length(subject_name)
    k=1;
    
    data = [];
    data.Feature = [];
    data.RemovedFeature = [];
    data.Condition = [];
    data.Delay = [];
    data.Frequency = [];
    data.Data = [];
    data.noTrials = [];
    data.Subject = [];
    for fr = 1:size(freq_band,1)
        for c=1:length(condition)
            for d=1:length(delay)
                if(d==1)
                    dd = '0.0';
                elseif(d==11)
                    dd = '1.0';
                else
                    dd = num2str(delay(d));
                end
                
                freqA = load([data_path 'partialCoh-trailLen-' num2str(trailLen) '-removedFirst-' num2str(removedFirst) 's-condition-' ...
                    condition{c} '-delay-' dd '-fr-' num2str(freq_band(fr,1)) '-' num2str(freq_band(fr,2)) 'Hz-' subject_name{s} '.mat']);
                freqA.crsspctrm = freqA.crsspctrm';
                freqA.powspctrm = freqA.powspctrm';
                freqA.label = freqA.label';
                
                
                cfgg            = [];
                cfgg.method     = 'coh';
                cfgg.complex     = 'abs';
                
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
                        data.Frequency{k} = ['fr-' num2str(freq_band(fr,1)) '-' num2str(freq_band(fr,2)) 'Hz'];
                        data.Data{k} = partialCoh;
                        data.Subject{k} = subject_name{s};
                        k=k+1;
                    end
                end
            end
        end
    end
    freq=freqA.freq;
    save([data_path 'PartialCoherence-' subject_name{s} '.mat'],'data','freq');
    
end


%% combine them all to single file
data = [];
data.Feature = [];
data.RemovedFeature = [];
data.Condition = [];
data.Delay = [];
data.Frequency = [];
data.Data = [];
data.Subject = [];

for k = 1:length(subject_name)
    a = load([data_path 'PartialCoherence-' subject_name{k} '.mat']);
    
    data.Feature{k} = a.data.Feature;
    data.RemovedFeature{k} = a.data.RemovedFeature;
    data.Condition{k} = a.data.Condition;
    data.Delay{k} = a.data.Delay;
    data.Frequency{k} = a.data.Frequency;
    data.Data{k} = a.data.Data;
    data.Subject{k} = a.data.Subject;
end
data1 = data;
data = [];
data.Feature = horzcat(data1.Feature{:});
data.RemovedFeature = horzcat(data1.RemovedFeature{:});
data.Condition = horzcat(data1.Condition{:});
data.Delay = horzcat(data1.Delay{:});
data.Frequency = horzcat(data1.Frequency{:});
data.Data = horzcat(data1.Data{:});
data.Subject = horzcat(data1.Subject{:});
save([data_path 'PartialCoherence_' num2str(removedFirst) '.mat'],'data','freq');









