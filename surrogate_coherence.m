%% coherence between speech and all egg channels using diiferent metric value
clc; clear; close all;

data_path = '..\data\partialCoh\';
%data_path = '\\10.96.7.1\projects\current\listen_italian_motor_entrainment\analysis\python\data\partialCoh\';

trailLen = 2;
removedFirst = 0.5;
feature = {'envelop';'jawaopening';'lipaparature';'TTCD';'TBCD';'TMCD';'lipProtrusion'};
condition = {'Hyper','Normal','Hypo','All'};
condition = {'All'};

delay = 0:0.1:1;
subject_name = {'Alice','Lucrezia','Elena','Jonluca','Manu','Sara','Marco','Elisa','Pasquale','Linda','Leonardo','Gianluca1','Federica','Silvia','Andrea','Giorgia','Laura','Daniel','Giada','Pagani','Silvia2',...
    'Elenora','Martina','Tommaso','Francesca'};

freq_band = [1:40];

permute_no = 500;

for s = 1:length(subject_name)
    k=1;
    
    data = [];
    data.Feature = [];
    data.Condition = [];
    data.Delay = [];
    data.Data = [];
    data.Surrogate = [];
    data.noTrials = [];
    data.Subject = [];
    for c=1:length(condition)
        for d=1:length(delay)
            if(d==1)
                dd = '0.0';
            elseif(d==11)
                dd = '1.0';
            else
                dd = num2str(delay(d));
            end
            fiff_file = [data_path 'partialCoh-trailLen-' num2str(trailLen) '-removedFirst-' num2str(removedFirst) 's-condition-' condition{c} '-delay-' dd 's-' subject_name{s} '.fif'];
            
            
            cfg = [];
            cfg.dataset = fiff_file;
            A = ft_preprocessing(cfg);
            
            cfg=[];
            cfg.channel ={'EEG'};
            eeg = ft_selectdata(cfg,A);
            eeg  = cat(3,eeg.trial{:});
            
            cfg.channel =feature;
            speech = ft_selectdata(cfg,A);
            speech  = cat(3,speech.trial{:});
            
            cfg            = [];
            cfg.output     = 'powandcsd';
            cfg.method     = 'mtmfft';
            cfg.taper      =  'hanning';
            cfg.foi        = freq_band;
            cfg.pad        = 'nextpow2';
            cfg.channelcmb = {'EEG' feature{1};'EEG' feature{2};'EEG' feature{3};'EEG' feature{4};'EEG' feature{5};'EEG' feature{6};'EEG' feature{7}};
            cfgg            = [];
            cfgg.method     = 'coh';
            cfgg.complex     = 'abs';
            
            % original
            freqA          = ft_freqanalysis(cfg, A);
            X = ft_connectivityanalysis(cfgg, freqA);
            
            B = zeros(413,length(freq_band),permute_no);
            a = 1:length(A.trial);
            parfor p=1:permute_no
                pp = datasample(a,length(a));
                b = cat(1, eeg, speech(:,:,pp));
                b = squeeze(num2cell(b,[1 2]))';
                
                BB = A;
                BB.trial = b;
                
                freqA          = ft_freqanalysis(cfg, BB);
                x = ft_connectivityanalysis(cfgg, freqA);
                B(:,:,p) = x.cohspctrm;
                
%                 aa = mean(x.cohspctrm(1:59,1:3),2);
%                 plot(aa);hold on;
                disp(['------------' num2str(p)])
%                 pause(0.5)
            end
            
            b = reshape(1:413,59,7);
            %%
            for ff=1:length(feature)
                data.Feature{k} = feature{ff};
                data.Condition{k} = condition{c};
                data.Delay{k} = num2str(delay(d));
                data.Data{k} = X.cohspctrm(b(:,ff),:);
                data.Surrogate{k} = mean(B(b(:,ff),:,:),3);
                data.noTrials{k} = length(A.trial);
                data.Subject{k} = subject_name{s};
                k=k+1;
            end
%             figure;imagesc(data.Data{1, 1} -data.Surrogate{1, 1});
        end
    end
    
    freq=freq_band;
    save(['..\data\SurrogateCoherence\SurrogateCoherence-' subject_name{s} '.mat'],'data','freq');
end



data_path = '..\data\SurrogateCoherence\';

%% combine them all to single file
data = [];
data.Feature = [];
data.Condition = [];
data.Delay = [];
data.Data = [];
data.Surrogate = [];
data.noTrials = [];
data.Subject = [];

for k = 1:length(subject_name)
    a = load([data_path 'SurrogateCoherence-' subject_name{k} '.mat']);
    
    data.Feature{k} = a.data.Feature;
    data.Condition{k} = a.data.Condition;
    data.Delay{k} = a.data.Delay;
    data.Data{k} = a.data.Data;
    data.Surrogate{k} = a.data.Surrogate;
    data.noTrials{k} = a.data.noTrials;
    data.Subject{k} = a.data.Subject;
end
data1 = data;
data = [];
data.Feature = horzcat(data1.Feature{:});
data.Condition = horzcat(data1.Condition{:});
data.Delay = horzcat(data1.Delay{:});
data.Data = horzcat(data1.Data{:});
data.Surrogate = horzcat(data1.Surrogate{:});
data.noTrials = horzcat(data1.noTrials{:});
data.Subject = horzcat(data1.Subject{:});
save([data_path 'SurrogateCoherence-' num2str(removedFirst) '.mat'],'data','freq');







