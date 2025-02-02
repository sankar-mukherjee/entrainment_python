import pandas as pd
from mne.event import define_target_events
import mne
import numpy as np
from mne.connectivity import spectral_connectivity
from IPython.display import clear_output
from itertools import permutations,combinations
from scipy import stats
from mne.stats import bonferroni_correction, fdr_correction

def epoch(raw, mat,Tmin, Tmax):
	# ignore stimuli shorter than tmax ms
	events = mne.find_events(raw, stim_channel='Trigger')
	reference_id = 105  # speech onset
	target_id = 106  # speech offset
	sfreq = raw.info['sfreq']  # sampling rate
	tmin = 0  
	new_id = 99  # the new event id for a hit. If None, reference_id is used.
	fill_na = 105  # the fill value for misses
	events_, lag = define_target_events(events, reference_id, target_id,sfreq, tmin, Tmax, new_id, fill_na)
	events_  = np.where(events_[:,2] == 105)[0] +1


	#behaviour (remove the wrong answer trials and seperate the in three conditions)	
	condition= mat['behaviour']['condition'][0]
	response= mat['behaviour']['response'][0]
	a = np.hstack((condition[0],response[0]))

	df = pd.DataFrame({'condition':a[:,0],'response':a[:,1]})
	df.index = df.index + 1

	hyper  = df.loc[(df['condition'] == 1) & (df['response'] == 0)]
	normal = df.loc[(df['condition'] == 2) & (df['response'] == 0)]
	hypo   = df.loc[(df['condition'] == 3) & (df['response'] == 0)]

	events = mne.find_events(raw, stim_channel='trial_no')
	hyper = np.intersect1d(events_, hyper.index.values)
	normal = np.intersect1d(events_, normal.index.values)
	hypo = np.intersect1d(events_, hypo.index.values)

	a = np.concatenate((hyper,normal,hypo),axis=0)
	IDX = np.sort(a, axis=0) -1

	hyper = events[hyper-1]
	hyper[:,2] = 1
	normal = events[normal-1]
	normal[:,2] = 2
	hypo = events[hypo-1]
	hypo[:,2] = 3
	a = np.vstack((hyper,normal,hypo))
	events = np.sort(a, axis=0)

	# epoching
	event_id = {'hyper': 1,'normal': 2,'hypo': 3}
	# Set up indices of channels to include in analysis
	picks = mne.pick_types(raw.info, meg=False, eeg=True, stim=False, eog=False,misc=False)
	epochs_params = dict(events=events, picks=picks,event_id=event_id, tmin=Tmin, tmax=Tmax,preload=True)
	epochs = mne.Epochs(raw, **epochs_params)

	#downsampling
	epochs_resampled = epochs.copy().resample(400, npad='auto')

	#add EMA, envelop signal as extra channels
	extra_channels = ['envelop','jawaopening','lipaparature','lipProtrusion','TTCD','TMCD','TBCD']
	envelop = mat['behaviour']['envelop'][0][0]
	jawaopening = mat['behaviour']['jawaopening'][0][0]
	lipaparature = mat['behaviour']['lipaparature'][0][0]
	lipProtrusion = mat['behaviour']['lipProtrusion'][0][0]
	TTCD = mat['behaviour']['TTCD'][0][0]
	TMCD = mat['behaviour']['TMCD'][0][0]
	TBCD = mat['behaviour']['TBCD'][0][0]

	# Initialize an info structure
	ch_names = np.hstack((epochs.ch_names,extra_channels))
	ch_types = np.repeat('eeg', len(ch_names))
	info = mne.create_info(ch_names = ch_names.tolist(),ch_types = ch_types,sfreq = epochs_resampled.info['sfreq'])

	# create new epoches
	b=abs(Tmin)*epochs_resampled.info['sfreq']
	first_samp=int(b)
	last_samp=int(epochs_resampled.tmax *epochs_resampled.info['sfreq'])

	A = np.zeros((len(epochs_resampled), len(extra_channels), first_samp+last_samp+1))

	for i in range(0, len(epochs_resampled)):
		a = np.vstack((envelop[IDX[i]][0],jawaopening[IDX[i]][0].T,lipaparature[IDX[i]][0].T,lipProtrusion[IDX[i]][0].T,TTCD[IDX[i]][0].T,TMCD[IDX[i]][0].T,TBCD[IDX[i]][0].T))
		b = np.zeros((a.shape[0],first_samp))
		a = np.hstack((b,a))
		A[i,:,:] = a[:,0:first_samp+last_samp+1]

	b = epochs_resampled.get_data()
	c = np.concatenate((b,A),axis=1)

	tmin = Tmin
	epochs = mne.EpochsArray(c, info, events, tmin, event_id)
	return epochs

def coherence_measure(epochs,fmin, fmax,sfreq,indices):
	# calculate coherence
	con, freqs, times, n_epochs, n_tapers = spectral_connectivity(epochs, method='coh',mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax,indices=indices,faverage=True, tmin=0, mt_adaptive=False, block_size=1000, n_jobs=1,verbose='ERROR')
	return con, freqs, times, n_epochs, n_tapers
	
def coherence_freq(epochs,fmin, fmax,feature):
	# Compute connectivity for band containing the evoked response.	
	a=np.where(np.array(epochs.ch_names)==feature)[0][0]

	indices = (np.repeat(a,59),np.arange(0,59))   


	hyper = epochs['hyper'].crop(0,epochs.tmax)
	normal = epochs['normal'].crop(0,epochs.tmax)
	hypo = epochs['hypo'].crop(0,epochs.tmax)


	hyper, freqs, times, n_epochs, n_tapers = coherence_measure(hyper,fmin, fmax,indices)
	normal, freqs, times, n_epochs, n_tapers = coherence_measure(normal,fmin, fmax,indices)
	hypo, freqs, times, n_epochs, n_tapers = coherence_measure(hypo,fmin, fmax,indices)

	return hyper, normal, hypo
	

	
def coherence_preprocess(epochs,d,trial_len,extra_channels,eeg_channles,info,ch_names,event_id):	

	eeg = epochs.copy().drop_channels(extra_channels)
	speech = epochs.copy().drop_channels(eeg_channles)
	speech.drop_channels(ch_names)
	
	E = eeg.copy().crop(d+0.5,d+0.5+trial_len)
	S = speech.copy().crop(1,1+trial_len)
	c = np.concatenate((E.get_data(),S.get_data()),axis=1)
	epochs = mne.EpochsArray(c, info, E.events, 0,event_id)
	
	hyper = epochs['hyper'].crop(0,epochs.tmax)
	normal = epochs['normal'].crop(0,epochs.tmax)
	hypo = epochs['hypo'].crop(0,epochs.tmax)
	
	return hyper, normal, hypo	
	
def coherence_preprocess_delay(epochs,remove_first,d,trial_len,extra_channels,eeg_channles,info,ch_names,event_id,iter_freqs,indices,condition):	

    eeg = epochs.copy().drop_channels(extra_channels)
    speech = epochs.copy().drop_channels(eeg_channles)
    speech.drop_channels(ch_names)

    E = eeg.copy().crop(d+remove_first,d+remove_first+trial_len)
    S = speech.copy().crop(0.5+remove_first,0.5+remove_first+trial_len)

    events = E.events
    sfreq = E.info['sfreq']    
    c = np.concatenate((E.get_data(),S.get_data()),axis=1)
    epochs = mne.EpochsArray(c, info, E.events, 0,event_id)


    fmin = []
    fmax = []
    for fr in range(0,len(iter_freqs)):
        fmin.append(iter_freqs[fr][1])
        fmax.append(iter_freqs[fr][2])

    frames = np.zeros((413,len(iter_freqs),len(condition)))
    for i in range(0,len(condition)):
        if condition[i] != 'All':
            c = epochs[condition[i]].copy()
        else:
            c = epochs.copy()

        for fr in range(0,len(iter_freqs)):
            coh, freqs, times, n_epochs, n_tapers = coherence_measure(c,fmin, fmax,sfreq,indices)
            frames[:,fr,i] = coh[:,fr]

    return frames	
	


def partialCoherence_preprocess_delay(epochs,remove_first,d,trial_len,feat,keep_feat,condition):	

    if condition != 'All':
        E = epochs[condition].copy()
    else:
        E = epochs.copy()
    
    
    a = np.setdiff1d(feat, keep_feat)
    eeg = E.copy().drop_channels(a)

    E = eeg.copy().crop(d+remove_first,d+remove_first+trial_len)

    c = E.get_data()

    return c


 
	
def coherence_preprocess_delay_surrogate(epochs,remove_first,d,trial_len,features,eeg_channles,info,ch_names,event_id,iter_freqs,indices):	

    eeg = epochs.copy().drop_channels(features)
    speech = epochs.copy().drop_channels(eeg_channles)
    speech.drop_channels(ch_names)

    E = eeg.copy().crop(d+remove_first,d+remove_first+trial_len)
    S = speech.copy().crop(0.5+remove_first,0.5+remove_first+trial_len)
    events = E.events
    sfreq = E.info['sfreq']
    
    E = E.get_data()
    S = S.get_data()

    fmin = []
    fmax = []
    for fr in range(0,len(iter_freqs)):
        fmin.append(iter_freqs[fr][1])
        fmax.append(iter_freqs[fr][2])
    
    # all possible combination
    trial_length=S.shape[0]
    a = list(permutations(np.arange(0,trial_length), 2))
    a = np.asarray(a)
    X = np.arange(0,trial_length)

    no_surrogates = 500 #dummy value
    B=[]
    for j in range(no_surrogates):
        X = np.roll(X,1)
        while True:
            A,a = get_combinations(X,a)        
            if A.shape[0] == trial_length:
                B.append(A)
                break
            elif len(a)==0:
                break
            else:
                X = np.roll(X,1)
                print('.',end=' ')
    
    B = np.asarray(B)
    no_surrogates = len(B)
    ##
    frames = np.zeros((413,len(iter_freqs),no_surrogates))
    for i in range(no_surrogates):
        print('--------------------'+str(i))
        EE = E.copy()
        SS = S.copy()       
        
        c = np.concatenate((EE[B[i][:,0]],SS[B[i][:,1]]),axis=1)
        #epochs = mne.EpochsArray(c, info, events, 0,event_id)
        for fr in range(0,len(iter_freqs)):
            coh, freqs, times, n_epochs, n_tapers = coherence_measure(c,fmin, fmax,sfreq,indices)
            frames[:,fr,i] = coh[:,fr]
        clear_output()  
    return frames	
	
	
def get_partialCoherence(conXY,conXR,conRY,fr):
    partial_coh_XY_R=[]
    for i in range(59):
        a = (abs(conXY[i,fr]-conXR[i,fr]*conRY[fr])**2) / ((1-abs(conXR[i,fr])**2)*(1-abs(conRY[fr])**2))
        partial_coh_XY_R.append(a)

    partial_coh_XY_R = np.asarray(partial_coh_XY_R)
    
    return partial_coh_XY_R



def PartialCoherence_preprocess_delay_surrogate(epochs,remove_first,d,trial_len,feat,
                                                eeg_channles,keep_feat,condition,iter_freqs):	

    
    keep_feat = np.unique(keep_feat)
    ##############
    if condition != 'All':
        E = epochs[condition].copy()
    else:
        E = epochs.copy()
    
    eeg = E.copy().drop_channels(feat)
    speech = E.copy().drop_channels(eeg_channles)
    
    a = np.setdiff1d(feat, keep_feat)
    speech.drop_channels(a)

    E = eeg.copy().crop(d+remove_first,d+remove_first+trial_len)
    S = speech.copy().crop(0.5+remove_first,0.5+remove_first+trial_len)
    
    sfreq = E.info['sfreq']
    
    E = E.get_data()
    S = S.get_data()

    label = np.concatenate((eeg.ch_names,speech.ch_names))
    
    ##################### all possible combination
    trial_length=S.shape[0]
    a = list(permutations(np.arange(0,trial_length), 2))
    a = np.asarray(a)
    X = np.arange(0,trial_length)

    no_surrogates = 500 #dummy value
    B=[]
    for j in range(no_surrogates):
        X = np.roll(X,1)
        while True:
            A,a = get_combinations(X,a)        
            if A.shape[0] == trial_length:
                B.append(A)
                break
            elif len(a)==0:
                break
            else:
                X = np.roll(X,1)
                print('.',end=' ')
    
    B = np.asarray(B)
    no_surrogates = len(B)
    
    #######################################à
    fmin = []
    fmax = []
    for fr in range(0,len(iter_freqs)):
        fmin.append(iter_freqs[fr][1])
        fmax.append(iter_freqs[fr][2])  
    
    #######################################
    
    indices = []
    for idx in keep_feat:
        b=np.where(label==idx)[0][0]
        b = (np.repeat(b,59),np.arange(0,59))
        indices.append(b)
    
    a = list(permutations(keep_feat, 2))
    for idx in range(0,len(a)):
        b=np.where(label==a[idx][0])[0][0]
        c=np.where(label==a[idx][1])[0][0]
        d = (np.repeat(b,1),np.repeat(c,1))
        indices.append(d)
    
    INDEX = []
    b=0
    for idx in range(0,len(indices)):
        a = np.arange(b,b+len(indices[idx][0]))
        INDEX.append(a)
        b = b+len(a)
    
    indices = np.concatenate((indices),axis=1)
    indices = (indices[0],indices[1])

    
    #####################################################################    
    frames = np.zeros((2,59,len(iter_freqs),no_surrogates))
    for i in range(no_surrogates):
        print('--------------------'+str(i))
        EE = E.copy()
        SS = S.copy()
        c = np.concatenate((EE[B[i][:,0]],SS[B[i][:,1]]),axis=1)
        
        coh = get_coherence(c,sfreq,fmin,fmax,indices)
        for iii in range(0,1):
            if iii ==0:
                conXY = coh[INDEX[0],:]
                conXR = coh[INDEX[1],:]
            else:
                conXY = coh[INDEX[1],:]
                conXR = coh[INDEX[0],:]
            
            conRY = coh[INDEX[2],:][0]

            for f in range(0,len(iter_freqs)):
                frames[iii,:,f,i] = get_partialCoherence(conXY,conXR,conRY,f)
        clear_output()  
        
    return frames

def get_coherence(epochs,sfreq,fmin,fmax,indices):
    con, freqs, times, n_epochs, n_tapers = mne.connectivity.spectral_connectivity(epochs, method='coh',mode='multitaper', 
                                                                                   sfreq=sfreq, 
                                                              fmin=fmin, fmax=fmax,indices=indices,faverage=True, 
                                                              tmin=0, mt_adaptive=False, block_size=1000,verbose='ERROR')

    return con

def get_combinations(X,a):
    aa = a
    A=[]
    EEG = []
    Speech = []
    for i in range(0,len(X)):
        b = np.where(a[:,0]==X[i])
        if not len(b[0]) == 0:
            for k in range(len(b[0])):
                if not a[b[0][k],1] in Speech:
                    A.append(a[b[0][k],:])
                    EEG.append(a[b[0][k],0])
                    Speech.append(a[b[0][k],1])
                    a = np.delete(a, b[0][k], 0)
                    break
    if len(A) == len(X):                    
        return np.asarray(A),a
    else:
        return np.asarray(A),aa	
    
    
    
def run_permutation_test(x,y,numSamples):
    pooled = np.hstack([x,y])
    d0 = abs(x.mean() - y.mean())
    d = np.zeros((numSamples,))
    for k in range(numSamples):
        np.random.shuffle(pooled)
        starZ = pooled[:len(x)]
        starY = pooled[-len(y):]    
        d[k] = abs(starZ.mean() - starY.mean())

    p = len(np.where( d >= d0 )[0])
    p = p / float(numSamples)
    p = np.round(p,decimals=3)

    return p
    
def run_permutation_test2(x,y,numSamples):
    #keeping the subject association
    pooled = np.vstack([x,y]).T
    d0 = abs(x.mean() - y.mean())
    d = np.zeros((numSamples,))
    for k in range(numSamples):
        for i in range(len(x)):
            a = pooled[i]
            np.random.shuffle(a)
            pooled[i] = a
            
        starZ = pooled[:,0]
        starY = pooled[:,1]    
        d[k] = abs(starZ.mean() - starY.mean())

    p = len(np.where( d >= d0 )[0])
    p = p / float(numSamples)
    p = np.round(p,decimals=3)

    return p    

def get_P_value2(data,feat_comb,freq_band,condition,delay,subject_name):
    A = []    
    for c in condition:
        for fr in freq_band:        
            for d in delay:
                group1 = data.loc[feat_comb[0][0],d,feat_comb[0][1],d,fr,c,subject_name]['partialCoh']
                group2 = data.loc[feat_comb[1][0],d,feat_comb[1][1],d,fr,c,subject_name]['partialCoh']
                p = run_permutation_test2(group1.get_values(),group2.get_values(),10000)

                #reject_fdr, p = fdr_correction(p, alpha=0.05, method='indep')
                a = pd.DataFrame({'Delay':d,'Freq':fr,'Condition':c,'p-value':[p]})
                A.append(a)
                
    P = pd.concat((A),axis=0)
    return P    
    
    

def get_P_value(data,feat_comb,freq_band,condition,delay,subject_name):
    A = []    
    for c in condition:
        for fr in freq_band:        
            for d in delay:
                group1 = data.loc[feat_comb[0][0],d,feat_comb[0][1],d,fr,c,subject_name]['partialCoh']
                group2 = data.loc[feat_comb[1][0],d,feat_comb[1][1],d,fr,c,subject_name]['partialCoh']
                t,p = stats.ttest_rel(group1.get_values(),group2.get_values())
                #reject_fdr, p = fdr_correction(p, alpha=0.05, method='indep')
                a = pd.DataFrame({'Delay':d,'Freq':fr,'Condition':c,'p-value':[p]})
                A.append(a)
                
    P = pd.concat((A),axis=0)
    return P 

def get_P_value_coherence(data,extra_channels,freq_band,condition,delay,subject_name,numSamples):
    A = []
    
    for c in condition:
        for fr in freq_band:        
            for d in delay:
                group1 = data.loc[extra_channels[0],fr,c,d,subject_name]['Coherence']
                group2 = data.loc[extra_channels[1],fr,c,d,subject_name]['Coherence']
                a = run_permutation_test2(group1.get_values(),group2.get_values(),numSamples)
                p = pd.DataFrame({'Delay':d,'Freq':fr,'Condition':c,'FeatureTest':[extra_channels],'p-value':a})
                A.append(p)
                
    P = pd.concat((A),axis=0)
    return P   
    
    
    
    
    
    
