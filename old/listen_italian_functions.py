import pandas as pd
from mne.event import define_target_events
import mne
import numpy as np
from mne.connectivity import spectral_connectivity

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

def coherence_measure(epochs,fmin, fmax,indices):
	# calculate coherence
	sfreq = epochs.info['sfreq']  # the sampling frequency
	tmin = 0	# exclude the baseline period
	con, freqs, times, n_epochs, n_tapers = spectral_connectivity(epochs, method='coh',mode='multitaper', sfreq=sfreq, fmin=fmin, fmax=fmax,indices=indices,faverage=True, tmin=tmin, mt_adaptive=False, n_jobs=1)
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
	

	
def coherence_preprocess(epochs,feature,d,trial_len):
	extra_channels = ['envelop','jawaopening','lipaparature','lipProtrusion','TTCD','TMCD','TBCD']
	eeg_channles = np.setdiff1d(epochs.ch_names, extra_channels)
	ch_names = np.hstack((eeg_channles,feature))
	ch_types = np.repeat('eeg', len(ch_names))
	info = mne.create_info(ch_names = ch_names.tolist(),ch_types = ch_types,sfreq = epochs.info['sfreq'])
	event_id = {'hyper': 1,'normal': 2,'hypo': 3}

	eeg = epochs.copy().drop_channels(extra_channels)
	speech = epochs.copy().drop_channels(eeg_channles)
	speech.drop_channels(np.setdiff1d(extra_channels,feature))
	
	E = eeg.copy().crop(d+0.5,d+0.5+trial_len)
	S = speech.copy().crop(1,1+trial_len)
	c = np.concatenate((E.get_data(),S.get_data()),axis=1)
	epochs = mne.EpochsArray(c, info, E.events, 0,event_id)
	
	hyper = epochs['hyper'].crop(0,epochs.tmax)
	normal = epochs['normal'].crop(0,epochs.tmax)
	hypo = epochs['hypo'].crop(0,epochs.tmax)
	
	return hyper, normal, hypo	
	
	
	
	
	
	
	
