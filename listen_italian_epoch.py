

import pandas as pd
from mne.event import define_target_events
import mne
import numpy as np

def listen_italian_epoch(raw, mat,Tmin, Tmax):
# extract trials of tmax second and remove the wrong answer trials and seperate the in three conditions
# start and end of an epoch in sec.

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

	hyper = events[hyper-1]
	hyper[:,2] = 1
	normal = events[normal-1]
	normal[:,2] = 2
	hypo = events[hypo-1]
	hypo[:,2] = 3
	a = np.vstack((hyper,normal,hypo))
	events = np.sort(a, axis=0)

	# epoching
	reject = dict(eeg=180e-6)
	event_id = {'hyper': 1,'normal': 2,'hypo': 3}
	# Set up indices of channels to include in analysis
	picks = mne.pick_types(raw.info, meg=False, eeg=True, stim=False, eog=False,misc=False)
	epochs_params = dict(events=events, picks=picks,event_id=event_id, tmin=Tmin, tmax=Tmax,reject=reject,preload=True)
	epochs = mne.Epochs(raw, **epochs_params)
	return epochs