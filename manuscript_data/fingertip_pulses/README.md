# fingertip_pulses
Individual pulse waveforms (width and amplitude normalized) used in CNN model, with age and participant ID as labels. 

1) Fingertip_all_waveforms.csv : All Fingertip-waveforms extracted (which means they passed the pre-processing and signal quality checks). These are only width-normalized (into 200 samples) but pulse amplitudes are not normalized. 
  
3) Fingertip_waveforms_split folder : Waveforms used for the training, validation, and test sets in the CNN model. The data were split within folds according to the 10-fold cross-validation scheme described in the Methods.
