Matlab file `spike_detector_hilbert_vXX.m` contains implementation of detector in MATLAB. Using of a detector and file structures is described in short tutorial `TUTORIAL.m`.

# IED detector

Interictal epileptiform discharges (spikes, IEDs) are electrographic markers of epileptic tissue and their quantification is utilized 
in planning of surgical resection. Visual analysis of long-term multi-channel intracranial recordings is extremely laborious and prone 
to bias. 

Development of new and reliable techniques of automatic spike detection represents a crucial step towards increasing the information 
yield of intracranial recordings and to improve surgical outcome. 

In this study, we designed a novel and robust detection algorithm that adaptively models statistical distributions of signal envelopes 
and enables discrimination of signals containing IEDs from signals with background activity. 

This detector demonstrates performance superior both to human readers and to an established detector. It is even capable of identifying 
low-amplitude IEDs which are often missed by experts and which may represent an important source of clinical information. 

Application of the detector to non-epileptic intracranial data from patients with intractable facial pain revealed the existence of sharp 
transients with waveforms reminiscent of interictal discharges that can represent biological sources of false positive detections. 

Identification of these transients enabled us to develop and propose secondary processing steps, which may exclude these transients, 
improving the detector’s specificity and having important implications for future development of spike detectors in general.

## Publication
Janca, R., Jezdik, P., Cmejla, R. et al. Detection of Interictal Epileptiform Discharges Using Signal Envelope Distribution Modelling: 
Application to Epileptic and Non-Epileptic Intracranial Recordings. Brain Topogr 28, 172–183 (2015). doi: [10.1007/s10548-014-0379-1](https://doi.org/10.1007/s10548-014-0379-1)

## Test data
Files contain dataset of seven patients labelled by three experts. Detailed information of electrode placement, diagnoses etc. are included.

Files are available on [Google Drive](https://drive.google.com/drive/folders/1lAjGZ7cXkD0bxNfMybVrERvPgCLtakJ_?usp=sharing).

## Contact
Questions and error reports sent to:

**Radek Janča**: jancarad [at] fel.cvut.cz

**Přemysl Jiruška**: jiruskapremysl [at] gmail.com

**Petr Marusič**: Petr.Marusic [at] fnmotol.cz
