# EM-for-BAYOMA
Bayesian operational modal analysis based on the expectation-maximization algorithm.
An expectation-maximization is proposed to speed up the original Bayesian FFT method for operational modal analysis.
Fast computations for both the most probable value (MPV) and the posterior covariance matrix (PCM) are invloved.
Before running the algorithm, the frequency bands should be manually selected by conbining the information in the spectra.
The significant peak in the power spectrum density (PSD) tells where the structural mode possibly locates.
The singular value (SV) spectrum consists of eiqenvalues of the PSD matrix at each frequency point.
In a resonance band, the number of lines significantly above others reveals the number of modes, and other lines indicate the noise.
