# Asokan_et_al_2023_CellReports
Codes written by Meenakshi M. Asokan used for data analysis in the study by Asokan et al., "Potentiation of cholinergic and corticofugal inputs to the lateral amygdala in threat learning", Cell Reports (2023).

## Codes corresponding to analysis of data presented in the Figures:
Source data deposited at https://doi.org/10.17632/m2s6dry5s5.1

1) Figure1_Behavior
2) Figure2_Neurograms_fft
3) Figure2_NeuralTrajectories
4) Figure3_CAmy_CellCounting_Density
5) Figure3_Phototagging
6) Figure3_CSResp_AssymetryIndices
7) Figure4_LFP_LaserResp_GrowthSlope
8) Figure5_stLFP 
9) Figure6_ACh

## Other data pre-processing codes:
### Video and Pupil analysis:
1) DLC_Video_PreProcess_AnteriorPosteriorPupillabelsOnly - To read the DLC output file and obtain the pupil diameter from Anterior and Posterior labels
2) get_orofacial_ROIs_compute_Motion - To draw a region of interest over the mouse cheek area (posterior to whisker pad) and compute motion energy
3) Motion_energy_video - To create a sample video demonstrating changes in facial motion energy
### Spike-triggered LFP analysis:
1) Compute_stLFP_linear_deconvolution - To compute stLFPs via linear deconvolution

