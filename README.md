# Asokan_et_al_2023_CellReports
Original codes used for data analysis in the study by Asokan et al., "Potentiation of cholinergic and corticofugal inputs to the lateral amygdala in threat learning", Cell Reports (2023).

## Video_and_Pupil_analysis

1) DLC_Video_PreProcess_AnteriorPosteriorPupillabelsOnly - To read the DLC output file and obtain the pupil diameter from Anterior and Posterior labels
2) get_orofacial_ROIs_compute_Motion - draw a region of interest over the mouse cheek area (posterior to whisker pad) and compute motion energy
3) Motion_energy_video - code to create a sample video demonstrating changes in facial motion energy

## Cell counting

1) Use ImageJ to mark the cells, add to the ROI manager, and save the coordinates
2) Mark the pial surface as a line and run the macro - Macro_DistanceToPialSurface.ijm (function written by Michael Cammer, Microscopy core NYU Langone Medical Center)
to compute the distance of the cells from the pial surface
3) Counting_Cells - reads in the spreadsheet with neuron counts and distance from pia, computes and plots cell counts and coordinates as rain cloud plots 