# Ocular biofeedback
Here are the MATLAB codes to apply microbreaks to impede fatigue using a biofeedback system based on eye tracking in connection with the following paper:
Marandi RZ, Madeleine P, Omland Ø, Vuillerme N, Samani A. An oculometrics-based biofeedback system to impede fatigue development during computer work: A proof-of-concept study. PloS one. 2019;14(5).

The files are described as follows:

WAME2f.m and WAME2f.fig are the MATLAB files made using GUIDE v2.5 in MATLAB to make graphical user interface.
ppl_fltr.mat is the digital filter used in the pre-processing of pupil data.
Mdl_EnsTree.mat is a pre-trained classification model used to detect fatigue based on eye movements and pupillary responses.
ET_data_streaming.m is to test the streaming of the data from the computer running the eye-tracker's software (Eyetrac7) to a computer that receives the data.
ET_RTcom directory contains the m-files provided by the manufacturer (Applied Science Laboratories, Bedford, MA, USA) required to stream the data between the eye-tracker's computer (as host) and a second computer (as receiver). 
Day1_pnt, Day2_pnt, and Day3_pnt are respectively the data of the predefined randomized points required to execute the pattern generation of the computer task explained in the paper.

To run the program it is needed to first configure the host and receiver computers for a TCP-IP connection. The eye-tracker's software must be also configured to provide the required data in the data selection.
The selected data are as follows:
    EH_gaze_hcoord
    EH_gaze_vcoord
    EH_gaze_length
    PplDiam
    EH_gd_x
    EH_gd_y
    EH_gd_z.

The WAME2f.m includes necessary comments required to understand each part of the program.
Note: The codes works for a video-based head-mounted eye tracker and it is required to be modified for other models or types of eye trackers.
The codes are provided only as an assistive tool for research not for commercial applications.
 
