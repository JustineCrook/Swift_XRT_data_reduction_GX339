"""
NOTES FOR GX339
Schedule: https://www.swift.psu.edu/operations/obsSchedule.php?source_name=Gx+339-4&ra=255.7057818297&dec=-48.78974665404001

"""


######################################################################################################################

## FUNCTION TO RUN
analysis_type = "get_final_results"
 

######################################################################################################################

## DATA RETRIEVAL -- get_data

# Target coordinates
target_coords = [255.7057818297, -48.78974665404001]


## get_data parameters
#target_names = ['GX 339-4']
#target_ids = ['00032898']
#segments_lc=[[172, 174,175,176,177,178,179,180,181,182,183,184,185,186,187,188, 189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209, 211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231, 232,233,234,235,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255]]

#target_names = ['GX 339-4']
#target_ids = ['00014052']
#segments_lc=[[1,2,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,102,103,104,105,106,108,107,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,154,153,156,155,157,158,162,161,164,163,165,166,168,167,169,170,171,172,175,174,176,177]]

#target_names = ['GX 339-4']
#target_ids = ['00032490']
#segments_lc= [[17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56]]

#target_names = ['GX_339m4']
#target_ids = ['00089196']
#segments_lc= [[1,2,3,4]]


#target_names = ['GX 339-4']
#target_ids = ['00014218']
#segments_lc= [[1,2,3,4,5,6,7,8,9,10]]

#target_names = ['GX_339m4']
#target_ids = ['00088938']
#segments_lc= [[1]]



# SPEC:

#target_names = ['GX 339-4', 'GX 339-4']
#target_ids = ['00032898','00014052']
#segments_spec = [[182,183,184,185,186,187,188, 232,233,234,235,236,237,238], [14,15,16,17,18, 158,161,163,164,165,166,167]] # observation number

#target_names = ['GX 339-4']
#target_ids = ['00032898']
#segments_spec = [[189,190,191,192,193,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209]]
#segments_spec = [[211,212,213,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231]]


#target_names = ['GX 339-4']
#target_ids = ['00032490']
#segments_spec = [[32, 33,35,36,37,38,39,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56]]


#target_names = ['GX 339-4']
#target_ids = ['00032898']
#segments_spec = [[239,240,241,243,245,246,248,249,250,252,253,254,255]]

#target_names = ['GX 339-4']
#target_ids = ['00014052']
#segments_spec = [[5,6,7,8,9,56,57,58,59,60,61,62,63,94,95,96,97,98,99,100,102,103,104,105,106,107,108,109,110]]
#segments_spec = [[111,112,113,114,115,116,117,118,136,137,140,141,142,143,144,145,146,147,148,150,151,156,157]]

#target_names = ['GX_339m4']
#target_ids = ['00089196']
#segments_spec = [[1,2]]


# Quiescence:
#target_names = ['GX 339-4']
#target_ids = ['00032898']
#segments_spec = [[174,175,176,177,178,179,180,181]]




# New IDs for additional SS spectra:

#############
#target_names = ['GX 339-4']
#target_ids = ['00032490']
#segments_spec = [[17,19,20,21,22,23,24,25,26,27,28,29,30,31]]

#target_names = ['GX 339-4']
#target_ids = ['00014218']
#segments_spec = [[9]]


#target_names = ['GX 339-4']
#target_ids = ['00032490']
#segments_spec = [[40]]



#############

#target_names = ['GX 339-4']
#target_ids = ['00014052']
#segments_spec = [[11,12,13,14,15,16,17,18,19,21,22,23,25,26,27,31,36,40,43,46,47,49,50,52,55]]

#target_names = ['GX 339-4']
#target_ids = ['00014218']
#segments_spec = [[7]]

#target_names = ['GX_339m4']
#target_ids = ['00089196']
#segments_spec = [[4]]



# New IDs for spectra that Lilia didn't include:

#target_names = ['GX 339-4']
#target_ids = ['00032898']
#segments_spec = [[172,214,242,244,247,251]]

#target_names = ['GX 339-4']
#target_ids = ['00032490']
#segments_spec = [[18,34,40]]

#target_names = ['GX_339m4']
#target_ids = ['00088938']
#segments_spec = [ [1]]


#target_names = ['GX 339-4']
#target_ids = [ '00014052']
#segments_spec = [[1,2,10,20,24,28,29,30,32,33,34,35,37,38,39,41,42,44,45,48,51,53,54,64]]


#target_names = ['GX 339-4']
#target_ids = [ '00014052']
#segments_spec = [[65,66,67,68,69,70,71,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89]]


#target_names = ['GX 339-4']
#target_ids = [ '00014052']
#segments_spec = [[90,91,92,93,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134]]


#target_names = ['GX 339-4']
#target_ids = [ '00014052']
#segments_spec = [[135,138,139,149,152,154,153,155,162,168,169,170,171,172,175,174,176,177]]

#target_names = ['GX_339m4']
#target_ids = ['00089196']
#segments_spec = [[3,4]]

#target_names = ['GX 339-4']
#target_ids = ['00014218']
#segments_spec = [[1,2,3,4,5,6,8,10]]


target_names = []
target_ids =[]
segments_lc =[]
segments_spec=[]




######################################################################################################################

## GROUPING SPECTRA -- group_spectra

# Minimum number of counts for the chi^2 bins when ncounts>300
min_counts_chi = 20
# Minimum energy of the first bin, which is the energy from which we start the binning
min_E_keV = 0.5 


######################################################################################################################

## PLOT LIGHT CURVE AND HARNESS RATIO -- plot_lightcurve_and_hr
# No inputs required


## RESULTS -- FROM LIGHT CURVE ROUTINE
# After running the light curve routine, the code outputs which observations are upper limits.
# These should not be fitted, so are input to the spectral fitting function (steps 3 and 5).
# The code also outputs the count rate, used when converting count rate to flux (step 6).
# Note that the upper limits from the pipeline are at the 3-σ level, calculated using the Bayesian method (see Kraft, Burrows & Nousek, 1991, ApJ, 374, 344).
uplims_IDs = ['00032898241pc', '00032898242pc', '00032898244pc', '00032898246pc', '00032898247pc','00032898248pc','00014052068pc', '00014052070pc','00014052071pc','00014052082pc','00014052084pc','00014052086pc','00014052087pc','00014052090pc','00014052091pc','00014052165wt']
uplims_count_rate_ar = [0.01447 ,0.012858,0.010053,0.011865,0.11846 ,0.028886,0.015947,0.017756,0.018122,0.018175,0.017761, 0.010015,0.010583,0.030287,0.016985,0.60695 ]
uplims_MJDs=  [59104.7319654777,59111.3084145078,59125.1141200413,59139.8468777103,59146.0163518411, 59153.1956094499,59629.0283928001,59643.6846448047,59650.3728300848,59720.0412057758, 59734.7508624447,59748.2866895142,59755.2434734378,59776.6737704285,59783.591131851 , 60193.6120083359]
uplims_MJDs_er = [0.0118400231481481,0.0103019907407407,0.0100117824074074,0.0111725694444444,0.0017412037037037, 0.0095184490740741,0.0105921759259259,0.0105341435185185,0.1411085648148148,0.0099247222222222, 0.0072549305555556,0.0092572685185185,0.0095765046296296,0.0083867361111111,0.348943449074074 , 0.406027916666666 ]



######################################################################################################################

## INPUTS --- INITIAL MODEL FITTING -- unconstrained_fit
# Uses uplims_IDs defined above

# Define models to fit -- used for fit_spec
models_unconstrained = ['pegged_powerlaw', 'powerlaw+diskbb','diskbb'] 



######################################################################################################################

## GET PARAMETER VALUES (N_H, GAMMA, T_IN) -- get_param_averages

# Based on the spectral fit results and any prior knowledge regarding spectral state ranges, define which spectral model to use for each point in time.
# model_indexes is a multi-dimensional array. Each outer element corresponds to the models defined above. The inner elements are ranges for each model.
# The index ranges are inclusive (i.e. [start, end]), and correspond to those listed in the spectral_results.txt file.
# e.g. models_indexes=[ [[6, 17], [19, 21]] , [[0, 5]] , [[]] ] means the first model should be used for the observations corresponding to indexes 6-17 and 19-21, the second model should be used for observations 0-5, and the last model that was fit should not be used
# If none, use models_indexes = []


#models_indexes_av =[ [[0,81], [95,140], [187,292]], [0,292], [[82,94], [141, 186] ] ]
#models_indexes_av =[ [[0,81], [97,138], [188,292]], [ [90,96], [139,143], [186,187]], [[82,89], [144, 185] ] ]
models_indexes_av =[ [[0,81], [95,138], [188,292]], [ [90,94], [139,143], [186,187]], [[82,89], [144, 185] ] ]



######################################################################################################################

## INPUTS --- MODEL FITTING WITH SOME FIXED PARAMETERS -- constrained_fit
# Re-run spectral fitting, but this time we want to fit the parameters for some observations -- using the results of step 4 and/or any prior knowledge (e.g. nH values in papers).
# Uses uplims_IDs, models, low_energy_fit_bound_keV, and low_count_indexes -- defined above


models_constrained = ['pegged_powerlaw', 'powerlaw+diskbb','diskbb'] 


## FROM MODEL FITTING RESULTS
# The spectral_results.txt file lists the low_count_indexes, which correspond to observations that cannot initial be fit (due to low count rate), so require some fitting parameters to be fixed.
# All the low count spectra are in the HS!!!!!!!!!!!!!!!!!!!!!
low_count_indexes = [  0,  1,  2,  3,  4,  5,  6,  7, 67, 121,122,123,124,201,202,203,204,205,206,207,208,209,210,211, 212,239,240,241,242,243,244,245,246,247,248,249,250,251,253,254,256,257,258,259,260,261,262,263, 264,265,266,267,268,269,271,272,274,283,284,287]
# Divide the set in two, since they have different behaviour. 
# Set 1:
low_count_indexes_1 = [  0,  1,  2,  3,  4,  5,  6,  7,121,122,123,124,201,202,203,204,205,206,207,208,209,210,211, 212]
# Set 2: For this second set of low-count spectra, we can see from the trend in photon index for the powerlaw that the photon index is a bit higher
low_count_indexes_2 = [67, 239,240,241,242,243,244,245,246,247,248,249,250,251,253,254,256,257,258,259,260,261,262,263, 264,265,266,267,268,269,271,272,274,283,284,287]



# Define parameters to fix in the fitting
# The indices below correspond to those listed in the spectral_results.txt file
# If we want to fix parameters: fix = {"nh": {"indices": "all", "value": 0.1989}, "gamma": {"indices": [-4, -3, -2, -1], "value": 1.7}, "Tin": {"indices": [-4, -3, -2, -1], "value": 0.5} }
# We will want to fix parameters for low_count_indexes, defined above, e.g. fix = {"nh": {"indices": "all", "value": 0.2}, "gamma": {"indices": low_count_indexes, "value": 1.62}, "Tin": {"indices": low_count_indexes, "value": 0.5}}
# NOTE: For nH, I don't do "all" as it is giving issues for some of the SS epochs
## TODO: I may want to fix nH.... especially for the HS epochs with low counts??

nH = 0.6
gamma1 = 1.55
gamma2 = 1.7
Tin = 0.5

#fix = {"nh": {"indices": low_count_indexes, "value": nH}, "gamma": {"indices": low_count_indexes_1, "value": gamma1}, "gamma": {"indices": low_count_indexes_2, "value": gamma2}, "Tin": {"indices": low_count_indexes, "value": Tin}}
fix = {
    "nh": {"indices": low_count_indexes, "value": nH},
    "gamma": [  
        {"indices": low_count_indexes_1, "value": gamma1},
        {"indices": low_count_indexes_2, "value": gamma2}
    ],
    "Tin": {"indices": low_count_indexes, "value": Tin}
}






######################################################################################################################

## F-TEST

fixing = True
simple_model = 'pegged_powerlaw'
complex_model = 'powerlaw+diskbb'



######################################################################################################################

## INPUTS --- GET UPPER LIMIT FLUXES -- get_uplims
# No further inputs required
# This step uses target_coords, the nH and gamma values above, and the upper limit inputs from step 2.


## STEP 6 RESULTS
import numpy as np
gamma = 1.55
uplims_fluxes = [np.nan] * len(uplims_IDs) 


######################################################################################################################

## INPUTS --- FINAL RESULTS -- get_final_results

#models_indexes =[ [[0,81], [95,140], [187,292]], [], [], [[82,94], [141, 186] ] ]
#models_indexes =[ [[0,81], [97,138], [188,292]], [ [90,96], [139,143], [186,187]], [[82,89], [144, 185] ] ]
models_indexes =[ [[0,81], [95,138], [188,292]], [ [90,94], [139,143], [186,187]], [[82,89], [144, 185] ] ]


######################################################################################################################

## INPUTS FOR HELPER FUNCTIONS

## For get_index_range
MJD_range = []

######################################################################################################################

