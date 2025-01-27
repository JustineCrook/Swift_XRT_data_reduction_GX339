
"""
NOTES FOR MAXI J1820+070:
- Target IDs: 00010627,00010754,00010774,00014223,00015032,00088657,00088779,00813771,00814259,00815603
- 00010627132wt is an upper limit, but there is a pc observation on the same date
"""


######################################################################################################################

## FUNCTION TO RUN
analysis_type = "step4"
 

######################################################################################################################

## STEP 1 INPUTS --- DATA RETRIEVAL

# Target coordinates
target_coords = [275.0914, 7.1854]

# 00010627
# target_names = ['MAXI J1820+070']
# target_ids = ['00010627']
#segments_lc   = [[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,27,28,29,30,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,53,54,55,56,57,58,60,62,63,65,66,67,68,69,70,72,73,74,76,77,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,100,101,102,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,139,140,141,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188,189,190,191,192,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,224,225,226,227,228,229]] 
#NOT USED: segments_spec = [[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,26,27,28,29,30,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,53,54,55,56,57,58,60,62,63,65,66,67,68,69,70,72,73,74,76,77,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97]] 
#segments_spec = [[98,100,101,102,104,105,106,107,108,109,110,111,112]]
#segments_spec = [[ 113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,139,140,141,143,144,145,146,147,148,149,150]]
#segments_spec =  [[ 151,152,153,154,155,156,157,158,159,160,161,162,163,164,166,167,168,169,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,187,188]]
#segments_spec =  [[ 189,190,191,192,194,195,196,197,198,199,200,201,202,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,224,225,226,227,228,229]]

# 00014223,00015032,00088657
target_names = ['MAXI J1820+070', 'MAXI J1820+070','MAXI_J1820p070']
target_ids = ['00014223', '00015032', '00088657']
segments_lc = [[1,3,4,5,6,7,8,9,10,11,12,13,15,14,16,17,19,18,24,25,26,27,28,29,30,31,32,33,34,35,37,38,39,40,41,42,43,44,46,47,48,49,50,51,52,53], [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], [10,11]]
segments_spec = segments_lc

######################################################################################################################

## STEP 2 
# No inputs required


## STEP 2 RESULTS -- FROM LIGHT CURVE ROUTINE
# After running the light curve routine, the code outputs which observations are upper limits.
# These should not be fitted, so are input to the spectral fitting function (steps 3 and 5).
# The code also outputs the count rate, used when converting count rate to flux (step 6).
# Note that the upper limits from the pipeline are at the 3-σ level, calculated using the Bayesian method (see Kraft, Burrows & Nousek, 1991, ApJ, 374, 344).
uplims_IDs = ['00010627132wt', '00010627194wt', '00010627195wt','00010627199pc','00010627214wt','00010627215wt','00010627216wt','00010627217wt','00010627218wt','00010627219pc','00010627222pc', '00014223005pc', '00014223011pc', '00014223027pc', '00014223035pc', '00014223037wt', '00014223041pc', '00014223042pc', '00015032005pc', '00015032006pc', '00015032014pc']
uplims_count_rate_ar = [0.36302 ,0.030703,0.06003 ,0.1088  ,0.017791,0.022321,0.034858,0.04014 ,0.039144,0.016746,0.015532, 0.012266,0.010932,0.01034 ,0.019378,0.60859 ,0.022678,0.065961,0.031431,0.026654,0.043244]
uplims_MJDs=  [58520.4826323133,58761.2063995436,58768.0468569869,58789.0998317126,58951.0876508723, 58958.0575253364,58965.0298386431,58972.5327235123,58980.841630267 ,59050.9571998144, 59053.9397797577,59305.3464732519,59308.6539219533,59314.8322711651,59318.6799727473, 59319.0694681906,59321.6755077681,59331.7556949267,59735.5479818653,59743.712824216 , 59784.4907593387]
uplims_MJDs_er = [1.293086111111112e-01,1.024032407407408e-02,1.122268518518518e-02,1.102754629629630e-03, 1.156958333333334e-02,1.075951388888888e-02,1.197460648148148e-02,1.029657407407408e-02, 1.564031944444444e-01,7.167870370370380e-03,1.149178240740740e-02,9.198101851851860e-03, 1.096081018518518e-02,1.116064814814814e-02,1.575245370370370e-02,1.216898148148148e-04, 8.020023148148141e-03,1.828240740740740e-03,4.730208333333340e-03,3.972305092592600e-01, 3.859606481481480e-03]


######################################################################################################################

## STEP 3 INPUTS --- INITIAL MODEL FITTING
# Uses uplims_IDs defined above

# Define models to fit -- used for fit_spec
models = ['pegged_powerlaw', 'powerlaw+diskbb'] 
low_energy_fit_bound_keV = 0.5 # must be 0.5 <= x <= 1.0 


## STEP 3 RESULTS --- FROM MODEL FITTING RESULTS
# The spectral_results.txt file lists the low_count_indexes, which correspond to observations that cannot initial be fit (due to low count rate), so require some fitting parameters to be fixed.
low_count_indexes = [ 33, 34, 35, 36, 65, 71, 72, 73, 74, 75, 76, 77, 78, 79, 91, 97, 98, 99,100,101,113,114,115,116, 117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140, 141,142,143,144,145,146,147,149,150,152,153,154,156,157,158,159,160,161,162,163,164,165,166,167, 168,169,170,171,172]


######################################################################################################################

## STEP 4 INPUTS --- GET PARAMETER VALUES (N_H, GAMMA, T_IN)

# Based on the spectral fit results (step 3) and any prior knowledge regarding spectral state ranges, define which spectral model to use for each point in time.
# model_indexes is a multi-dimensional array. Each outer element corresponds to the models defined above. The inner elements are ranges for each model.
# The index ranges are inclusive (i.e. [start, end]), and correspond to those listed in the spectral_results.txt file.
# e.g. models_indexes=[ [[6, 17], [19, 21]] , [[0, 5]] , [[]] ] means the first model should be used for the observations corresponding to indexes 6-17 and 19-21, the second model should be used for observations 0-5, and the last model that was fit should not be used
# If none, use models_indexes = []
models_indexes=[ [[8,172]], [[0,7]]]


######################################################################################################################

## STEP 5 INPUTS --- MODEL FITTING WITH SOME FIXED PARAMETERS
# Re-run spectral fitting, but this time we want to fit the parameters for some observations -- using the results of step 4 and/or any prior knowledge (e.g. nH values in papers).
# Uses uplims_IDs, models, low_energy_fit_bound_keV, and low_count_indexes -- defined above

# Define parameters to fix in the fitting
# The indices below correspond to those listed in the spectral_results.txt file
# If we want to fix parameters: fix = {"nh": {"indices": "all", "value": 0.1989}, "gamma": {"indices": [-4, -3, -2, -1], "value": 1.7}, "Tin": {"indices": [-4, -3, -2, -1], "value": 0.5} }
# We will want to fix parameters for low_count_indexes, defined above, e.g. fix = {"nh": {"indices": "all", "value": 0.2}, "gamma": {"indices": low_count_indexes, "value": 1.62}, "Tin": {"indices": low_count_indexes, "value": 0.5}}
nH = 0.09100 # 0.2175 for the 1-10keV fit
gamma = 1.6153 # 1.6880 for the 1-10keV fit
Tin = 0.5
fix = {"nh": {"indices": "all", "value": nH}, "gamma": {"indices": low_count_indexes, "value": gamma}, "Tin": {"indices": low_count_indexes, "value": Tin}}


######################################################################################################################

## STEP 6 INPUTS --- GET UPPER LIMIT FLUXES
# No further inputs required
# This step uses target_coords, the nH and gamma values above, and the upper limit inputs from step 2.


## STEP 6 RESULTS
uplims_fluxes = [5.3173796344917639e-11, 2.0000729759577985e-12, 2.3756583459663299e-12, 4.1494539494675587e-12, 1.8295031357803069e-12 ,1.8524819073723543e-12, 1.8778650910737478e-12, 1.9718023476986085e-12, 1.8917034869406333e-12 ,7.9991742280550490e-13, 4.7276899557394712e-13 ,7.5151876043236059e-13, 1.0176161157365480e-12 ,6.9557222471669509e-13 ,1.4958790969538238e-12 ,3.1588648376955921e-10, 9.8239681080999278e-13, 2.9703151345188235e-12, 1.2038759554012382e-12 ,6.9462315127741660e-13, 1.6929224952518906e-12]


######################################################################################################################

## STEP 7 INPUTS --- FINAL RESULTS
# No further inputs required


######################################################################################################################

## INPUTS FOR HELPER FUNCTIONS

## For get_index_range
MJD_range = []

######################################################################################################################

