PARAMETER ESTIMATION

1) In order to estimate the parameters of the pharmacokinetics and pharmacodynamics, run the program “Estimate_BG_finger_tapping_random.m”. It is worth noting that the results can vary significantly from one trial to the next, due to noise variability.

Estimate_BG_finger_tapping_random.m: this program performs the estimation of the parameters using one initial guess for the pharmacokinetics and Np (at present  Np = 10, line 59) random combinations of the parameters as the initial guess for pharmacodynamics. The synapse parameters are loaded from the data file "W_tot_new_W0e5_D1e0.mat" file (line 12), which contains the synapses obtained after various periods of training, starting from a naïve network (as shown in Figs. 2 and 3 of Ursino and Baston, EJN 2018). At present the program uses synapses after 100 training periods (lines 13-16), considered as typical of a PD subject with moderate skill performance. 
Then, the program asks for the name of the patient data file (lines 22-24).
Subsequently the file "curve_tapping_3.mat" is loaded, which contains the curve "finger tapping rate vs. concentration of levodopa" obtained previously with the same synapses of lines12-16. The inversion of this curve is used to find the first initial guess of parameter D0.
To estimate the parameters of the pharmacokinetics, the function "Cost_levodopa.m" is minimized (row 82): that function calculates the sum of the squares of the difference between plasma concentration in the model and in the patient.
Then Np estimates of the parameters of the pharmacodynamics are carried out (using the cycle “for” starting at line 123), with Np different first attempt values ​​of the parameters generated randomly (lines 124-128) . The estimation works by minimizing the function "Cost_tapping.m". This function calculates a least squares cost function + a term that takes into account the maximum error (line 55). In turn, this function calls the function "BG_model_function_tapping_mauro_3 .m" which performs the simulation of finger tapping with the neurocomputational model.
The estimated parameters are saved at the end in a file "parametri_paziente.mat" file (line 197). It is recommended to change the name of this file each time you move from one patient to the next.

2) To graph the results of the previous estimates, you can call the program  “Plot_many_random_estimates.m". The data of patients used in the paper are contained in files named ‘pazienteGxSx.mat’ and the corresponding optimal parameters are contained in the files named ‘parametri_pazienteGxSxx_max.mat’. The program looks at the 10 different cost functions, and perform all graphs.


Plot_many_random_estimates.m: the program uploads the synapse file "W_tot_new_W0e5_D1e0.mat" (lines 12-187; asks for the name of the patient's file and the name of the file where estimated parameters are stored (lines 22-27); calls the function "BG_model_function_tapping_mauro_3.m" to perform the finger tapping simulation and creates 11 graphs: in the first there is the comparison between patient and simulated plasma concentration; in the subsequent 10 graphs there is the comparison between the patient's finger tapping frequency and that simulated with the model for each of the 10 different estimates of pharmacodynamics parameters.







GRAPHICS OF OPTIMAL RESULTS

To graph the results obtained with the optimal values of the parameters (these parameters include the Np different possibilities estimated above) you can use the file "Generic_patient.m".
The data of patients used in the paper are contained in files named ‘pazienteGxSx’ and the corresponding optimal parameters are contained in the files named ‘parametri_pazienteGxSxx_max’. The program looks at the 10 different cost functions, and choose the combination with minimal cost-function to make the graphs. 

Generic_patient.m: the program asks for the patient's name and the file name with the estimated values of the parameters (lines 22-27), and looks at the combination of parameters with minimal cost function (line 107). Then it simulates the patient by calling the function "BG_model_function_tapping_mauro_3.m" to perform the finger tapping simulation. Then it creates 2 graphs: in the first there is the comparison between plasma concentration in the patient and that simulated with the model; in the second there is the comparison between the patient's finger tapping frequency and that simulated with the model  with the optimal parameter values.


TRAINING THE PATIENT

To train a patient at different minutes after LD administration, run the program  “Train_the_patient.m”. It is worth noting that the results can vary significantly from one trial to the next, due to noise variability.

Train_the_patient.m: this program reads the patient's parameters
(lines 25-55 in which 5 cases are presented), and performs  first a finger tapping test without training, and then repeat the test, by carrying out training at a predetermined time. The training instant is assigned in line 109.
The program calls the program "calculate_levodopa.m" to perform the simulation of the levodopa plasma concentration.
To perform the training, at the assigned moment the program calls another program named "Synapse_training.m" (line 120). This program, in turn, calls the function "BG_model_function_Ach.m"
to perform a single training test (i.e. a single finger movement followed by reward or punishment)
repeated for N epochs (N_epochs is assigned in line 58).
To perform finger tapping, the function "BG_model_function_tapping_mauro_3.m" is called, which reads the value of the dopaminergic DA input.
