# Code for 'An evaluation of influenza season onset prediction methods'
## William Shepard

The code written is for the poster 'An evaluation of influenza season onset prediction methods', as part of the 2025/2026 Mathematics and Statistics Vacation Scholarship Program at the University of Melbourne. There are three R files
- code_final_results: runs all the code necessary to reproduce the results found in the poster
- preprocessing and model_generation: sepearte functions to run the described method in the poster for new data

In order to run the files:
1. Download current fluID and fluNet data from https://www.who.int/teams/global-influenza-programme/surveillance-and-monitoring/influenza-surveillance-outputs and place in the 'Data' folder
2. Run code_final_results for replaicated resutls. Note: this will take around 15-20 minutes due the ammount of models needed to be generated
3. Run preprocessing and then model_genearation to test on current data. Modify the country name and start and end dates depending on what you are analysing.
