# Bachelor End Project on assessing forecast comparison tests under proper scoring rules

The files that are used in the simulations is adapted from the publicly accessible repository of Lerch et al., 2020, https://doi.org/10.5194/npg-2019-62. (Github repository: https://github.com/slerch/multiv_pp). The copula extensions and artificial data from M.G.J. Flos et al. (2025) and M. G. J. Flos (2022) are accessible at https://github.com/mgjfl/BEP-BAM and https://github.com/elisaperrone/COBASE.

## References:
- Lerch, S., Baran, S., Möller, A., Groß, J., Schefzik, R., Hemri, S., & Gräter, M. (2020). Simulation-based comparison of multivariate ensemble post-processing methods. Nonlinear Processes in Geophysics Discussions. https://doi.org/10.5194/npg-2019-62
- Flos, M.G.J., Fran¸cois, B., Schicker, I., Whan, K., & Perrone, E. (2025). COBASE: A new copula-based shuffling method for ensemble weather forecast postprocessing. http://arxiv.org/abs/2510.25610
- Flos, M. G. J. (2022). Copula-based statistical postprocessing for weather data (tech. rep.). Eindhoven University of Technology. https://research.tue.nl/en/studentTheses/copula-based-statistical-postprocessing-for-weather-data/

# Code Structure
```python

BAM_BEP_DM_test/
├── Artificial data expiriments/           # R Project file
│    ├── data provided by paper COBASE/    # artificial dataset of Flos, M.G.J. et al. (2025)
|    │   └── Mock_data.csv
│    ├── data files/                       # data files can be stored in this map
|    │   
│    ├── run_setting_mock.R                  # forecast generation of artificial data
|    ├── run_setting_mock_not_timewindow.R    # forecast generation of artificial data without time window
|    │   
│    ├── Multivariate DM/                     # multivariate DM statistic values calculations
|    │   └── multivariateDM_mock.R  
│    └── univariate DM/                      # univariate DM statistic values calculations
|        └── DM_bootstrapping_mock_data.R
├── Simulation code/      # Simulation setting files  
|    ├── run_setting_archimedean_shuffle_bep.R
│    ├── Simulation data/
|    │   ├── Rdata_dir simulation
|    │   └── Rout_dir simulation
│    ├── source code/
|    │   ├── evaluation_functions.R
|    │   ├── generate_ensfc.R
|    │   ├── generate_observations.R
|    │   ├── mvpp_arch_shuffle.R
|    │   ├── mvpp_mock.R
|    │   ├── postprocess_ensfc_arch.R
|    │   └── postprocess_ensfc_mock.R
│    ├── Multivariate DM code/
|    │   ├── nieuwe MDM selection.R
|    |   └── MDM_scoretest.R
│    └── Univariate DM code/
|        └── DM_TestStatistic_computation.R
|
├── Analysis files/             # Code to generate figures for answering the Research Questions (RQ)
│    ├── Plot generation/      
|    │   └── Boxplots for RQ 1.R
│    ├── RQ1/
|    │   ├── DM_q1_onderzoek.R
|    │   └── Q1_for_mock_research.R
│    ├── RQ2/
|    │   └── MDM_and_agreement_analyses.R
│    └── RQ3/
|        └── power_and_size_calculations_for_MDM.R
|
└── README.md
```
