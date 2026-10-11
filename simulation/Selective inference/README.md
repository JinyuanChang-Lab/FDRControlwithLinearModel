## To get the results in the paper:

* This section provides results of Figures F31-F36 in the Supplementary Material.

* Note for Selective inference: The simulation scripts for this section are located in their respective method subfolders (e.g., knockoff_recycle, knockoff_split, OLS_BH, Our_methods, PSI_BH, and tian_BH), while the plotting scripts (SI_plot_S1.R and SI_plot_S2.R) and their corresponding output folders are kept in the main directory.

* Run `S1_simulation_knockoff_recycle.R`, `S1_simulation_knockoff_split.R`, `S1_simulation_OLS_BH.R`, `S1_simulation_SI.R`, `S1_simulation_PSI_BH.R`, and `S1_simulation_Tian_BH.R` to get the results of the knockoff recycle method using the square-root Lasso, the knockoff split method using the square-root Lasso, the OLS method combined with the BH procedure, our proposed method, the post-selection inference method combined with the BH procedure, and Tian's method based on square-root Lasso combined with the BH procedure for the global null setting.

* Run `S2_simulation_knockoff_recycle.R`, `S2_simulation_knockoff_split.R`, `S2_simulation_OLS_BH.R`, `S2_simulation_SI.R`, `S2_simulation_PSI_BH.R`, and `S2_simulation_Tian_BH.R` to get the results of the corresponding methods for the power comparison setting.

* Put simulation results of all the S1 scripts above in the folder `SI_plot_S1`.

* Put simulation results of all the S2 scripts above in the folder `SI_plot_S2`.

* Run `SI_plot_S1.R` to generate all subfigures of Figures F31 and F32.

* Run `SI_plot_S2.R` to generate all subfigures of Figures F33-F36.

* The folder `SI_plot_S1` contains all subfigures of Figures F31 and F32. 

* The folder `SI_plot_S2` contains all subfigures of Figures F33-F36.