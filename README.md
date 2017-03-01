Here are some of the scripts I use during my PhD. Most of them are written in Python 2.7 and non numerically optimized.

Among other things, there are PDF calculations, statistical analyses of coordination networks in ZIFs, conversion scripts for CRYSTAL14, scripts that clean the output of CP2K (xyz, frc, vel, ener)

- CNH2Im : xyz (C,N,H,Zn) --> xyz (Zn, Im)

- CNH2Im_2250 : CNH2Im, in case of imidazolate breakage

- CNH2mIm : almost identical to CNH2Im but made for methylimidazolate instead of imidazolate

- Cij2elas : using ELATE (progs.coudert.name/elate) to convert a Cij matrix (second order elastic constants) to sensible mechanical properties

- Imorient : Computing mean square angular deviation of Im or mIm (along N-N or along C-(middle of C-C))

- PDF : Computing pair distribution functions with atomic factors at 59 Kev --> to compare with X-ray experiments (see http://scripts.iucr.org/cgi-bin/paper?th0051)

- PDFtoFQ : Convert PDF to a F(Q) (see http://scripts.iucr.org/cgi-bin/paper?th0051)

- PMF : Computing pair correlation function (g(r)) ; PMF stands for potential of mean force given by -kTlog(g(r))

- RMSD : Computing mean square displacement for any given element (used for Zn and X (representing Im or mIm))

- SuperCell : xyz --> xyz expanded periodically along the 3 axes by integer factors

- Tang_avg : Tetrahedral angle (N-Zn-N))

- Tang_traj : Tetrahedral angle during exchanges

- ZnImang_avg : Zn-Im-Zn angle

- ZnImang_traj : Zn-Im-Zn angle during exchanges

- avg_gr : Averages output of PDF or PMF 

- center_xyz : Removes the drift from molecular simulation xyz outputs

- cif_coord : cif --> xyz

- cif_d12 : cif --> d12 (CRYSTAL14 input)

- clean_cp2k : for ener, frc, vel CP2K outputs --> removes successive identical steps that might have been put there during restarts ; for xyz, as for the others and wraps in the unit cell

- crystal2Cij : CRYSTAL14 output to Cij matrix (result can be further analyzed by ELATE at progs.coudert.name/elate)

- crystal2cif : CRYSTAL14 output to cif file with several structures (opt steps) written by FX Coudert (coudert.name)

- distance_analysis : xyz --> dist : distance file containing the distance between a given pair of elements (mainly used for N-Zn pairs)

- energies_Si : Energy per Si from CRYSTAL14 output

- Cv : heat capacity

- Udist : Energy distribution

- U2Cv : Cv from U

- cleaform_traj : all cleavages and formations averaged

- coord_dist : coordination numbers distribution

- coord_stats : coordination numbers

- ex_traj : exchange trajectories

- timing_stats : cleavage frequencies...

- wtimes : waiting times before succesful exchange distribution
