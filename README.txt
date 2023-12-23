This is the readme file for the model associated with the following paper:
      Kim SY, Lim W (2021) Effect of Diverse Recoding of Granule Cells on Optokinetic Response in A 
      Cerebellar Ring Network with Synaptic Plasticity. Neural Networks 134: 173-204

Implementation of simulation work in the above paper was done by Kim and Lim, to whom questions 
should be addressed; sykim@icn.re.kr or wclim@icn.re.kr

Our cerebellar spiking neural network is well shown in Fig. 2 of the above paper. The following source 
codes for simulation in the above paper were written in the C programming language:
      (1) raster_plot.c: code to get the raster plots of spikes for the constituent cells [Figs. 4(a), 5, 7(a), 
          8(a), 9, 11(a), 12(a), 13(b), 14(a), 17(a)]
      (2) IPSR.c: code to get the IPSR (instantaneous population spike rate) for the constituent cells 
          [Figs. 4(b), 5, 11(b)-(d2), 14(b), 15(c)-(d), 17(b)]
      (3) CI.c: code to get the conjunction index (CI) [Figs. 6, 14(c)-(e), 17(c)]

Refer to Tables A.1 – A.4 for the system parameter values. There was one variable: connection  
probability p_c from the Golgi cells to the granule cells; we mainly considered the optimal case of 
p_c=0.06 where the spiking patterns of the granule cells are the most diverse. 

Then, using the GCC complier we ran the above source codes on personal computers with Intel i5-
10210U CPUs (1.6 GHz) and 8 GB of RAM. The number of used personal computers vary (from 1 to 70) 
depending on the type of jobs. 

Using the OriginPro 2016 (Graphing & Analysis Software), we handled the output data, and all the figures 
in the above paper were made.  
