In the Data_preprocessing folder: 
==
The data doesn't work very well when we set the tubal-rank too low during data preprocessing. <br>
When tubal-rank is 15, we can see that the rank approximation result has some new noise.  <br>
When tubal-rank is 32, we can get a good approximation result. <br>


**The CDF of singular values for the seismic data, compared with processed data.**
![](https://github.com/hust512/Seismic_Sensory_Data_Analysis/blob/master/Data_preprocessing/原始数据与预处理数据的CDF图对比.png)


**Original data.**
![](https://github.com/hust512/Seismic_Sensory_Data_Analysis/blob/master/Data_preprocessing/original_data.png)

**Processed data.**
![](https://github.com/hust512/Seismic_Sensory_Data_Analysis/blob/master/Data_preprocessing/tubal_rank_32_approximate_data.png)

These figures are obtained by running the test_r_error.m in the Data_preprocessing folder.


Result:
==

The result of TNN algorithm and Tubal-Alt-Min algorithm:
--
**data size: t * m * n: 300 * 120 * 80, tubal-rank:15.**<br>
**RSE of TNN:6.1e-03. RSE of Tubal-Alt-Min:1.56e-02.**<br>
![](https://github.com/hust512/Seismic_Sensory_Data_Analysis/blob/master/Result/slicemissingrecovery.png)

(a) Slice view of the original complete data. <br>
(b) Slice view of reconstructed data using TNN. <br>
(c) Slice view of reconstructed data using Tubal-Alt-Min. <br>

run Tubal_Alt_Min.m to get the result.

The reconstruction error of TNN algorithm and Tubal-Alt-Min algorithm:
--
![](https://github.com/hust512/Seismic_Sensory_Data_Analysis/blob/master/Result/error_v3.png)

The figure shows the RSE of the Tubal-Alt-Min algorithm and TNN algorithm for varying frontal slice sampling rates from 80% to 98%.<br>

run Tubal_TNN.m in the error folder to get the figure

Two m_files in the error folder:
==
1. Tubal_Alt_Min_error.m is used to calculate RSE in the case of different slice missing. <br>
2. Tubal_TNN.m compared the recovery accuracy and convergence speed between Tubal_Alt_Min algorithm and TNN algorithm. <br>
3. The TNN_solver folder is some functions of the TNN algorithm. <br>

In the Plot_CDF folder:
==
the function of the test low-tubal-rank character of the seismic data in three different dimensions. <br>


**T_synthetic_tubal_rank_2.mat is a synthetic seismic data set of size: 64 * 64 * 256. Three dimensions are inline,crossline,time respectively.**<br>

**The experimental data of size 326 * 431 * 531. The three dimensions are: time,inline and crossline. tubal-rank: 45.**<br>

**The SeisPlot.zip is the toolbox to plot the seismic data. You should add path in MATLAB after decompression** <br>


Tubal-sampling result
==
**data size:64 * 64 * 256, tubalrank:2.**<br>
**tubal-sampling : 50%.**<br>

Tubal_Alt_Min algorithm result and TNN algorithm result:
--
**Tubal_Alt_Min RSE = 1.289423e-15.TNN RSE = 1.475729e-04．**<br>

corrupted data<br>
![](https://github.com/hust512/Seismic_Sensory_Data_Analysis/blob/master/Result/Corrupted_data.png)

Tubal_Alt_Min recovery data<br>
![](https://github.com/hust512/Seismic_Sensory_Data_Analysis/blob/master/Result/Tubal_Recovery_data.png)

TNN recovery data<br>
![](https://github.com/hust512/Seismic_Sensory_Data_Analysis/blob/master/Result/TNN_recovery_data.png)

run Tubal_Alt_Min_TNN_tubal_sampling.m to get the result.<br>

reconstruction error (RSE) under different tubal sampling rates:<br>
![](https://github.com/hust512/Seismic_Sensory_Data_Analysis/blob/master/Result/tubal_sampling_error.png)

run tubal_sampling_error.m to get the result.<br>


GCG algorithm result：
--
**GCG RSE:4.907527e-03.**<br>
GCG recovery data<br>
![](https://github.com/hust512/Seismic_Sensory_Data_Analysis/blob/master/GCGalgorithm/Result/Recovery_data.png)

run gcgexample.m in the GCGalgorithm folder to get the result.<br>

conclusion:
--
tubal sampling: 70%. Tubal_Alt_Min RSE:5.504735e-16. TNN RSE:7.84263e-05.  GCG RSE: 2.165258e-03.<br>
tubal sampling: 30%. Tubal_Alt_Min RSE:4.928588e-07. TNN RSE:3.27454e-03.  GCG RSE: 3.460236e-02.<br>
