data CytokinesT1;
input Raw_P;
datalines;
0.0063
0.0046
0.1097
0.0017
0.0190
0.0025
0.8031
;

ods graphics on;
proc multtest inpvalues=cytokinesT1 fdr pfdr afdr plots=all;
run;

data CytokinesT2;
input Raw_P;
datalines;
0.3913
0.0124
0.2928
0.1349
0.2515
0.0839
0.7543
;
ods graphics on;
proc multtest inpvalues=cytokinesT2 fdr pfdr afdr plots=all;
run;

data VTROI;
input Raw_P;
datalines;
0.0030
0.0104
0.9732
0.4175
0.0604
0.8024
0.0749
0.1635
0.5736
0.3671
0.5553
0.3952
0.0810
0.8423
;
ods graphics on;
proc multtest inpvalues=VTROI fdr pfdr afdr plots=all;
run;

data K1ROI;
input Raw_P;
datalines;
0.0008
0.0013
0.0215
0.0246
0.0003
0.0005
0.0162
0.0422
0.0343
0.0105
0.0026
0.0010
0.0350
0.0699
;
ods graphics on;
proc multtest inpvalues=K1ROI fdr pfdr afdr plots=all;
run;

data SCFAs;
input Raw_P;
datalines;
0.6683
0.8466
0.1787
0.2235
;
ods graphics on;
proc multtest inpvalues=SCFAs fdr pfdr afdr plots=all;
run;