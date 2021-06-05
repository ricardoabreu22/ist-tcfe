%Lab5                                                  Academic Year 2020/2021
%                                                            TCFE
%Authors:
%Joana Matos, 95799
%Ricardo Abreu, 95842
%Vasco Emídio, 95856
%MEAer
close all 
clear all

format long;

%----Data---
C1=220e-6;
C2=220e-6;
R1=1e3;
R2=10e3;
R3=100e3;
R4a=1e3;
R4b=1e3;

 printf("Data_TAB \n");
 printf("$R_{1}$ = %f Ohm\n", R1); 
 printf("$R_{2}$ = %f Ohm \n", R2);
 printf("$R_{3}$ = %f Ohm \n", R3);
 printf("$R_{4a}$ = %f Ohm \n", R4a);
 printf("$R_{4b}$ = %e Ohm \n", R4b);
 printf("$C_{1}$ = %e F \n", C1);
 printf("$C_{2}$ = %e F \n", C2);
 printf("Data_END \n \n");
%------------
%--Replace R4a and R4b by an equivalent R4
R4=1/((1/R4a)+(1/Rab));
%------------------------------------------

%----1. Gain, Zi, Zout---------------------
fo=1e3;
wL = 1/(R1*C1);
wH = 1/(R2*C2);
wO = sqrt(wL*wH);

%---------------------------------------------

%----2. --------------------------------------


%--------------------------------------------

