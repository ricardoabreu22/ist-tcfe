%Lab4                                                  Academic Year 2020/2021
%                                                            TCFE
%Authors:
%Joana Matos, 95799
%Ricardo Abreu, 95842
%Vasco Emídio, 95856
%MEAer
close all 
clear all

format long;

f = 100*1e3;
w = 2*pi*f;

%--------------------------------Gain Stage-------------------------------------
%-----Data-----
VCC=12;
RS=100;
RB1=80000;
RB2=20000;
RC1=1000;
RE1=100;
CS=10*1e-3;
Cb=10*1e-3;

 printf("GSData_TAB \n"); 
 printf("$V_{CC}$ = %e V\n", VCC);
 printf("$R_{in}$ = %e Ohm\n", RS); 
 printf("$R_{1}$ = %e Ohm \n", RB1);
 printf("$R_{2}$ = %e Ohm \n", RB2);
 printf("$R_{c}$ = %e Ohm \n", RC1);
 printf("$R_{e}$ = %e Ohm \n", RE1);
 printf("$C_{in}$ = %e F \n", CS);
 printf("$C_{b}$ = %e F \n", Cb);
 printf("GSData_END \n \n");

VT=25e-3;           
BFN=178.7;         %Bijector
VAFN=69.7;
VBEON=0.7;

 printf("GSBijector_TAB \n"); 
 printf("$V_{T}$ = %f V\n", VT);
 printf("$beta$ = %f \n", BFN); 
 printf("$V_{A}$ = %f V\n", VAFN);
 printf("$V_{BEON}$ = %f V \n", VBEON);
 printf("GSBijector_END \n \n");

%-----Operating Point Analysis-----
RB=1/(1/RB1+1/RB2);                 %Req   <---- Thevenin's equivalent
VEQ=RB2/(RB1+RB2)*VCC;              %Veq   <----
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1);   %IB
IC1=BFN*IB1;                        %IC
IE1=(1+BFN)*IB1;                    %IE
VE1=RE1*IE1;                        %VE
VO1=VCC-RC1*IC1;                    %VO
VCE=VO1-VE1;                        %VCE
 printf("GSOP_TAB \n"); 
 printf("$R_{B}$ = %e Ohm\n", RB);
 printf("$V_{eq}$ = %e V\n", VEQ); 
 printf("$I_{B1}$ = %e A \n", IB1);
 printf("$I_{C1}$ = %e A \n", IC1);
 printf("$I_{E1}$ = %e A \n", IE1);
 printf("$V_{E1}$ = %e V \n", VE1);
 printf("$V_{O1}$ = %e V \n", VO1);
 printf("$V_{CE}$ = %e V \n", VCE);
 printf("GSOP_END \n \n");


%-----Incremental Circuit-----            (for medium frequencies)

gm1=IC1/VT                         %gm    --Incremental Parameters--
rpi1=BFN/gm1                       %rpi
ro1=VAFN/IC1                       %ro

RSB=RB*RS/(RB+RS)                  %RS||RB

%AC Analysis (Incremantal circuit gain)
AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)   %vo/vi
AVI_DB = 20*log10(abs(AV1))                                                                     %vo/vi(dB)
AV1simple = RB/(RB+RS) * gm1*RC1/(1+gm1*RE1)
AVIsimple_DB = 20*log10(abs(AV1simple))

RE1=0
AV1 = RSB/RS * RC1*(RE1-gm1*rpi1*ro1)/((ro1+RC1+RE1)*(RSB+rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)
AVI_DB = 20*log10(abs(AV1))
AV1simple =  - RSB/RS * gm1*RC1/(1+gm1*RE1)
AVIsimple_DB = 20*log10(abs(AV1simple))

RE1=100
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZX = ro1*((RSB+rpi1)*RE1/(RSB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RSB)+1/RE1+gm1*rpi1/(rpi1+RSB)))
ZX = ro1*(   1/RE1+1/(rpi1+RSB)+1/ro1+gm1*rpi1/(rpi1+RSB)  )/(   1/RE1+1/(rpi1+RSB) ) 
ZO1 = 1/(1/ZX+1/RC1)

RE1=0
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZO1 = 1/(1/ro1+1/RC1)
%-------------------------------------------------------------------------------

%---------------------------Ouput Stage-----------------------------------------
BFP = 227.3
VAFP = 37.2
RE2 = 100
VEBON = 0.7
VI2 = VO1
IE2 = (VCC-VEBON-VI2)/RE2
IC2 = BFP/(BFP+1)*IE2
VO2 = VCC - RE2*IE2

gm2 = IC2/VT
go2 = IC2/VAFP
gpi2 = gm2/BFP
ge2 = 1/RE2

AV2 = gm2/(gm2+gpi2+go2+ge2)
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)
ZO2 = 1/(gm2+gpi2+go2+ge2)
%-------------------------------------------------------------------------------

%------------------------------Total--------------------------------------------
gB = 1/(1/gpi2+ZO1)
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1
AV_DB = 20*log10(abs(AV))
ZI=ZI1
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)
