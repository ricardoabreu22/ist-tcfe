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

%--------------------------------Gain Stage-------------------------------------
%-----Data-----
VCC=12;
RS=100;
RB1=80000;
RB2=20000;
RC1=940;
RE1=775;
RE2 = 2335;
CS=690*1e-6;
Cb=4180*1e-6;
Co=2250 *1e-6;
RL=8;
Vinput=1;

 printf("GSData_TAB \n");
 printf("$R_{S}$ = %e Ohm\n", RS); 
 printf("$R_{B1}$ = %e Ohm \n", RB1);
 printf("$R_{B2}$ = %e Ohm \n", RB2);
 printf("$R_{C1}$ = %e Ohm \n", RC1);
 printf("$R_{E1}$ = %e Ohm \n", RE1);
 printf("$R_{E2}$ = %e Ohm \n", RE2);
 printf("$R_{L}$ = %e Ohm \n", RL);
 printf("$C_{I}$ = %e F \n", CS);
 printf("$C_{b}$ = %e F \n", Cb);
 printf("$C_{O}$ = %e F \n", Co);
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

RE1=775
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZX = ro1*((RSB+rpi1)*RE1/(RSB+rpi1+RE1))/(1/(1/ro1+1/(rpi1+RSB)+1/RE1+gm1*rpi1/(rpi1+RSB)))
ZX = ro1*(   1/RE1+1/(rpi1+RSB)+1/ro1+gm1*rpi1/(rpi1+RSB)  )/(   1/RE1+1/(rpi1+RSB) ) 
ZO1 = 1/(1/ZX+1/RC1)

RE1=0
ZI1 = 1/(1/RB+1/(((ro1+RC1+RE1)*(rpi1+RE1)+gm1*RE1*ro1*rpi1 - RE1^2)/(ro1+RC1+RE1)))
ZO1 = 1/(1/ro1+1/RC1)

 printf("GSAV_TAB \n"); 
 printf("$A_{V1}$ = %e dB\n", AVI_DB);
 %printf("$A_{VIsimple}$ = %e dB\n", AVIsimple_DB); 
 printf("GSAV_END \n \n"); 
 
 printf("GSZ_TAB \n"); 
 printf("$Z_{I1}$ = %e Ohm\n", ZI1);
 printf("$Z_{O1}$ = %e Ohm\n", ZO1); 
 printf("GSZ_END \n \n");

%-------------------------------------------------------------------------------

%---------------------------Ouput Stage-----------------------------------------
%-----Data-----
BFP = 227.3
VAFP = 37.2
VEBON = 0.7
 printf("OSBijector_TAB \n"); 
 printf("$beta$ = %f \n", BFP); 
 printf("$V_{AFP}$ = %f V\n", VAFP);
 printf("$V_{BEON}$ = %f V \n", VEBON);
 printf("OSBijector_END \n \n");

%-----OP Analysis-----
VI2 = VO1
IE2 = (VCC-VEBON-VI2)/RE2
IC2 = BFP/(BFP+1)*IE2
VO2 = VCC - RE2*IE2
 printf("OSOP_TAB \n"); 
 printf("$V_{I2}$ = %e V\n", VI2);
 printf("$I_{E2}$ = %e A \n", IE2);
 printf("$I_{C2}$ = %e A \n", IC2);
 printf("$V_{O2}$ = %e V \n", VO2);
 printf("OSOP_END \n \n");
 
 
%-----Incremental Circuit-----
gm2 = IC2/VT                               %--Incremental Parameters--
go2 = IC2/VAFP
gpi2 = gm2/BFP
ge2 = 1/RE2

%AC Analysis (Voltage gain)
AV2 = gm2/(gm2+gpi2+go2+ge2)
ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2)
ZO2 = 1/(gm2+gpi2+go2+ge2)

 printf("OSAV_TAB \n"); 
 printf("$A_{V2}$ = %e dB\n", AV2);
 printf("OSAV_END \n \n"); 
 
 printf("OSZ_TAB \n"); 
 printf("$Z_{I2}$ = %e Ohm\n", ZI2);
 printf("$Z_{O2}$ = %e Ohm\n", ZO2); 
 printf("OSZ_END \n \n");

%-------------------------------------------------------------------------------

%------------------------------Total--------------------------------------------
gB = 1/(1/gpi2+ZO1)
AV = (gB+gm2/gpi2*gB)/(gB+ge2+go2+gm2/gpi2*gB)*AV1
AV_DB = 20*log10(abs(AV))
ZI=ZI1
ZO=1/(go2+gm2/gpi2*gB+ge2+gB)

 printf("Final_TAB \n"); 
 printf("$A_{V}$ = %e dB\n", AV_DB);
 printf("$Z_{I}$ = %e Ohm\n", ZI);
 printf("$Z_{O}$ = %e Ohm\n", ZO); 
 printf("Final_END \n \n");
 
%------------------------------Teoria Ponto 3------------------------------------
f = 100*1e3;
w = 2*pi*f;

f_H=f;
%f_L= (1/(3*CS)+ 1/((ZO+RL)*Co) +1/((ZI+RS)*Cb))/(2*pi);
f_L= (1/((ZI+RS)*CS) + 1/(((rpi1+RSB)/(rpi1*gm1))*Cb) + 1/((ZO+RL)*Co))/(2*pi)

band=f_H - f_L

 printf("LC_TAB \n"); 
 printf("$Lower CO freq$ = %e Hz\n", f_L);
 printf("LC_END \n \n");
 
freq=logspace(1,8);

for i=1:50
 wfreq=2*pi*freq(i)
ZCS   = 1/(wfreq*j*CS);
ZCB   = 1/(wfreq*j*Cb);
ZEB = 1/(1/RE1 + 1/ZCB);
zpi2  = 1/gpi2;
zo2   = 1/go2;
ZE2 = 1/(1/zo2 + 1/RE2);
ZCO   = 1/(wfreq*j*Co);

 HOT = [  RS+ZCS+RB , -RB, 0 , 0 , 0 , 0 , 0; ...
 	0 , -ZEB , -ro1, ZEB + ro1 + RC1, -RC1, 0, 0 ; ...
        0 , rpi1*gm1 , 1  , 0 , 0, 0 , 0  ; ...
        0 , 0 , 0  , -RC1, zpi2+RC1, ZE2 , -ZE2; ...
        0 , 0  , 0  , 0  , -1-zpi2*gm2 , 1  , 0 ; ...
        0, 0 , 0  , 0 , 0 , -ZE2, ZE2+ZCO+RL ; ...
        -RB  , RB + rpi1 + ZEB , 0 , -ZEB  , 0 , 0 , 0 ];
        
DOG = [Vinput; 0; 0 ; 0 ; 0; 0; 0];

Gain= HOT\DOG;

GaindB(i) = 20*log10(abs(Gain(7)*RL/Vinput))
endfor

fig1=figure(1)
plot(log10(freq),GaindB,"r");
xlabel("log(f) (Hz)");
ylabel("Gain (dB)");
grid on;
print (fig1, "Gain", "-depsc");
close(fig1)

