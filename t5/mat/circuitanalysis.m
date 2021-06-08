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

%-----Data-----
C1=220e-9;
C2a=220e-9;
C2b=220e-9;
R1a=1e3;
R1b=10e3;
R2=1e3;
R3a=100e3;
R3b=100e3;
R3c=100e3;
R4a=1e3;
R4b=10e3;
%------------
%--Replace R_a and R_b by an equivalent R_---
C2=(C2a*C2b)/(C2a+C2b);
R1=(R1a*R1b)/(R1a+R1b);
R3=R3a+(R3b*R3c)/(R3b+R3c);
R4=(R4a*R4b)/(R4a+R4b);
%R4=R4a

 printf("Data_TAB \n");
 printf("$R_{1a}$ = %e Ohm \n", R1a);
 printf("$R_{1b}$ = %e Ohm \n", R1b);
 printf("$R_{1}$ = %e Ohm\n", R1);
 printf("$R_{2}$ = %e Ohm \n", R2);
 printf("$R_{3a}$ = %e Ohm \n", R3a);
 printf("$R_{3b}$ = %e Ohm \n", R3b);
 printf("$R_{3c}$ = %e Ohm \n", R3c);
 printf("$R_{3}$ = %e Ohm \n", R3);
 printf("$R_{4a}$ = %e Ohm \n", R4a);
 printf("$R_{4b}$ = %e Ohm \n", R4b);
 printf("$R_{4}$ = %e Ohm \n", R4);
 printf("$C_{1}$ = %e F \n", C1);
 printf("$C_{2a}$ = %e Ohm \n", C2a);
 printf("$C_{2b}$ = %e Ohm \n", C2b);
 printf("$C_{2}$ = %e F \n", C2);
 printf("Data_END \n \n");
%------------------------------------------

%----1. Gain, Zi, Zout---------------------
fo=1e3;
woo=2*pi*fo;
wL = 1/(R1*C1);
wH = 1/(R2*C2);
wO = sqrt(wL*wH);
ZC1=1/(j*C1*wO);
ZC2=1/(j*C2*wO);

AVHP1=(R1*C1*j*wO/(1+R1*C1*j*wO));
AVOP1=(1+R3/R4);
AVLP1=(1/(1+R2*C2*j*wO));
T_1=AVHP1*AVOP1*AVLP1;
%dB
T_h1dB=20*log10(abs(AVHP1));
T_op1dB=20*log10(abs(AVOP1));
T_l1dB=20*log10(abs(AVLP1));
T_1dB=20*log10(abs(T_1));

wCdB=wO/(2*pi);

AVHP0=(R1*C1*j*woo/(1+R1*C1*j*woo));
AVOP0=(1+R3/R4);
AVLP0=(1/(1+R2*C2*j*woo));
T_fo=AVHP0*AVOP0*AVLP0;
%dB
T_h0dB=20*log10(abs(AVHP0));
T_op0dB=20*log10(abs(AVOP0));
T_l0dB=20*log10(abs(AVLP0));
T_fodB=20*log10(abs(T_fo));


Zin=R1+ZC1
Zout=R2*ZC2/(R2+ZC2)
Zinmod=sqrt(real(Zin)^2+imag(Zin)^2)
Zoutmod=sqrt(real(Zout)^2+imag(Zout)^2)

 printf("Freq_TAB \n");
 printf("$w_{L}$ = %e rad/s \n", wL);
 printf("$w_{H}$ = %e rad/s \n", wH);
 printf("$w_{O}$ = %e rad/s \n", wO);
 printf("$f_{O}$ = %e Hz \n", wCdB);
 printf("Freq_END \n \n");
 
 printf("Results1_TAB \n");
 printf("$Z_{in}$ = %e + j%e Ohm \n", real(Zin), imag(Zin));
 printf("$Z_{out}$ = %e + j%eOhm \n", real(Zout), imag(Zout));
 printf("$|Z_{in}|$ = %e Ohm \n", Zinmod);
 printf("$|Z_{out}|$ = %e Ohm \n", Zoutmod);
 printf("Results1_END \n \n");
 
 printf("Results2_TAB \n");
 printf("$AV_{HP}$ = %e dB \n",T_h1dB);
 printf("$AV_{OPAMP}$ = %e dB \n",T_op1dB);
 printf("$AV_{LP}$ = %e dB \n",T_l1dB);
 printf("$AV$ = %e dB \n",T_1dB);
 printf("Results2_END \n \n");

 printf("Results3_TAB \n");
 printf("$AV_{HP}$ = %e dB \n",T_h0dB);
 printf("$AV_{OPAMP}$ = %e dB \n",T_op0dB);
 printf("$AV_{LP}$ = %e dB \n",T_l0dB);
 printf("$AV$ = %e dB \n",T_fodB);
 printf("Results3_END \n \n");
 
%----2.--------------------------------------
f=logspace(1,8,70);
s=j*2*pi*f;
for n=1:length(s)
	T(n)=(R1*C1*s(n))./(1+R1*C1*s(n)).*(1+R3/R4).*(1/(1+R2*C2*s(n)));
	T_dB(n)=20*log10(abs(T(n)));
endfor

figure1=figure();
semilogx(f,180*arg(T)/pi);
xlabel("Frequency [Hz]");
ylabel("Phase [Deg]");
title("Phase/Frequency response");
print(figure1,"phase_response.eps","-depsc");
figure2=figure();
semilogx(f,T_dB);
xlabel("Frequency [Hz]");
ylabel("Gain [dB]");
title("Gain/Frequency response");
print(figure2,"gain_response.eps","-depsc");

%--------------------------------------------