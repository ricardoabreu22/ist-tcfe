%Lab2-Circuit analysis                             Academic Year 2020/2021
%                                                            TCFE
%Authors:
%Joana Matos, 95799
%Ricardo Abreu, 95842
%Vasco Emídio, 95856
%MEAer
close all
clear all

format long e

%---Open the file and save the data---
data=fopen("../data.txt",'r');
val = fscanf(data,'\n\nPlease enter the lowest student number in your group: \n\nUnits for the values: V, mA, kOhm, mS and uF\n\n\nValues:  R1 = %e \nR2 = %e \nR3 = %e \nR4 = %e \nR5 = %e \nR6 = %e \nR7 = %e \nVs = %e \nC = %e \nKb = %e \nKd = %e\n\n');
fclose (data);
 %---Data generated by Python 2.7---
 printf("Data:  V, A, Ohm, S and F\n")
 R1 = val(1)*10^(3)
 R2 = val(2)*10^(3)
 R3 = val(3)*10^(3)
 R4 = val(4)*10^(3)
 R5 = val(5)*10^(3)
 R6 = val(6)*10^(3)
 R7 = val(7)*10^(3)
 Vs = val(8)
 C = val(9)*10^(-6)
 Kb = val(10)*10^(-3)
 Kd = val(11)*10^(3)
 
 %----2. Theoretical Analysis----
 %--Data table--
 printf("Data_NTAB \n");                             %nodal table
 printf("$R_{1}$ = %e \n", R1);
 printf("$R_{2}$ = %e \n", R2);
 printf("$R_{3}$ = %e \n", R3);
 printf("$R_{4}$ = %e \n", R4);
 printf("$R_{5}$ = %e \n", R5);
 printf("$R_{6}$ = %e \n", R6);
 printf("$R_{7}$ = %e \n", R7);
 printf("$V_{s}$ = %e \n", Vs);
 printf("$C$ = %e \n", C);
 printf("$K_{b}$ = %e \n", Kb);
 printf("$K_{d}$ = %e \n", Kd);
 printf("Data_NEND \n \n");
 %--Ngspice--
 ng1=fopen("../sim/ng1.cir", "w");
 fprintf(ng1, ".OP\n");
 fprintf(ng1, "Vs v1 GND %e\n", Vs);
 fprintf(ng1, "R1 v1 v2 %e\n", R1);
 fprintf(ng1, "R2 v3 v2 %e\n", R2);
 fprintf(ng1, "R3 v2 v5 %e\n", R3);
 fprintf(ng1, "R4 v5 GND %e\n", R4);
 fprintf(ng1, "R5 v5 v6 %e\n", R5);
 fprintf(ng1, "R6 GND v4 %e\n", R6);
 fprintf(ng1, "V6 v4 v7 %e\n",0); %we add this voltage source to sense the current through R6
 fprintf(ng1, "R7 v7 v8 %e\n", R7);
 fprintf(ng1, "C v6 v8 %e\n", C);
 fprintf(ng1, "Gb v6 v3 v2 v5 %e\n", Kb);
 fprintf(ng1, "Hd v5 v8 V6 %e\n", Kd);
 fprintf(ng1, ".END\n");
 fclose(ng1);
 %--------------------------------------------
 %--Point 1--
 printf("\n--2.1) Node Analysis--\n");
 G1=1/R1; G2=1/R2; G3=1/R3; G4=1/R4; G5=1/R5; G6=1/R6; G7=1/R7; %calculate 
                                                                %conductances
 An=[1 ,0        ,0  ,0        ,0  ,0     ,0  ;
     G1,-G1-G2-G3,G2 ,G3       ,0  ,0     ,0  ;
     0 ,G2+Kb    ,-G2,-Kb      ,0  ,0     ,0  ;
     0 ,G3       ,0  ,-G3-G4-G5,G5 ,G7    ,-G7;
     0 ,-Kb      ,0  ,G5+Kb    ,-G5,0     ,0  ;
     0 ,0        ,0  ,0        ,0  ,-G6-G7,G7 ;
     0 ,0        ,0  ,1        ,0  ,G6*Kd ,-1 ];          %matrix form
 Bn=[Vs;0;0;0;0;0;0];

 Vn=An\Bn;                                           %nodal voltages array

 printf("Nodal_NTAB \n");                             %nodal table
 printf("$V_{1}$ = %e \n", Vn(1));                   
 printf("$V_{2}$ = %e \n", Vn(2));
 printf("$V_{3}$ = %e \n", Vn(3));
 printf("$V_{5}$ = %e \n", Vn(4));
 printf("$V_{6}$ = %e \n", Vn(5));
 printf("$V_{7}$ = %e \n", Vn(6));
 printf("$V_{8}$ = %e \n", Vn(7));
 printf("Nodal_NEND \n \n");

 printf("It is now possible to determine:\n");
 i1=(Vn(1)-Vn(2))*G1;                                %branch currents
 i2=(Vn(3)-Vn(2))*G2;
 i3=(Vn(2)-Vn(4))*G3;
 i4=(Vn(4)-0)*G4;
 i5=(Vn(4)-Vn(5))*G5;
 i6=(0-Vn(6))*G6;
 i7=(Vn(6)-Vn(7))*G7;
 is=-i1;
 ib=Kb*(Vn(2)-Vn(4));
 ivd=-i7;
 ic=0;
 printf("ValN_MNTAB\n");                             
 printf("$I_{R1}$ = %e \n", i1);
 printf("$I_{R2}$ = %e \n", i2);
 printf("$I_{R3}$ = %e \n", i3);
 printf("$I_{R4}$ = %e \n", i4);
 printf("$I_{R5}$ = %e \n", i5);
 printf("$I_{R6}$ = %e \n", i6);
 printf("$I_{R7}$ = %e \n", i7);
 printf("$I_{s}$ = %e \n", is);
 printf("$I_{b}$ = %e \n", ib);
 printf("$I_{Vd}$ = %e \n", ivd);
 printf("$I_{c}$ = %e \n", ic);
 printf("ValN_MNEND\n");
 
  %--Point 2--
  Vx = Vn(5)-Vn(7)

printf("\n----------------2.2)-------------\n");

D = [1    ,0              ,0    ,0   ,0,0         ,0    ;
     1/R1 ,-1/R1-1/R3-1/R2,1/R2 ,1/R3,0,0         ,0    ;
     0    ,Kb+1/R2        ,-1/R2,-Kb ,0,0         ,0    ;
     -1/R1,1/R1           ,0    ,1/R4,0,1/R6      ,0    ;
     0    ,0              ,0    ,1   ,0,Kd/R6     ,-1   ;
     0    ,0              ,0    ,0   ,1,0         ,-1   ;
     0    ,0              ,0    ,0   ,0,-1/R6-1/R7,1/R7];

e = [0;0;0;0;0;Vx;0];

f = D\e;
V6o=f(5)
Ix=-(f(4)-f(5))/R5 - Kb*(f(2)-f(4))

printf("Ix = %f\n", Ix*1000);

 i12=(f(1)-f(2))*G1;                                %branch currents
 i22=(f(3)-f(2))*G2;
 i32=(f(2)-f(4))*G3;
 i42=(f(4)-0)*G4;
 i52=(f(4)-f(5))*G5;
 i62=(0-f(6))*G6;
 i72=(f(6)-f(7))*G7;
 is2=-i12;
 ib2=Kb*(f(2)-f(4));
 ivd2=-i72;

Req=Vx/Ix
T=Req*C

 printf("Val21_NTAB21 \n");                          %table 1
 printf("$I_{R1}$ = %e \n", i12);
 printf("$I_{R2}$ = %e \n", i22);
 printf("$I_{R3}$ = %e \n", i32);
 printf("$I_{R4}$ = %e \n", i42);
 printf("$I_{R5}$ = %e \n", i52);
 printf("$I_{R6}$ = %e \n", i62);
 printf("$I_{R7}$ = %e \n", i72);
 printf("$I_{s}$ = %e \n", is2);
 printf("$I_{b}$ = %e \n", ib2);
 printf("$I_{Vd}$ = %e \n", ivd2);
 printf("$V_{1}$ = %e \n", f(1));                   
 printf("$V_{2}$ = %e \n", f(2));
 printf("$V_{3}$ = %e \n", f(3));
 printf("$V_{5}$ = %e \n", f(4));
 printf("$V_{6}$ = %e \n", f(5));
 printf("$V_{7}$ = %e \n", f(6));
 printf("$V_{8}$ = %e \n", f(7));
 printf("Val21_NEND21 \n \n");                        
 printf("Val22_NTAB22 \n");                             %table 2
 printf("$V_{x}$ = %e \n", Vx);
 printf("$I_{x}$ = %e \n", Ix);
 printf("$R_{eq}$ = %e \n", Req);
 printf(" tau = %e\n", T);
 printf("Val22_NEND22 \n \n");   
%--Ngspice--
 ng2=fopen("../sim/ng2.cir", "w");
 fprintf(ng2, ".OP\n");
 fprintf(ng2, "Vs v1 GND %e\n", 0);
 fprintf(ng2, "R1 v1 v2 %e\n", R1);
 fprintf(ng2, "R2 v3 v2 %e\n", R2);
 fprintf(ng2, "R3 v2 v5 %e\n", R3);
 fprintf(ng2, "R4 v5 GND %e\n", R4);
 fprintf(ng2, "R5 v5 v6 %e\n", R5);
 fprintf(ng2, "R6 GND v4 %e\n", R6);
 fprintf(ng2, "V6 v4 v7 %e\n",0); %we add this voltage source to sense the current through R6
 fprintf(ng2, "R7 v7 v8 %e\n", R7);
 fprintf(ng2, "Gb v6 v3 v2 v5 %e\n", Kb);
 fprintf(ng2, "Hd v5 v8 V6 %e\n", Kd);
 fprintf(ng2, "Vx v6 v8 %e\n", Vx);
 fprintf(ng2, ".END\n");
 fclose(ng2);
 %-------------------------------------------- 
 
 %--Point 3--
printf("----------------Ponto 3------------\n");
t = 0:1e-6:20e-3;
V6n = Vx * exp(-t/T);

y = figure(1);
plot(t*1e3, V6n, "r");
xlabel("t [ms]");
ylabel("V6n [V]");
legend("vn6","location","northeast");
grid on;
print(y, "theoretical.eps", "-depsc");
close(y);

%--Ngspice 3--
 ng3=fopen("../sim/ng3.cir", "w");
 fprintf(ng3, ".OP\n");
 fprintf(ng3, "Vs v1 GND %e\n", 0);
 fprintf(ng3, "R1 v1 v2 %e\n", R1);
 fprintf(ng3, "R2 v3 v2 %e\n", R2);
 fprintf(ng3, "R3 v2 v5 %e\n", R3);
 fprintf(ng3, "R4 v5 GND %e\n", R4);
 fprintf(ng3, "R5 v5 v6 %e\n", R5);
 fprintf(ng3, "R6 GND v4 %e\n", R6);
 fprintf(ng3, "V6 v4 v7 %e\n",0); %we add this voltage source to sense the current through R6
 fprintf(ng3, "R7 v7 v8 %e\n", R7);
 fprintf(ng3, "Gb v6 v3 v2 v5 %e\n", Kb);
 fprintf(ng3, "Hd v5 v8 V6 %e\n", Kd);
 fprintf(ng3, "C v6 v8 %e\n", C);
 fprintf(ng3, ".ic v(v6)=%e v(v8)=%e", f(5), f(7));
 fprintf(ng3, ".END\n");
 fclose(ng3);
 %-------------------------------------------- 
 %--Point 4--
printf("----------------Ponto 4------------\n");
fr=1000;
w=2*pi*fr;
Zc=1/j*w*C;
cVs=exp(-j*pi/2);

cMA = [1,0,0,0,0,0,0;
      -1/R1,1/R1,0,1/R4,0,1/R6,0;
      1/R1,-1/R1-1/R3-1/R2,1/R2,1/R3,0,0,0;
      0,1/R2+Kb,-1/R2,-Kb,0,0,0;
      0,-Kb,0,Kb+1/R5,-1/R5-1/Zc,0,1/Zc;
      0,0,0,0,0,-1/R6-1/R7,1/R7;
      0,0,0,1,0,Kd/R6,-1];

cMb = [cVs;0;0;0;0;0;0];

cMc = cMA\cMb;

for i=1:7
	if(i>3)
		printf("V%if=%f + j*%f \n",i+1, real(cMc(i)),imag(cMc(i)));
	else	
		printf("V%if=%f + j*%f \n",i, real(cMc(i)),imag(cMc(i)));
	endif
endfor 

polar= [abs(cMc(1)), angle(cMc(1));
 	abs(cMc(2)), angle(cMc(2)); 
 	abs(cMc(3)), angle(cMc(3)); 
 	abs(cMc(4)), angle(cMc(4)); 
 	abs(cMc(5)), angle(cMc(5)); 
 	abs(cMc(6)), angle(cMc(6)); 
 	abs(cMc(7)), angle(cMc(7))];
 for i=1:7
 	k=1;
 	if(i>3)
		printf("Amplitude and Phase of V%if = %f %f \n",i+1, polar(i,k),polar(i,k+1));
	else	
		printf("Amplitude and Phase of V%if = %f %f \n",i, polar(i,k),polar(i,k+1));
	endif
endfor
  printf("Val41_TAB41 \n");                             %table Amplitude
  printf("$V_{1}$ = %e \n", polar(1,1));
  printf("$V_{2}$ = %e \n", polar(2,1));
  printf("$V_{3}$ = %e \n", polar(3,1));
  printf("$V_{5}$ = %e \n", polar(4,1));
  printf("$V_{6}$ = %e \n", polar(5,1));
  printf("$V_{7}$ = %e \n", polar(6,1));
  printf("$V_{8}$ = %e \n", polar(7,1));
  printf("Val41_END41 \n \n");
   
  printf("Val42_TAB42 \n");                             %table Phase
  printf("$V_{1}$ = %e \n", polar(1,2));
  printf("$V_{2}$ = %e \n", polar(2,2));
  printf("$V_{3}$ = %e \n", polar(3,2));
  printf("$V_{5}$ = %e \n", polar(4,2));
  printf("$V_{6}$ = %e \n", polar(5,2));
  printf("$V_{7}$ = %e \n", polar(6,2));
  printf("$V_{8}$ = %e \n", polar(7,2));
  printf("Val42_END42 \n \n");
  
 temp=0:1e-6:20e-3;
 V6f=polar(5,1)*cos(w*temp-polar(5,2));
 
 yf=figure(2);
 plot(temp*1000,V6f,"r");
 xlabel("t [ms]");
 ylabel("V6f [V]");
 legend("vf6","location","northeast");
 grid on;
 print(yf, "Forced4.eps", "-depsc");
 close(yf);

printf("----------------Ponto 5------------\n");

t = 0:1e-6:20e-3;
V6t=f(5)*exp(-t/T) + polar(5,1)*cos(w*t -polar(5,2));
VSF=sin(w*t);

y2 = figure(3);
r=-5e-3:1e-6:0
line([-5 0],[Vn(5) Vn(5)],"linestyle","-","color","r");
hold on;
line([-5 0],[Vs  Vs],"linestyle","-","color","b");

plot(t*1000, V6t, "r");
hold on;
plot(t*1000, VSF, "b");
xlabel("t [ms]");
ylabel("[V]");
grid on;
legend("v6","vs","location","northeast");
print(y2, "Theoretical5.eps", "-depsc");
close(y2);


printf("----------------Ponto 6------------\n");

fre=logspace(-1,6);
w2=2*pi*fre;
ZI2=j*w2*C;
for o=1:50
cMD = [1,0,0,0,0,0,0;-1/R1,1/R1,0,1/R4,0,1/R6,0;1/R1,-1/R1-1/R3-1/R2,1/R2,1/R3,0,0,0;0,1/R2+Kb,-1/R2,-Kb,0,0,0;0,-Kb,0,Kb+1/R5,-1/R5-ZI2(o),0,ZI2(o);0,0,0,0,0,-1/R6-1/R7,1/R7;0,0,0,1,0,Kd/R6,-1];
cMh = [cVs;0;0;0;0;0;0];
cMl = cMD\cMh;

listaAMP6(o)=20*log10(abs(cMl(5)/cVs));
listaAMPs(o)=20*log10(abs(cMl(1)/cVs));
listaAMPc(o)=20*log10(abs((cMl(5)-cMl(7))/cVs));
listaPH6(o)=arg(cMl(5)/cVs);
listaPHs(o)=arg(cMl(1)/cVs);
listaPHc(o)=arg((cMl(5)-cMl(7))/cVs);
endfor


Amp=figure(4);
semilogx(fre,listaAMP6,"r");
hold on;
semilogx(fre,listaAMPs,"b");
hold on;
semilogx(fre,listaAMPc,"g");
hold on;
xlabel("log(f) (N/A)");
ylabel("Magnitude (dB)");
grid on;
legend("v6","vs","vc","location","northeast");
print(Amp,"freq6SC.eps","-depsc");
close(Amp);

Phase=figure(5);
semilogx(fre,listaPH6*180/pi,"r");
hold on;
semilogx(fre,listaPHs*180/pi,"b");
hold on;
semilogx(fre,listaPHc*180/pi,"g");
hold on;
xlabel("log(f) (N/A)");
ylabel("Phase (Degrees)");
grid on;
legend("v6","vs","vc","location","northeast");
print(Phase,"ph6SC.eps","-depsc");
close(Phase);
%--Ngspice4--
 ng4=fopen("../sim/ng4.cir", "w");
 fprintf(ng4, ".OP\n");
 fprintf(ng4, "Vs v1 GND 0 ac 1 sin(0 1 1k)\n");
 fprintf(ng4, "R1 v1 v2 %e\n", R1);
 fprintf(ng4, "R2 v3 v2 %e\n", R2);
 fprintf(ng4, "R3 v2 v5 %e\n", R3);
 fprintf(ng4, "R4 v5 GND %e\n", R4);
 fprintf(ng4, "R5 v5 v6 %e\n", R5);
 fprintf(ng4, "R6 GND v4 %e\n", R6);
 fprintf(ng4, "V6 v4 v7 %e\n",0); %we add this voltage source to sense the current through R6
 fprintf(ng4, "R7 v7 v8 %e\n", R7);
 fprintf(ng4, "Gb v6 v3 v2 v5 %e\n", Kb);
 fprintf(ng4, "Hd v5 v8 V6 %e\n", Kd);
 fprintf(ng4, "C v6 v8 %e\n", C);
 fprintf(ng4, ".ic v(v6)=%e v(v8)=%e",f(5),f(7));
 fprintf(ng4, ".END\n");
 fclose(ng4);
 
 %----------------------------------------------------
