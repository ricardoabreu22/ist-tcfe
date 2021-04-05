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
 Kd = val(11)*10^(-3)
 
 %----2. Theoretical Analysis----

 %--2.1)Node Analysis--
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
 i1=(Vn(2)-Vn(1))*G1;                                %branch currents
 i2=(Vn(3)-Vn(2))*G2;
 i3=(Vn(2)-Vn(4))*G3;
 i4=(Vn(4)-0)*G4;
 i5=(Vn(4)-Vn(5))*G5;
 i6=(0-Vn(6))*G6;
 i7=(Vn(6)-Vn(7))*G7;
 is=i1;
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
 
  %--2.2)Node Analysis--
  Vx = Vn(5)-Vn(7)

printf("\n----------------PONTO 2-------------\n");

D = [1,0,0,0,0,0,0;1/R1,-1/R1-1/R3-1/R2,1/R2,1/R3,0,0,0;0,Kb+1/R2,-1/R2,-Kb,0,0,0;-1/R1,1/R1,0,1/R4,0,1/R6,0;0,0,0,1,0,Kd/R6,-1;0,0,0,0,1,0,-1;0,0,0,0,0,-1/R6-1/R7,1/R7];

e = [0;0;0;0;0;Vx;0];

f = D\e;
V6o=f(5)
Ix=-(f(4)-f(5))/R5 - Kb*(f(2)-f(4))

printf("Ix = %f\n", Ix*1000);

Req=Vx/Ix
T=Req*C
printf("----------------Ponto 3------------\n");
t = 0:1e-6:20e-3;
V6n = Vx * exp(-t/T);

y = figure(1);
plot(t*1e3, V6n, "r");
xlabel("t [ms]");
ylabel("V6n [V]");
print(y, "Theoretical3.eps", "-depsc");
close(y);

printf("----------------Ponto 4------------\n");
fr=1000;
w=2*pi*fr;
Zc=1/j*w*C;
cVs=exp(-j*pi/2);

cMA = [1,0,0,0,0,0,0;-1/R1,1/R1,0,1/R4,0,1/R6,0;1/R1,-1/R1-1/R3-1/R2,1/R2,1/R3,0,0,0;0,1/R2+Kb,-1/R2,-Kb,0,0,0;0,-Kb,0,Kb+1/R5,-1/R5-1/Zc,0,1/Zc;0,0,0,0,0,-1/R6-1/R7,1/R7;0,0,0,1,0,Kd/R6,-1];

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
printf("----------------Ponto 5------------\n");

function fu = pieceWise(t, amplitude,phase, w, Volt, T)
	fu= Volt* ones(size(t));
	ind = t >= 0;
		fu(ind)=Volt*exp(-t(ind)/T) + amplitude*cos(w*t(ind) -phase);
endfunction
function fuVS = pieceWise2(t, w,VS)
	fuVS= VS* ones(size(t));
	ind = t >= 0;
		fuVS(ind)=sin(w*t(ind));
endfunction

t = -5e-3:1e-6:20e-3;

y2 = figure(2);
plot(t*1000, pieceWise(t, polar(5,1), polar(5,2), w, Vn(5), T), "r");
hold on;
plot(t*1000, pieceWise2(t,w,Vs), "b");
xlabel("t [ms]");
ylabel("[V]");
print(y2, "Theoretical5.eps", "-deps");