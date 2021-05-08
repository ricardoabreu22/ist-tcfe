%Lab3                                                  Academic Year 2020/2021
%                                                            TCFE
%Authors:
%Joana Matos, 95799
%Ricardo Abreu, 95842
%Vasco Emídio, 95856
%MEAer
close all 
clear all

format long;

%---------------------Primary Circuit (given data)------------------------------
f=50;
T=1/f;
Vin = 230;

%--------------------Secondary Circuit (chosen data)----------------------------
R1 = 60000;
R2 = 60000;
C = 100e-6;
%-----------------------Transformer---------------------------------------------
n = 11.203488888888;
A = Vin/n;

%Data table
 printf("Data_TAB \n"); 
 printf("$V_{in}$ = %e V\n", Vin);
 printf("$f$ = %e Hz\n", f); 
 printf("$n$ = %e \n", n);
 printf("$R_{1}$ = %e Ohm \n", R1);
 printf("$R_{2}$ = %e Ohm \n", R2);
 printf("$C$ = %e F \n", C);
 printf("Data_END \n \n");

%-----------------------Envelope Detector---------------------------------------
t=linspace(0, 10*(T/2), 1000);
w=2*pi*f;

vs = A*cos(w*t);

vfr = zeros(1, length(t));     %Full-wave rectifier

  for i=1:length(t)
    vfr(i) = abs(vs(i));
  endfor

venv = zeros(1, length(t));   %Envelope Voltage

tOFF = (1/w) * atan(1/(w*R1*C));               %tOFF
vOexp = A*cos(w*tOFF)*exp(-(t-tOFF)/(R1*C));   %exponencial


for i=1:length(t)                               %Envelope Voltage
	if t(i) < tOFF
	  venv(i) = vfr(i);

	elseif vOexp(i) > vfr(i)
    venv(i) = vOexp(i);
    
	else
	  tOFF = tOFF + 1/(f*2);
    vOexp = A*abs(cos(w*tOFF))*exp(-(t-tOFF)/(R1*C));
    venv(i) = vfr(i);
	 endif
endfor
	
average = mean(venv)
ripple_env = max(venv) - min(venv)
average_env = ripple_env/2 + min(venv)

printf ("Envelope_TAB\n");
printf ("Ripple Envelope = %e \n", ripple_env);
printf ("Average Envelope = %e \n", average_env);
printf ("Envelope_END\n\n");

%-----------------------Voltage Regulator---------------------------------------
ndiode = 20;
VON = 0.6;

vreg = zeros(1, length(t));
dc_vreg = 0;
ac_vreg = zeros(1, length(t));
	
if average_env >= VON*ndiode       %DC regulator
  dc_vreg = VON*ndiode;
else
  dc_vreg = average_env;
endif
	
 
VT = 0.025;                   %AC component regulator
Is = 1e-14;
eta = 1;

RD = eta*VT/(Is*exp(VON/(eta*VT)));

ac_vreg = ((ndiode*RD)/(ndiode*RD +R2))*(venv-average_env);

vreg = dc_vreg+ac_vreg;

%plots of the values
	
%output voltages at rectifier, envelope detector and regulator

%figure 1
%fig1 = figure(1);
%title("Regulator and envelope output voltage v_o(t)")
%plot (t*1000,venv, ";vo_{envelope}(t);");
%xlabel ("t[ms]")
%ylabel ("v_O [Volts]")
%legend('Location','northeast');
%print (fig1, "all_vout.eps", "-depsc");


fig1 = figure(1);
title("Envelope and regulator output voltages")
plot (t*1000, vfr, ";vo_{rectified}(t);", t*1000,venv, ";vo_{envelope}(t);", t*1000,vreg, ";vo_{regulator}(t);");
xlabel ("t[ms]")
ylabel ("v_O [Volts]")
legend('Location','northeast');
print (fig1, "Output", "-depsc");
close(fig1)
	
%Deviations (vOenv - 12) 
fig2 = figure(2);
title("Deviation from the DC voltage")
plot (t*1000,venv-12, ";vo-12 (t);",t*1000,vreg-12,"color","r",";vOreg-12 (t);");
xlabel ("t[ms]")
ylabel ("v_O [Volts]")
legend('Location','northeast');
print (fig2, "Deviation.eps", "-depsc");
close(fig2)
	
average_reg = mean(vreg);
ripple_reg = max(vreg)-min(vreg);

printf ("V_TAB\n");
printf ("$V_{max}$= %e \n", max(vreg));
printf ("$V_{min}$= %e \n", min(vreg));
printf ("Average Voltage = %e \n", average_reg);
printf ("V_END\n\n");

printf ("Ripple_TAB\n");
printf ("Ripple Voltage= %e \n", ripple_reg);
printf ("Average Voltage Difference = %e \n", average_reg - 12);
printf ("Ripple_END\n\n");

total_cost = ((R1+R2)/1000)+C*1000000+((ndiode+4)*0.1);
merit=1/(total_cost*(ripple_reg+abs(average_reg-12)+10e-6));

printf ("Total cost of the components = %e \n", total_cost);
printf ("Merit= %e \n", merit);
