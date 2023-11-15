clear 

[B1,S1]=textread('C:\Users\juanp\Downloads\Técnicas A\CALCITA ISLANDIA 5.dat','%*d %*d %n %n %*d'); 
B1=B1/10;
S1=2980000*(S1)-42;

Sys.Nucs = '55Mn';
c = 2.998;
Sys.S = 5/2;
Sys.lw = 1;

gperp = 2.00123;   
gtang = 2.00131;
Sys.g = [gperp gtang];
Aperp = -88.23*c;
Atang = -87.6*c;
Sys.A = [Aperp Atang];
%
% ZFS terms
%
Sys.B4 = [0 1.12 0 0 -0.047 0 0 0 0]* c; %B34 and B04 terms only  
Sys.D = -76.0 * c; %consider adding Dstrain
%Sys.DStrain =2;
Sys.E = 0; % consider adding Estrain
%
% Nuclear Quadrupole
%
PD = -0.176*c;
PE = 0.13*c;
Qxx = PE - PD/3;
Qyy = -PE - PD/3;
Qzz =2*PD/3;
Sys.Q = [Qxx Qyy Qzz];
%
% Nuclear zeeman?

%
% Spectrometer settings
%
Exp.mwFreq = 9.43;
Exp.Range = [220 460];
Exp.nPoints = 2048;
[B,S] = pepper(Sys,Exp);

%Vary.g = [0.2 0.2];
%Vary.A = [200 200];
%Vary.D=80;
%Vary.E=15;
%SimOpt.Method = 'matrix';

%FitOpt.Method = 'simplex diff'; % simplex algorithm, integrals of spectra

%esfit('pepper',S1,Sys,Vary,Exp,SimOpt,FitOpt);

plot(B1,S1, 'r', B,S, 'b');