clear

[B1,S1]=textread('C:\Users\juanp\Downloads\Técnicas A\CALCITA ISLANDIA 5.dat','%*d %*d %n %n %*d'); 
B1=B1/10;
S1=2980000*S1-42;
Sys.Nucs = '55Mn';
Sys.S = 5/2;
Sys.lw = 1;
g1 = 1.9998;   
g2 = 1.9996;
g3 = 2.0007;
Sys.g = [g1 g2 g3];
Sys.gStrain= [2 2 2];
A1 = 20.00;
A2 = 264.18;
A3 = 264.09;
Sys.A = [A1];
%Sys.AStrain = [20];
%Poner un A isotrópico de en torno a 30 variación de A en torno a 60
Sys.D = [124.4 -112.7 -11.7]; 
%Sys.DStrain = [200 200 200]; 
Sys.E = [10 -20 -5]
Sys.EStrain = [200 200 500]
Exp.mwFreq = 9.43;
Exp.Range = [220 460]; 
Exp.nPoints = 2048;

[B,S] = pepper(Sys,Exp);

Vary.g = [0.5 0.5 0.5];
Vary.A = [20];
Vary.D=[1000 1000 1000];
SimOpt.Method = 'matrix';

FitOpt.Method = 'simplex diff'; 
esfit('pepper',S1,Sys,Vary,Exp,SimOpt,FitOpt);
%plot(B1,S1, 'r', B,S, 'b');