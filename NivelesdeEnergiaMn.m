clear, clf

n1 = 500; n2 = 1500; n = n2 - n1 + 1;
H1 = 2200; H2 = 4600; N = 2048; % H1 es donde inicia la medición, H2 es donde finaliza, N es el número de puntos (Exp.nPoints)
dH = (H2-H1)/(N-1); h1 = (H1 + (n1-1)*dH)/10; h2 = (H1 + (n2-1)*dH)/10;
%data = textread('C:\Users\juanp\Downloads\Técnicas A\CALCITA ISLANDIA 5.dat','%*d %*d %*n %n %*d'); 
Xexp = textread('C:\Users\juanp\Downloads\Técnicas A\CALCITA ISLANDIA 5.dat','%*d %*d %n %*n %*d'); 
Yexp = textread('C:\Users\juanp\Downloads\Técnicas A\CALCITA ISLANDIA 5.dat','%*d %*d %*n %n %*d'); 
Exp.CenterSweep=[(h1+h2)/2 h2-h1]; 
%Exp.CenterSweep=[340 580];  
% Remove baseline of Yexp
a = 0;
Yexp = Yexp + a; 

Exp = struct('mwFreq',9.43,'Range',[h1 h2],'nPoints',n);
Exp.Harmonic = 1;


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
Sys.E = 40; % consider adding Estrain
%Sys.EStrain =1;
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
Exp.ModAmp = 4;
[B,S] = pepper(Sys,Exp);

% ================= Powder spectrum - transitions ===================

c = 0.06;
b= 0.14;
Opt.Transitions = [1 2]; % 'all'; 
[B,S] = pepper(Sys,Exp,Opt);

% ===================== Plot =======================================

% B = 10*B;
plot(Xexp/10,Yexp+2.5,'k',B,S-1,'r');
xlim([150 500]); 
%ylim([-50 250]);
B = B';
%,B,Fit5-10,'b',B,Fit6-12,'b'
% ============ Common parameters ====================================

Opt.Transitions = 'all'; % for resfields
Opt.Threshold = 0.0;  % for resfields
Par.PlotThreshold = 0.0; % for levelsplot

% ============ Energy levels - 1 ====================================

Ori = [0 0 0]; % Euler angles
figure; levelsplot(Sys,Ori,[0 800],Exp.mwFreq,Par); xlim([150 500]);
xlabel('magnetic field (mT)');

% ============ Crystal rotation - 1 ==================================

% Sys.gFrame = [10 20 30]*pi/180;
% Exp.CrystalSymmetry = 'P212121';

% Generate orientations in a single rotation plane
rotN = [0 1 0];

% rotation axis
% Output of rotplane is [phi,theta].
% theta : angle between vector and z-axis.
% phi : angle between projection of v in xy-plane and x-axis.
% x, y, z : laboratory axes.

N = 91;
[phi,theta] = rotplane(rotN,[pi 0],N);
chi = zeros(N,1);
Exp.CrystalOrientation = [phi(:) theta(:) chi];

% Simulate spectra
Opt.Output = 'separate'; % make sure spectra are not added up
Bres1 = resfields(Sys,Exp,Opt);
Bres1 = Bres1';
[Pos,Amp,Wid,Trans,Grad] = resfields(Sys,Exp,Opt);
Pos = Pos'; Amp = Amp';

% plotting
figure; plot(Bres1,theta*180/pi); xlim([150 500]);
xlabel('magnetic field (mT)');
ylabel('theta (°)');

% ============ Energy levels - 2 ====================================

Ori = [0 pi/2 0]; % Euler angles
figure; levelsplot(Sys,Ori,[0 800],Exp.mwFreq,Par); xlim([150 500]);

% ============ Crystal rotation - 2 ==================================

% Sys.gFrame = [10 20 30]*pi/180;
% Exp.CrystalSymmetry = 'P212121';

% Generate orientations in a single rotation plane
rotN = [0 0 1];

% rotation axis
N = 91;
[gamma,beta] = rotplane(rotN,[pi 0],N);
alpha = zeros(N,1);
Exp.CrystalOrientation = [alpha beta(:) gamma(:)];

% Simulate spectra
Opt.Output = 'separate'; % make sure spectra are not added up
Bres2 = resfields(Sys,Exp,Opt);
Bres2 = Bres2';

% plotting
figure; plot(Bres2,theta*180/pi); xlim(Exp.Range);
xlabel('magnetic field (mT)');
ylabel('theta (°)'); 

% ===============================================================  
