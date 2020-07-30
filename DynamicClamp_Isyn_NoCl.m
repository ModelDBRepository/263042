
%% This is the main function to do dynamic clamp to differnt SK mutant fly photoreceptors, 
%% during the naturalistic stimuli. 

%% This program calls DynamicClampIter to calculate Isyn
%% and all the K+ current. Then it calculates the ATP that is consumed during the naturalistic 
%% stimuli pattern

%% Input: datasource saves the photoreceptor's intracellular's recordings; param is the parameter for the 
%% Hodgkin Huxley model 

%datasource = 'SKSlo4Vol_All'
%datasource = 'SKVol';
%datasource = 'WTVol';
%datasource = 'sloVol'

%% The parameters used in roger's clamp data

%WT
% param =  [RestPos{1,2}+correct_restvol 20  -57.1  0  -85  4*0.085e-3   0 0 -5 0   2.4e-3   5e-3     0.11e-3 ]; 

% Skslo
% param =  [RestPos{4,2}+correct_restvol 20  -57.1  0  -85  4*0.085e-3  0 0 -5 0   2.4e-3   0.85*5e-3     0.11e-3 ]; 

% Sk
% param =  [-65  20  -57.1  0  -85  4*0.085e-3   0 0 -5 0   0.6*2.4e-3    0.6*5e-3     0.11e-3 ];

% slo
% %param =  [-65  20  -57.1  0  -85  4*0.085e-3   0 0 -5 0   0.8*2.4e-3   0.65*5e-3     0.11e-3 ];



%% the parameters used for Table 1. 
% WT
%param = [-65   1  -57.1  0  -85    0.85e-3  0  0   -5   0    2.4e-3   5e-3   0.11e-3 0];

% SKSlo
% param = [-65   1  -57.1  0  -85    0.85e-3  0  0   -5   0    2.4e-3  0.85*5e-3   0.11e-3 0];  

% SK
% param =  [-65  1  -57.1   0  -85  0.85e-3      0 0 -5 0   0.6*2.4e-3    0.6*5e-3     0.11e-3 ];

% Slo
% param =  [-65  1  -57.1  0  -85  0.85e-3      0 0 -5 0   0.8*2.4e-3   0.65*5e-3     0.11e-3 ];


function  DynamicClamp_Isyn_NoCl(datasource,param)

load(datasource)

% WTVol are the selected NS mean responses for different WT photoreceptors.
% the slected data are in the folder /Users/zhuoyisong/MATLAB/SK-Zhuoyi/DynamicConductance_Feedback/newOct2017/DynamicClamp_Diana_NS/NS_EXP_Diana/WT
wt_rest = -57; sk_rest = -47.45; slo18_rest = -46; SKslo18_rest = -41;
RestPos = cell(4,2);
RestPos{1,1} = 'wt';RestPos{2,1} = 'sk';RestPos{3,1} = 'slo18';RestPos{4,1} = 'SKslo18';
RestPos{1,2} = wt_rest; RestPos{2,2} = sk_rest; RestPos{3,2} = slo18_rest; RestPos{4,2} = SKslo18_rest;

eval(['VolEXP = ' datasource ';']);
%VolEXP = SKSlo4Vol_All;
%VolEXP = WTVol;
%VolEXP = SKVol;
%VolEXP = sloVol;

correct_restvol = 0; % a parameter to correct the resting potential for the experimental data
samprate        = 1000;

VolSim = []; IK =[]; IShaker=[];  IShab=[];  INew=[];  ICleak=[]; 
ICl = []; ILeak=[]; ATP_m = [];ATPk_m = []; ATPcl_m = [];ATP_simon_m = [];
ISyn = [];

for gg = 1:size(VolEXP,2)


TT = 980;

Vol_exp = VolEXP(1:TT,gg);

load NS_xf_BG105_MacroC
II = mean(NS_xf_BG105_MacroC,2)*0.001;
I = II(length(II)-TT+1:end);

Isyn = zeros(size(I));
gsyn = zeros(size(I));
gsynm = 0;

samprate = 1000;

[t, y, Isyn,Ishaker,Ishab,Ileak,Inew,Icleak,Icl]= DynamicClampGoodIter(I, param, samprate,Vol_exp);

%figure(1);hold all; plot(t,y(:,1));plot(t,Vol_exp); title('Voltage');xlabel('ms'); ylabel('mV');
%figure(2); plot(Isyn); title('LIC'); xlabel('ms'); ylabel('nA');

Vol = y(:,1);
Ik = Ishaker + Ishab + Inew + Ileak;

VolSim = [VolSim, Vol]; IK =[IK,Ik]; IShaker=[IShaker,Ishaker];  
IShab=[IShab,Ishab];  INew=[INew,Inew];  ICleak=[ICleak,Icleak]; 
ICl = [ICl,Icl]; ILeak=[ILeak,Ileak]; ISyn = [ISyn,Isyn];

timelength = length(I)/samprate; 
% here outward current is positive, Icl is positive outward currrent, Cl
% in, K out
Ip = 0.5*(Ishaker + Ishab + Inew + Ileak - I*0.0086*0.001) - 0.25*(Icl+Icleak); 
Ipk = 0.5*(Ishaker + Ishab + Inew + Ileak - I*0.0086*0.001);
Ipcl = -0.25*(Icl+Icleak); % the K extruded out by Na-K-Cl cotranspoter, does not require ATP

%To increase stability, 10^(-12) is multiplied to NA
NA = 6.02*10^(11);   % avacado constant, 6.02*10^23/mol
F = 96485;           % farady constant, C/mol
Ip_s = 201; Ip_e = length(I);
ATP = sum(Ip)*NA/F/(timelength)

Isimon = I+Isyn+Icl;
ATP_simon = (1/3)*(sum(Isimon))*NA/F/(timelength)

ATP_m = [ATP_m,ATP];
ATP_simon_m = [ATP_simon_m,ATP_simon];

figure;
subplot(2,2,1);hold all;
plot(t,Vol_exp);plot(t,y(:,1));
title('Voltage');xlabel('ms'); ylabel('mV');axis([1,TT,-110,-0]);
legend('experiment','simulation');
text(200,-80,['ATP ' num2str(ATP, '%10.5e\n')])
text(200,-90,['ATP simon ' num2str(ATP_simon, '%10.5e\n')])
 
subplot(2,2,2);hold all;
hold all;plot(Isyn,'k');plot(I,'b');
title('Isyn vs. LIC');xlabel('Time (ms)'); ylabel('Current (nA)');
legend('Isyn','LIC');

subplot(2,2,3);hold all;
plot(Ishaker,'b');plot(Ishab,'r');plot(Inew,'k');
title('K+ currents');xlabel('Time (ms)'); ylabel('Current (nA)');
legend('shaker','shab','new');

subplot(2,2,4);hold all;
plot(Ileak,'b');plot(Icleak);plot(Icl);
title('leak currents');xlabel('Time (ms)'); ylabel('Current (nA)');
legend('K leak','Cl leak','Icl','Location','east');

end

% SaveDataFile = [datasource '_ATP_NSDynamicClamp_0Cl_wtC'];
% save(SaveDataFile,'param','VolEXP','VolSim', 'ISyn', 'ICleak', 'ICl', 'IShaker',...
%     'IShab', 'INew', 'ILeak', 'ATP_m', 'ATP_simon_m');

