close all
clear all
%System Config
sys = mr.opts('MaxGrad', 32, 'GradUnit', 'mT/m', ...
 'MaxSlew', 130, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 30e-6, ...
 'rfDeadTime', 100e-6);
%% Define the hyperparameters
seq=mr.Sequence(sys); % Create a new sequence object
%roDuration=sys.gradRasterTime*160; % ADC duration
roDuration=3.2e-3; % ADC duration

fov=256e-3; Nx=256; Ny=Nx; % Define FOV and resolution
thickness=6e-3; % slice thickness
rfDuration=0.7e-3;

balanced=1; %enable b-ssfp or nb-ssfp (refocused FLASH)
inversion=0; %enable inversion pulses
TI1 = 3000e-3; % inversion time 1 in ms
TI2 = 2000e-3; % inversion time 2 in ms

TR=10e-3;
TE=TR/2; 
alpha=60;
%% Generating the necessary pulses
%rf180 = mr.makeBlockPulse(pi,'Duration',rfDuration,'system',sys, 'use','refocusing'); %Non-selective Inversion Pulse
rf180 = mr.makeAdiabaticPulse('hypsec',sys,'Duration',rfDuration,'dwell',1e-5);

% Create alpha/2-degree slice selection pulse and gradient for pre-steady
% part of the pulse. After applying alpha/2-degree slice selection pulse
% we will not sample for some time and then, we will assume it is in
% steady-state. 
[rf_preSteady, gz_preSteady] = mr.makeSincPulse(alpha/2*pi/180,'Duration',rfDuration,...
 'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'system',sys);


% Create alpha-degree slice selection pulse and gradient for SSFP
[rf, gz] = mr.makeSincPulse(alpha*pi/180,'Duration',rfDuration,... 
'SliceThickness',thickness,'apodization',0.5,'timeBwProduct',4,'system',sys);

% Define other gradients and ADC events
deltak=1/fov;
gx = mr.makeTrapezoid('x','FlatArea',Nx*deltak,'FlatTime',roDuration,'system',sys);
adc = mr.makeAdc(Nx,'Duration',gx.flatTime,'Delay',gx.riseTime,'system',sys);
gxPre = mr.makeTrapezoid('x','Area',-gx.area/2,'Duration',0.7e-3,'system',sys);
gzReph = mr.makeTrapezoid('z','Area',-gz.area/2,'Duration',0.6e-3,'system',sys);
gzReph_preSteady = mr.makeTrapezoid('z','Area',-gz_preSteady.area/2,'Duration',0.6e-3,'system',sys);
phaseAreas = -((0:Ny-1)-Ny/2+0.5)*deltak; 

if balanced == 1
    delayTE=ceil((TE - mr.calcDuration(gz) - mr.calcDuration(gxPre) - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
    delayTR_preSteady=ceil((TR/2 - mr.calcDuration(gz_preSteady) -mr.calcDuration(gzReph_preSteady)*2)/seq.gradRasterTime)*seq.gradRasterTime;
    delayTR=ceil((TR - TE - mr.calcDuration(gx)/2 -mr.calcDuration(gxPre)-mr.calcDuration(gzReph))/seq.gradRasterTime)*seq.gradRasterTime;
else
    delayTE=ceil((TE - mr.calcDuration(gz) - mr.calcDuration(gxPre) - mr.calcDuration(gx)/2)/seq.gradRasterTime)*seq.gradRasterTime;
    delayTR_preSteady=ceil((TR/2 - mr.calcDuration(gz_preSteady) -mr.calcDuration(gzReph_preSteady))/seq.gradRasterTime)*seq.gradRasterTime;
    delayTR=ceil((TR - TE - mr.calcDuration(gx)/2 -mr.calcDuration(gxPre))/seq.gradRasterTime)*seq.gradRasterTime;    
end

TI1delay=round((TI1-(mr.calcDuration(rf180)-mr.calcRfCenter(rf180)-rf180.delay))/seq.gradRasterTime)*seq.gradRasterTime;
TI2delay=round((TI2-(mr.calcDuration(rf180)-mr.calcRfCenter(rf180)-rf180.delay))/seq.gradRasterTime)*seq.gradRasterTime;



assert(all(delayTE>=0));
assert(all(delayTR>=0));
assert(all(delayTR_preSteady>=0));
assert(all(TI1delay>=0));
assert(all(TI2delay>=0));

%% Pulse Sequence
seq.addBlock(mr.makeLabel('SET','REV', 1));
if inversion==1 
    % Inversion
    seq.addBlock(rf180);
    seq.addBlock(mr.makeDelay(TI1delay));
    seq.addBlock(rf180);
    seq.addBlock(mr.makeDelay(TI2delay));
end
if balanced == 1
    % Initialization: 
    seq.addBlock(rf_preSteady,gz_preSteady);
    seq.addBlock(gzReph_preSteady);
    seq.addBlock(mr.makeDelay(delayTR_preSteady));
    seq.addBlock(gzReph_preSteady);
    for initial=1:20 %Cover 20 Scans until Steady-State
     rf.signal=-rf.signal; %Change Flip Direction
     seq.addBlock(rf,gz);
     gyPre =mr.makeTrapezoid('y','Area',phaseAreas(initial),'Duration',mr.calcDuration(gxPre),'system',sys);
     seq.addBlock(gxPre,gyPre,gzReph);
     seq.addBlock(mr.makeDelay(delayTE));
     seq.addBlock(gx);
     gyPre.amplitude=-gyPre.amplitude;
     seq.addBlock(gxPre,gyPre);
     seq.addBlock(gzReph);
     seq.addBlock(mr.makeDelay(delayTR))
    end

    for i=1:Ny % loop over phase encodes and define sequence blocks
         rf.signal=-rf.signal; %Change Flip Direction
         seq.addBlock(rf,gz);
         gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',mr.calcDuration(gxPre),'system',sys);
         seq.addBlock(gxPre,gyPre,gzReph);
         seq.addBlock(mr.makeDelay(delayTE));
         seq.addBlock(gx,adc);
         gyPre.amplitude=-gyPre.amplitude;
         seq.addBlock(gxPre,gyPre, mr.makeLabel('SET','LIN', i));
         seq.addBlock(mr.makeDelay(delayTR)); 
         seq.addBlock(gzReph);
    end
else 
    % Initialization: 
    seq.addBlock(rf_preSteady,gz_preSteady);
    seq.addBlock(gzReph_preSteady);
    seq.addBlock(mr.makeDelay(delayTR_preSteady));
    for initial=1:20 %Cover 20 Scans until Steady-State
     rf.signal=-rf.signal; %Change Flip Direction
     seq.addBlock(rf,gz);
     gyPre =mr.makeTrapezoid('y','Area',phaseAreas(initial),'Duration',mr.calcDuration(gxPre),'system',sys);
     seq.addBlock(gxPre,gyPre,gzReph);
     seq.addBlock(mr.makeDelay(delayTE));
     seq.addBlock(gx);
     gyPre.amplitude=-gyPre.amplitude;
     gxPre.amplitude=-gxPre.amplitude;
     seq.addBlock(gxPre,gyPre);
     seq.addBlock(mr.makeDelay(delayTR))
     gxPre.amplitude=-gxPre.amplitude;
    end

    for i=1:Ny % loop over phase encodes and define sequence blocks
         rf.signal=-rf.signal; %Change Flip Direction
         seq.addBlock(rf,gz);
         gyPre = mr.makeTrapezoid('y','Area',phaseAreas(i),'Duration',mr.calcDuration(gxPre),'system',sys);
         seq.addBlock(gxPre,gyPre,gzReph);
         seq.addBlock(mr.makeDelay(delayTE));
         seq.addBlock(gx,adc);
         gyPre.amplitude=-gyPre.amplitude;
         gxPre.amplitude=-gxPre.amplitude;
         seq.addBlock(gxPre,gyPre, mr.makeLabel('SET','LIN', i));
         seq.addBlock(mr.makeDelay(delayTR)); 
         gxPre.amplitude=-gxPre.amplitude;
    end    
end
%% Check whether the timing of the sequence is correct
[ok, error_report]=seq.checkTiming;
if (ok)
 fprintf('Timing check passed successfully\n');
else
 fprintf('Timing check failed! Error listing follows:\n');
 fprintf([error_report{:}]);
 fprintf('\n');
end

%% Plot sequence and k-space diagrams
%seq.plot('timeRange', TI1+TI2+[0 4]*(TR), 'TimeDisp', 'ms', 'label', 'lin');
% % k-space trajectory calculation
% [ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] = 
%[ktraj_adc, t_adc, ktraj, t_ktraj, t_excitation, t_refocusing] =seq.calculateKspacePP();
% % plot k-spaces
% figure; plot(t_ktraj, ktraj'); % plot the entire k-space trajectory
% hold; plot(t_adc,ktraj_adc(1,:),'.'); % and sampling points on the kxaxis
% title('k-space components as functions of time');
%figure; plot(ktraj(1,:),ktraj(2,:),'b'); % a 2D plot
%axis('equal'); % enforce aspect ratio for the correct trajectory display
%hold;plot(ktraj_adc(1,:),ktraj_adc(2,:),'r.'); % plot the sampling points
%title('2D k-space');
% 

rep = seq.testReport;
%fprintf([rep{:}]);

%% Prepare sequence export
seq.setDefinition('Name', 'writebSSFP_draft');
% 
if (ok)
  if balanced == 1
    if inversion == 1
        seq.write(['TI1_',num2str(round(TI1*1e3)),'_TI2_',num2str(round(TI2*1e3)),'_flip',num2str(round(alpha)),'_TR',num2str(round(TR*1e3)),'_FOV',num2str(fov*1e3),'_N',num2str(Nx),'_AdiabaticDIbSSFP.seq']) % Write to pulseq file
    else
        seq.write(['flip',num2str(round(alpha)),'_TR',num2str(round(TR*1e3)),'_TE',num2str(round(TE*1e3)),'_FOV',num2str(fov*1e3),'_N',num2str(Nx),'_AdiabaticDIbSSFP.seq']) % Write to pulseq file
    end
  else
    if inversion == 1
        seq.write(['nb_TI1_',num2str(round(TI1*1e3)),'_TI2_',num2str(round(TI2*1e3)),'_flip',num2str(round(alpha)),'_TR',num2str(round(TR*1e3)),'_FOV',num2str(fov*1e3),'_N',num2str(Nx),'_AdiabaticDIbSSFP.seq']) % Write to pulseq file
    else
        seq.write(['nb_flip',num2str(round(alpha)),'_TR',num2str(round(TR*1e3)),'_TE',num2str(round(TE*1e3)),'_FOV',num2str(fov*1e3),'_N',num2str(Nx),'_AdiabaticDIbSSFP.seq']) % Write to pulseq file
    end      
  end
else
    fprintf('Not Saved')
end

