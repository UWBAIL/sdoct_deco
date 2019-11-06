%% Benchmark SD_OCE system performance
% Test axial resolution
% Output Phase Stability
% Modify to analyze decorrelation statistics
%
% Written by Mitchell A. Kirby
% November 2019

clear all; %close all;


% input prefix
rawname='t';

% crop z-pixel at a specified depth
depth=500;
dz = 0.004593E-3;

% # of A-scans to skip when calculating phase difference (FIR
% differentiator)
frameshift=5;

% # of time-points
fs=46500;
linenum=750;

% # of x-positions
framenum=3;
xloc=2; %location to analyze variance data

% # of active pixels on line-scan camera
pixel=1024;

% A-line scan rate
Fs=46500;

% Display range for log-compressed oct intensity data
imgrg=[2.4 4.5];

% OCT calibration coefficients (dispersion compensation)
%coefs=[4.2271E+1	6.2963E-1	4.3917E-4	-6.6336E-8]; % UW 04/23/2019
coefs=[1.9580E+1	6.3021E-1	6.5232E-4	-2.5845E-7];
coefs=fliplr(coefs);

% Modified Blue-Red colormap to display +_ phase variations
brmap=[0 0.749019622802734 0.749019622802734;0 0.749019622802734 0.749019622802734;0 0.749019622802734 0.749019622802734;0 0.749019622802734 0.749019622802734;0 0.714973270893097 0.760427832603455;0 0.680926918983459 0.771836042404175;0 0.646880567073822 0.78324419260025;0 0.612834215164185 0.79465240240097;0 0.578787863254547 0.806060612201691;0 0.544741570949554 0.817468822002411;0 0.510695219039917 0.828877031803131;0 0.476648837327957 0.840285241603851;0 0.442602515220642 0.851693391799927;0 0.408556163311005 0.863101601600647;0 0.374509811401367 0.874509811401367;0 0.34046345949173 0.885918021202087;0 0.306417107582092 0.897326231002808;0 0.272370785474777 0.908734381198883;0 0.238324418663979 0.920142590999603;0 0.204278081655502 0.931550800800323;0 0.170231729745865 0.942959010601044;0 0.136185392737389 0.954367220401764;0 0.102139040827751 0.965775430202484;0 0.0680926963686943 0.97718358039856;0 0.0340463481843472 0.98859179019928;0 0 1;0 0 0.857142865657806;0 0 0.714285731315613;0 0 0.571428596973419;0 0 0.428571432828903;0 0 0.28571429848671;0 0 0;0 0 0;0 0 0;0.242016807198524 0.0459383763372898 0;0.363025218248367 0.0689075663685799 0;0.484033614397049 0.0918767526745796 0;0.605042040348053 0.114845938980579 0;0.726050436496735 0.13781513273716 0;0.847058832645416 0.160784319043159 0;0.847127020359039 0.161739140748978 0;0.847195208072662 0.162693947553635 0;0.84726345539093 0.163648769259453 0;0.847331643104553 0.164603590965271 0;0.847399830818176 0.165558397769928 0;0.847468018531799 0.166513219475746 0;0.847536206245422 0.167468041181564 0;0.84760445356369 0.168422847986221 0;0.847672641277313 0.169377669692039 0;0.847740828990936 0.170332491397858 0;0.847809016704559 0.171287298202515 0;0.847877204418182 0.172242119908333 0;0.84794545173645 0.173196941614151 0;0.848013639450073 0.174151748418808 0;0.848081827163696 0.175106570124626 0;0.85089510679245 0.214492753148079 0;0.853708446025848 0.253878951072693 0;0.856521725654602 0.293265134096146 0;0.859335064888 0.332651317119598 0;0.862148344516754 0.372037500143051 0;0.864961624145508 0.411423712968826 0;0.867774963378906 0.450809895992279 0;0.87058824300766 0.490196079015732 0;0.87058824300766 0.490196079015732 0];

fn=1;

while exist([rawname,num2str(fn),'noise.oct'],'file'),
    fn=fn+1;
end;

fn=fn-1;


for fileloop=[8:9]
    
    %% Load data
    clear Frame phframe phframe_c
    close all;
    clc;
    
    useref=1;
    refname=['ref_data'];
    lpcontrol=fileloop;
    
    filename=[rawname,num2str(lpcontrol),'noise.oct'];
    disp(['Extracting raw data....',rawname,num2str(lpcontrol)])
    
    % extract .oct data to complex array
    [Frame]=frameextractv5(pixel,coefs,filename,useref,refname,linenum,framenum);
    
    [nz,nt,nx]=size(Frame);
    %% Calculate intensity, phase, and complex arrays
    %save raw data
    %save([rawname,num2str(lpcontrol),'.mat'],'Frame','img','-v7.3');
    
    % log compress complex OCT data to generate OCT intensity array
    img=log10(abs(Frame));
    img_xz=squeeze(img(:,10,:));
    
    % test the intensity signal (SNR,resolution,etc)
    % Note:Need to do some sort of FWHM measurement, calibrate for dz, and
    % qantify noise floor/Signal (SNR)
    ax=img_xz(:,2);
    z=(1:nz)*dz;
   
    %% function to calculate surface for variance measurement
    windowlength=10;
    maxjump=5;,
    minseg=5;
    surface_z = sd_detect_surface(img, img_xz,windowlength, maxjump, minseg);
    
    % Plot the surface detection
    figure;
    subplot(211),imagesc(img_xz),hold on, plot(surface_z,'r.','MarkerSize',18)
    xlabel('x (pixels)')
    ylabel('z (pixels)')
    subplot(234),plot(img_xz(:,1),z*10^3),set(gca,'Ydir','reverse'),
    ylabel('depth (mm)')
    subplot(235),plot(img_xz(:,2),z*10^3),set(gca,'Ydir','reverse'),
    xlabel('log compressed intensity (a.u.)')
    subplot(236),plot(img_xz(:,3),z*10^3),set(gca,'Ydir','reverse'),
    saveas(gcf,[filename(1:end-4),'_fig1.png'])
    
    %% Compute Intensity Variation
    img_surf=img(surface_z(xloc),:,xloc);
    imgdif = img - circshift(img, -frameshift, 2);
    imgdifs=imgdif(surface_z(xloc),:,xloc);
    
    stdevimg=std(imgdifs);
    avgimg = moment(imgdifs,1);
    varianceimg= moment(imgdifs,2);
    
    % Compute phase data
    phraw = angle(Frame);
    
    nt = size(Frame, 2);
    ph = angle(Frame) - circshift(angle(Frame), -frameshift, 2);% Compute phase difference from OCT complex data
    ph(ph > pi) = ph(ph > pi) - 2*pi;
    ph(ph < -pi) = ph(ph < -pi) + 2*pi;
    ph(:, (nt-frameshift+1):nt, :) = 0;
    
    % Compute complex data
    Frame_surf=Frame(surface_z(xloc),:,xloc);
    comp = Frame - circshift(Frame, -frameshift, 2);
    comp_surf=comp(surface_z(xloc),:,xloc);
    comp_surfi=imag(comp_surf);
    comp_surfr=real(comp_surf);
    
    %% Anaylyze phase variance at the surface
    phraw_surf=phraw(surface_z(xloc),:,xloc);
    ph_surf=ph(surface_z(xloc),:,xloc);
    
    stdev=std(ph_surf);
    avg = moment(ph_surf,1);
    variance= moment(ph_surf,2);
    
    time=(1:nt)*1/fs;
    
    %% Make Plots
    fig1=figure;
    set(gcf,'Position',[100 100 1400 700])
    
    subplot(231)
    yyaxis left
    plot(time*10^3,img_surf*1000,'b'),hold on,
    yyaxis right
    plot(time*10^3,imgdifs*1000,'r')
    %ylim([-7 7])
    xlabel('time (ms)')
    ylabel('log comp intensity (a.u.)')
    legend('intensity signal','intensity difference')
    
    subplot(234)
    histogram(imgdifs*1000)
    xlabel('log comp intensity (a.u.)')
    ylabel('counts')
    title(['[mean, stdev, variance] = [',num2str(avgimg*10^3),' , ',num2str(stdevimg*10^3),' , ',num2str(varianceimg*10^3),'] (a.u.)'])
    
    % Plot Phase Signal
    subplot(232)
    plot(time*10^3,phraw_surf,'b.'),hold on,plot(time*10^3,ph_surf,'r.')
    legend('phase signal','phase difference')
    xlabel('time (ms)')
    ylabel('phase (rad)')
    
    subplot(235)
    histogram(ph_surf),
    xlabel('phase (rad)')
    ylabel('counts')
    title(['[mean, stdev, variance] = [',num2str(avg*10^3),' , ',num2str(stdev*10^3),' , ',num2str(variance*10^3),'] mrad'])
    
    %Plot Complex Signal
    subplot(233)
    plot(real(Frame_surf),  imag(Frame_surf),'b.'),hold on,plot(real(comp_surf),  imag(comp_surf),'r.'),
    legend('complex signal','complex difference')
    xlabel('real part')
    ylabel('imaginary part')
    
    subplot(236),
    histogram2(comp_surfr,comp_surfi,'DisplayStyle','tile'),
    xlabel('real part')
    ylabel('imaginary part')
    cm=colorbar
    colormap(jet)
    title(cm,'counts')
    saveas(gcf,[filename(1:end-4),'_fig2.png'])
    %% Calculate decorrelation statistics
    %time- dependent complex-valued autocorrelation of the signal (check
    %this)
    [acf_img,lags_img,bounds_img] = autocorr(imgdifs,nt-1);
    [acf,lags,bounds] = autocorr(comp_surf,nt-1);
    [acf_phase,lags_phase,bounds_phase] = autocorr(ph_surf,nt-1);
    [acf,lags,bounds] = autocorr(comp_surf,nt-1);
    
    figure;
    set(gcf,'Position',[100 100 1000 600])
    plot(time*10^3,abs(acf_img),'r.--','MarkerSize',18),hold on
    plot(time*10^3,abs(acf_phase),'b.--','MarkerSize',18),
    plot(time*10^3,abs(acf),'k.--','MarkerSize',18),
    xlim([0 .5])
    xlabel('lags (time (ms))')
    ylabel('autocorrelation (r)')
    title('autocorrelation function')
    legend('real', 'phase','complex')
    
    saveas(gcf,[filename(1:end-4),'_fig3.png'])
    %%
    
    close all
end
