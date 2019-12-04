%% Analyze data from SD OCT system designed to measure autocorrelation statistics
% Load raw data
% Detect Surface
% Simple example to analyze decorrelation statistics
%
% Written by Mitchell A. Kirby
% December 2019

clear all; 

%% System input parameters (from LabView)
% input prefix
rawname='t';
% A-Line scan rate (set on camera/LabView)
fs=90000;
% # of X-positions 
linenum=800;
% # of TIME positions
framenum=800;
% # of active pixels on line-scan camera
pixel=2048;
% OCT calibration coefficients (dispersion compensation)
%coefs=[4.2271E+1	6.2963E-1	4.3917E-4	-6.6336E-8]; % UW 04/23/2019
coefs=[1E+1	6.892E-1	4.521E-4	-1.7E-7]; % SDOCE system 11/14/2019
coefs=fliplr(coefs);

%% Display metrics
% input physical parameter related to axial pixel spacing
dz = 0.004593E-3; %unit of m
% Display range for log-compressed oct intensity data (leave for now)
imgrg=[2.4 4.5];

%x-location to analyze variance data
xloc=50; 
%t-location used to detect surface
tloc=20;
% # of A-scans to skip when calculating phase difference (FIR
% differentiator)
frameshift=1; %frameshift/fs will determine the maximum time lag (unit s)


%% Load data (looped for batch processing)
% Initialize loop to load data in current directory
fn=1;

while exist([rawname,num2str(fn),'.oct'],'file'),
    fn=fn+1;
end;

fn=fn-1;

% Begin main loop
for fileloop=[7]
    clear Frame phframe phframe_c
    close all;
    clc;
    
    %use reference file ('shade sample arm...')
    useref=1;
    refname=['ref_data'];
    lpcontrol=fileloop;
      
    %load data
    filename=[rawname,num2str(lpcontrol),'.oct'];
    disp(['Extracting raw data....',rawname,num2str(lpcontrol)])
    
    % extract .oct data to complex array
    [Frame]=frameextractv5(pixel,coefs,filename,useref,refname,linenum,framenum);
    
    % nx- set by 'linenum' nt- set by 'framenum'
    [nz,nx,nt]=size(Frame);
    
    %save raw data
    %save([rawname,num2str(lpcontrol),'.mat'],'Frame','img','-v7.3');
   
  
    %% function to calculate surface location
    % log compress complex OCT data to generate OCT intensity array
    img=20*log10(abs(Frame)); %structural image!
    
    % save a b-scan image for display
    img_xz=squeeze(img(:,:,tloc));

    % array of z-pixels with unit m
    z=(1:nz)*dz;
    
    % median filter on structural data to smoothen surface signal
    img_xz_filt=medfilt2(img_xz,[1 10]); % 1D filter on each axes preserves edge
    img_xz_filt=medfilt2(img_xz_filt,[20 1]); % 1D filter on each axes preserves edge
    
    windowlength=10; % set length of gaussian derivative to calculate surface signal (10 pixels works well)
    maxjump=5; % max corresponding jump (used only when detecting 2D surface)
    minseg=5;% minimum distance to correct (used only when detecting 2D surface)
    
    % function to locate the surface of the material
    surface_z = sd_detect_surface(img,img_xz_filt,windowlength, maxjump, minseg);
    
    % This looks at a location 10 pixels below the surface
    surface=round(surface_z(xloc))+10; 
    
    % Plot the surface detection
    figure;
    subplot(211),imagesc(img_xz),hold on, plot(surface_z,'r.','MarkerSize',10)
    xlabel('x (pixels)')
    ylabel('z (pixels)')
    subplot(234),plot(img_xz(:,1),z*10^3),set(gca,'Ydir','reverse'),
    ylabel('depth (mm)')
    subplot(235),plot(img_xz(:,2),z*10^3),set(gca,'Ydir','reverse'),
    xlabel('log compressed intensity (a.u.)')
    subplot(236),plot(img_xz(:,3),z*10^3),set(gca,'Ydir','reverse'),
    suptitle('test surface detection')
    %saveas(gcf,[filename(1:end-4),'_fig1.png'])
    
    %% Compute Intensity Variation
    %chose surface location. Will probably do some averaging here later
    img_surf=squeeze(img(surface,xloc,:));
    
    %FIR differentiator
    imgdif = img_surf - circshift(img_surf, -frameshift);
    imgdifs=imgdif;
    %chose surface location. Will probably do some averaging here later
    %imgdifs=imgdif(surface,:,xloc);
    
%     %FIR differentiator
%     imgdif = img - circshift(img, -frameshift, 2);
%     %chose surface location. Will probably do some averaging here later
%     imgdifs=imgdif(surface,:,xloc);
    
    % calculate variance metrics
    stdevimg=std(imgdifs);
    avgimg = moment(imgdifs,1);
    varianceimg= moment(imgdifs,2);
      
    %% Compute Phase Variation
    %raw phase signal
    %phraw = angle(Frame);
    phraw = angle(squeeze(Frame(surface,xloc,:)));
    phraw(phraw > pi) = phraw(phraw > pi) - 2*pi;
    phraw(phraw < -pi) = phraw(phraw < -pi) + 2*pi;
    %ph(:, (nt-frameshift+1):nt, :) = 0;
    phraw((nt-frameshift+1):nt, :) = 0;
    
    % phase difference with basic unwrapping
    %ph = angle(Frame) - circshift(angle(Frame), -frameshift, 2);% Compute phase difference from OCT complex data
    ph = (phraw) - circshift((phraw), -frameshift);% Compute phase difference from OCT complex data
    ph(ph > pi) = ph(ph > pi) - 2*pi;
    ph(ph < -pi) = ph(ph < -pi) + 2*pi;
    %ph(:, (nt-frameshift+1):nt, :) = 0;
    ph((nt-frameshift+1):nt, :) = 0;

    
    %chose surface location. Will probably do some averaging here later
%    phraw_surf=phraw(surface,:,xloc);
    phraw_surf=phraw;
    %ph_surf=ph(surface,:,xloc);
    ph_surf=ph;
    
    %calculate variance metrics
    stdev=std(ph_surf);
    avg = moment(ph_surf,1);
    variance= moment(ph_surf,2);
    
    
    %% Compute complex data (need to re-visit/expand this section)
    Frame_surf=Frame(surface,:,xloc);
    %comp = Frame - circshift(Frame, -frameshift, 2);
    comp = Frame_surf - circshift(Frame_surf, -frameshift);
    comp_surf=comp;%(surface,:,xloc);
    comp_surfi=imag(comp_surf);
    comp_surfr=real(comp_surf);
    
    %% Make Plots of signal histogram distributions
    time=(1:nt)*1/fs;

    fig1=figure;
    set(gcf,'Position',[100 100 1400 700])
    
    subplot(231)
    yyaxis left
    plot(time*10^3,img_surf*1000,'b.'),hold on,
    yyaxis right
    plot(time*10^3,imgdifs*1000,'r.')
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
    ylim([-3.14 3.14])
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
    cm=colorbar;
    colormap(jet)
    title(cm,'counts')
    saveas(gcf,[filename(1:end-4),'_fig2.png'])
    
     %% Simple example to calculate decorrelation statistics
    %we will need to greatly expand the processing in this section. This
    %simply gives us an idea of where to start.....
    
    %time- dependent complex-valued autocorrelation of the signal
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
    
    
    %close all
end
