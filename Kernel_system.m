%% OPM with all recons- Kernal OPM Phantom data
% data available OSF: https://osf.io/teygz/

clear
%% constant variables 
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80
center1= [-0.00350699, 0.01138051, 0.05947857]; 
center2= [-0.00433911, 0.04081329, 0.05194245]; 
%adjuct to device coordinate system
center1 = center1 - [0,0,0.05];
center2 = center2 - [0,0,0.05];


%% Kernel opm data: AKCLEE_110 updated April 2024 correct sensor positions
coordsys = 'device'; 
filename= '~/audio_ERF_portal_raw.fif';
info = fiff_read_meas_info(filename);
nchan=info.nchan;
ch_types=ones(nchan,1);
[raw] = fiff_setup_read_raw(filename);
[data,times] = fiff_read_raw_segment(raw);
t_start=206022; %200sec
t_end=260420; %260
phi_raw=data(:,t_start:t_end);
for i=1:nchan
    opm_matrix(:,i)=info.chs(i).loc(1:3,:);
    theta_hat(:,i)=info.chs(i).loc(4:6,:);
    phi_hat(:,i)=info.chs(i).loc(7:9,:);
    R_hat(:,i)=info.chs(i).loc(10:12,:);
end


%read evoked
file= '~/audio_ERF_portal_evoked.fif';
[evoked] = fiff_read_evoked(file);
evoked_data=evoked.evoked.epochs;
phi_0p=evoked_data;
evoked_times = evoked.evoked.times;
time=evoked.evoked.times;

%% SSS expansions
%calculate single in single out
[Sin_p,SNin_p] = Sin_vsh_vv([0,0,0]',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lin);
[Sout_p,SNout_p] = Sout_vsh_vv([0,0,0]',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lout);

%calculate multi-vsh in and single-vsh out
thresh = 0.005;
[SNin_tot_p, SNout] = multi_sss(center1,center2,opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lin, Lout, thresh);

%calculate spheroidal in/out
[semi_major,semi_minor,origin]=find_ellipse_axis(opm_matrix');
[Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(opm_matrix',R_hat',origin,semi_major,semi_minor,Lin,Lout);
for j = 1:size(Sin_spm_p,2)
  SNin_spm_p(:,j) = Sin_spm_p(:,j)/norm(Sin_spm_p(:,j));
end

for j = 1:size(Sout_spm_p,2)
  SNout_spm_p(:,j) = Sout_spm_p(:,j)/norm(Sout_spm_p(:,j));
end


%% reconstrct internal data
%single in, single out
pS_p=pinv([SNin_p SNout_p]);
XN_p=pS_p*phi_0p;
data_rec_vsh_p=real(SNin_p*XN_p(1:size(SNin_p,2),:));
%multi in, vsh out
pS_multi_vsh_p=pinv([SNin_tot_p SNout_p]);   
XN_multi_vsh_p=pS_multi_vsh_p*phi_0p;
data_rec_multi_vsh_p=real(SNin_tot_p*XN_multi_vsh_p(1:size(SNin_tot_p,2),:)); 
%spheroidal in, spheroidal out
pS_sph_sph_p=pinv([SNin_spm_p,SNout_spm_p]);  
XN_sph_sph_p=pS_sph_sph_p*phi_0p;
data_rec_sph_sph_p=real(SNin_spm_p*XN_sph_sph_p(1:size(SNin_spm_p,2),:)); 
%spheroidal in, single vsh out
pS_sph_vsh_p=pinv([SNin_spm_p SNout_p]);   
XN_sph_vsh_p=pS_sph_vsh_p*phi_0p;
data_rec_sph_vsh_p=real(SNin_spm_p*XN_sph_vsh_p(1:size(SNin_spm_p,2),:));


%% check condition numbers
cond_vsh_vsh_p=cond([SNin_p SNout_p]);
cond_SNin_p=cond(SNin_p);
cond_SNin_tot_p = cond(SNin_tot_p);
cond_SNout_p= cond(SNout_p);
condition_multi_vsh_p = cond([SNin_tot_p SNout_p]);
cond_SNin_spm_p=cond(SNin_spm_p);
cond_SNout_spm_p= cond(SNout_spm_p);
condition_sph_sph_p = cond([SNin_spm_p SNout_spm_p]);
condition_sph_vsh_p = cond([SNin_spm_p SNout_p]);


%% plot data to check
ymin=-1.5e-12;
ymax=1.5e-12;
chan_num=3;

figure(2);
hold on;
plot(time, phi_0p(chan_num,:))
plot(time,data_rec_vsh_p(chan_num,:))
plot(time,data_rec_multi_vsh_p(chan_num,:))
plot(time,data_rec_sph_sph_p(chan_num,:))
plot(time,data_rec_sph_vsh_p(chan_num,:))
title('All SSS Methods, Kernel OPM Auidio Evoked, R-hat Sensing')
xlabel('Time')
ylabel('Mag 1 (T)')
ylim([ymin, ymax])
%legend({'Raw Data','VSH/VSH','Spm/Spm'},'location','northwest')
legend({'Raw Data','VSH/VSH','Multi/VSH','Spm/Spm','Spm/VSH'},'location','northwest')
hold off

% Create block plots
figure(6);
hold on
t = tiledlayout(2,3); %3
ax1=nexttile;
plot(ax1,time,phi_0p)
ylim(ax1,[ymin, ymax])
title(ax1, 'Evoked Data')
ax2 = nexttile;
plot(ax2,time,data_rec_vsh_p)
ylim(ax2,[ymin, ymax])
title(ax2,'sVSH/sVSH')
ax3=nexttile;
plot(ax3,time,data_rec_multi_vsh_p)
ylim(ax3,[ymin, ymax])
title(ax3,'mVSH/sVSH')
ax5=nexttile;
plot(ax5,time,data_rec_sph_vsh_p)
ylim(ax5,[ymin, ymax])
title(ax5, 'Spheroid/sVSH')
ax6=nexttile;
plot(ax6,time,data_rec_sph_sph_p)
ylim(ax6,[ymin, ymax])
title(ax6,'Spheroid in/out')
title(t, 'SSS Processed Kernel OPM Auidio Evoked, R-hat Sensing')
xlabel(t,'Time (sec)')
ylabel(t,'Evoked Signal (T)')
% Move plots closer together
t.TileSpacing = 'compact';
hold off


%% SNR calculations (signal/noise)
%noise from time 0.4 to 0.5, indicies (601,701), 30
%peak from time 0.11 to 0.21, indicies (311,411), 20
%If A is a matrix whose columns are random variables and whose rows are observations 
% then S is a row vector containing the standard deviation corresponding to each column
std_noise_raw=mean(std(phi_0p(:,601:701))); %each entry is the std of all the channels at 1 time, then average over all time
std_peak_raw=mean(std(phi_0p(:,311:411)));

SNR_raw = std_peak_raw/std_noise_raw;
SNR_vsh_p = mean(std(data_rec_vsh_p(:,311:411)))/mean(std(data_rec_vsh_p(:,601:701)));
SNR_mvsh_p = mean(std(data_rec_multi_vsh_p(:,311:411)))/mean(std(data_rec_multi_vsh_p(:,601:701)));
SNR_sphsph_p = mean(std(data_rec_sph_sph_p(:,311:411)))/mean(std(data_rec_sph_sph_p(:,601:701)));
SNR_sphvsh_p = mean(std(data_rec_sph_vsh_p(:,311:411)))/mean(std(data_rec_sph_vsh_p(:,601:701)));


%% compare data reconstructions
sVSH_sVSH_p=[SNin_p SNout_p];
mVSH_sVSH_p=[SNin_tot_p SNout_p];
oid_oid_p=[SNin_spm_p,SNout_spm_p];
oid_sVSH_p=[SNin_spm_p,SNout_p];

for i=(1:size(time,2))
    check_data_vsh_vsh_p(i) = subspace(phi_0p(:,i), sVSH_sVSH_p)*180/pi;
    check_data_mvsh_vsh_p(i) = subspace(phi_0p(:,i), mVSH_sVSH_p)*180/pi;
    check_data_oid_oid_p(i) = subspace(phi_0p(:,i), oid_oid_p)*180/pi;
    check_data_oid_vsh_p(i) = subspace(phi_0p(:,i), oid_sVSH_p)*180/pi;
end

check_data_vsh_vsh_pmin = min(check_data_vsh_vsh_p);
check_data_vsh_vsh_pmax = max(check_data_vsh_vsh_p);
check_data_vsh_vsh_pav = mean(check_data_vsh_vsh_p);

check_data_mvsh_vsh_pmin = min(check_data_mvsh_vsh_p);
check_data_mvsh_vsh_pmax = max(check_data_mvsh_vsh_p);
check_data_mvsh_vsh_pav = mean(check_data_mvsh_vsh_p);

check_data_oid_oid_pmin = min(check_data_oid_oid_p);
check_data_oid_oid_pmax = max(check_data_oid_oid_p);
check_data_oid_oid_pav = mean(check_data_oid_oid_p);

check_data_oid_vsh_pmin = min(check_data_oid_vsh_p);
check_data_oid_vsh_pmax = max(check_data_oid_vsh_p);
check_data_oid_vsh_pac = mean(check_data_oid_vsh_p);







