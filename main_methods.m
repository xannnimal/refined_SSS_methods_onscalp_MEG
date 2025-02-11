%% Main code for refined SSS methods paper
% Alexandria McPherson Feb 2025
% see README for details on sensor systems/data sets from Nottingham, etc
% designed to comment certain section out depending on the desired
% combination of sensor systems (SQUID, Kernel OPM, QuSpin, etc) and type of
% data (simulated mag dipole, current dipole, combo, or raw data)
% function 'multi_sss.m' is implementation of novel mSSS method
% Spheroidal harmonic functions can be found using README info
clear

%% constant variables 
coordsys = 'device';
magscale = 100;
Lin = 8; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80


%% Sensor systems
%306cnah SQUID system 
% rawfile = 'sample_audvis_raw.fif';
% [R,EX,EY,EZ] = fiff_getpos(rawfile, coordsys);
% point_mags = 0; %0: account for grads and mags, %1: only mags

%102chan SQUID mag system 
% rawfile = 'sample_audvis_raw.fif';
% [Rf,EXf,EYf,EZf] = fiff_getpos(rawfile, coordsys);
% point_mags = 1; %0: account for grads and mags, %1: only mags
% k=1;
% for i=(1:size(Rf,2))
%     if mod(i,3)==0
%         R(:,k)=Rf(:,i);
%         EX(:,k)=EXf(:,i);
%         EY(:,k)=EYf(:,i);
%         EZ(:,k)=EZf(:,i);
%         k=k+1;
%     end
% end


% % Kernel opm data
% filename= '~/audio_ERF_portal_raw.fif';
% info = fiff_read_meas_info(filename);
% nchan=info.nchan;
% for i=1:nchan
%     R(:,i)=info.chs(i).loc(1:3,:);
%     EX(:,i)=info.chs(i).loc(4:6,:);
%     EY(:,i)=info.chs(i).loc(7:9,:);
%     EZ(:,i)=info.chs(i).loc(10:12,:);
% end
% point_mags = 1;

% % QuSpin system from Nottingham
% nroom=2; % nroom=1 is our 5-layer MSR, nroom=2 is the lightweight MSR
% F.OPM_data.path=[pwd,'\']; % path to files
% nroom =2;
% if nroom==1
%     bad_chans=[7:9,34,35,109:111,130:132,178:180]; % main room bad channels
%     F.OPM_data.name1='N1_noise_main.lvm'; % filename
% elseif nroom==2
%     bad_chans=[61:63,175:180]; % light room bad channels
%     F.OPM_data.name1='N1_noise_light.lvm';
% end
% d=get_N1_data_1(F); % load data and extract trials
% pos=d.Sens_pos; % channel positions (64 triaxial sensors each of the three channels per sensor has the same position)
% pos=pos-mean(pos,1); % make origin centre of mass of helmet
% ors=d.Sens_ors; % channel orientations
% pos(bad_chans,:)=[];
% ors(bad_chans,:)=[];
% R = pos';
% EX = ors'; EY=ors'; EZ=ors';
% point_mags = 1;

%% plot sensor helmet
figure(1)
hold on
%quiver3(R(:,1),R(:,2),R(:,3),theta_hat(:,1),theta_hat(:,2),theta_hat(:,3))
scatter3(R(1,:),R(2,:),R(3,:),12,'filled')
% scatter3(center1(:,1),center1(:,2),center1(:,2),'*','g')
% scatter3(center2(:,1),center2(:,2),center2(:,2),'*','g')
grid on
rotate3d
view(135, 20);
hold off


%% specify point mags or no
k=1;
if point_mags ==0 
    for i=(1:size(EZ,2))
        if mod(i,3)==0 %every third is a magnetometer
            ch_types(i)=1;
            mags(k)=i;
            k=k+1;
        else
            ch_types(i)=0;
            k=k;
        end
    end
else
    for i=(1:size(EZ,2))
        ch_types(i)=1;
        mags(i)=i;
    end
end


%% SSS methods
%VSH SSS
[Sin,SNin] = Sin_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lout);

% spheroidal harmonics
[semi_major,semi_minor,origin]=find_ellipse_axis(R');
[Sin_spm_p,Sout_spm_p] = spheroidIN_spheroidOUT(R',EZ',origin,semi_major,semi_minor,Lin,Lout);
for j = 1:size(Sin_spm_p,2)
  SNin_spm(:,j) = Sin_spm_p(:,j)/norm(Sin_spm_p(:,j));
end
for j = 1:size(Sout_spm_p,2)
  SNout_spm(:,j) = Sout_spm_p(:,j)/norm(Sout_spm_p(:,j));
end


%calculate multi-vsh in and single-vsh out
thresh = 0.005;
% % centers for Nottingham (x-y transposed)
% center1= [0.01138051,-0.00350699, 0.05947857] - [0,0,0.05]; 
% center2= [0.04081329,-0.00433911, 0.05194245] - [0,0,0.05];
% % centers for else
center1= [-0.00350699, 0.01138051, 0.05947857] - [0,0,0.05]; 
center2= [-0.00433911, 0.04081329, 0.05194245] - [0,0,0.05];
[SNin_tot, SNout] = multi_sss(center1,center2,R,EX,EY,EZ,ch_types,Lin, Lout, thresh);

%create SSS method matricies
sVSH_sVSH=[SNin SNout];
mVSH_sVSH=[SNin_tot SNout];
oid_oid=[SNin_spm SNout_spm];
oid_sVSH=[SNin_spm SNout];

%% generate time dependent dipoles
rs=[0,0,0];
q=[0,1,0];
r0=[0.05,0,0]; %5cm along x axis
dip_mom_out=[1,0,0];
dip_pos_out = [0,0,1.5]; %1.5 meters
%add time dependence to dipole moment
f_start = 8; % start frequency
f_end = 3; % end frequency
f_start_out = 5; % start frequency
f_end_out = 3; % end frequency
timestep = 0.001;
T = 1.00;
rate_of_change = (f_start - f_end)/T;
rate_of_change_out=(f_start_out-f_end_out)/T;
times = timestep:timestep:T;
for i=(1:3)
    q_t(i,:) = q(i)*sin(2*pi*(f_start*times - times.^2*rate_of_change/2));
    dip_mom_t_out(i,:) = dip_mom_out(i)*sin(2*pi*(f_start_out*times - times.^2*rate_of_change_out/2))*5e8;
end

for i=(1:size(times,2))
    phi_out(:,i) = magneticDipole(R,EX,EY,EZ,dip_pos_out',dip_mom_t_out(:,i),ch_types)'*1e-12;
    phi_in(:,i) = dipole_field_sarvas(rs',q_t(:,i),r0',R,EX,EY,EZ,mags)'*1e-12;

end
rng("default") %set seed
noise = randn(size(phi_in,1),size(phi_in,2));
amplitude = 0.15 * phi_in;
%%%% modify this line to do only internal, in+ext, or add noise %%%
phi_0 = phi_in + phi_out + amplitude.* noise; %
% apply magscale to simulated data
if point_mags ==0
    for i=(1:size(phi_0,1))
        if mod(i,3)==0 %every third is a magnetometer
            phi_0(i,:)=phi_0(i,:)*magscale;
            phi_in(i,:)=phi_in(i,:)*magscale;
            phi_out(i,:)=phi_out(i,:)*magscale;
        else
            phi_0(i,:)=phi_0(i,:);
            phi_in(i,:)=phi_in(i,:);
            phi_out(i,:)=phi_out(i,:);
        end
    end
else %if point_mags=1, all are mags, need to scale all 
    for i=(1:size(phi_0,1))
        phi_0(i,:)=phi_0(i,:)*magscale;
        phi_in(i,:)=phi_in(i,:)*magscale;
        phi_out(i,:)=phi_out(i,:)*magscale;
    end
end


%% Data reconstruction
%single in, single out
pS=pinv(sVSH_sVSH);
XN=pS*phi_0;
data_rec_vsh=real(SNin*XN(1:size(SNin,2),:));
%multi in, vsh out
pS_multi_vsh=pinv(mVSH_sVSH);   
XN_multi_vsh=pS_multi_vsh*phi_0;
data_rec_multi_vsh=real(SNin_tot*XN_multi_vsh(1:size(SNin_tot,2),:)); 
%spheroidal in, spheroidal out
pS_sph_sph=pinv(oid_oid);   
XN_sph_sph=pS_sph_sph*phi_0;
data_rec_sph_sph=real(SNin_spm*XN_sph_sph(1:size(SNin_spm,2),:)); 
%spheroidal in, single vsh out
pS_sph_vsh=pinv(oid_sVSH);   
XN_sph_vsh=pS_sph_vsh*phi_0;
data_rec_sph_vsh=real(SNin_spm*XN_sph_vsh(1:size(SNin_spm,2),:));

%% angles between in and out basis
a_sVSHin_sVSHout = subspace(SNin,SNout)*180/pi;
a_mVSHin_sVSHout = subspace(SNin_tot, SNout)*180/pi;
a_spherin_spherout = subspace(SNin_spm, SNout_spm)*180/pi;
a_spherin_sVSHout = subspace(SNin_spm, SNout)*180/pi;


%% pairwise angles between all collumns in each basis
for i=1:size(SNin,2)
    for j=1:size(SNin,2)
        pa_SNin(i,j)=subspace(SNin(:,i),SNin(:,j));
        pa_SNin_spm(i,j)=subspace(SNin_spm(:,i),SNin_spm(:,j));
    end
end
for i=1:size(SNin_tot,2)
    for j=1:size(SNin_tot,2)
        pa_SNin_tot(i,j)=subspace(SNin_tot(:,i),SNin_tot(:,j));
    end
end
%calculate frobenius norm
pa_SNin_n = norm(pa_SNin,"fro");
pa_SNin_spm_n = norm(pa_SNin_spm,"fro");
pa_SNin_tot_n = norm(pa_SNin_tot,"fro");

%vizualize: for Appendix figures
% figure(2)
% hold on
% imagesc(pa_SNin)
% axis([0 inf 0 inf])
% title("SNin pairwise subspace angles")
% colorbar
% hold off
% 
% figure(3)
% hold on
% imagesc(pa_SNin_tot)
% axis([0 inf 0 inf])
% title("SNin,tot pairwise subspace angles")
% colorbar
% hold off
% 
% figure(4)
% hold on
% imagesc(pa_SNin_spm)
% axis([0 inf 0 inf]);
% title("SNin,spm pairwise subspace angles")
% colorbar
% hold off


%% check condition numbers
cond_vsh_vsh=cond(sVSH_sVSH);
cond_SNin=cond(SNin);
cond_SNin_tot = cond(SNin_tot);
cond_SNout= cond(SNout);
cond_multi_vsh = cond(mVSH_sVSH);
cond_SNin_spm=cond(SNin_spm);
cond_SNout_spm= cond(SNout_spm);
cond_sph_sph = cond(oid_oid);
cond_sph_vsh = cond(oid_sVSH);



%% subsapce angles
for i=(1:size(times,2))
    check_data_vsh_vsh_d(i) = subspace(phi_0(:,i), sVSH_sVSH)*180/pi;
    check_data_mvsh_vsh_d(i) = subspace(phi_0(:,i), mVSH_sVSH)*180/pi;
    check_data_oid_oid_d(i) = subspace(phi_0(:,i), oid_oid)*180/pi;
    check_data_oid_vsh_d(i) = subspace(phi_0(:,i), oid_sVSH)*180/pi;
end
check_data_vsh_vsh_dmin = min(check_data_vsh_vsh_d);
check_data_vsh_vsh_dmax = max(check_data_vsh_vsh_d);
check_data_vsh_vsh_dav = mean(check_data_vsh_vsh_d);

check_data_mvsh_vsh_dmin = min(check_data_mvsh_vsh_d);
check_data_mvsh_vsh_dmax = max(check_data_mvsh_vsh_d);
check_data_mvsh_vsh_dav = mean(check_data_mvsh_vsh_d);

check_data_oid_oid_dmin = min(check_data_oid_oid_d);
check_data_oid_oid_dmax = max(check_data_oid_oid_d);
check_data_oid_oid_dav = mean(check_data_oid_oid_d);

check_data_oid_vsh_dmin = min(check_data_oid_vsh_d);
check_data_oid_vsh_dmax = max(check_data_oid_vsh_d);
check_data_oid_vsh_dav = mean(check_data_oid_vsh_d);


% plot data - simulated results
% chan_num =6;
% t_max = 1000;
% figure(5);
% hold on;
% plot(times(:,:),phi_in(chan_num,:),'LineWidth',1)
% plot(times(:,:),phi_out(chan_num,:),'--','LineWidth',1)
% plot(times(:,:),data_rec_vsh(chan_num,:),'LineWidth',1)
% plot(times(:,:),data_rec_multi_vsh(chan_num,:), 'LineWidth',1)
% title('183-Chan QuSpin system, current dipole 5cm x, exterior dipole 1.5m z')
% xlabel('Time (sec)')
% ylabel('Dipole Signal (T) chan 6')
% %ylim([-8e-12 8e-12])
% legend({'B-in', 'B-out','sVSH/sVSH','mVSH/sVSH'},'location','northwest')
% hold off