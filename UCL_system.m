%% OPM with all recons- UCL OPM auditory data
% data/tutorial can be found MNE-Python: https://mne.tools/stable/auto_tutorials/preprocessing/80_opm_processing.html
clear
%% constant variables 
Lin = 6; % Truncation order of the internal VSH basis
Lout = 3; % Truncation order of the external VSH basis
dim_in = (Lin+1)^2 - 1; % Dimension of the internal SSS basis, should be 80
center1= [-0.00350699, 0.01138051, 0.05947857]; 
center2= [-0.00433911, 0.04081329, 0.05194245]; 
%adjuct to device coordinate system
center1 = center1 - [0,0,0.05];
center2 = center2 - [0,0,0.05];
ymin=-1e-10;
ymax=1e-10;


%% read in data
coordsys = 'device'; 
fileraw= '~/UCL_OPM_audtitory_raw.fif';
info = fiff_read_meas_info(fileraw);
nchan=info.nchan;
[raw] = fiff_setup_read_raw(fileraw);
[datan,time] = fiff_read_raw_segment(raw);
for i=1:nchan
    opm_matrixn(:,i)=info.chs(i).loc(1:3,:);
    theta_hatn(:,i)=info.chs(i).loc(4:6,:);
    phi_hatn(:,i)=info.chs(i).loc(7:9,:);
    R_hatn(:,i)=info.chs(i).loc(10:12,:);
end

%% remove all bad NaN channels
k=1;
for i=(1:size(R_hatn,2))
    if isnan(R_hatn(1,i))
        k=k;
    else
        R_hat(:,k)=R_hatn(:,i);
        phi_hat(:,k)=phi_hatn(:,i);
        theta_hat(:,k)=theta_hatn(:,i);
        opm_matrix(:,k)=opm_matrixn(:,i);
        data(k,:)=datan(i,:);
        k=k+1;
    end
end


%% SSS expansions
%calculate single in single out
[Sin,SNin]= Sin_vsh_vv([0,0,0]',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lout);

%calculate multi-vsh in and single-vsh out
%previously 0.005
thresh = 0.05;
%calculate single VSH expansions from two optimized origins
[~,SNin1] = Sin_vsh_vv(center1',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lin);
[~,SNin2] = Sin_vsh_vv(center2',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lin); 
[U,sigma,~] = svd([SNin1 SNin2]);
sig_num = diag(sigma)';
for i=1:size(sig_num,2)
    ratio(i) = sig_num(i)/sig_num(1);
    if ratio(i) >= thresh
        SNin_tot(:,i) = U(:,i);
    end
end

%calculate spheroidal in/out
[semi_major,semi_minor,origin]=find_ellipse_axis(opm_matrix');
[Sin_spm,Sout_spm] = spheroidIN_spheroidOUT(opm_matrix',R_hat',origin,semi_major,semi_minor,Lin,Lout);
for j = 1:size(Sin_spm,2)
  SNin_spm_p(:,j) = Sin_spm(:,j)/norm(Sin_spm(:,j));
end

for j = 1:size(Sout_spm,2)
  SNout_spm_p(:,j) = Sout_spm(:,j)/norm(Sout_spm(:,j));
end


%% reconstrct internal data
%single in, single out
pS_p=pinv([SNin SNout]);
XN_p=pS_p*phi_0p;
data_rec_vsh_p=real(SNin*XN_p(1:size(SNin,2),:));
%multi in, vsh out
pS_multi_vsh_p=pinv([SNin_tot SNout]);   
XN_multi_vsh_p=pS_multi_vsh_p*phi_0p;
data_rec_multi_vsh_p=real(SNin_tot*XN_multi_vsh_p(1:size(SNin_tot,2),:)); 
%spheroidal in, spheroidal out
pS_sph_sph_p=pinv([SNin_spm_p,SNout_spm_p]);  
XN_sph_sph_p=pS_sph_sph_p*phi_0p;
data_rec_sph_sph_p=real(SNin_spm_p*XN_sph_sph_p(1:size(SNin_spm_p,2),:)); 
%spheroidal in, single vsh out
pS_sph_vsh_p=pinv([SNin_spm_p SNout]);   
XN_sph_vsh_p=pS_sph_vsh_p*phi_0p;
data_rec_sph_vsh_p=real(SNin_spm_p*XN_sph_vsh_p(1:size(SNin_spm_p,2),:));

%% check condition numbers
cond_vsh_vsh_p=cond([SNin SNout]);
cond_SNin_p=cond(SNin);
cond_SNin_tot_p = cond(SNin_tot);
cond_SNout_p= cond(SNout);
condition_multi_vsh_p = cond([SNin_tot SNout]);
cond_SNin_spm_p=cond(SNin_spm_p);
cond_SNout_spm_p= cond(SNout_spm_p);
condition_sph_sph_p = cond([SNin_spm_p SNout_spm_p]);
condition_sph_vsh_p = cond([SNin_spm_p SNout]);

%% compare data reconstructions
sVSH_sVSH_p=[SNin SNout];
mVSH_sVSH_p=[SNin_tot SNout];
oid_oid_p=[SNin_spm_p,SNout_spm_p];
oid_sVSH_p=[SNin_spm_p,SNout];

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
