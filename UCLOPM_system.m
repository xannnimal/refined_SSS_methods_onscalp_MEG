%% Load UCL OPM dataset from MNE-Python, then process with each SSS method
% reconstructed data is then saved to ".fif" format to be loaded back into
% Python for topo/evoked plotting using MNE-Python functions
%% constant variables 
clear
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

origin = [0,0,0];
coordsys = 'device'; 
infile = 'C:/Users/xanmc/RESEARCH/UCL_OPM_example/UCL_OPM_audtitory_fix_raw.fif';
outfile = 'C:/Users/xanmc/OneDrive/Documents/MATLAB/SSS/refined_SSS_methods/fileio/UCL_OPM_audtitory_oidvsh_check.fif';

info = fiff_read_meas_info(infile);
nchan=info.nchan;
[raw] = fiff_setup_read_raw(infile);
[data,time] = fiff_read_raw_segment(raw);
for i=1:nchan
    opm_matrix(:,i)=info.chs(i).loc(1:3,:);
    theta_hat(:,i)=info.chs(i).loc(4:6,:);
    phi_hat(:,i)=info.chs(i).loc(7:9,:);
    R_hat(:,i)=info.chs(i).loc(10:12,:);
end

keep_37 = data(37,:);
keep_95 = data(95,:);

%% remove all bad NaN channels
%mark bad chans
bad_chans = [37:40,43,44,55,56,95];
raw.info.bads =[]; %list of bad channel names
k=1;
for i =1:info.nchan
    if ismember(i, bad_chans) 
        raw.info.bads{1, k} = info.ch_names{1,i};
        k=k+1;
    end
end

data(bad_chans,:)=[];
opm_matrix(:,bad_chans)=[];
theta_hat(:,bad_chans)=[];
phi_hat(:,bad_chans)=[];
R_hat(:,bad_chans)=[];

ch_types=ones(size(R_hat,2),1);
phi_0p=data;

%% calculate SSS basis variants
%calculate single in single out
[Sin,SNin]= Sin_vsh_vv([0,0,0]',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lin);
[Sout,SNout] = Sout_vsh_vv([0,0,0]',opm_matrix,theta_hat,phi_hat,R_hat,ch_types,Lout);

%calculate multi-vsh in and single-vsh out
%was 0.005
thresh = 0.005;
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


%% save files
% data rec is now 86xntimes because we removed bad channels, add them back
% as zeros?
data_rec_vsh = zeros(raw.info.nchan,size(data_rec_vsh_p,2));
data_rec_mvsh = zeros(raw.info.nchan,size(data_rec_multi_vsh_p,2));
data_rec_oid = zeros(raw.info.nchan,size(data_rec_sph_sph_p,2));
data_rec_oidvsh = zeros(raw.info.nchan,size(data_rec_sph_vsh_p,2));
k=1;
for i=1:95
    if ismember(i, bad_chans)
        data_rec_vsh(i,:) = data_rec_vsh(i,:);
        data_rec_mvsh(i,:) = data_rec_mvsh(i,:);
        data_rec_oid(i,:) = data_rec_oid(i,:);
        data_rec_oidvsh(i,:) = data_rec_oidvsh(i,:);
    else
        data_rec_vsh(i,:)=data_rec_vsh_p(k,:);
        data_rec_mvsh(i,:)=data_rec_multi_vsh_p(k,:);
        data_rec_oid(i,:)=data_rec_sph_sph_p(k,:);
        data_rec_oidvsh(i,:)=data_rec_sph_vsh_p(k,:);
        k=k+1;
    end

end
data_rec_oidvsh(37,:)=keep_37;
data_rec_oidvsh(95,:)=keep_95;

[outfid,cals] = fiff_start_writing_raw(outfile,raw.info);
fiff_write_raw_buffer(outfid,data_rec_oidvsh,cals);
fiff_finish_writing_raw(outfid);
