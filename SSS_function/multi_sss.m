function [SNin_tot, SNout] = multi_sss(center1,center2,R,EX,EY,EZ,ch_types,Lin, Lout, thresh)
% calculate two-origin mSSS basis
% Xan McPherson, 2024
% give two centers, code calculates two SSS expansions, then combines them
% based on eSSS and SVD concatenation
% INPUT
%   center1, center2: 3x1 (x,y,z) locations of expansion centers
%   R: 3xnchan matrix of sensor locations
%   EX,EY,EZ: three 3xnchan normal vectors for coilori, EZ=sensing
%       direction
%   ch_types: 1xnchan vector of 1's for magnetometers, 0's gradiometers
%   Lin,Lout: vsh truncation order, typically (8,3)
% OUTPUT 
%   SNin_tot: nchan x order matrix of norm combined interior expansion
%   SNout: nchan x order matrix norm exterior expansion

%calculate single VSH expansions from two optimized origins
[~,SNin1] = Sin_vsh_vv(center1',R,EX,EY,EZ,ch_types,Lin);
[~,SNin2] = Sin_vsh_vv(center2',R,EX,EY,EZ,ch_types,Lin); 

%comebine VSH expansions into one interior basis
[U,sigma,~] = svd([SNin1 SNin2]);
sig_num = diag(sigma)';
%keep vectors over a significance value > thresh = 0.005
for i=1:size(sig_num,2)
    ratio(i) = sig_num(i)/sig_num(1);
    if ratio(i) >= thresh
        SNin_tot(:,i) = U(:,i);
    end
end

%calculate exterior VSH basis at origin of system
[~,SNout] = Sout_vsh_vv([0,0,0]',R,EX,EY,EZ,ch_types,Lout);

end