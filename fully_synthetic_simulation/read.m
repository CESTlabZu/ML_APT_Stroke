%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulation Codes for Machine Learning based APT imaging using fully
% synthetic data
% Please uncomment required lines of text
%
% Authors: Malvika Viswanathan, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
% tt mean saturation power (uT) along with B1 inhomogeneities -0.1:0.05:0.1uT
tt_7pT = [0.9, 0.95, 1, 1.05, 1.1];

% sequence parameters
pulseduration=5;
exciteduration=0.002;
exciteangle=90;
gauss=100;


% saturation offsets 
offppm= 300; 
sep1_7pT=3.6*offppm;
sep2_7pT=3*offppm;
sep3_7pT=2*offppm
sep4_7pT=-1.6*offppm;
sep5_7pT=-3.3*offppm;


% relaxations
R1S=1/1.5;
R2S1=1/0.002;
R2S2=1/0.01;
R2S3=1/0.01;
R2S4=1/0.001;
R2S5=1/0.0005;
R1M=1/1.5;
R2M=1/0.00005;

% creating offset frequency array
max=1500;
step=50;
offset= -max:step:max;
k=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000];
k_7pT = [k-40;k-20;k;k+20;k+40]';
satangle=tt_7pT*42.6*360*pulseduration;

% constant sample parameters
fs3=0.0003;
fs4=0.003;
fs5=0.01;
ksw2 =5000;
ksw3=500;
ksw4=50;
ksw5=20;
kmw=25;

% simulation sample parameters
num_T1W=5;
num_T2W=5;
num_fs1=5;
num_fs2=3;
num_fs3=1;
num_fs5=3;
num_fm=3;
num_ksw1=5;
num_T2S=3;
num_tt = 5;
num_B0 = 5;
 

T1W_matrix=[1.6, 1.8, 2.0, 2.2, 2.4]; %s
T2W_matrix=[40, 60, 80, 100, 120]*0.001; %ms
T2S_matrix=[0.002, 0.004, 0.006]; %s

fs2_matrix=[0.5,1,1.5]*0.003;
% uncomment for type 2 simulation
fs3_matrix =0.0003; %[0.0001,  0.0003,  0.0005];
fm_matrix=[0.08, 0.1, 0.12];

i=1;

for ii_T1W=1:num_T1W
    ii_T1W
    R1W=1./T1W_matrix(ii_T1W);
  for ii_T2W=1:num_T2W
      ii_T2W
      R2W=1./T2W_matrix(ii_T2W);
      for ii_T2S=1:num_T2S
          R2S_cal = 1./T2S_matrix(ii_T2S);
              for ii_fs1=1:num_fs1
                fs1=0.0009+0.0004*(ii_fs1-1);
                for ii_fs2=1:num_fs2
                    fs2=fs2_matrix(ii_fs2);
                    for ii_fs3=1:num_fs3
                        fs3=fs3_matrix(ii_fs3);
                           for ii_fs5=1:num_fs5
                               fs5=0.008+0.002*(ii_fs5-1);
                               for ii_fm=1:num_fm
                                   fm=fm_matrix(ii_fm);
                                   for ii_ksw1=1:num_ksw1
                                        ksw1=40+30*(ii_ksw1-1);
                                        for ii_tt=1:num_tt
                                            tt_shift = satangle(ii_tt); 
                                            for ii_kk = 1:num_B0
                                                B0_shift = k_7pT(:,ii_kk); 
                                              
                                 
R1W_cal_obs=(R1W+(fm*R1M))./(1+fm); 
R1W_cal_matrix(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs3, ii_fs5,ii_fm,ii_ksw1,ii_tt,ii_kk)=R1W_cal_obs;  

% spectra with B0 and B1 inhomogeneites
a25mspulse = runsteadysimgauss(ksw1, ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, fs3,fs4, fs5, 1, fm, R1S, R2S_cal, R2S2, R2S3,R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M,sep1_7pT*2*pi,sep2_7pT*2*pi,sep3_7pT*2*pi,sep4_7pT*2*pi, sep5_7pT*2*pi, pulseduration, gauss, tt_shift, 1, 2, 1, .00, 1, 1, B0_shift*2*pi, 1);
a25mspulse_ref = runsteadysimgauss(ksw1,ksw2, ksw3, ksw4, ksw5, kmw, 0, fs2, fs3,fs4, fs5, 1, fm, R1S, R2S_cal, R2S2, R2S3,R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M,sep1_7pT*2*pi,sep2_7pT*2*pi,sep3_7pT*2*pi,sep4_7pT*2*pi, sep5_7pT*2*pi, pulseduration, gauss, tt_shift, 1, 2, 1, .00, 1, 1, B0_shift*2*pi, 1);
disp("step1")
aa_7pT(:,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5,ii_fm,ii_ksw1,ii_tt,ii_kk)=a25mspulse(:,6);
aa_7pT_ref(:,ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5,ii_fm,ii_ksw1,ii_tt,ii_kk)=a25mspulse_ref(:,6);
disp("Zspec")
Slab = a25mspulse(:,6);
Sref = a25mspulse_ref(:,6);

% spectra without B0 and B1 inhomogeneites
a25mspulse_ns = runsteadysimgauss(ksw1, ksw2, ksw3, ksw4, ksw5, kmw, fs1, fs2, fs3,fs4, fs5, 1, fm, R1S, R2S_cal, R2S2, R2S3,R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M,sep1_7pT*2*pi,sep2_7pT*2*pi,sep3_7pT*2*pi,sep4_7pT*2*pi, sep5_7pT*2*pi, pulseduration, gauss, satangle(3), 1, 2, 1, .00, 1, 1, k_7pT(:,3)*2*pi, 1);
a25mspulse_ref_ns = runsteadysimgauss(ksw1,ksw2, ksw3, ksw4, ksw5, kmw, 0, fs2, fs3,fs4, fs5, 1, fm, R1S, R2S_cal, R2S2, R2S3,R2S4, R2S5, R1W_cal_obs, R2W, R1M, R2M,sep1_7pT*2*pi,sep2_7pT*2*pi,sep3_7pT*2*pi,sep4_7pT*2*pi, sep5_7pT*2*pi, pulseduration, gauss, satangle(3), 1, 2, 1, .00, 1, 1, k_7pT(:,3)*2*pi, 1);
Slab_ns = a25mspulse_ns(:,6);
Sref_ns = a25mspulse_ref_ns(:,6);
S0 = 1;

 % save essential parameters 
fm_cal_matrix(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs3,ii_fs5,ii_fm,ii_ksw1,ii_tt,ii_kk)=fm;
params_matrix(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs3,ii_fs5,ii_fm,ii_ksw1,ii_tt,ii_kk) = [R1W_cal_obs; R2W;R2S_cal; fs1; fs2; fs3; fs5; fm; ksw1; tt_shift; B0_shift(35)];

% CESTR
mtr = reshape(((1-Slab_ns)-(1-Sref_ns)), [69 1]);
mtr_amp(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5,ii_fm,ii_ksw1,ii_tt,ii_kk) =mtr(14,:);
mtr_width(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs3,ii_fs5,ii_fm,ii_ksw1,ii_tt,ii_kk)= fwhm(mtr,k_7pT(:,3));

% AREX
arex=((1./Slab_ns) - (1./Sref_ns)).*R1W_cal_obs*(1+fm);
arex_amp(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5,ii_fm,ii_ksw1,ii_tt,ii_kk) = arex(14);
arex_width(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2,ii_fs3, ii_fs5,ii_fm,ii_ksw1,ii_tt,ii_kk)= fwhm(arex(1:30),k_7pT(:,3));

X = sprintf("-------------------------------------%d",i);
disp(X)
i = i+1;
% 
                                      end
                                end
                            end
                        end
                     end
                  end
              end
           end
        end
    end
end

% make matrices for required data for training!
matrix_input_all=reshape(aa_7pT(:,:),  [69, num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_tt*num_B0 ]);   
sz = size(matrix_input_all(:,1));
matrix_input_all_noise = matrix_input_all + 0.01*(randn(sz));
R1W_cal_matrix_output_all=reshape(R1W_cal_matrix,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_tt*num_B0  1]);
fm_cal_matrix_output_all=reshape(fm_cal_matrix,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_tt*num_B0  1]);

matrix_output1(:,1)=reshape(mtr_amp,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_tt*num_B0 1]);     
matrix_output1(:,2)=reshape(mtr_width,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_tt*num_B0  1]);

matrix_output2(:,1)=reshape(arex_amp,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_tt*num_B0  1]);     
matrix_output2(:,2)=reshape(arex_width,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs3*num_fs5*num_fm*num_ksw1*num_tt*num_B0  1]);

