%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulation Codes for Machine Learning based APT imaging using partially
% synthetic data from tissue mimicking data
% Please uncomment required lines of text
%
% Authors: Malvika Viswanathan, Zhongliang Zu
%
% Correspondance: zhongliang.zu@vumc.org 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
clc;

% Load tissue mimicking data here!! Data should be a numerical array of size [1,69]
% load('...');
% amine_mean = ...; % create AREX quantified amine-CEST spectrum 
% MT_mean =... % create AREX quantified MT spectrum 
% fm =...; % tissue mimicking data

% saturation power along with B1 inhomogeneities -0.1:0.05:0.1uT
tt=1.0; % B1 loop
tt =[0.9, 0.95, 1, 1.05, 1.1];

% sequence parameters
pulseduration=5;
i_SNR=5.0; 
exciteduration=0.002;
exciteangle=90;
TR=0.02;

% saturation offsets 
offppm = 300; % Hz
sep1=3.6*offppm;
sep2=3*offppm;
sep3=2*offppm
sep4=-1.6*offppm;
sep5=-3.3*offppm;

% creating offset frequency array
max=1500;
step=50;
offset= -max:step:max;
k=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000]; %Hz

% adding B0 shifts -40:20:40 Hz
k_7pT = [k-40;k-20;k;k+20;k+40];


% constant sample parameters
R1M=1/1.5;
ksw5 = 20;
R2S5=1/0.0005;

% simulation sample parameters
num_T1W=5;
num_T2W=5;
num_fs1=5;
num_fs2=3;
num_fs5=3;
num_fm=3;
num_ksw1=5;
num_T2S=3;
num_tt = 5;
num_B0 = 5;
 
T1W_matrix=[1.6, 1.8, 2.0, 2.2, 2.4]; %s
T2W_matrix=[40, 60, 80, 100, 120]*0.001; %ms
T2S_matrix=[0.002, 0.004, 0.006]; %s

fs2_matrix=[0.5, 1.0, 1.5]; %r_amine
fm_matrix=[0.8, 1.0, 1.2]; %r_MT

fNOE_matrix=[0.008, 0.01, 0.012];
 
i=1;
 
for ii_T1W=1:num_T1W
    ii_T1W
    R1W_cal=1./T1W_matrix(ii_T1W);
  for ii_T2W=1:num_T2W
      ii_T2W
      R2W_cal=1./T2W_matrix(ii_T2W);
      for ii_T2S=1:num_T2S
          R2S_cal = 1./T2S_matrix(ii_T2S);
              for ii_fs1=1:num_fs1
                fs1_cal=0.0006+0.0002*(ii_fs1-1);
                for ii_fs2=1:num_fs2
                    fs2_cal=fs2_matrix(ii_fs2);
                       for ii_fs5=1:num_fs5
                           fs5_cal=fNOE_matrix(num_fs5);
                           for ii_fm=1:num_fm
                               fm_cal=fm_matrix(ii_fm);
                               for ii_ksw1=1:num_ksw1
                                    ksw1_cal=40+30*(ii_ksw1-1);
                                    for ii_tt = 1:num_tt
                                        tt_shift = tt(ii_tt); 
                                        for ii_kk = 1:num_B0
                                            B0_shift = k_7pT(ii_kk,:);


 R1W_cal_obs=(R1W_cal+fm_cal*fm*R1M)./(1+fm_cal*fm);

 % spectra with B0 and B1 inhomogeneites
 cal_Lorentzian1_cal=(fs1_cal.*ksw1_cal.*(tt_shift.*42.6*2*pi).^2./((tt_shift.*42.6*2*pi).^2+(R2S_cal+ksw1_cal)*ksw1_cal+ksw1_cal./(R2S_cal+ksw1_cal).*((B0_shift+sep1)*2*pi).^2)); %
 cal_Lorentzian2_cal=fs2_cal*interp1(k_7pT(3,:),amine_mean,B0_shift)*((tt_shift)^2); 
 cal_Lorentzian5_cal = (fs5_cal.*ksw5.*(tt_shift.*42.6*2*pi).^2./((tt_shift.*42.6*2*pi).^2+(R2S5+ksw5)*ksw5+ksw5./(R2S5+ksw5).*((B0_shift+sep5)*2*pi).^2));
 cal_Lorentzian6_cal=fm_cal*interp1(k_7pT(3,:),MT_mean,B0_shift)*((tt_shift)^2);
 cal_eff_cal=R1W_cal_obs.*((B0_shift)*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2)+R2W_cal.*(tt_shift.*42.6*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2);
 
 % spectra without B0 and B1 inhomogeneites
 cal_Lorentzian1_cal_1=(fs1_cal.*ksw1_cal.*(tt(3).*42.6*2*pi).^2./((tt(3).*42.6*2*pi).^2+(R2S_cal+ksw1_cal)*ksw1_cal+ksw1_cal./(R2S_cal+ksw1_cal).*((k_7pT(3,:)+sep1)*2*pi).^2)); %
 cal_Lorentzian2_cal_1=fs2_cal*amine_mean';
 cal_Lorentzian5_cal_1 = (fs5_cal.*ksw5.*(tt(3).*42.6*2*pi).^2./((tt(3).*42.6*2*pi).^2+(R2S5+ksw5)*ksw5+ksw5./(R2S5+ksw5).*((k_7pT(3,:)+sep5)*2*pi).^2));
 cal_Lorentzian6_cal_1=fm_cal*MT_mean';
 cal_eff_cal_1=R1W_cal_obs.*((k_7pT(3,:))*2*pi).^2./((tt(3).*42.6*2*pi).^2+((k_7pT(3,:))*2*pi).^2)+R2W_cal.*(tt(3).*42.6*2*pi).^2./((tt(3).*42.6*2*pi).^2+((k_7pT(3,:))*2*pi).^2);

 % inverse summation to create Z spectra with inhomogeneities
 sscal = R1W_cal_obs./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm)+cal_Lorentzian2_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_cal).*(((B0_shift)*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2));
 SS_cal(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1,ii_tt,ii_kk)=R1W_cal_obs./(cal_eff_cal+cal_Lorentzian1_cal./(1+fm_cal*fm)+cal_Lorentzian2_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_cal).*(((B0_shift)*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2));
 SS_cal_value_ref=R1W_cal_obs./(cal_eff_cal+0./(1+fm_cal*fm)+cal_Lorentzian2_cal./(1+fm_cal*fm)+cal_Lorentzian5_cal./(1+fm_cal*fm)+cal_Lorentzian6_cal).*(((B0_shift)*2*pi).^2./((tt_shift.*42.6*2*pi).^2+((B0_shift)*2*pi).^2));


 % inverse summation to create Z spectra with inhomogeneities
 sscal_1 = R1W_cal_obs./(cal_eff_cal_1+cal_Lorentzian1_cal_1./(1+fm_cal*fm)+cal_Lorentzian2_cal_1./(1+fm_cal*fm)+cal_Lorentzian5_cal_1./(1+fm_cal*fm)+cal_Lorentzian6_cal_1).*(((k_7pT(3,:))*2*pi).^2./((tt(3).*42.6*2*pi).^2+((k_7pT(3,:))*2*pi).^2));
 SS_cal_1(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1,ii_tt,ii_kk)=R1W_cal_obs./(cal_eff_cal_1+cal_Lorentzian1_cal_1./(1+fm_cal*fm)+cal_Lorentzian2_cal_1./(1+fm_cal*fm)+cal_Lorentzian5_cal_1./(1+fm_cal*fm)+cal_Lorentzian6_cal_1).*(((k_7pT(3,:))*2*pi).^2./((tt(3).*42.6*2*pi).^2+((k_7pT(3,:))*2*pi).^2));
 SS_cal_value_ref_1=R1W_cal_obs./(cal_eff_cal_1+0./(1+fm_cal*fm)+cal_Lorentzian2_cal_1./(1+fm_cal*fm)+cal_Lorentzian5_cal_1./(1+fm_cal*fm)+cal_Lorentzian6_cal_1).*(((k_7pT(3,:))*2*pi).^2./((tt(3).*42.6*2*pi).^2+((k_7pT(3,:))*2*pi).^2));

 % save essential parameters 
 R1W_cal_matrix(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1,ii_tt,ii_kk)=R1W_cal_obs;
 fm_cal_matrix(ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1,ii_tt,ii_kk)=fm_cal*fm;
 params_matrix(:, ii_T1W,ii_T2W,ii_T2S, ii_fs1, ii_fs2, ii_fs5, ii_fm, ii_ksw1,ii_tt,ii_kk) = [R1W_cal_obs; R2W_cal;R2S_cal; fs1_cal; fs2_cal; fs5_cal; fm_cal; ksw1_cal; tt_shift; B0_shift(35)];

 %CESTR
 mtr =(1-sscal_1) - (1-SS_cal_value_ref_1);
 mtr_amp(1,i) = mtr(1,14);
 mtr_width(1,i)= fwhm(mtr,k_7pT(3,:));

 % AREX
 arex=((1./sscal_1) - (1./SS_cal_value_ref_1)).*(R1W_cal_obs).*(1+(fm_cal*fm));
 arex_amp(1,i) = arex(1,14); % amplitude
 arex_width(1,i)= fwhm(arex,k_7pT(3,:)); %width

 i = i+1;


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

matrix_input_all(:,:)=reshape(SS_cal,  [length(k') num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1*num_tt*num_B0]); 
sz = size(matrix_input_all(:,1));
matrix_input_all_noise= matrix_input_all+ 0.01*(randn(sz));

R1W_cal_matrix_output_all=reshape(R1W_cal_matrix,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1*num_tt*num_B0 1]);
fm_cal_matrix_output_all=reshape(fm_cal_matrix,  [num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1*num_tt*num_B0 1]);

params_matrix = reshape(params_matrix, [10 num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1*num_tt*num_B0]);
matrix_MTR_output(:,1)=reshape(mtr_amp,  [1 num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1*num_tt*num_B0]);   
matrix_MTR_output(:,2)=reshape(mtr_width,  [1 num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1*num_tt*num_B0]);
matrix_AREX_output(:,1)=reshape(arex_amp,  [1 num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1*num_tt*num_B0]);   
matrix_AREX_output(:,2)=reshape(arex_width,  [1 num_T1W*num_T2W*num_T2S*num_fs1*num_fs2*num_fs5*num_fm*num_ksw1*num_tt*num_B0]);