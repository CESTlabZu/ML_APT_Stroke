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
% load CEST Z spectrum, R1W and fm values from simulations or measured
% data.
% load('...')

% required initial parameters
max=1500;
step=50;
offset= -max:step:max;
k_7pT=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000];
k_7pT=k_7pT';

 for i=1:length(matrix_input_all)

    sig=(1-matrix_input_all(:,i)); 
    R1W_AREX=R1W_cal_matrix_output_all(i); % loaded R1W values
    fm_AREX=fm_cal_matrix_output_all(i); % loaded PSR values
   
    x =k_7pT(index);
    beta0= [0.9, 0, 420,           0.025, -1050, 150,       0.01, -600, 450,         0.001, 450, 300,         0.02, 1050, 900,       0.1, 0, 7500]; % initial test
    lb=[  0.02, -300, 30,          0, -1200,120,       0, -900, 150,            0, 300, 0,           0, 750, 300,           0, -1200, 3000]; % lower bound
    ub=[ 1, 300,   3000,      0.2, -900, 900,         0.2,-300, 1500,          0.2, 600.001, 450,       1, 1350, 1500,         1, 1200, 30000]; % upper bound
    
    Delta=[1]; 
    options=optimset('lsqcurvefit') ; 
    options = optimset('Display','off','TolFun',1e-8,'TolX',1e-8,'MaxFunEvals',5e4*length(x),'MaxIter',2e5) ;
    
    [beta,resnorm,residual,exitflag,output,lambda,jacobian] = ...
        lsqcurvefit(@matsolv, beta0, x, sig, lb, ub, options, Delta) ;
    
    
    
    % amide
    beta_amide=beta;
    sig_simur_amide=matsolv(beta_amide,x,Delta);
    mor_AREX_amide2(:,i) = spectrum(beta_amide(4), beta_amide(6), -abs(beta_amide(5)));

    beta_amide(4)=0;
    sig_simur_ref_amide=matsolv(beta_amide,x,Delta);
    % save required values!
    mor_MTR_amide(:,i)=(sig_simur_amide-sig_simur_ref_amide);
    mor_AREX_amide(:,i)=(1./(1-sig_simur_amide)-1./(1-sig_simur_ref_amide))*R1W_AREX*(1+fm_AREX);
    
    % amine
    beta_amine=beta;
    sig_simur_amine=matsolv(beta_amine,x,Delta);
    beta_amine(7)=0;
    sig_simur_ref_amine=matsolv(beta_amine,x,Delta);

    mor_AREX_amine(:,i)=(1./(1-sig_simur_amine)-1./(1-sig_simur_ref_amine))*R1W_AREX*(1+fm_AREX);
    mor_MTR_amine(:,i)=(sig_simur_amine-sig_simur_ref_amine);
    
    % MT
    beta_MT=beta;
    sig_simur_MT=matsolv(beta_MT,x,Delta);
    beta_MT(16)=0;
    sig_simur_ref_MT=matsolv(beta_MT,x,Delta);
    
    mor_MT=(sig_simur_MT-sig_simur_ref_MT);
    mor_AREX_MT(:,i) = (mor_MT./(1-mor_MT))*R1W_AREX;
    mor_dir_MT(:,i)=(sig_simur_MT-sig_simur_ref_MT);

    % NOE3p5
    beta_NOE3p5=beta;
    sig_simur_NOE3p5=matsolv(beta_NOE3p5,x,Delta);
    beta_NOE3p5(13)=0;
    sig_simur_ref_NOE3p5=matsolv(beta_NOE3p5,x,Delta);
    
    mor_NOE3p5(:,i)=(1./(1-sig_simur_NOE3p5)-1./(1-sig_simur_ref_NOE3p5))*R1W_AREX*(1+fm_AREX);
    mor_dir_NOE3p5(:,i)=(sig_simur_NOE3p5-sig_simur_ref_NOE3p5);

    % NOE1p6
    beta_NOE1p6=beta;
    sig_simur_NOE1p6=matsolv(beta_NOE1p6,x,Delta);
    beta_NOE1p6(10)=0;
    sig_simur_ref_NOE1p6=matsolv(beta_NOE1p6,x,Delta);
    
    mor_NOE1p6=(1./(1-sig_simur_NOE1p6)-1./(1-sig_simur_ref_NOE1p6))*R1W_AREX*(1+fm_AREX);
    mor_dir_NOE1p6=(sig_simur_NOE1p6-sig_simur_ref_NOE1p6);
   
    sprintf("----------------------- %d",i)
        
 end

 