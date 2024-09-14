clear;
load('testing_data_tm3.mat','matrix_input_all','R1W_cal_matrix_output_all','fm_cal_matrix_output_all')
tt_7pT=1.0;
tt_test = [0.9, 0.95, 1, 1.05, 1.1];
pulseduration=5;

exciteduration=0.002;
exciteangle=90;

gauss=100;

TR=0.02;
offppm= 300; 

max=1500;
step=50;
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


offset= -max:step:max;
k_7pT=[-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000];
k_7pT=k_7pT';
k = [-4000, -3500, -3000, -2500, offset, 2500, 3000,3500,4000];
k_7pT_test = [k-40;k-20;k;k+20;k+40]';
satangle=tt_7pT*42.6*360*pulseduration;
satangle_test = tt_test*42.6*360*pulseduration;

%index = [1,10,12,13,14,25,29,33,36,38,44,51,54,64]+1;
% index = [1,6,8,9,10,11,12,13,14,15,20,22,38,48,51,64,65,66]+1;
% index = [3 ,8, 11, 12, 13, 14, 15, 21, 22, 49, 52, 65]+1;
index = [1, 11, 12, 13, 14, 22, 23, 29, 33, 38, 53, 54, 55]+1;
 for i=1:length(matrix_input_all)

    sig=(1-matrix_input_all(index,i));
    R1W_AREX=R1W_cal_matrix_output_all(i);
    fm_AREX=fm_cal_matrix_output_all(i);
   
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
