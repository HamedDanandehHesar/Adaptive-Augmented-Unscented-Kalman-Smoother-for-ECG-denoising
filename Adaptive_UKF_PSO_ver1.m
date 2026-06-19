% cubature Kalman Smoother with non-additive noise (special case of ukf with alpha = 1 beta and kappa =0)
close all

tic
clear
 
%% Load ECG data from MAT file
[file, path] = uigetfile('*.mat','Select ECG mat file');
data = load(file);


fs = data.fs;                % Sampling frequency of ECG signal (Hz)

x = data.x;               % Extract stored signal matrix
ecg = x(1,:);             % Use first channel as ECG signal
length_sig = length(ecg); % Total number of ECG samples

ecg_bins = round(fs/2);           % Number of phase bins for ECG mean calculation

dt = 1/fs; % time step

%% adding white gaussian Noise
SNR = 6;   % you can change the Noise SNR level
x_noisy = zeros(3,length(x(1,:)));
x_noisy(1,:) = awgn(x(1,:),SNR,'measured');
%% -------- R-peak detection using Pan–Tompkins algorithm
[qrs_positions] = pantompkins_qrs(x_noisy(1,:),fs);
figure(1),plot(x_noisy(1,:),'b'),hold on,plot(qrs_positions,x_noisy(1,qrs_positions),'*r'),hold off
legend({'Noisy ECG Signal','R Peaks'})
title([file '   at SNR = ' num2str(SNR)])
axis('tight')
%% -------- Phase calculation
% Linear phase based on RR intervals
[Linearphase,~] = calculate_linear_phase_ver2(qrs_positions,length_sig,fs);
x_noisy(2,:) = Linearphase; 



[ECGsd,ECGmean,meanphase] = ecgsd_extractor_ver1(x_noisy(1,:),Linearphase,ecg_bins);

% further smoothing of ECG mean using wavelet
ECGmean = wdenoise(ECGmean,5,Wavelet="bior4.4",DenoisingMethod="BlockJS");



%% -------- ECG parameter extraction using Gaussian mixture model

MaxNumGaussian = 50;   % Maximum number of Strongest Gaussian components

% ========================building new myfun based on L Gaussians
L_num_of_Gaussian_kernels = 50;
ecg_mean_temp = 0;
ai = [];
bi = [];
tetai  = [];
for i=1:L_num_of_Gaussian_kernels
% disp(num2str(i))
ecg_mean_temp1 = ECGmean - ecg_mean_temp;
lb = [-1.5*max(ecg_mean_temp1).*ones(1,1)   0.000001*ones(1,1)   (-pi+.014)*ones(1,1)  ];
ub = [(1.5*max(ecg_mean_temp1)).*ones(1,1)  5*ones(1,1)  (pi-.014)*ones(1,1)  ];  
myfun1 = @(params)  norm(ecg_mean_temp1'-sum((repmat(params(1:1),ecg_bins,1).*exp(-(rem(repmat(meanphase,1,1)'-repmat(params(3),ecg_bins,1)+pi,2*pi)-pi) .^2 ./ (2*(repmat(params(2),ecg_bins,1)) .^ 2))),2));


% options = optimoptions('particleswarm','SwarmSize',30,'HybridFcn',@fmincon,'MaxIter',1000);
options = optimoptions('particleswarm','SwarmSize',50,'MaxIter',100,'Display','off');

OptimumParams = particleswarm(myfun1,3*1,lb,ub,options);

% L = (length(OptimumParams)/3);

ai_1 = OptimumParams(1);
bi_1 = OptimumParams(2);
tetai_1 = OptimumParams(3);
ai = [ai ai_1];
bi = [bi bi_1];
tetai  = [tetai tetai_1];
dtetai_1 = rem(meanphase - tetai_1 + pi,2*pi)-pi;
ecg_mean_temp = ecg_mean_temp + ai_1 .* exp(-dtetai_1 .^2 ./ (2*bi_1 .^ 2));
figure(41),plot(ecg_mean_temp,'b'),hold on,plot(ECGmean,'r')
legend({'Synthetic ECG','ECG Mean'}),hold off
title([num2str(i) 'th' '  Gaussian Kernel found'])
axis tight
% pause(3)
end

%% selection of the Strongest Peaks
[~,indx_strongest_peaks] = sort(abs(ai),'descend');

ai = ai(indx_strongest_peaks(1:MaxNumGaussian));
bi = bi(indx_strongest_peaks(1:MaxNumGaussian));
tetai = tetai(indx_strongest_peaks(1:MaxNumGaussian));


Alpha_i = ai;   % Gaussian amplitudes
Beta_i  = bi;   % Gaussian widths
Theta_i = tetai;   % Gaussian centers


% Sorting of parameters from based on tetai from -pi to pi

[Theta_i,idx] = sort(Theta_i,'ascend');
Alpha_i = Alpha_i(idx);
Beta_i = Beta_i(idx);
OptimumParams = [Alpha_i Beta_i Theta_i];
params = OptimumParams;
size_params = length(params);


%%  ANGULAR FREQUENCY MEASUREMENT

ind=1*qrs_positions;
ind2 = ind-[0 ind(1:end-1)];
RR = mean(ind2(1,2:end)); % mean of RR Intervals
w=2*pi*fs/RR; % angular frequency
RR_var=std(2*pi*fs./ind2(1,2:end));% standard deviation of RR Intervals







RR = mean(diff(ind(2:end-1)));

stepteta=2*pi/RR;
w_1 = fs*stepteta;
for j=ind(1,1):-1:1
    
     x_noisy(3,j) = w_1;

end   
for i=1:length(ind)-1



bins = ind(1,i+1)-ind(1,i);

 stepteta = 2*pi/(bins);
w_1 = fs*stepteta;
for j=ind(1,i)+1:ind(1,i+1)
 
    x_noisy(3,j) = w_1;

   
end

end
stepteta=2*pi/RR;
teta= 0;
w_1 = fs*stepteta;

for j=min(ind(end,end)+1,size(x_noisy,2)):size(x_noisy,2)

        x_noisy(3,j) = w_1;

 
end



%% CKF3 Initializations
ai = params(1,1:size_params/3);
bi = params(1,size_params/3+1:2*size_params/3);
tetai = params(1,2*size_params/3+1:size_params);
 
alpha= 1; 
kappa_state=  0;
kappa_measure = kappa_state;
beta= 2;2.5; % The parameter beta is an added degree of freedom to include a priori knowledge on the 
         %original PDF. In the case of a Gaussian PDF, beta = 2 is optimal.

RR_var = var(x_noisy(3,:));

R = diag([ 0.1*mean(ECGsd)^2 (.001) RR_var ]);

% Q = diag( [   .01*R(1,1)  (.1) (.1*RR_var)^2  (.05.*ai.*ones(1,size_params/3)).^2 (.05*bi.*ones(1,size_params/3)).^2 (.005*tetai.*ones(1,size_params/3)).^2  ] );

Q = diag( [   .01*R(1,1)  (.1) (.1*RR_var)^2  (.05.*1.*ones(1,size_params/3)).^2 (.0005.*ones(1,size_params/3)).^2 (.00005.*ones(1,size_params/3)).^2  ] );


y = [x_noisy(1,:);x_noisy(2,:);x_noisy(3,:)]; % noisy measurements 



X0=[x_noisy(1,1);x_noisy(2,1);x_noisy(3,1);0;0;0;ai';bi';tetai']; %initial state; 

P0 = diag([(1*max(abs(x(1,:)'))).^2   (2*pi).^2  w.^2 ]);
Pp0  = repmat(P0,[1,1,size(x_noisy,2)]);
Pp =Pp0;








 
Xukf_update = zeros(3,length(x_noisy(1,:)));
Xukf_predict = zeros(size(Xukf_update));
% UKF parameters
N_state=6+size_params; %space dimension
N_measure = 3; % measurement Dimension
Num_sigma_points = 2*N_state+1; %number of sigma points

lambda_state= ((alpha^2)*(N_state+kappa_state))-N_state; % lambda is a scaling parameter
lambda_measure= ((alpha^2)*(N_measure+kappa_state))-N_measure; % lambda is a scaling parameter

%The parameter ? determines the spread of the sigma points around the centre, and
%usually takes a small positive value equal or less than one. The parenthesis term
% (N + kappa) is usually equal to 3


%% weights of sigma points
% w0m = ?/(N + ?) ; w0c = w0m + (1 ? ?^2 + ?)  
% wi = 1/2(N + ?) , i = 1, ..., 2N
aab=(1-(alpha^2)+beta);
lN_state=lambda_state+N_state; 
LaN_state=lambda_state/lN_state; 
w0c_state =aab+LaN_state;

lN_measure=lambda_measure+N_measure; 
LaN_measure=lambda_measure/lN_measure; 


w0m_state = lambda_state/(lambda_state+N_state);
wm_state = 1/(2*lambda_state+2*N_state);
w0c_state = lambda_state/(lambda_state+N_state)+(1-alpha^2+beta);
wc_state = 1/(2*lambda_state+2*N_state);

w0m_measure = lambda_measure/(lambda_measure+N_measure);
wm_measure = 1/(2*lambda_measure+2*N_measure);
w0c_measure = lambda_measure/(lambda_measure+N_measure)+(1-alpha^2+beta);
wc_measure = 1/(2*lambda_measure+2*N_measure);




X_initial=[x_noisy(1,1);x_noisy(2,1);x_noisy(3,1);0;0;0;ai';bi';tetai']; %initial state


X00=X_initial; % initial value of filter state
X_intermediate=X0; %initial intermediate state

X_sigma_state=zeros(length(X_initial),2*N_state+1); %space for sigma points for state space model
X_sigma_measure=zeros(3,2*N_measure); %space for sigma points for measurement model

X_propagate_sigma=zeros(length(X0),2*N_state+1); %space for propagated sigma points
y_propagate_sigma=zeros(3,2*N_measure); %%space for propagated measurement equation on sigma points
P_update_matrix(:,:,1)= diag([max(abs(x_noisy(1,:)))   (.2*pi).^2  .2*w.^2]); 
P0 = P_update_matrix(:,:,1); 

P00 = diag([P_update_matrix(1,1,1) P_update_matrix(2,2,1) P_update_matrix(3,3,1) .01*R(1,1) (.1) (.1*RR_var)^2  ...
    (.05.*ones(1,size_params/3)).^2 (.005.*ones(1,size_params/3)).^2 (.005.*ones(1,size_params/3)).^2  ] );



X_propagate_sigma_matrix = zeros(length(X_initial),2*N_state+1,length(y(1,:)));
X_propagate_predict_sigma_matrix = zeros(length(X0),2*N_state+1,length(y(1,:)));

P_predict_matrix = P_update_matrix(:,:,1);
P_cross_covariance_xk_xk_plus_1 = P_update_matrix(:,:,1); 

%% Covariance Adaptation matrix parameters
Memory_R = [];
window_size = round(fs/3);
forgetting_factor = 0.99;
%% UKF filtering
for ii=1:size(Xukf_update,2)
    disp(['Forward UKF- Processing sample number:' num2str(ii) '/' num2str(size(Xukf_update,2))])
    
    
     X0 = X0(1:3);
    Xukf_predict(:,ii)= X0;
    P_predict_matrix(:,:,ii) = P0;
    M(:,:,ii) = P0;
    ym=y(:,ii); %measurement
    %propagation of sigma points (measurement)
    
%sigma points
Square_root_P0=(chol((lN_measure*P0))); %matrix square root
X_sigma_measure(:,7)=X0;
X_sigma_measure(:,1)=X0+Square_root_P0(1,:)'; 
X_sigma_measure(:,2)=X0+Square_root_P0(2,:)';
X_sigma_measure(:,3)=X0+Square_root_P0(3,:)';
X_sigma_measure(:,4)=X0-Square_root_P0(1,:)'; 
X_sigma_measure(:,5)=X0-Square_root_P0(2,:)';
X_sigma_measure(:,6)=X0-Square_root_P0(3,:)';    
mm=1:(2*N_measure+1);
y_propagate_sigma(:,mm)= X_sigma_measure(:,mm);

%measurement mean
y_propogate_sigma_mean=0;
for mm=1:2*N_measure
y_propogate_sigma_mean= y_propogate_sigma_mean+wm_measure*y_propagate_sigma(:,mm);
end;
% y_propogate_sigma_mean = y_propogate_sigma_mean/(2*lN_measure);
% y_propogate_sigma_mean = y_propogate_sigma_mean/(2*N+1);
y_propogate_sigma_mean = y_propogate_sigma_mean+(w0m_measure*y_propagate_sigma(:,7));
%measurement cov.
P_xy=0;
for mm=1:6,
P_xy=P_xy+wc_measure*((y_propagate_sigma(:,mm)-y_propogate_sigma_mean)*(y_propagate_sigma(:,mm)-y_propogate_sigma_mean)');

end;
% P_xy=P_xy/(2*lN_measure);
% P_xy=P_xy+(w0c_state*((y_propagate_sigma(:,7)-y_propogate_sigma_mean)*(y_propagate_sigma(:,7)-y_propogate_sigma_mean)'));
P_xy=P_xy+(w0c_measure*((y_propagate_sigma(:,7)-y_propogate_sigma_mean)*(y_propagate_sigma(:,7)-y_propogate_sigma_mean)'));

Syy=P_xy+R;

%cross cov
P_xy=0;
for mm=1:6,
% P_xy=P_xy+((X_propagate_sigma(:,mm)-X_intermediate(:))*(y_propagate_sigma(:,mm)-y_propogate_sigma_mean)');
P_xy=P_xy+wc_measure*((X_sigma_measure(:,mm)-X0)*(y_propagate_sigma(:,mm)-y_propogate_sigma_mean)');

end;
% P_xy=P_xy/(2*lN_measure);
P_xy=P_xy+(w0c_measure*((X_sigma_measure(:,end)-X0)*(y_propagate_sigma(:,end)-y_propogate_sigma_mean)'));
% Sxy=0.5*(P_xy+P_xy');
Sxy= P_xy;

%Kalman gain, etc.
K =Sxy*inv(Syy);
P_update_matrix(:,:,ii)=M(:,:,ii)-(K*Syy*K');
X0 = X_intermediate(1:3)+(K*(ym-y_propogate_sigma_mean)); %estimated (a posteriori) state

Xukf_update(:,ii) = X0;
Memory_R = [Memory_R (ym-y_propogate_sigma_mean)];
error = (ym-y_propogate_sigma_mean);
if ii>(window_size+length(X0)+1)
    R(1,1) = (1-forgetting_factor)*mean(Memory_R(1,end-window_size+1:end-1).^2)+forgetting_factor*R(1,1);

    temp = K*(error*error')*K';
    Q(1,1)= (1-forgetting_factor)*Q(1,1)+(forgetting_factor)*temp(1,1);  % it is working
    Memory_R(:,1) = [];
 
end


%Prediction
%sigma points


P00 = diag([P_update_matrix(1,1,ii) P_update_matrix(2,2,ii) P_update_matrix(3,3,ii) Q(1,1)  (.1) (.1*RR_var)^2 ...
    (.05.*ones(1,size_params/3)).^2 (.0005*bi.*ones(1,size_params/3)).^2 (.00005*tetai.*ones(1,size_params/3)).^2]);



P00(1:3,1:3) = P_update_matrix(1:3,1:3,ii);


ai = params(1,1:size_params/3);
bi = params(1,size_params/3+1:2*size_params/3);
tetai = params(1,2*size_params/3+1:size_params);
X00 = [X0;0;0;0;ai';bi';tetai'];
Square_root_P00= (chol((lN_state*P00)));     
    
    
    
    


X_sigma_state(:,end)=X00;
for LL =1:N_state
X_sigma_state(:,LL)=X00+Square_root_P00(LL,:)'; 
X_sigma_state(:,LL+N_state) =X00-Square_root_P00(LL,:)';
end
for nn=1:(2*N_state+1)
X_propagate_sigma_matrix(:,nn,ii) = X_sigma_state(:,nn);
end
%a priori state
%propagation of sigma points (state transition)
 nn=1:(2*N_state+1);
   
   X_propagate_sigma_matrix(:,nn,ii) = X_sigma_state(:,nn);

   ai_matrix = X_sigma_state(7:7+size_params/3-1,nn)';
   bi_matrix = X_sigma_state(7+size_params/3:7+2*size_params/3-1,nn)';
   tetai_matrix = X_sigma_state(7+2*size_params/3:6+3*size_params/3,nn)';

   dtetai = rem(repmat(X_sigma_state(2,nn),size_params/3,1)'-tetai_matrix+pi,2*pi)-pi;
   X_propagate_sigma(1,nn) = X_sigma_state(1,nn)-(dt*sum(repmat(X_sigma_state(3,nn),size_params/3,1)'.*ai_matrix./(bi_matrix.^2).*dtetai.*exp(-dtetai.^2./(2*bi_matrix.^2)),2))'+X_sigma_state(4,nn); 
   X_propagate_sigma(2,nn) = rem(X_sigma_state(2,nn)+X_sigma_state(3,nn)*dt+X_sigma_state(5,nn)+pi,2*pi)-pi;

   
   X_propagate_sigma(3,nn)= X_sigma_state(3,nn)+X_sigma_state(6,nn);
   X_propagate_predict_sigma_matrix(1:3,nn,ii) = X_propagate_sigma (1:3,nn);

%a priori state mean (a weighted sum)
X_intermediate=0;
for nn=1:2*N_state
X_intermediate=X_intermediate+wm_state*(X_propagate_sigma(:,nn));
end;
% X_intermediate=X_intermediate/(2*lN_state);
X_intermediate=X_intermediate+(w0m_state*X_propagate_sigma(:,end));
% Xukf_predict(:,ii+1);
X0= X_intermediate(1:3);
%a priori cov.
P_xx=0; 
P_xx1=0;
for nn=1:2*N_state
P_xx=P_xx+wc_state*((X_propagate_sigma(1:3,nn)-X0)*(X_propagate_sigma(1:3,nn)-X0)');
end;
% P_xx=P_xx/(2*lN_state);
% P_xx1=((X_propagate_sigma(1:3,end)-X_intermediate(1:3))*(X_propagate_sigma(1:3,end)-X_intermediate(1:3))');
P_xx=P_xx+(w0c_state*((X_propagate_sigma(1:3,end)-X0)*(X_propagate_sigma(1:3,end)-X0)'));
M(:,:,ii)=P_xx(1:3,1:3);

P0 = M(:,:,ii);

% P_cross_covariance_xk_xk_plus_1(:,:,ii+1) = P_xx2; 
end



Xuks = Xukf_update;
Puks = P_update_matrix; 

for ii=size(Xukf_update,2)-1:-1:1
 
D = 0;
for nn=1:2*N_state    
D = D + wc_state*(X_propagate_predict_sigma_matrix(1:3,nn,ii+1)-Xukf_predict(1:3,ii+1))*(X_propagate_sigma_matrix(1:3,nn,ii)-Xukf_update(1:3,ii))';
end
% D = D/(2*lN_state);
D = D + w0c_state*(X_propagate_predict_sigma_matrix(1:3,end,ii+1)-Xukf_predict(1:3,ii+1))*(X_propagate_sigma_matrix(1:3,end,ii)-Xukf_update(1:3,ii))';

% D = 0;
% for nn=1:6    
% D = D + (X_propagate_predict_sigma_matrix(:,nn,ii)-Xukf_predict(:,ii+1))*(X_propagate_sigma_matrix(:,nn,ii)-Xukf_update(:,ii))';
% end
% D = D/(2*lN);
% D = D + LaN*(X_propagate_predict_sigma_matrix(:,7,ii)-Xukf_predict(:,ii+1))*(X_propagate_sigma_matrix(:,7,ii)-Xukf_update(:,ii))';


L = D*inv(P_predict_matrix(:,:,ii+1));
Xuks(:,ii) = Xukf_update(:,ii) + L*(Xuks(:,ii+1)-Xukf_predict(:,ii+1));
Puks(:,:,ii) = P_update_matrix(:,:,ii) + L*(Puks(:,:,ii+1)-P_predict_matrix(:,:,ii+1))*L';


disp(['Backward UKS- Processing sample number:' num2str(ii) '/' num2str(size(Xukf_update,2))])


end
% 



AUKF_SNR = 10*log10(mean((x(1,:)-y(1,:)).^2)/mean((x(1,:)-Xukf_update(1,:)).^2))
AUKS_SNR = 10*log10(mean((x(1,:)-y(1,:)).^2)/mean((x(1,:)-Xuks(1,:)).^2))


figure(2),
subplot(3,1,1)
        plot(1:length(x),x(1,:),'k')
                legend({'Original'})
                        axis('tight')
        title([file])
        subplot(3,1,3)
        plot(1:length(x),Xukf_update(1,:),'r')
                legend({'Adaptive UKF'})
                        axis('tight')

        subplot(3,1,2),plot(1:length(x),x_noisy,'b')
        legend({'Noisy'})
        axis('tight')


figure(3),
subplot(3,1,1)
        plot(1:length(x),x(1,:),'k')
                legend({'Original'})
                        axis('tight')
        title([file])
        subplot(3,1,3)
        plot(1:length(x),Xuks(1,:),'r')
                legend({'Adaptive UKS'})
                        axis('tight')

        subplot(3,1,2),plot(1:length(x),x_noisy,'b')
        legend({'Noisy'})
        axis('tight')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
function [Phase,Omega] = calculate_linear_phase_ver2(locs,length_sig,fs)

% locs       : indices of detected R-peaks
% length_sig : total number of ECG samples
% fs         : sampling frequency

ind = locs(:)';                % Convert R‑peak indices to row vector

Phase = zeros(1,length_sig);   % Phase of each ECG sample
Omega = zeros(1,length_sig);   % Instantaneous angular frequency

RR = mean(diff(ind));          % Mean RR interval (samples)

%% -------- Phase before the first R‑peak

stepTheta = 2*pi/RR;           % Average phase increment per sample
omega_val = fs*stepTheta;      % Instantaneous angular frequency

theta = 0;                     % Initialize phase

for j = ind(1)-1:-1:1          % Move backward from first R‑peak
    theta = theta - stepTheta; % Decrease phase
    theta = mod(theta+pi,2*pi)-pi; % Wrap phase into [-pi , pi]

    Phase(j) = theta;          % Store phase
    Omega(j) = omega_val;      % Store frequency
end

%% -------- Phase between consecutive R‑peaks

for k = 1:length(ind)-1

    bins = ind(k+1)-ind(k);    % Number of samples between R-peaks

    stepTheta = 2*pi/bins;     % Phase increment so phase spans one cycle
    omega_val = fs*stepTheta;  % Corresponding angular frequency

    theta = 0;
    Phase(ind(k)) = 0;         % Define phase at R‑peak as zero

    for j = ind(k)+1 : ind(k+1)-1
        theta = theta + stepTheta; % Linear phase progression
        if theta>pi
            theta = -pi;
        end
        Phase(j) = theta;
        Omega(j) = omega_val;
    end

    Phase(ind(k+1)) = 0;       % Next R‑peak also set to zero phase
end

%% -------- Phase after the last R‑peak

stepTheta = 2*pi/RR;           % Use mean RR again
omega_val = fs*stepTheta;

theta = 0;

for j = ind(end)+1:length_sig
    theta = theta + stepTheta; % Continue phase linearly
    theta = mod(theta+pi,2*pi)-pi;

    Phase(j) = theta;
    Omega(j) = omega_val;
end

end



function [ecgsd,ecg_mean,phase_mean] = ecgsd_extractor_ver1(ecg,phase,bins)

x1 = ecg;                        % ECG signal
meanPhase = zeros(1,bins);       % Mean phase per bin
ECGmean = zeros(1,bins);         % Mean ECG per bin
ECGsd = zeros(1,bins);           % ECG standard deviation per bin

% Handle wrap-around phase bin near -pi / +pi
I = find( phase >= (pi-pi/bins) | phase < (-pi+pi/bins) );

if(~isempty(I))
    meanPhase(1) = -pi;
    ECGmean(1) = mean(x1(I));
    ECGsd(1) = std(x1(I));
else
    ECGsd(1) = -1;               % Mark empty bins
end

% Loop over phase bins
for i = 1 : bins-1
    I = find( phase >= 2*pi*(i-0.5)/bins - pi & ...
              phase <  2*pi*(i+0.5)/bins - pi );

    if(~isempty(I))
        meanPhase(i+1) = mean(phase(I));
        ECGmean(i+1) = mean(x1(I));
        ECGsd(i+1) = std(x1(I));
    else
        ECGsd(i+1) = -1;
    end
end

% Interpolate missing bins
K = find(ECGsd==-1);

for i = 1:length(K)
    switch K(i)
        case 1
            meanPhase(1) = -pi;
            ECGmean(1) = ECGmean(2);
            ECGsd(1) = ECGsd(2);
        case bins
            meanPhase(bins) = pi;
            ECGmean(bins) = ECGmean(bins-1);
            ECGsd(bins) = ECGsd(bins-1);
        otherwise
            meanPhase(K(i)) = mean(meanPhase([K(i)-1 K(i)+1]));
            ECGmean(K(i))   = mean(ECGmean([K(i)-1 K(i)+1]));
            ECGsd(K(i))     = mean(ECGsd([K(i)-1 K(i)+1]));
    end
end

phase_mean = meanPhase;
ecg_mean   = ECGmean;
ecgsd      = ECGsd;

end





