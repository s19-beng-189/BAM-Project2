% Solve delay-differential equations for absence seizures
% For the paper 'Generalized Seizures in a Neural Field Model with
% Bursting Dynamics'
% NB: for different states the population connection strengths must be changed

% by Xuelong Zhao

function absence_seizure()

clear all; clc

% Parameters from Table 1
gamma = 100; % s-1
t0 = 80e-3; % s
Qmax = 250; % s-1
theta = 15.0e-3; % V
sigma = 3.3e-3; % V
Ic = 0.176; % A m-2
taux = 0.15; % s
tauh = 0.60; % s
a_x = 0.003; % s
mu = 16.0; % S m-2
Qmax_r = 80; % s-1
sigma_r = 3.3e-3; % V
Vk = -95e-3; % V
Vx = 140e-3; % V
Veff = -66e-3; % V
phin_0 = 1.0; % s-1

% State-dependent parameters for seizures (Table 3)
nu_ee = 0.0010; % V s
nu_ei = -0.0018;
nu_es = 0.0032;
nu_ie = 0.0010;
nu_ii = -0.0018;
nu_is = 0.0032;
nu_re = 0.0016;
nu_rs = 0.0006;
nu_se = 0.0017; 
nu_sr = -0.0008;
nu_sn = 0.002;
alpha = 50; % s-1
beta = 200; % s-1
gX = 3.0; % S m-2
gH = gX*2.3678; % S m-2
Xinit = 0.0595;
Hinit = 0.1788;
Ib = -gX*(Veff-Vx)/3;
Ia = Ib-gH*(Veff-Vk);
phin = 5e-5;

% Time (s)
start = 0.0;
stop = 100.0 + t0/2; 
dt = 1e-4;
k0 = t0/2/dt; % the delay in timesteps
time_array = [start:dt:stop];
vec_len = length(time_array);

% Ramp nuse variable
t1 = 20; t2 = 40; t3 = 60; t4 = 80; 
nu_se1 = 0.0017; nu_se2 = 0.0033; 
upslope = [nu_se1:(nu_se2-nu_se1)/(t2/dt-t1/dt):nu_se2];
downslope = [nu_se2:(nu_se1-nu_se2)/(t4/dt-t3/dt):nu_se1];
horizontal1 = nu_se1*ones(1,(t0/dt+t1/dt));
horizontal2 = nu_se2*ones(1,(t3/dt-t2/dt));
horizontal3 = nu_se1*ones(1,((stop-t0/2)/dt-t4/dt));
nu_se = [horizontal1 upslope(1:end-1) horizontal2 downslope(1:end-1) horizontal3]; % ramp

% Noise
rng(0,'twister'); 
noise = alpha*beta*nu_sn*sqrt(phin*dt)*randn(1,vec_len); 

% Outputs
phie = zeros(1,vec_len);
Ve = zeros(1,vec_len);
Vr = zeros(1,vec_len);
Vs = zeros(1,vec_len);
phiedot = zeros(1,vec_len);
Vedot = zeros(1,vec_len);
Vrdot = zeros(1,vec_len);
Vsdot = zeros(1,vec_len);
X = zeros(1,vec_len);
H = zeros(1,vec_len);
modtheta = zeros(1,vec_len);

% Initialize outputs
phie(1:k0+1) = 3.175;
Ve(1:k0+1) = 0.0006344; 
Vr(1:k0+1) = 0.005676;
Vs(1:k0+1) = -0.003234;
X(1:k0+1) = Xinit;
H(1:k0+1) = Hinit;
modtheta(1:k0+1) = Ic/mu;

% Delay-differential equations
% Solved using Euler-Maruyama scheme
for i = (k0+2):(stop/dt+1) 
    Ve(i) = Ve(i-1) + Vedot(i-1)*dt;
    Vedot(i) = Vedot(i-1) + dt*( alpha*beta * ( nu_ee*phie(i-1) + nu_ei*sig(Ve(i-1),Qmax,theta,sigma) + nu_es*sig(Vs(i-1-k0),Qmax,theta,sigma) - (1/alpha + 1/beta)*Vedot(i-1) - Ve(i-1) ) );
       
    Vr(i) = Vr(i-1) + Vrdot(i-1)*dt;
    Vrdot(i) = Vrdot(i-1) + dt*( alpha*beta * ( nu_re*phie(i-1-k0) + nu_rs*sig(Vs(i-1),Qmax,theta,sigma) - (1/alpha + 1/beta)*Vrdot(i-1) - Vr(i-1) ));
    
    Vs(i) = Vs(i-1) + Vsdot(i-1)*dt;
    Vsdot(i) = Vsdot(i-1) + dt*( alpha*beta * ( nu_se(i-1)*phie(i-1-k0) + nu_sr*sigr(Vr(i-1),Qmax_r,modtheta(i-1),sigma_r) - (1/alpha + 1/beta)*Vsdot(i-1) - Vs(i-1) + nu_sn*phin_0 )) + noise(i-1);
    
    phie(i) = phie(i-1) + phiedot(i-1)*dt;
    phiedot(i) = phiedot(i-1) + dt*( gamma^2 * ( sig(Ve(i-1),Qmax,theta,sigma) - 2/gamma*phiedot(i-1) - phie(i-1) ));
    
    modtheta(i) = ( Ic - 3*Ib*X(i-1) + (Ib-Ia)*H(i-1) )/mu;
    X(i) = X(i-1) + dt*( -1/taux * ( X(i-1) - a_x*sigr(Vr(i-1),Qmax_r,modtheta(i-1),sigma_r) )); % X_inf = a_x*sigr(Vr(i)); 
    H(i) = H(i-1) + dt*( -1/tauh * ( H(i-1) - 3*X(i-1) ));
end

% Remove initial t0-length from outputs
phie = phie((k0+1):end);
Ve = Ve(k0+1:end);
Vr = Vr(k0+1:end);
Vs = Vs(k0+1:end);
X = X(k0+1:end);
H = H(k0+1:end);
time_array = time_array(k0+1:end);

% Time series
figure;
plot(time_array,phie,'k'), title('phie')

% Power pectra
figure; X=downsample(phie(1:end-1),100);
[S,F,T,P] = spectrogram(X,256, 250, 256, 1E2); 
surf(T,F,10*log10(P),'edgecolor','none'); axis tight; 
view(0,90); xlabel('time (s)'); ylabel('Frequency (Hz)');

% Functions
function firing_rate = sig(v,Qmax,theta,sigma)
	% sigmoid for voltage-rate relationship
	firing_rate = Qmax / (1 + exp(-(v-theta) / sigma));
 
function firing_rate = sigr(v,Qmax_r,modtheta,sigma_r)
     % sigmoid for voltage-rate relationship in reticular nucleus
     firing_rate = Qmax_r / (1 + exp(-(v-modtheta) / sigma_r));