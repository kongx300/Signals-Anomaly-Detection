% function anamoly_dectector()
% 
global params;
params.NUM_PTS = 1024;
params.NUM_FRAMES = 100;
params.MAX_RANK = 20;
params.ZEROTHRESH = 1.0e-16;

NUM_PTS = params.NUM_PTS;
NUM_FRAMES = params.NUM_FRAMES;
MAX_RANK = params.MAX_RANK;
ZEROTHRESH = params.ZEROTHRESH;

fprintf("\n\n\n")

% % Create a sample singal mutually
% Set the total length of the artificial signal 
Nx = NUM_PTS + NUM_FRAMES - 1;
% Set the boundary and build the time series vector
Tmin = 0;
Tmax = 13;
t = linspace(Tmin,Tmax*(Nx-1)/Nx, Nx); 
t = t'; % t is col vector
% Set step
dt = t(2) - t(1);
fsamp = 1/dt;
% Set the frequency of the sine wave and build the wave
f0 = 2.0;
v_sample = sin(2*pi*f0*t); 
% Plot the signal
figure(1)
plot(t,v_sample);
hold on;
plot(t,v_sample,'ro');

% % Training Phase 
% Initialize the matrix Yf for svd decompo
Yf = [];
% Loop over the artificial signal
for i =1:NUM_FRAMES
    % Grab the i-th 1024 samples of sound from A/D each time,
    v = v_sample(i:NUM_PTS+i-1);
    % Compute mean
    v_mean = mean(v);
    % Substruct mean
    y = v - v_mean;
    
    % Do the FFT
    yspec = fft(y);
    % Compute the power spectral density
    ypsd = abs(yspec);
    ypsd = ypsd.^2;
    % Normalize and threshold spectrum
    tmp1 = sum(ypsd);
    pow_learn_thresh = tmp1/1000.0;
   
    ypsd = ypsd/tmp1;
    idx = find(ypsd < ZEROTHRESH);
    for j=1:length(idx)
        ypsd(idx(j)) = 0.0;
    end
    % Create the matrix Yf
    Yf = horzcat(Yf,ypsd);
end


% Do svd
[U,S,V] = svd(Yf);
% Get effective rank of M0, then pull out vectors
pr = sum(sum(S));
X = S/pr;
H = 0;
for j=1:NUM_FRAMES
    if(X(j,j)>0)
        entry = X(j,j);
        H = H + entry*log(entry);
    end
end
r = floor(exp(-H)); 
% fprintf("r = %d\n",r);

% If the dimensionality beyond the maximum, set dimen = r
if(r > MAX_RANK)
    r = MAX_RANK;
end
% Extract
Ur = U(:,r);
% Create projection
Nr = eye(size(Ur)) - Ur*((Ur'*Ur)\Ur');

% Compute the norm
yn = norm(Nr*ypsd)/norm(ypsd);
% Use K = 2.8 for example
K = 2.8;
% Threshold
Tr = K*yn;
% fprintf("Tr = %f\n",Tr)

fprintf("----------------------------\n")
fprintf("Now do the operating phase\n")
% Take a segment of the sine wave in the training phase
operating_phase(v,tmp1,Nr,Tr)

fprintf("----------------------------\n")
fprintf("Now try a different signal\n")
% Set the frequency of a signal
f1 = 3;
% Create a new signal wave
% different from the one in the training phase
v_sample2 = cos(2*pi*f1*t); 
vd = v_sample2(10: NUM_PTS+10-1);
% Plot the new wave
figure(3);
plot(t(10: NUM_PTS+10-1),vd)
hold on;
plot(t(10: NUM_PTS+10-1),vd,"ro")
% Do the operating phase
operating_phase(vd,tmp1,Nr,Tr)
   



    
    

    


    
















% end