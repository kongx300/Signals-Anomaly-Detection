function operating_phase(v,tmp1,Nr,Tr)

    global params;
    NUM_PTS = params.NUM_PTS;
    NUM_FRAMES = params.NUM_FRAMES;
    MAX_RANK = params.MAX_RANK;
    ZEROTHRESH = params.ZEROTHRESH;
    
    x = v - mean(v);
    % Compute fft
    xspec = fft(x);
    xpsd = abs(xspec).^2;
    % Normalize and threshold spectrum
    xpsd = xpsd/tmp1;
    idx = find(xpsd < ZEROTHRESH);
    for j=1:length(idx)
        xpsd(idx(j)) = 0.0;
    end

    % Project onto nullspace
    xproj = Nr*xpsd;
    % Get the norm of xproj
    xprojNorm = norm(xproj);
    % Get the norm of xpsd
    xpsdNorm = norm(xpsd);
    fprintf("The ratio of two norms is %f, ", xprojNorm/xpsdNorm)
    fprintf("and Tr = %f\n", Tr)
    if(xprojNorm/xpsdNorm > Tr)
        fprintf("Anomaly detected\n");
    else
        fprintf("Seems like it is normal\n");
    end    
end