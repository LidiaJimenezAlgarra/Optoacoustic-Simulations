function S = timeFiltS(S, fLim, dt)

filtOrder = 4;
Ns        = size(S,1);
Nt        = size(S,2);

f = linspace(-1/2/dt,1/2/dt,Nt);
Sf = fftshift(fft(ifftshift(S,2),[],2),2);
filt = (exp(-(f/fLim(2)).^filtOrder)).*(1-(exp(-(f/fLim(1)).^filtOrder)));
filt = repmat(filt,Ns,1);
S = real(fftshift(ifft(ifftshift(Sf.*filt,2),[],2),2));
end
