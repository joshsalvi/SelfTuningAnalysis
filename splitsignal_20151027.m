load('/Users/joshsalvi/Documents/Lab/Lab/Clamp Data/2015-10-22.01/Ear 1/Cell 11/Sens-ind1  2  3-data.mat')

Fs = 5e3;   % input sample rate (Hz)
ksf = 200;  % input fiber stiffness (µN/m)
dt = 1/Fs;
tvec = 0:dt:length(Xd_pulse{1,ind(1)})*dt-dt;
NFFT = (2^4)*2^nextpow2(numel(tvec));   % CHOOSE number of FFT points (2^n)*(n_p)

% Split the data
numsplit = 5;
sizeX = size(Xd_split);
for j = 1:sizeX(1)
    for k = 1:sizeX(2)
        for l = 1:sizeX(3)
            for m = 1:numsplit
                Xo_split2{j,k,l,m} = Xo_split{j,k,l}(1+length(Xo_split{j,k,l})/numsplit*(m-1):length(Xo_split{j,k,l})/numsplit*(m));
                Xd_split2{j,k,l,m} = Xd_split{j,k,l}(1+length(Xd_split{j,k,l})/numsplit*(m-1):length(Xd_split{j,k,l})/numsplit*(m));
            end
        end
    end
end

stimlength=length(Xd_split2{1,1,1});

% Calculate FFT
q=0;
for j = 1:sizeX(1)
    for k = 1:sizeX(2)
        for l = 1:sizeX(3)
            for m = 1:numsplit
                Xd_fft2{j,k,l,m} = fft(Xd_split2{j,k,l,m},NFFT)./stimlength; Xd_fft2{j,k,l,m} = Xd_fft2{j,k,l,m}(1:NFFT/2+1);
                Xo_fft2{j,k,l,m} = fft(Xo_split2{j,k,l,m},NFFT)./stimlength; Xo_fft2{j,k,l,m} = Xo_fft2{j,k,l,m}(1:NFFT/2+1);
                freq0 = findnearest(f,freqrange(l)); freq0 = freq0(1);
                Xopk2(j,k,l,m) = Xo_fft2{j,k,l,m}(freq0);
                Xdpk2(j,k,l,m) = Xd_fft2{j,k,l,m}(freq0);
                q = q+1;
            end            
            disp([num2str(round(q/(sizeX(1)*sizeX(2)*sizeX(3)*numsplit)*100)) '% complete']);
        end
    end
end

%% Calculate means and SEMs
close all
analyzedind = 1:5;  % choose indices to analyze


for j = 1:sizeX(1)
    for k = 1:sizeX(2)
        for l = 1:sizeX(3)
            Xdpkmean2(j,k,l) = abs(mean(Xdpk2(j,k,l,analyzedind))) * 1e-9;
            Fepkmean2(j,k,l) = abs(mean(Xopk2(j,k,l,analyzedind))) * ksf * 1e-15;
            sensmean2(j,k,l) = Xdpkmean2(j,k,l) / Fepkmean2(j,k,l);
            Xdpksem2(j,k,l) = abs(std(Xdpk2(j,k,l,analyzedind)))*1e-9/sqrt(length(analyzedind));
            Fepksem2(j,k,l) = abs(std(Xopk2(j,k,l,analyzedind)))*ksf*1e-15/sqrt(length(analyzedind));
            senssem2(j,k,l) = sensmean2(j,k,l) * sqrt( (Xdpksem2(j,k,l)/Xdpkmean2(j,k,l))^2 + (Fepksem2(j,k,l)/Fepkmean2(j,k,l))^2 );
        end
    end
end
%{
for j = 1:sizeX(1)
    for k = 1:sizeX(2)
        figure(k)
        errorbar(freqrange,Xdpkmean2(j,k,:),Xdpksem2(j,k,:)); hold all;
    end
end
xlabel('Frequency (Hz)'); ylabel('Response amplitude (m)');
set(3,'WindowStyle','docked')

for j = 1:sizeX(1)
    for k = 1:sizeX(2)
        figure(k+sizeX(2))
        errorbar(freqrange,sensmean2(j,k,:),senssem2(j,k,:)); hold all;
    end
end
xlabel('Frequency (Hz)'); ylabel('Sensitivity (m/N)');
set(2,'WindowStyle','docked')
%}

for j = 1:sizeX(1)
    for l = 1:sizeX(3)
        meansensmean2(j,l) = mean(sensmean2(j,:,l));
        semsensmean2(j,l) = std(sensmean2(j,:,l))/sqrt(sizeX(1));
        meanXdpkmean2(j,l) = mean(Xdpkmean2(j,:,l));
        semXdpkmean2(j,l) = std(Xdpkmean2(j,:,l))/sqrt(sizeX(1));
    end
end

figure
setfiguredefaults(sizeX(1));
for j = 1:sizeX(1)
    errorbar(freqrange,squeeze(meansensmean2(j,:)).*1e-3,squeeze(semsensmean2(j,:)).*1e-3); hold all;
end
xlabel('Frequency (Hz)'); ylabel('Sensitivity (km/N)');
