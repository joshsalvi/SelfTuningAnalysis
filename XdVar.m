function XdVar(timeseries,window)
%Running variance calculation
endnum = length(timeseries)-window;

for k = 1:endnum
    MovVar(k) = var(timeseries(k:k+window-1));
end

figure
plot(MovVar)

end

