fs = 1000;
t = 0:0.001:1-0.001;
x = cos(2*pi*100*t)+randn(size(t));
periodogram(x,[],length(x),fs,'centered')