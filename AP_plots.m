x = [20:3:50]
peak = [42	42.1001	41.6433	40.9273	40.0049	38.5542	36.3282	33.1411	30.1315	25.0022	18.919];
width = [15.9 12.1 9 6.4 4.3 3.4 2.6 2 1.6 1.2 0.9];

figure(1)
plot(x,peak)
title('Temperature vs AP Peak')
xlabel('Temperature (°C)')
ylabel('AP Peak (mV)')

figure(2)
plot(x,width)
title('Temperature vs AP Width')
xlabel('Temperature (°C)')
ylabel('AP Width (ms)')
