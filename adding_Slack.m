clc
clear all
close all

%%

t = -2:0.1:2;
thresh = 2;
n = 1/(log(thresh^2-1));
n = 0.8
thresh = n*thresh;

o1 = t.^2;    % normal square

for i=1:length(t)
    if abs(t(i)) > thresh
        o2(i) = (abs(t(i))-thresh).^2;
    else
        o2(i) = 0;
    end
end

o3 = exp(abs(t)/thresh) -1;

plot(t,o1,'-b');
hold on
plot(t,o3.^2,'-k');
axis([min(t) max(t) -1 max(o1)*2])
grid on
