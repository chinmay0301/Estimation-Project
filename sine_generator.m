function [ x ] = sine_generator(size, freq)
t = linspace(0,1,size);
x = sin(2*pi*freq*t)';
end

