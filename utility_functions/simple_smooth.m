function ysmooth = simple_smooth(y)
%simple smoothing based on three element moving average

ysmooth = y;
ysmooth(2:end-1) = (y(1:end-2)+y(2:end-1)+y(3:end))/3;

end

