function width = fwhm_dakis(x)
% return the Full Width at Half Maximum of the pulse x
nsize = size(x,1);
if nsize <=1
    nsize = size(x,2);
end 

% half_peak = max(x)/2;
[peak ind] = max(x);
half_peak = peak/2; % 3 dB
for iii=1:nsize-1
    if (x(iii)<=half_peak) && (x(iii+1)>half_peak)
        break;
    end
end

for jjj=ind:nsize-1
    if (x(jjj)>=half_peak) && (x(jjj+1)<half_peak)
        break;
    end
end
x1a=iii;
x2a=iii+1;
y1a=x(x1a);
y2a=x(x2a);
f1=x1a + (x2a-x1a)*(half_peak-y1a)/(y2a-y1a);

x1b=jjj;
x2b=jjj+1;
y1b=x(x1b);
y2b=x(x2b);
f2=x1b + (x2b-x1b)*(half_peak-y1b)/(y2b-y1b);

width = f2-f1;
end
