function amp = amplitude_generator(k0,k1,type,ph)
if type=='HSim'
    amp=k0+k1*abs(ph)/max(abs(ph(:)));
elseif type=='LSim'
    amp=k0+k1*abs(cos(15*ph));
end
end

