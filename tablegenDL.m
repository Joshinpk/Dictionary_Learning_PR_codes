function [Tphase Tamp]= tablegenDL(RMSE_phi,RMSE_ampl,ALGparam)

if ALGparam.Pois_Gaus_sel==1
    noise_level_sets=ALGparam.KAPPA_set;
else
    noise_level_sets=ALGparam.SNR_Gaus_set;
end
noise_level_len=length(noise_level_sets);
SIG_set=ALGparam.SIG_set;
num_sig=length(SIG_set);
for i=1:noise_level_len
    Tphase_temp=table(RMSE_phi(:,i));
    Tamp_temp=table(RMSE_ampl(:,i));
    kappa_str=num2str((i));
    cname=strcat('Kappa',kappa_str);
    Tphase_temp.Properties.VariableNames{1}=cname;
    Tamp_temp.Properties.VariableNames{1}=cname;

    if i==1
        Tphase=Tphase_temp;
        Tamp=Tamp_temp;
    else
        Tphase=[Tphase Tphase_temp];
        Tamp  =[Tamp Tamp_temp];
    end
end

for i=1:num_sig
    Tphase.Properties.RowNames{i}= strcat('signal:',num2str(SIG_set(i)));
    Tamp.Properties.RowNames{i}= strcat('signal:',num2str(SIG_set(i)));
end


end
