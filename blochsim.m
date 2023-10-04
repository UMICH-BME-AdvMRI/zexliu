function [M_mtx] = blochsim(rf_abs, rf_angle, g, dT, df, zlocs, b1_range)
% rf (mT), g (mT/mm), dT: time step (ms), df (Hz), zlocs (m)

nt = length(rf_abs);
nf = length(df);
nb1 = length(b1_range);

M_init = zeros(3,nf,nb1);
M_init(3,:,:) = 1;

M_mtx = M_init;

for b1 = 1:nb1
    for f = 1:nf
        M_fb1 = M_init(:,f,b1);
        for t = 1:nt
            Rgrad = zrot_new((42.58*zlocs*g(t)+df(f)/1000)*dT*360);
            alpha = b1_range(b1)*rf_abs(t)*dT*42.58*360;
            phi = rf_angle(t)/pi*180;
            Rrf = throt_new(alpha, phi);
            M_fb1 = Rrf*Rgrad*M_fb1;
        end
        M_mtx(:,f,b1) = M_fb1;
    end
end
end