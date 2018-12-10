% transonic_nozzle.m
% Author: Andrew Foster
%
% Flow analysis in the transonic region of a convergent-divergent nozzle
% with a circular convergent contour, based to the work of Kliegel and
% Levine, 1969

clear;
clc

gm = 1.4;
rcsub = 2;
M_spec = 1.005;

nu_w = 0.5*log((1+((1+2*rcsub)^0.5)/(1+rcsub))/(1-((1+2*rcsub)^0.5)/(1+rcsub)));
nu = 0:0.01*nu_w:nu_w;

xi = zeros(1,length(nu));
xi(1) = 0.2;
xi_delta = pi/4;

for i=1:length(nu)
    
    h = 0;
    
    xi_low = xi(i) - xi_delta;
    xi_hi = xi(i) + xi_delta;
    
    [M_low, ~, ~] = trans_nozzle(nu(i), xi_low, gm, rcsub);
    [M_hi, ~, ~] = trans_nozzle(nu(i), xi_hi, gm, rcsub);
            
    xi(i) = (xi_hi - xi_low)*((M_spec - M_low)/(M_hi - M_low)) + xi_low;
    
    while h ~= 1
        
        [M, z(i), r(i)] = trans_nozzle(nu(i), xi(i), gm, rcsub);
        
        if M < (1+1e-6)*M_spec && M > (1-1e-6)*M_spec
            if i<length(nu)
                xi(i+1) = xi(i);
                h = 1;
            else
                h = 1;
            end
            
        else
            xi_delta = 0.7*xi_delta;
            xi_low = xi(i) - xi_delta;
            xi_hi = xi(i) + xi_delta;
            
            [M_low, ~, ~] = trans_nozzle(nu(i), xi_low, gm, rcsub);
            [M_hi, ~, ~] = trans_nozzle(nu(i), xi_hi, gm, rcsub);
            
            xi(i) = (xi_hi - xi_low)*((M_spec - M_low)/(M_hi - M_low)) + xi_low;
            
        end
    end
end

xi_w = min(xi):(max(xi)-min(xi))*0.01:max(xi);

for i = 1:length(xi_w)
    
    r_w(i) = ((1+2*rcsub)^0.5)*sinh(nu_w)/(cosh(nu_w)+cos(xi_w(i)));
    z_w(i) = ((1+2*rcsub)^0.5)*sin(xi_w(i))/(cosh(nu_w)+cos(xi_w(i)));
    
end

hold on
plot(z,r)
plot(z_w,r_w,'r')

ax_range = max(r_w);

axis([-0.5 0.5 0 1]*ax_range)
