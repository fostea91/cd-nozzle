% transonic_nozzle.m
% 
% Copyright 2018 Andrew Foster <fost22@protonmail.com>
%
% Flow analysis in the transonic region of a convergent-divergent nozzle
% with a circular, parabolic, or rectangular hyperbolic convergent contour,
% based to the work of Dutton and Addy, 1981
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.

%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <https://www.gnu.org/licenses/>.

clear;
clc

gm = 1.4;
rcsub = 0.625;
M_spec = 1;
contour_type = 1; % 1 for circular, 0 for parabolic, -1 for hyperbolic
nu = 1;

expan = 1/(rcsub + nu);

x = 0.3;
x_delta = 0.2;

y = 0;
y_delta = 0.01;

a = 0;
i = 1;

while a~=1
    
    h = 0;
    
    x_low = x(i) - x_delta;
    x_hi = x(i) + x_delta;
    
    M_low = trans_nozzle(nu, x_low, y(i), gm, expan);
    M_hi = trans_nozzle(nu, x_hi, y(i), gm, expan);
            
    x(i) = (x_hi - x_low)*((M_spec - M_low)/(M_hi - M_low)) + x_low;
    
    while h ~= 1
        
        M(i) = trans_nozzle(nu, x(i), y(i), gm, expan);
        y_w(i) = 1 + 0.5*expan*x(i)^2 + (0.5*nu*x(i)^2)*expan^2 + ...
            (0.5*(nu^2)*(x(i)^2)+0.125*contour_type*x(i)^4)*expan^3;
        
        if M(i) < (1+1e-6)*M_spec && M(i) > (1-1e-6)*M_spec
            if y(i)>=y_w
                h = 1;
                a = 1;
            else
                x(i+1) = x(i);
                x_delta = 0.075;
                y(i+1) = y(i) + y_delta;
                h = 1;
            end
            
        else
            x_delta = 0.7*x_delta;
            x_low = x(i) - x_delta;
            x_hi = x(i) + x_delta;
            
            M_low = trans_nozzle(nu, x_low, y(i), gm, expan);
            M_hi = trans_nozzle(nu, x_hi, y(i), gm, expan);
            
            x(i) = (x_hi - x_low)*((M_spec - M_low)/(M_hi - M_low)) + x_low;
            
        end
    end
    
    i = i+1;
    
end

cd2 = -(8*gm+21-48*nu)/2304;
cd3 = (754*gm^2+(1971-2880*nu)*gm+2007-7560*nu+8640*nu^2)/276480;
cd = 1 - (gm+1)*(expan^2)*(1/96 + cd2*expan + cd3*expan^2);

hold on
plot(x,y)
plot(x,y_w,'r')

ax_range = max(y_w);

axis([-0.5 0.5 0 1]*ax_range)
