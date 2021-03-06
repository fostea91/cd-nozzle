% sonic_line.m
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

%   This function provides the line at which sonic conditions are achieved
%   in the throat of a nozzle.
%   gm = specific heat ratio
%   rcsub = the wall radius of curvature in the subsonic region
%   M_spec = the Mach value for the calculated line
%   contour_type =  1 for circular arc
%                   0 for parabolic
%                   -1 for hyperbolic
%   eta = should be set to 1 for most solutions

function sonic_line_calc = sonic_line(gm, rcsub, M_spec, contour_type, eta)

expan = 1/(rcsub + eta);

x = 0.3;
x_delta = 0.2;

y = 0;
y_delta = 0.01;

a = 0;
i = 1;

while a~=1
    
    h = 0;
    
    while h ~= 1
        
        M(i) = trans_nozzle_mach(eta, x(i), y(i), gm, expan);
        y_w(i) = 1 + 0.5*expan*x(i)^2 + (0.5*eta*x(i)^2)*expan^2 + ...
            (0.5*(eta^2)*(x(i)^2)+0.125*contour_type*x(i)^4)*expan^3;
        
        if M(i) >= (1+1e-6)*M_spec
            x_delta = 0.7*x_delta;
            x2 = x(i) - x_delta;
            
            M_low = trans_nozzle_mach(eta, x2, y(i), gm, expan);
            
            x(i) = (x(i) - x2)*((M_spec - M_low)/(M(i) - M_low)) + x2;
        elseif M(i) <= (1-1e-6)*M_spec
            x_delta = 0.7*x_delta;
            x2 = x(i) + x_delta;
            
            M_low = trans_nozzle_mach(eta, x2, y(i), gm, expan);
            
            x(i) = (x(i) - x2)*((M_spec - M_low)/(M(i) - M_low)) + x2;
        else
            if y(i)>=y_w
                h = 1;
                a = 1;
            else
                x(i+1) = x(i);
                x_delta = 0.075;
                y(i+1) = y(i) + y_delta;
                h = 1;
            end
        end
    end
    
    i = i+1;
    
end

x_delta = 0.2;
    
xw1 = (y(length(y)) - y_w(length(y_w)))/...
    ((y_w(length(y_w)) - y_w(length(y_w)-1))/(x(length(x)) - x(length(x)-1)) - ...
    (y(length(y)) - y(length(y)-1))/(x(length(x)) - x(length(x)-1)));
    
yw1 = 1 + 0.5*expan*xw1^2 + (0.5*eta*xw1^2)*expan^2 + ...
    (0.5*(eta^2)*(xw1^2)+0.125*contour_type*xw1^4)*expan^3;
    
h = 0;
    
while h ~= 1
        
    Mw1 = trans_nozzle_mach(eta, xw1, yw1, gm, expan);
    
    if Mw1 > (1+1e-6)*M_spec
        x_delta = 0.7*x_delta;
        xw2 = xw1 - x_delta;
        
        M_low = trans_nozzle_mach(eta, xw2, yw1, gm, expan);
        
        xw1 = (xw1 - xw2)*((M_spec - M_low)/(Mw1 - M_low)) + xw2;
    elseif Mw1 < (1-1e-6)*M_spec
        x_delta = 0.7*x_delta;
        xw2 = xw1 + x_delta;
        
        M_low = trans_nozzle_mach(eta, xw2, yw1, gm, expan);
            
        xw1 = (xw1 - xw2)*((M_spec - M_low)/(Mw1 - M_low)) + xw2;
    else
        yw2 = 1 + 0.5*expan*xw1^2 + (0.5*eta*xw1^2)*expan^2 + ...
            (0.5*(eta^2)*(xw1^2)+0.125*contour_type*xw1^4)*expan^3;
            
        if yw2>yw1*(1+1e-8)
            x_delta = 0.7*x_delta;
            xw1 = xw1 - x_delta;
            yw1 = 1 + 0.5*expan*xw1^2 + (0.5*eta*xw1^2)*expan^2 + ...
                (0.5*(eta^2)*(xw1^2)+0.125*contour_type*xw1^4)*expan^3;
        elseif yw2<yw1*(1-1e-8)
            x_delta = 0.7*x_delta;
            xw1 = xw1 + x_delta;
            yw1 = 1 + 0.5*expan*xw1^2 + (0.5*eta*xw1^2)*expan^2 + ...
                (0.5*(eta^2)*(xw1^2)+0.125*contour_type*xw1^4)*expan^3;
        else
            h = 1;
        end
    end
end

sonic_line_calc.x = [x(1:(i-3)), xw1];
sonic_line_calc.y = [y(1:(i-3)), yw1];
sonic_line_calc.y_w = [y_w(1:(i-3)), yw2];
sonic_line_calc.M = [M(1:(i-3)), Mw1];
[sonic_line_calc.p, sonic_line_calc.S] = polyfit(sonic_line_calc.y, sonic_line_calc.x, 2);

end


function M = trans_nozzle_mach(eta, x, y, gm, expan)

z = (((gm+1)*expan*0.5)^(-0.5))*x;

u1 = 0.5*y^2 - 0.25 + z;
v1 = 0.25*y^3 - 0.25*y + y*z;
    
u2 = ((2*gm+9)/24)*y^4 - ((4*gm+15-12*eta)/24)*y^2 + ...
     ((10*gm+57-72*eta)/288) + z*(y^2 + (4*eta-5)/8) - ((2*gm-3)/6)*z^2;
    
v2 = ((gm+3)/9)*y^5 - ((20*gm+63-36*eta)/96)*y^3 + ...
     ((28*gm+93-108*eta)/288)*y + z*(((2*gm+9)/6)*y^3 - ...
     ((4*gm+15-12*eta)/12)*y) + y*z^2;
    
u3 = ((556*gm^2+1737*gm+3069)/10368)*y^6 - ...
     ((388*gm^2+(1161-384*eta)*gm+(1881-1728*eta))/2304)*y^4 + ...
     ((304*gm^2+(831-576*eta)*gm+(1242-2160*eta+864*eta^2))/1728)*y^2 - ...
     ((2708*gm^2+(7839-5760*eta)*gm+(14211-32832*eta+20736*eta^2))/82944) + ...
     (((52*gm^2+51*gm+327)/384)*y^4 - ...
     ((52*gm^2+75*gm+279-288*eta)/192)*y^2 + ...
     ((92*gm^2+180*gm+639-1080*eta+432*eta^2)/1152))*z + ...
     ((-(7*gm-3)/8)*y^2 + (((13-16*eta)*gm-(27-24*eta))/48))*z^2 + ...
     ((4*gm^2-57*gm+27)/144)*z^3;
    
v3 = ((6836*gm^2+23031*gm+30627)/82944)*y^7 - ...
     ((3380*gm^2+(11391-3840*eta)*gm+15291-11520*eta)/13824)*y^5 + ...
     ((3424*gm^2+(11271-7200*eta)*gm+15228-22680*eta+6480*eta^2)/13824)*y^3 - ...
     ((7100*gm^2+(22311-20160*eta)*gm+30249-66960*eta+38880*eta^2)/82944)*y + ...
     (((556*gm^2+1737*gm+3069)/1728)*y^5 - ...
     ((388*gm^2+(1161-384*eta)*gm+1881-1728*eta)/576)*y^3 + ...
     ((304*gm^2+(831-576*eta)*gm+1242-2160*eta+864*eta^2)/864)*y)*z + ...
     (((52*gm^2+51*gm+327)/192)*y^3 - ...
     ((52*gm^2+75*gm+279-288*eta)/192)*y)*z^2 - ...
     ((7*gm-3)/12)*y*z^3;
    
u = 1 + u1*expan + u2*expan^2 + u3*expan^3;
v = sqrt((gm+1)*expan*0.5)*(v1*expan + v2*expan^2 + v3*expan^3);

M = 1 + 0.5*(gm+1)*(u1*expan + (u2 + 0.75*(gm-1)*u1^2)*expan^2 + ...
    (u3 + 0.25*(gm+1)*v1^2 + 1.5*(gm-1)*u1*u2 + 0.125*(5*gm^2 - 8*gm + 3)*u1^3)*expan^3);
end

