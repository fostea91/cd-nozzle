% pressure_dist.m
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

function p_x = pressure_dist(q, l, p_t, p_e, x)

p_x = exp(((q*l + 2*(log(p_t)-log(p_e)))/(l^3))*(x^3) - ((q*l*2 + 3*(log(p_t)-log(p_e)))/(l^2))*(x^2) + q*x + log(p_t));

end
