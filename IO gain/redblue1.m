function c = redblue1(cormin,cormax)
%REDBLUE    Shades of red and blue color map
%   REDBLUE(M), is an M-by-3 matrix that defines a colormap.
%   The colors begin with bright blue, range through shades of
%   blue to white, and then through shades of red to bright red.
%   REDBLUE, by itself, is the same length as the current figure's
%   colormap. If no figure exists, MATLAB creates one.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(redblue)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

%   Adam Auton, 9th October 2009
n0 = 500;
m1 = floor((1-cormin)*n0);
m2 = floor((cormax-1)*n0);

% From [0 0 1] to [1 1 1], then [1 1 1] to [1 0 0];
r0 = (0:m1-1)'/max(m1-1,1);
b0 = flipud((0:m2-1)'/max(m2-1,1));
r = [r0; ones(m2,1)];
g = [r0; b0];
b = [ones(m1,1);b0];

r = r(2:end-2);
g = g(2:end-2);
b = b(2:end-2);

c = [r g b]; 


