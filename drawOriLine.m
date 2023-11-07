% drawOriLine.m
%
%      usage: drawOriLine(x, y, ph, factor, linewidth, color);
%         by: eli merriam
%       date: 06/10/15
%    purpose: 
%
function h = drawOriLine(x, y, ph, factor, linewidth, color);
% ph: 0 = horizontal, then goes counterclockwise

% check arguments
if ~any(nargin == [6])
  help drawOriLine
  return
end



[xx,yy] = pol2cart(ph, factor);

h = line([x-xx x+xx], [y-yy y+yy]);

set(h, 'linewidth', linewidth, 'color', color);
