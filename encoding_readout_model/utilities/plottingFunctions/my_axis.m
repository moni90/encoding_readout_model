
function my_axis(ax, squared, offset_x, offset_y)

if nargin<3
    offset_x=0;
    offset_y=0;
end
    
if squared
    pos=[offset_x+2.5 offset_y+2.5 6 6];
else
    pos=[offset_x+2.5 offset_y+2.5 4.5 6];
end

set(ax, 'Fontsize', 13, 'LineWidth', 1.5);
set(ax,'Units','Centimeters','OuterPosition',[0 0 7 7],'Position',pos)

% 