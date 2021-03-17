function colorMap = myColorMapRedBlue_log(this_color)

red_top=[255,0,0]/255;
blue_top=[0,0,255]/255;

tmp=log(linspace(1,5,100))/log(5);


blue_scale = [(tmp*(blue_top(1)-.5)+.5)' (tmp*(blue_top(2)-.5)+.5)' (tmp*(blue_top(3)-.5)+.5)'];
red_scale = [(tmp*(red_top(1)-.5)+.5)' (tmp*(red_top(2)-.5)+.5)' (tmp*(red_top(3)-.5)+.5)'];

colorMap = [flipud(blue_scale); red_scale];

% colorMap = [linspace(green_top(1),.5,128)' linspace(green_top(2),.5,128)' linspace(green_top(3),.5,128)'; linspace(.5,purple_top(1),128)' linspace(.5,purple_top(2),128)' linspace(.5,purple_top(3),128)'];
% colorMap = [(tmp*(.5-red_top(1))+red_top(1))' (tmp*(.5-red_top(2))+red_top(2))' (tmp*(.5-red_top(3))+red_top(3))'; ];
% colorMap = flipud(colorMap);

% colorMap = [[linspace(1,1,64) linspace(1,.5,64)]' [linspace(0,.5,64) linspace(.5,.5,64)]'  [linspace(0,0,64) linspace(0,.5,64)]'; [linspace(.5,0,64) linspace(0,0,64)]' [linspace(.5,.5,64) linspace(.5,0,64)]' [linspace(.5,1,64) linspace(1,1,64)]']
% colorMap = flipud(colorMap);

% colorMap = [linspace(0,1,30)' linspace(0,0,30)' linspace(1,0,30)'];

if nargin>0
colorMap=colorMap(this_color,:);
end


% figure;
% hold on;
% plot([1:40]/40-1/(40),log(1:40)/log(40),'b');
% plot([1:20]/20-1/(20),log(1:20)/log(20),'r');
% plot([1:10]/10-1/(10),log(1:10)/log(10),'r');
% plot([1:5]/5-1/(5),log(1:5)/log(5),'r');