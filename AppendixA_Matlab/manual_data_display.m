% Demonstrates how to enter and display data

data = [
%   X       deltaX  Y       deltaY
    0.9     0.2     2.1     0.2
    1.4     0.3     3.4     0.3
    2.1     0.2     3.1     0.6
    2.7     0.2     4.8     0.5
    3.4     0.2     5.1     0.2
    ];

dxneg = data(:,2);      % leftside x error bars
dxpos = data(:,2);      % rightside x error bars
dyneg = data(:,4);      % lower y error bars
dypos = data(:,4);      % upper y error bars

errorbar(data(:,1),... % continues on next line
		 data(:,3),dyneg,dypos,dxneg,dxpos,'o');
axis([0 4 0 7]);		% set axes limits
grid('on');

hold('on');				% allow overplotting
plot([0.5,3.75],[1.5, 5.8],...
     'r-');             % "by eye" fit
hold('off');			

xlabel('{\Delta}L   ( nm )');
ylabel('T   ( {^\circ}C )');
shg;                    % bring the plot forward
