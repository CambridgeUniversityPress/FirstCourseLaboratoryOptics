% -------------------------------------------------------------------------------------------
% DATA FITTING DEMONSTRATION 
% Language: Matlab (see Line 158 for modifications needed to run under Octave)
% -------------------------------------------------------------------------------------------
% This script makes some simulated data for the intensity at the output of a michelson as a 
% function of mirror position with added noise.  It then fits the data to a straight line
% using one of Matlab's built-in curve-fitting tools: lsqnonlin. It handles non-uniform
% uncertainties correctly and propagates the data uncertainties into the best fit parameters.
% -------------------------------------------------------------------------------------------

%-----------------------------------------------------------------
% Make the simulated data
%-----------------------------------------------------------------

d=[
   -5.0419    0.4463    0.0089
   -3.6250    0.4798    0.0096
   -2.4634    0.4917    0.0098
   -1.4560    0.4921    0.0098
    0.1898    0.5638    0.0113
    1.3114    0.5110    0.0102
    2.5270    0.5328    0.0107
    3.8530    0.5403    0.0108
    5.0523    0.5406    0.0108
    ];
x = d(:,1);
y = d(:,2);
yerr = d(:,3);
%-----------------------------------------------------------------
% Display the data
%-----------------------------------------------------------------

format compact
format short
disp([x,y,yerr]);                           % The experimenters lab book table
figure(1);                                  % Open a figure in which to plot
h=errorbar(x,y,yerr,'d');                   % Plot the data with the unc. estimates
set(h,'linewidth',1.5);                     % Make the plot lines a bit thicker
set(gca,'xlim',[-10 10]);                   % Choose the x axis limits of the plot
set(gca,'ylim',[0.4 0.65]);                 % Choose the y axis limits of the plot]
xticks(-10:2:10);                           % Make the x axis have ticks on the even integers
grid on;                                    % draw the grid
set(gca,'fontsize',14);                     % Bigger default font for this plot
box off                                     % No box around the plot
xlabel('Mirror motion (nm)');               % Label the x axis
ylabel('P_{out} / P_{one arm}');            % Label the y axis
hold on;                                    % Allow the next plot to share the same axes


%-----------------------------------------------------------------
% Fit the data to "fitfunc"                 % see end of script for functions
%-----------------------------------------------------------------

a0 = [0;1];                                 % Your initial guess at the best fit values.
[a,~,res,~,~,~,jac] = ...                   % "a" are the best fit values, res and jac 
    lsqnonlin(...                           % are used to find the uncertainty in a(1),a(2),... 
        @(a)weighted_residuals(a,@fitfunc,x,y,yerr),... % This "anonymous function" is minimized.
        a0...                               % Note: @(a)weigh... defines the anonymous function for
        );                                  % the weighted_residual with the (x,y,yerr) data given.
xfit = linspace(-8,8,100);                  % The x values for plotting the fit function.
plot(xfit,fitfunc(a,xfit),'--',...          % Plot the best fit (i.e.n use the best-fit function to
    'linewidth',2);                         % generate the y-values of the plot).
hold off;                                   % Allow the next plot to clear the current one.
legend('Data','Best fit line',...           % Put a legend in the lower right corner,
    'location','southeast');                % a.k.a the "southeast" corner.
title('Data and fit');
residual = y-fitfunc(a,x);                  % residual is data minus the fit
X2red = 1/(length(x)-length(a))*sum(...     % The reduced chi-squared of the best fit...
    (residual).^2./yerr.^2);          % should be close to 1.
disp(['X2red = ',num2str(X2red)]);          % In the Matlab command window


%-----------------------------------------------------------------
% Direct estimate of the uncertainties from X2 curvature
%-----------------------------------------------------------------
J=zeros(length(a),1); da=zeros(length(a),1); % Set up the variables
for r=1:length(a)                           
    Jsqr(r)=sum( (jac(:,r)).^2 );           % curvature of the a_r chi-square cut at the minimum
    da(r) =rmss(residual./yerr)/sqrt(Jsqr(r)); % is approx C_r =  2*sum(jac(:,r).^2. Gives parabola.
end
disp('Solution +/- uncertaintes:')          
disp([a,da]);                               % display best-fit values and uncertainties 


%-----------------------------------------------------------------
% OPTIONAL:
%
% Make and display chi-square cuts for each fit variable. Also,
% draw the line corresponding to chi-square increasing by 1.
%-----------------------------------------------------------------

a_1 = linspace(0.4,0.6,200);            % The domain of a_1-axis chi-square cut
a_2 = linspace(0.005,0.015,200);              % The domain of a_2-axis chi-square cut
X2_a1cut = zeros(size(a_1));                % will hold the y-values of the a_1 cut
X2_a2cut = zeros(size(a_2));                % will hold the y-values of the a_2 cut
for s = 1:length(a_1)
    X2_a1cut(s) = ChiSqr([a_1(s),a(2)],@fitfunc,x,y,yerr); % the a_1 cut chi-square values 
end
for s = 1:length(a_2)
    X2_a2cut(s) = ChiSqr([a(1),a_2(s)],@fitfunc,x,y,yerr); % the a_2 cut chi-square values
end

figure(2);                                              % this figure 2 will hold two plot windows
subplot(2,1,1);                                         % this sets up the first of two plot windows
h=plot(a_1,X2_a1cut,'-',a(1),min(X2_a1cut),'o',...      % plot X^2 cut in the a_1 dir, the min, ...
    [min(a_1),max(a_1)],[min(X2_a1cut),min(X2_a1cut)]+1,'--');  % and min+1 line
set(h,'linewidth',2);                                   % use bolder lines
grid on;                                                % draw a grid on the plots
xlabel('a_1','fontname','Times New Roman','fontangle','italic');    % x-axis label
ylabel('\chi^2');                                       % y-axis label
set(gca,'fontsize',16);                                 % make the font size bigger
title('\chi^2 cuts');                                   % add a title to the graph
subplot(2,1,2);                                         % set up the second plot window
h=plot(a_2,X2_a2cut,'-',a(2),min(X2_a2cut),'o',...      % plot X^2 cut in the a_2 dir, the min, ..
    [min(a_2),max(a_2)],[min(X2_a2cut),min(X2_a2cut)]+1,'--'); % and min+1 line
set(h,'linewidth',2);                                   
grid on;
xlabel('a_2','fontname','Times New Roman','fontangle','italic');
ylabel('\chi^2');
set(gca,'fontsize',16);


%-----------------------------------------------------------------
% OPTIONAL:
%
% Generate & Display the chi-squared surface
% (This section for illustrative purposes and can be omitted.)
%-----------------------------------------------------------------

[a1,a2] = meshgrid(...                      % Set up the domain. a1 and a2 are 
    linspace(0.47,0.55,100),...             % matrixes of coordinates (parameters)
    linspace(-0.01,0.03,100)...             % at which to calculate chi-squared
    );

X2 = zeros(size(a1));                       % This is the surface we will be finding
for k=1:size(a1,1)                          % step through all the values of a1 and a2
    for s=1:size(a2,2)                      % in the desired range
        X2(k,s) = ChiSqr([a1(k,s),a2(k,s)],@fitfunc,x,y,yerr); % Uses ChiSqr function def. below
    end                                     % Formula for the chi-squared
end

figure(3);                                  % Open a figure to hold the chi-sqr plot
pcolor(a1,a2,X2);                           % plot the chi-square surface.
colormap(gray(20));                         % White is high, black is low
shading interp;                             % Makes it a bit smoother
caxis(min(min(X2))+[0,100])                 % Set the color range
cbar=colorbar;                              % Color key
set(gca,'fontsize',16);                     % Bigger fonts are more visible
ylabel(cbar,'   \chi^2','rotation',0,'fontsize',16); % y-axis text label
daspect([1,1,1]);                           % Make the axes equally spaced
hold on;
contour(a1,a2,X2,min(min(X2))+[1,1],'w:','linewidth',2); % 1 contour at min(X2) + 1
plot(a(1),a(2),'w.');                       % Best fit values of a1 and a2
xlabel('a_1','fontname','Times New Roman','fontangle','italic');
ylabel('a_2','fontname','Times New Roman','fontangle','italic');
title('\chi^2 surface');
hold off


%-------------------------------------------------------------------
% Functions used in script (these must be at the end of the script
% file in Matlab but at the beginning of the script file in Octave.)
%-------------------------------------------------------------------

1;                                          % (Octave only) script can't start with function def.
% The function which we fit to the data
function y=fitfunc(a,x)
y=a(1)+a(2).*x;                             % straight line (can be changed to anything)
end

% This is the quantity whos least square is to be minimized
function r = weighted_residuals(a,fhandle,x,y,yerr)
r=(feval(fhandle,a,x)-y)./yerr;             % gets the fit function via its "handle"
end

% Convenience function for calculating the chi-square
function C = ChiSqr(a,fhandle,x,y,yerr) 
wr = weighted_residuals(a,fhandle,x,y,yerr);% chisqr is just the quadrature sum of ...
C = sum(wr.^2);                             % the weighted residuals 
end

% Function for to calculate the root mean square of a vector
function result=rmss(v)
result = sqrt( mean( (v).^2 ) );
end