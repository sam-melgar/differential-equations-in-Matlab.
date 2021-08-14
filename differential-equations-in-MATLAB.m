%LA3

%________________FUNCTION 1_____________
subplot(2,1,1);

%To see the graphs I shared, make the step value 0.01 and final value 2
step = input('For the FIRST equation: \nPlease enter a step value:\n');

final_x = input('For the FIRST equation: \nPlease enter a final x value:\n');

x= [0:step:final_x]; %initialize x vector with user input

y(1)=1; %at index 1 place 

final_i = (final_x/step);

for i = [1:final_i]  % 10
previousx = x(i);
previousy = y(i);
y(i+1) = find_next_y(previousx, previousy, step);


end
% // Analytical solutions using dot operators 
u=x;
v = 1 + 0.5*exp(-4 .* u) - 0.5*exp(-2.* u);

% // Matlab built-in ODE solver
tspan = x;
y0 = 1;
[t,s] = ode45(@(t,s) 2-exp(-4*t)-2*s, tspan, y0);




plot(x,y, '-mo', 'linewidth', 2, 'markersize', 3, 'markeredgecolor', 'g', 'markerfacecolor', 'y')
text(0.1, 0.95,'\leftarrow Euler''s y_i_+_1 = y_i + h*(2-e^-^4^x-2y)', 'Color', 'm', 'FontSize', 10)
text(1, 0.87,'Graph: step value <-> 0.01 final x value <-> 2.', 'Color', 'black', 'FontSize', 10)
hold on
plot(u,v,'r','linewidth', 2)
text(0.13,0.92,'\leftarrow Analytical solution: y(x) = 1 + 0.5*e^-^4^x - 0.5*e^-^2^x','Color','red','FontSize',10)
plot(t,s,'-bo')
text(1,0.97,'MATLAB45 solution \rightarrow','Color','blue','FontSize',10)
hold off 


xlabel('\bf x values');
ylabel('\bf y values');
title('\bf LA3: \it Euler''s Method \rm Function 1 I DID IT!! :D');
 
legend('Euler''s approximation answer', 'Analytical Solution', 'MATLAB ODE Solver', 'Location','best')

%___________________FUNCTION 2_____________________________________

subplot(2,1,2);

%To see the graphs I shared, make the step value 0.1 and final value 3
step = input('For the SECOND equation: \nPlease enter a step value:\n');
 
final_x = input('For the SECOND equation: \nPlease enter a final x value:\n');
 
x_two= [0:step:final_x]; %11 values in this array for xf=1
 
y_two(1)=0; %at index 1 place 
 
final_i = (final_x/step);
 
for i = [1:final_i]  % 10
previousx = x_two(i);
previousy = y_two(i);
y_two(i+1) = find_next_y2(previousx, previousy, step);
 
 
end

% // Analytical solutions using dot operators 
u_two=x_two;
v_two = exp(u_two./ 2).*sin(5.*u_two);

% // Matlab built-in ODE solver
tspan2 = x_two;
y2_0 = 0;
[t_two,s_two] = ode45(@(t_two,s_two) s_two - 0.5*exp(t_two/2)*sin(5*t_two) + 5*exp(t_two/2)*cos(5*t_two), tspan2, y2_0);


plot(x_two,y_two, '-mo', 'linewidth', 2, 'markersize', 3, 'markeredgecolor', 'g', 'markerfacecolor', 'y')
text(1.8, 5,'y_i+1 = y_i + h*(y - 0.5*e^(x/2)*sin(5*x)+5*e^(x/2)*cos(5*x)) \rightarrow', 'Color', 'm')
text(0,7.5,'Graph: step value <-> 0.1 final x value <-> 3.', 'Color', 'black', 'FontSize', 10)
hold on
plot(u_two,v_two,'r', 'linewidth', 2)
text(2.3,-3,'\leftarrow Analytical solution: y(x) = e^x^/^2 * sin(5x)','Color','red','FontSize',10)
plot(t_two,s_two,'-bo')
text(2.5,0,'\leftarrow MATLAB45 solution','Color','blue','FontSize',10)
hold off 


xlabel('\bf x values');
ylabel('\bf y values');
title('\bf LA3: \it Euler''s Method \rm Function 2')
legend('Euler''s approximation answer','Analytical Solution', 'MATLAB ODE Solver', 'Location','best')


%comment paragraph:
% Using the analyitical solutions as an "answer reference" these graphs
% show that Euler's method is only an approximation. While yes, if you make
% the step value smaller you will APPROACH a more accurate answer which can
% be seen by the Euler's approximation answer line getting closer to the
% analytical solution line. This phenomena is seen in both examples.
% Another important finding is that the MATLAB ODE solver outputs a very
% accurate answer, the most accurate. If you zoom in super super close you
% can see a divergence between the MATLAB ODE solver and the analytical
% solution however between Euler's approximation and MATLAB, Matlab is the
% most accurate. 
