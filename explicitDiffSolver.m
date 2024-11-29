
%%Debugger
end_value=10;
step_size=0.01;
y_values_exp=expl_euler(2,step_size,end_value,@(x)x);
%y_values_heun=Heun_Solver(2,step_size,end_value,@(x)x);
y_values_Runge=Runge_Kutta(2,step_size,end_value,@(x)x);
step_count=floor(end_value/step_size);
t_values=linspace(0,end_value,step_count);
reference=@(x) x.^2 +2;
hold on
plot(t_values,y_values_exp);
plot(t_values,y_values_heun);
plot(t_values,y_values_Runge);
plot(t_values,reference(t_values));
hold off
%

function y = expl_euler(y_0,dt,t_end, f)
step_count=floor(t_end/dt);
y=zeros(step_count,1);
y(1)=y_0;
for i=2:(step_count)
y(i)= y(i-1)+f(y(i-1))*dt;
end
end


function y = Heun_Solver(y_0,dt,t_end, f)
step_count=floor(t_end/dt);
y=zeros(step_count,1);
y(1)=y_0;
for i=2:(step_count)
next_prediction= y(i-1)+f(y(i-1))*dt;
y(i)=y(i-1)+((f(y(i-1))+f(next_prediction))*dt)/2;
end
end


function y = Runge_Kutta(y_0,dt,t_end, f)
step_count=floor(t_end/dt);
y=zeros(step_count,1);
y(1)=y_0;
for i=2:(step_count)
    k1=f(y(i-1));
    k2=f(y(i-1)+(k1*dt)/2);
    k3=f(y(i-1)+(k2*dt)/2);
    k4=f(y(i-1)+(k3*dt));
y(i)=y(i-1)+dt*((k1+2*k2+2*k3+k4)/6);
end
end