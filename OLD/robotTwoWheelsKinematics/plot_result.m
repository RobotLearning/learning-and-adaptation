function plot_result(t,x_iter,u_app,iter)

figure(2);
for i = 1:size(x_iter)
    subplot(2,3,i); plot(t, x_iter(i,:)); 
end

for j = 1:size(u_app,1)
    subplot(2,3,i+j); plot(t, u_app(j,:,iter)); 
end