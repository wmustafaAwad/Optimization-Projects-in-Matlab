%To Draw Required Figures:

function [] = DrawIPHistory(history,problem,model_name)

h=history;
p=problem;

if(size(h.x_hist,1)<2)
    error("The dimensionality of the problem must be at least 2 to plot the results.");
end


x1=linspace(min(h.x_hist(1,:)),max(h.x_hist(1,:)),4);

x1pts=[];
x2pts=[];
for i=1:size(p.A,1)
    x1pts=[x1pts,x1];
    x2pts=[x2pts, (p.b(i)-p.A(i,1).*x1)./p.A(i,2)];
end

x1pts=reshape(x1pts,[numel(x1pts),1]);
x2pts=reshape(x2pts,[numel(x2pts),1]);
%fill(x1pts,x2pts,'r')


%Plot 1: Function value vs iteration
figure(1)
f_vals=p.c(1:2)'*h.x_hist(1:2,:); %Apply transformation (c) to all history points to calculate function
scatter(1:h.iters,f_vals,'*r');
hold on
plot(1:h.iters,f_vals,'--b');
if(p.maximize)
    title(strcat(model_name,' : Function value vs iteration (Maximizing)'))
    ylim([  min(f_vals) max(f_vals)+0.05*max(f_vals)  ])
else
    title(strcat(model_name,' : Function value vs iteration (Minimizing)'))
    ylim([  min(f_vals)-0.05*min(f_vals) max(f_vals)  ])
end
ylabel('Function Value')
xlabel('Iterations')

if(problem.save_figures)
    saveas(figure(1),strcat('./results/',model_name,'_function_val_vs_iter.jpg'))
end


%Plot2: Central Path
figure(2)
P=[x1pts,x2pts];
try
    [k,av] = convhull(P);
    fill(P(k,1),P(k,2),'g')
    alpha(0.15)
    hold on
catch
    for i=1:size(p.A,1)
        pt1=( p.b(i)-p.A(i,1).*min(h.x_hist(1,:)) )./ p.A(i,2);
        pt2=( p.b(i)-p.A(i,1).*max(h.x_hist(1,:)) )./ p.A(i,2);
        plot( [min(h.x_hist(1,:)) max(h.x_hist(1,:))], [pt1 pt2],'black');
        hold on;
    end
end
    

%plotregion(-plottingA,b,[2.5,1.45],[3.3,1.8],[0.7,0.2,0.3]);
plot(h.x_hist(1,:),h.x_hist(2,:),'+r') %Plot first two dimensions only
hold on
plot(h.x_hist(1,:),h.x_hist(2,:),'--b') %Plot first two dimensions only
title(strcat(model_name,' : Centeral Path'))
xlabel('X1')
ylabel('X2')
if(problem.save_figures)
    saveas(figure(2),strcat('./results/',model_name,'central_path.jpg'));
end




%Plot3: Complementarity Condition
figure(3)
x1s1=diag(h.x_hist(1,:))*diag(h.s_hist(1,:))*ones(h.iters,1);
x2s2=diag(h.x_hist(2,:))*diag(h.s_hist(2,:))*ones(h.iters,1);
scatter(x1s1,x2s2,'or');
hold on
plot(x1s1,x2s2,'--b');
title(strcat(model_name,' : Complementarity Condition'))
xlabel('X1S1')
ylabel('X2S2')

if(problem.save_figures)
    saveas(figure(3),strcat('./results/',model_name,'complementarity.jpg'))
end

end





