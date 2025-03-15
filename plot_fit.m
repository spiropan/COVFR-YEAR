function [b,stats]=plot_fit(x,y,robust,dots,line,labels)
% usage: [b,stats]=plot_fit(x,y,0,'ko','r',labels)
% labels is OPTIONAL

if nargin < 6
    labels = [];
end

if robust==1
    [b,stats]=robustfit(x,y);
else
    [b,~,stats]=glmfit(x,y);
end

yFitted=b(1)+b(2)*x;

plot(x(:,end),y,dots,x(:,end),yFitted,line,'MarkerSize',14,'LineWidth',1.5);
%xlim([min(x)-0.25,max(x)+0.25])
padx=0.05*(max(x)-min(x));
pady=0.05*(max(y)-min(y));

xlim([min(x)-padx,max(x)+padx])
ylim([min(y)-pady,max(y)+pady])

if ~ isempty(labels)
    labelpoints(x,y,labels,'NE',0.1);
end

