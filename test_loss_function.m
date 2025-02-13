x = -1:0.01:1;
y = x;
[X,Y] = meshgrid(x,y);
c = 0.001;

% f=(abs(alpha-2)/alpha).*(((x./c).^2 /abs(alpha-2)+1).^(alpha/2)-1);


% for alpha=-inf
%     f   = (abs(alpha-2)/alpha).*(((x./c).^2 /abs(alpha-2)+1).^(alpha/2)-1);
    f1  = (1-exp(-0.5.*(x./c).^2));
    df1 = (x./c^2).*exp(-0.5.*(x./c).^2);
    %
    f2 = log(0.5.*(x./c).^2+1);
    df2 = (2*x)./(x.^2+2*c^2);
    %
    f3 = (x/1).^2;
    %
    f4 = 20.*exp(-(0.5/1000).*(x./c).^2);
    %
    f5 = log(0.5.*(x./c).^2+1) + 20.*exp(-(0.5/1000).*(x./c).^2);
%     fg = 
    
% end
figure;
subplot(1,2,1); hold on; grid on;
plot(x,f1,'r',x,f2,'b',x,f3,'g',x,f4,'y');
subplot(1,2,2); hold on; grid on;
plot(x,f5)