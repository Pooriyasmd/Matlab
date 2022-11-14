x = -10:0.01:10;
plot(x,Rect(x+7),'LineWidth',1.5)

function y = Rect(x)

    y=zeros(size(x)); 
    ind=find(x>0 & x<1);
    y(ind)=1;
    ind=find(x==0 | x==1);
    y(ind)=0.5; % to add up correctly at boundary

end