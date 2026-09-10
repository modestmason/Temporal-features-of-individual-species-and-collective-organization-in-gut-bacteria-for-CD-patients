function [tv] = meanfluc(n)
str=['CD',num2str(n),'.csv'];
data=readmatrix(str, 'NumHeaderLines', 1);
data(:, 1) = [];
[N,T]=size(data);


    for j=1:N
        tv(j)=0;
        for i=1:T-1
            tv(j)=tv(j)+(data(j,i+1)-data(j,i))^2;
        end
        tv(j)=tv(j)/(T-1);
    end
end

