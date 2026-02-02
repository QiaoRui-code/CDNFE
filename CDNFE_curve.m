clear;
clc;
close all;
[data_entropy,empty] = xlsread('E:\CNFE\CESC_entropy_matrix08.csv');

pnum=[1,3,2,1];
data_entropy_size=size(data_entropy);

n=1;
for t=1:length(pnum)
    count=0;
    for c=n:sum(pnum(1:t))   %psize(3)
        count=count+1;
        tmp_entropy=sort(data_entropy(:,c),'descend');
        aver_entropy(t,count)=mean(tmp_entropy(1:floor(data_entropy_size(1)*0.5)));
    end
    case_result(t)=mean(aver_entropy(t,1:pnum(t)));
    n=1+sum(pnum(1:t));
end

tmp_entropy2=sort(abs(data_entropy),1,'descend');
% result(1)=mean(case_result(1:2));
% result(2:6)=case_result(3:7);
figure('Position', [100, 100,400,400]);

t=[1 2 3 4 ];
plot(t,case_result,'r-*','LineWidth',0.8);
set(gca,'XTick',1:4);
B={'I'  'II'  'III' 'IV' };
set(gca,'XTickLabel',B);
% ylim([33,44])
%set(gca,'Ylim',[33 44]);
set(gca,'YTick');
xlabel('Stages','FontName', 'Arial', 'FontSize', 10);
ylabel('CDNFE scores','FontName', 'Arial', 'FontSize', 10);
%plot(t,aver_comidx,'r','LineWidth',3);
title('CDNFE scores of DNBs for CESC ','FontName', 'Arial', 'FontSize', 10);



% 创建图形窗口
figure('NumberTitle', 'off', 'Name', 'The landscapes of DNB for data','Position', [100, 100,400,400]);

% 绘制热图
surf([1:size(abs(data_entropy),2)], [1:size(abs(data_entropy),1)], abs(data_entropy));

% 设置坐标轴和标签
set(gca, 'XTick', 1:4,'LineWidth',0.8);  % 假设有4个时间点或条件
B = {'I'  'II'  'III'  'IV'};  % 时间点或条件的标签
set(gca, 'XTickLabel', B);
xlabel('Time','FontName', 'Arial', 'FontSize', 7);
ylabel('Gene','FontName', 'Arial', 'FontSize', 7);
zlabel('Score','FontName', 'Arial', 'FontSize', 7);
%screenposition = get(gcf,'Position');
%set(gcf,...
   % 'PaperPosition',[0 0 screenposition(3:4)],...
   % 'PaperSize',[screenposition(3:4)]);
%print -dpdf -painters test.pdf


% 添加标题和着色
title('The landscapes of DNB for CESC','FontName', 'Arial', 'FontSize', 7);
%exportgraphics(gcf,'peaks.pdf','ContentType','vector');
shading interp;



