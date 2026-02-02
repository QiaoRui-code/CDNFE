
clear;
clc;
close all;
%pipi储存第一列
[profile,pipi]= xlsread('Luohe_cesc.xlsx');
[empty,network]= xlsread('net.xlsx');

normal=profile(:,1:63);
case_mprofile=profile(:,1:720);

edge_source=network(:,1);
edge_target=network(:,2);
edge_source1=network(:,1);
edge_target1=network(:,2);
%unique 筛选重复值，结果按照升序排列  node1中心基因
node1=unique(network(:,1));
reference_sample_num=3;
case_mprofile=fillmissing(case_mprofile,'constant',0);
es=0.00000001;
patients_num=[63,140,219,298];

tempcase(:,1,1:patients_num(1))=case_mprofile(:,1:63);     % Stage IA 
tempcase(:,2,1:patients_num(2))=case_mprofile(:,64:203);   % Stage IB
tempcase(:,3,1:patients_num(3))=case_mprofile(:,204:422);  
tempcase(:,4,1:patients_num(4))=case_mprofile(:,423:720);


psize=size(tempcase);
Land_entropy=zeros(length(node1),63,4);
Land_entropy_mix=zeros(length(node1),63,4);
Land_entropy_1=zeros(length(node1),63,4);
Land_entropy_1_mix=zeros(length(node1),63,4);
stage=4;

%%%确定局部网络邻域基因
for l=1:psize(2)
    for  s=1:patients_num(l)
         for na=1:length(node1)
             
             edges_list1=[]; 
             %判定pipi中的元素是否在node1(na)里面出现，
             %cencter_idd 若出现length(node2)中元素为1，否为0； cencter_liang中的元素则为在node1(na)出现的索引
             %find 找非零元素的线性索引值
             [cencter_liang1,cencter_idd1]= find(ismember(pipi,node1(na)));
             %判定是否为空
             if isempty(cencter_liang1)
                 continue;
             else
                 [liang1,nei_idd1]= find(ismember(edge_target1,node1(na)));
                 nei_gene1=edge_source1(liang1);
                 e_=0;
                 for n=1:length(nei_gene1)
                     [nei_liang1,nie_idd1]= find(ismember(pipi,nei_gene1(n)));
                     if ~isempty(nei_liang1)
                         e_=e_+1;
                         edges_list1(e_,:)=[cencter_liang1 nei_liang1];
                     end
                 end
                 
                 [liang1_1,nei_idd1_1]= find(ismember(edge_source1,node1(na)));
                 nei_gene1_1=edge_target1(liang1_1);
                 e1_1=0;
                 for n=1:length(nei_gene1_1)
                     [nei_liang1_1,nie_idd1_1]= find(ismember(pipi,nei_gene1_1(n)));
                     if ~isempty(nei_liang1_1)
                         e1_1=e1_1+1;
                         edges_list1_1(e1_1,:)=[cencter_liang1 nei_liang1_1];
                     end
                 end
             end
            if e_<2
                continue;
            end
            for i=1:e_
                curr_pcc1=abs(corr(normal(edges_list1(i,1),:)',normal(edges_list1(i,2),:)'));
                temp_add_onecase1=[normal(edges_list1(i,1),:),reshape(tempcase(edges_list1(i,1),l,s),1,1)];
                temp_add_onecase2=[normal(edges_list1(i,2),:),reshape(tempcase(edges_list1(i,2),l,s),1,1)];
                curr_pcc_add_onecase=abs(corr(temp_add_onecase1',temp_add_onecase2'));
                delt1=abs(curr_pcc_add_onecase-curr_pcc1);
                p_val1=SNN(delt1,curr_pcc1,reference_sample_num);
            end  
%             if e1_1<2
%                 continue;
%             end
            for i=1:e1_1
                curr_pcc1_1=abs(corr(normal(edges_list1_1(i,1),:)',normal(edges_list1_1(i,2),:)'));
                temp_add_onecase_1=[normal(edges_list1_1(i,1),:),reshape(tempcase(edges_list1_1(i,1),l,s),1,1)];
                temp_add_onecase_2=[normal(edges_list1_1(i,2),:),reshape(tempcase(edges_list1_1(i,2),l,s),1,1)];
                curr_pcc_add_onecase_2=abs(corr(temp_add_onecase_1',temp_add_onecase_2'));
                delt1_1=abs(curr_pcc_add_onecase_2-curr_pcc1_1);
                p_val1_1=SNN(delt1_1,curr_pcc1_1,reference_sample_num);
            end
              JQ=1/5*p_val1+4/5*p_val1_1;
              JQ1(na,2)=JQ;
         end 
              
    end
     l
end
[JQ2,idx]=sort(JQ1,1,'descend');
%快捷缩进 crtl+i
JQ3=idx(1:997,2);
edge_target3=node1(JQ3(1:997));
[liang2,idx1]= find(ismember(edge_target,edge_target3));
network2=network(liang2,1:2);
% 假设 network2 是一个矩阵
csvwrite('network.csv', network2);
