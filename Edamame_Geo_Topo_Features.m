% This is the Matlab code for extracting geometric and topological features
% used in the study described in:
% Kshitiz Dhakal, Qian Zhu, Bo Zhang, Mao Li, Song Li. (2020). 
% Analysis of shoot architecture traits in edamame reveals potential strategies to improve harvest efficiency" 

load Scale  %This is the image resolution information
load Secondary %First column records the sample ID which has secondary branches; 
               %Second column records the seconday branch ID and its
               %connected primary branch ID (third column)
file=dir('IMG*.csv'); %These data has the branch label information


for idx=1:length(file)
    % only consider the primary branches
    X=csvread(file(idx).name,1,2);
    w=find(Secondary(:,1)==idx);
    if ~isempty(w)
        a=ismember(X(:,1),Secondary(w,2));
        a=find(a==1);
        X(a,:)=[];
    end
    clear V
    V(:,1)=X(:,2);V(:,2)=-X(:,3);
    MB=find(X(:,1)==0); %Main branch was labelled as 0
    E=[];
    s=unique(X(:,1));s(1)=[];
    for j=1:length(s)
        a=find(X(:,1)==s(j));
        E=[E;a(1:end-1) a(2:end)];
        MB=[MB;a(1)];
    end
    [m1 m2]=sort(V(MB,2));
    E=[E;MB(m2(1:end-1)) MB(m2(2:end))];
    
    W=sqrt((V(E(:,1),1)-V(E(:,2),1)).^2+(V(E(:,1),2)-V(E(:,2),2)).^2);
    G=graph(E(:,1),E(:,2),W);
    fun=distances(G,MB(1));
    fun=fun/Scale(idx); % compute the geodesic distance
    s=unique(X(:,1));
    
    %Extract topological features (persistent homology)
    for j=1:length(s)
        a=find(X(:,1)==s(j));
        Diagram(j,:)=[-fun(a(end)),-fun(a(1))]; % persitence barcode could be directly computed for this case
    end
    if idx<10
        dlmwrite(['diagram000',num2str(idx),'.txt'],Diagram,'delimiter',' ');
    elseif idx>9&&idx<100
        dlmwrite(['diagram00',num2str(idx),'.txt'],Diagram,'delimiter',' ');
    elseif idx>99&&idx<1000
        dlmwrite(['diagram0',num2str(idx),'.txt'],Diagram,'delimiter',' ');
    end
    
    %Extract geometric features
    clear Length
    for j=1:length(s)
        a=find(X(:,1)==s(j));
        Length(j)=fun(a(end))-fun(a(1));
        Node(j)=fun(a(1));
    end
    TotalLength(idx)=sum(Length);
    [m1 m2]=sort(Node);
    MBLength(idx)=Length(m2(1));
    NumberPB(idx)=length(s)-1;
    AvgLength(idx)=(TotalLength(idx)-MBLength(idx))./NumberPB(idx);
    w=m1(2:end)-m1(1:end-1);
    a=find(w>=0.5);
    InternodeLength(idx,1:4)=w(a(1:4));
    if isempty(find(w(1:4)<0.5))
        InternodeLength(idx,5)=0;
    else
        InternodeLength(idx,5)=1;
    end
end
%Bottleneck distance was computed on Danforth Center server.
%The instruction can be found https://github.com/danforthcenter/persistent_homology.git
%The output distance matrix is 'analysis.matrix.csv');
BD=csvread('analysis.matrix.csv');
BD=BD+BD';
dlmwrite(['BottleneckDistance_PB.txt'],BD,'delimiter','\t');

Geometry=[TotalLength' MBLength' NumberPB' AvgLength' InternodeLength];
dlmwrite(['Geometry.txt'],Geometry,'delimiter','\t');
