% generate a Full Graph Model
% and save it for database lookup

clear all;

dim=3;

[OPMLString OPMLMat OQSMLStrings] =QoperatorsIon1(dim);


numOp=size(OPMLString,2);

maxDepth=4;
Id=eye(2^dim);
maxNodes=1;%5000;
maxEdges=1;%12000;
Nodes=zeros(maxNodes,2^dim,2^dim); %resulting Unitary
Nodes(1,:,:)=Id;
%HashNodes=zeros(maxNodes,16);
%HashNodes(1,:) = DataHash(Id, Opt);

ResultDepth(1)=0;

AllEdges=zeros(maxEdges,3);
currentNode=1;
currentEdge=0;
ops=1:numOp;
FreeOperations{1}=ops;  %to be optimized later
FreeOps(1)=size(ops,2);
B=repmat(ops,maxNodes,1);
FreeOperations=num2cell(B,2);
%temporary:
GModel = digraph();
GModel = addnode(GModel,1);
Dists=distances(GModel,1,'Method','unweighted');
%maxNodes=10000;
indN1=find(Dists<(maxDepth));
indN2=find(FreeOps>0);
indN=intersect(indN1,indN2);
while(size(indN,2)>0)
%for iterNT=1:50
    %sel=indN(1);
   % FreeOperations
    R(:,:)=Nodes(indN(1),:,:);
    ops=FreeOperations{indN(1)};
    sOp=randi(size(ops,2));
    code=ops(sOp);
    [RM allMstr]=MatQCode(code, OPMLString,dim);
    %R=R*RM; %right multiplied? 
    R=RM*R;
  
    %result up to global phase

    pos=0;
    for k=1:size(Nodes,1)
        MT(:,:)=Nodes(k,:,:);
        er=compare(R,MT);
        if(er<0.0001)
            pos=k;
            %pos
            break;
        end
    end


    % new node if it is allowed
    if((pos==0))%&&(Dists(indN(1))<maxDepth))
        currentNode=currentNode+1;
        Nodes(currentNode,:,:)=R;
       % HashNodes(currentNode,:) = DataHash(R, Opt);

        currentEdge=currentEdge+1;
        AllEdges(currentEdge,1:3)=[indN(1),currentNode,ops(sOp)];
        ResultDepth(currentNode)=ResultDepth(indN(1))+1;

        FreeOperations{currentNode}=1:numOp;
        FreeOps(currentNode)=numOp;%size(ops,2);

        %  tic
        GModel = addnode(GModel,1);
        GModel = addedge(GModel, indN(1), currentNode, ops(sOp));
        %toc
        % Dists=distances(GModel,1,'Method','unweighted');
       

    else
       % if(pos>0)
        currentEdge=currentEdge+1;
        AllEdges(currentEdge,1:3)=[indN(1),pos,ops(sOp)];
        GModel = addedge(GModel, indN(1), pos, ops(sOp));
      %  end
    end

    ops=ops([1:sOp-1,sOp+1:end]);
    FreeOperations{indN(1)}=ops;
    FreeOps(indN(1))=size(ops,2);

    %Compute new distance transform
    % tic

    if(rand()>0.7)
        Dists=distances(GModel,1,'Method','unweighted');
        ResultDepth=Dists;
    end

    %new index
    indN1=find(Dists<(maxDepth));
    indN2=find(FreeOps>0);
    indN=intersect(indN1,indN2);
 if(rand()>0.98)
       GModel
       Unitaries=[];
OptStrings={};
OptOrder={};
OptLength=[];
count=1;
for k=1:size(Dists,2)
    if(Dists(k)>0)
        [P,d,E]=shortestpath(GModel,1,k,'Method','unweighted');
        vals=GModel.Edges.Weight(E);
        allStr=[];
        for i=1:size(vals,1)
            allStr=[allStr ' ' OPMLString{vals(i)}];%OPSString{vals(i)}];
        end
        allStr;
        Unitaries(count,:,:)=Nodes(k,:,:);
        OptStrings{count}=allStr;
        OptOrder{count}=vals';
        %vals'
        OptLength(count)=size(vals,1);
        count=count+1;
    end
end
       save(['OptUnitDim3_4Temp.mat'],'OptLength','Unitaries','OptOrder','GModel','Nodes');

 end

end

Dists=distances(GModel,1,'Method','unweighted');

GModel
%plot(GModel,'MarkerSize',6,'EdgeColor','k','Layout','force')

%save all shortest paths
Unitaries=[];
OptStrings={};
OptOrder={};
OptLength=[];
count=1;
for k=1:size(Dists,2)
    if(Dists(k)>0)
        [P,d,E]=shortestpath(GModel,1,k,'Method','unweighted');
        vals=GModel.Edges.Weight(E);
        allStr=[];
        for i=1:size(vals,1)
            allStr=[allStr ' ' OPMLString{vals(i)}];%OPSString{vals(i)}];
        end
        allStr;
        Unitaries(count,:,:)=Nodes(k,:,:);
        OptStrings{count}=allStr;
        OptOrder{count}=vals';
        %vals'
        OptLength(count)=size(vals,1);
        count=count+1;
    end
end

save('OptUnitDim3_4.mat','OptLength','Unitaries','OptOrder','GModel','Nodes');

%Just for testing
alers=0;
for i=1:size(Unitaries,1)
 M(:,:)=Unitaries(i,:,:);
 code=OptOrder{i};
  [RM allMstr]=MatQCode(code, OPMLString,dim);
 RM;

%compare(M,RM);
if(( compare(M,RM)>0.001)|| (compare(RM,M)>0.001))
    alers=alers+1;
    code
    break;
end
end
alers