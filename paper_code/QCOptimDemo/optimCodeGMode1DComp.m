%simple code optimizer


clear all;
%close all;

plotit=0

dim=5;
dimD=3;
[OPMLString OPMLMat OQSMLStrings] =QoperatorsIon1(dim);
%dataset information
[OPMLStringDS OPMLMatDS OQSMLStringsDS] =QoperatorsIon1M(dimD);

startCode=randi(size(OPMLString,2),[1,300]);
%load('startCode.mat');

[allMstr]=Code2str(startCode, OPMLString,dim);
startgates=eval(allMstr);

load OptUnitDim3_4.mat % up to 5

maxUnit=size(Unitaries,1);

%Test which parts are commutative
ComCode=zeros(size(startCode));
ComCount=1;
for i=1:size(startCode,2)-1
    gates=startgates(i:i+1);
    gates2=startgates(i+1:-1:i);
    [W WN] =ExtractWireMap(gates);
    gatesN=WireMapF(gates,W,WN);
    gatesN2=WireMapF(gates2,W,WN);

    R = quantumCircuit(gatesN);
    RM=getMatrix(R);
    R = quantumCircuit(gatesN2);
    RM2=getMatrix(R);

    er=compare(RM,RM2);
    if(er<0.001)
        ComCode(i)=ComCount;
        ComCount=ComCount+1;
    end
end

%identify basis operations in substrings
dgates=startgates;
if(plotit==1)
    
    gates=dgates;
    c2 = quantumCircuit(gates,dim);
    figure(1)
    plot(c2);
    drawnow();
end

while(true)%for i=1:999000000
  
    if(rand()>0.3)
        shuffeled=zeros(size(ComCode));
        for j=1:max(ComCode)
            ind=find(ComCode==j);
            if((shuffeled(ind)==0)&&(rand()<0.6))
                dgates(ind:ind+1)=dgates(ind+1:-1:ind);
                shuffeled(ind+1)=1;
            end
        end


        %update ComCode after each shuffle
        ComCode=zeros(1,size(dgates,1));
        ComCount=1;
        for i=1:size(dgates,1)-1
            gates=dgates(i:i+1);
            gates2=dgates(i+1:-1:i);
            [W WN] =ExtractWireMap(gates);
            gatesN=WireMapF(gates,W,WN);
            gatesN2=WireMapF(gates2,W,WN);
            R = quantumCircuit(gatesN);
            RM=getMatrix(R);
            R = quantumCircuit(gatesN2);
            RM2=getMatrix(R);
            er=compare(RM,RM2);
            if(er<0.001)
                ComCode(i)=ComCount;
                ComCount=ComCount+1;
            end
        end
    end
    %toc

    Found=1;
    while(Found==1)
        s1=randi(size(dgates,1)-1);
        e1=min(randi(size(dgates,1)-s1),7); %test blocks up to X

        testgates=dgates(s1:s1+e1);
        [W WN] =ExtractWireMap(testgates);
        gatesN=WireMapF(testgates,W,WN);
        if(size(W,2)<(dimD+1))
            Found=0;
        end
    end
    R = quantumCircuit(gatesN,dimD);
    RM=getMatrix(R);
    %compare to pairs of sampled operators
    foundit=0;
    optPos=0;
    ind=find(OptLength<size(testgates,1));
    size(ind);
    for j=1:size(ind,2)
        RT(:,:)=Unitaries(ind(j),:,:);
        er=compare(RM,RT);
        if(er<0.001)
            foundit=1;
            optPos=ind(j);
            break;
        end
    end

    if(foundit==1)
        % dCodeOld=dCode;
        [allMStr] =Code2str(OptOrder{optPos}, OPMLStringDS,dimD);
        igates=eval(allMStr);
        gatesN2=WireMapF(igates,WN,W);

        dgates=[dgates(1:s1-1); gatesN2; dgates(s1+e1+1:end)];
        dgates=simplereduce(dgates);
        size(dgates)
        % figure(2)
        % subplot(1,4,1)
        % plot(quantumCircuit(testgates,dim))
        % subplot(1,4,2)
        % plot(quantumCircuit(gatesN,dimD))
        % subplot(1,4,3)
        % plot(quantumCircuit(igates,dimD))
        % subplot(1,4,4)
        % plot(quantumCircuit(gatesN2,dim))
        % drawnow();
        generateCodeNr(dgates, startgates,dim,10);
      
        if(plotit==1)
            c2 = quantumCircuit(dgates,dim);
            figure(3)
            plot(c2);
            drawnow();
        end
        %update ComCode after each shuffle
        ComCode=zeros(1,size(dgates,1));
        ComCount=1;
        for i=1:size(dgates,1)-1
            gates=dgates(i:i+1);
            gates2=dgates(i+1:-1:i);
            [W WN] =ExtractWireMap(gates);
            gatesN=WireMapF(gates,W,WN);
            gatesN2=WireMapF(gates2,W,WN);
            R = quantumCircuit(gatesN);
            RM=getMatrix(R);
            R = quantumCircuit(gatesN2);
            RM2=getMatrix(R);
            er=compare(RM,RM2);
            if(er<0.001)
                ComCode(i)=ComCount;
                ComCount=ComCount+1;
            end
        end
    end
    if(size(dgates,1)<5)
        break;
    end
end
dgates

%toc





