function [gatesN] =simplereduce(startgates)
code=[];
gatesN=startgates;
i=1;
%for i=1:size(gatesN,1)-1
while(i<size(gatesN,1)-1)
    gates=gatesN(i:i+1);
    [W WN] =ExtractWireMap(gates);
    gatesNT=WireMapF(gates,W,WN);
    R = quantumCircuit(gatesNT);
    RM=getMatrix(R);
    RM2=eye(size(RM));
    er=compare(RM,RM2);
    if(er<0.001)
        gatesN=gatesN([1:i-1,i+2:end]);
    end
    i=i+1;
end


end
