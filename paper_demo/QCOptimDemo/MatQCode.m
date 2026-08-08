function [RM allMStr] =MatQCode(order, OPMLString,dim)

allStr=[];
RM=[];
 CurrentOP={};  
 for i=1:size(order,2)
     CurrentOP{end+1}=OPMLString{order(i)};
 end


%tic
 x0=[];

allMStr=['['];
for i=1:size(order,2)
    if(i==1)
      allMStr=[allMStr ' ' CurrentOP{i}];
    else
        allMStr=[allMStr ' ; ' CurrentOP{i}];
    end
end
allMStr=[allMStr ']' ];

gates=eval(allMStr);
if(size(gates,1)>0)
R = quantumCircuit(gates,dim);
RM=getMatrix(R);
end


end
