function [gates] =WireMapF(gates, W, WN)

  %generate new wire set
  nstr=['['];
  for i=1:size(gates,1)
      cq=gates(i).ControlQubits;
      tq=gates(i).TargetQubits;
      aq=gates(i).Angles;
      if(size(cq,2)==0)
          cqs='';
      else
          cqs=[num2str(WN(find(W==cq))) ','];
      end
      if(size(aq,2)==0)
          aqs='';
      else
          aqs=[',' string(aq)];
      end

      if(size(tq,2)==2)
        
          cqs=[num2str(WN(find(W==(tq(1))))) ','];
          tqs=[num2str(WN(find(W==(tq(2))))) ];
          else
       tqs=num2str(WN(find(W==tq)));
      end
     


     
   nstr=strcat(nstr, strjoin([gates(i).Type 'Gate(' cqs tqs aqs ');'],''));
   nstr=strjoin(nstr);
  end
  nstr=strcat(nstr,']');
  gates=eval(nstr);
end