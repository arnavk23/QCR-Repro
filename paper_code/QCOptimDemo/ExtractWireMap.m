function [W WN] =ExtractWireMap(gates)

%Find the involved wires
 wires=[];
 for i=1:size(gates,1)
wires=[wires gates(i).ControlQubits gates(i).TargetQubits ];
 end
  W=unique(wires);
  WN=1:size(W,2);
  

end
