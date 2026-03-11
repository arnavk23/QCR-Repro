function val=compare(Res,FB)
%example for global phase:
% compare(RXF(1,pi/4,1),[1 -i ; -i 1])
% compare(RZF(1,pi/8,1),TF(1,1))
% compare(RZF(1,pi/4,dim),[1 0 ; 0 i])


Rv=reshape(Res,[size(Res,1)*size(Res,2),1]);
Fv=reshape(FB,[size(FB,1)*size(FB,2),1]);

%R*x=F ->
scale=  Rv'*Fv /(Rv'*Rv);
Res2=Res*scale;

%val=norm(Res-FB,2);

% % up to global phase
% %check if topology is the same
% f1=find(Res==0);
% t1=sum(abs(FB(f1)));
% 
% if(t1==0)
% 
% out1 = Res./FB;      %// Perform the elementwise division
% out1(FB==0)=0 ; % shift zeros back
% 
% out2 = FB./Res;      %// Perform the elementwise division
% out2(Res==0)=0 ; % shift zeros back
% 
% 
% 
% %prüfe ob Einträge gleich sind)
% mO1=out1(find(abs(out1)~=0,1));
% eqC=abs(sum(sum(out1))-out1(find(abs(out1)~=0,1))*size(find(abs(out1)~=0),1));
% 
% Res2=Res;
% if(eqC<0.00001)
%   Res2=Res2./mO1;
% end
% 
% else
%     Res2=Res;
% end

val=norm(Res2-FB,2);

end
