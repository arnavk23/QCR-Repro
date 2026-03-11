function [] =generateCode( dgates, startgates, dim, Nr)

%generate code for qiskit


%adapt dimension !
allStr=['OPENQASM 2.0;' sprintf('\n') 'include "qelib1.inc";' sprintf('\n') ... 
    ' qreg q[' num2str(dim) ']; ' sprintf('\n') ' creg c[1]; \n ' sprintf('\n') ''];
allMStr=allStr;

%generate new wire set
  nstr=[''];
  for i=1:size(dgates,1)
      cq=dgates(i).ControlQubits;
      tq=dgates(i).TargetQubits;
      aq=dgates(i).Angles;
      if(size(cq,2)==0)
          cqs='';
      else
          cqs=[num2str(cq-1) ];
      end
      if(size(aq,2)==0)
          aqs='';
      else
          aqs=[ string(aq)];
      end
      tqs=num2str(tq-1);
      
   if(strcmp(dgates(i).Type,'rxx'))
          %   nstr=strcat(nstr, strjoin([dgates(i).Type  ' q[' cqs '], ' 'q[' tqs '];\n'  ],''));
     nstr=strcat(nstr, strjoin([dgates(i).Type '(' aqs ')'  ' q[' num2str(tq(1)-1) '], ' 'q[' num2str(tq(2)-1) '];\n'  ],''));
  
   else
         nstr=strcat(nstr, strjoin([dgates(i).Type '(' aqs ')'   ' q[' tqs '];\n'  ],''));
   end

   nstr=strjoin(nstr);

  end
allStr=strcat(allStr, nstr);
FNstr=strcat('shortcode', num2str(Nr),'.txt');
fileID = fopen(FNstr,'w');
fprintf(fileID,allStr);
fclose(fileID);


%adapt dimension !
allStr=['OPENQASM 2.0;' sprintf('\n') 'include "qelib1.inc";' sprintf('\n') ... 
    ' qreg q[' num2str(dim) ']; ' sprintf('\n') ' creg c[1]; \n ' sprintf('\n') ''];
allMStr=allStr;

%generate new wire set
  nstr=[''];
  for i=1:size(startgates,1)
      cq=startgates(i).ControlQubits;
      tq=startgates(i).TargetQubits;
      aq=startgates(i).Angles;
      if(size(cq,2)==0)
          cqs='';
      else
          cqs=[num2str(cq-1) ];
      end
      if(size(aq,2)==0)
          aqs='';
      else
          aqs=[ string(aq)];
      end
      tqs=num2str(tq-1);
      
   if(strcmp(startgates(i).Type,'rxx'))
            % nstr=strcat(nstr, strjoin([startgates(i).Type  ' q[' cqs '], ' 'q[' tqs '];\n'  ],''));
     
     nstr=strcat(nstr, strjoin([startgates(i).Type '(' aqs ')'  ' q[' num2str(tq(1)-1) '], ' 'q[' num2str(tq(2)-1) '];\n'  ],''));
  else
         nstr=strcat(nstr, strjoin([startgates(i).Type '(' aqs ')'   ' q[' tqs '];\n'  ],''));
   end

   nstr=strjoin(nstr);

  end
allStr=strcat(allStr, nstr);
FNstr=strcat('longcode', num2str(Nr),'.txt');
fileID = fopen(FNstr,'w');
%fileID = fopen('longcode4_1.txt','w');
fprintf(fileID,allStr);
fclose(fileID);



end
