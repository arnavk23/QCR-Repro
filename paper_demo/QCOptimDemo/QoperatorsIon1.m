function [OPMLString OPMLMat OQSMLStrings] =Qoperators(dim)
%Get the family of operators to parse here just the rxx Gates

OPMLString={}; %strings for matlab code
OPMLMat={};  %matrices for fast comparison
OQSMLStrings={};

OP={};
OP={};
cN=1;

%OPMLString{end+1}='';
stepsS=pi/2;
%OP=cat(3,OP,XO);
for i=1:dim
    for j=i+1:dim
        if(not(i==j))

            for theta=pi/2:stepsS:pi/2
                str=['rxxGate(' num2str(i) ',' num2str(j) ',' num2str(theta)   ')'];
                strQ=['rxx(' num2str(theta) ') q[' num2str(i-1) '], q[' num2str(j-1) ']; \n'; ];
                %remove identity matrices
                cC=rxxGate(i,j,theta);
                c = quantumCircuit(cC,dim);

                if(cN~=0)
                    OPMLString{end+1}=str;
                    OQSMLStrings{end+1}=strQ;
                end
            end
        end
    end
end



for i=1:dim
    for theta=-pi/2:stepsS:pi/2
        str=['rxGate(' num2str(i) ','  num2str(theta)  ')'];
        %remove identity matrices
        % strQ=['rx(' num2str(theta) ') q[' num2str(i-1)  ']; \n'; ];
        if(theta<0)
            strQ=['rx(-pi/2' ') q[' num2str(i-1)  ']; \n'; ];
        else
            strQ=['rx(pi/2' ') q[' num2str(i-1)  ']; \n'; ];
        end

        cC=rxGate(i,theta);
        c = quantumCircuit(cC,dim);

        if(theta~=0)
            OPMLString{end+1}=str;
            %  OPMLMat{end+1}=targetC;
            OQSMLStrings{end+1}=strQ;
        end
        %OPMLString{end+1}=str;
    end
end


for i=1:dim
    for theta=-pi/2:stepsS:pi/2

        %remove identity matrices
        str=['ryGate(' num2str(i) ','  num2str(theta)  ')'];
        strQ=['ry(' num2str(theta) ') q[' num2str(i-1)  ']; \n'; ];

        cC=ryGate(i,theta);
        c = quantumCircuit(cC,dim);

        if(theta~=0)
            OPMLString{end+1}=str;
            % OPMLMat{end+1}=targetC;
            OQSMLStrings{end+1}=strQ;
        end
    end
end
for i=1:dim
    for theta=-pi/2:stepsS:pi/2

        %remove identity matrices
        str=['rzGate(' num2str(i) ','  num2str(theta)  ')'];
        strQ=['rz(' num2str(theta) ') q[' num2str(i-1)  ']; \n'; ];

        cC=rzGate(i,theta);
        c = quantumCircuit(cC,dim);

        if(theta~=0)
            OPMLString{end+1}=str;
            % OPMLMat{end+1}=targetC;
            OQSMLStrings{end+1}=strQ;
        end
    end
end


end



