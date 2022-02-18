%FEM code to solve the transient 2D Laplace equation
%10.4167*dT/dt = (d^2)T/dx^2 +(d^2)T/dy^2
%boundary conditions
%T(0,y,t)=25, T(50,y,t)=40, dT/dy (x,0,t)=0 ,dT/dy (x,6,t)=q
%Initial condition: T(x,y,0)=25 over the domain.
%using bilinear rectangular elements and crank nicholson method
%Variable descriptions
%K=element matrix for spatial term (T,xx +T,yy)
%M=element matrix for time-dependent term (T,t)
%F - element forcing vector
%KK -system matrix of K
%MM- system matrix of M
%FF- system forcing vector
%FN - effective system vector
%KN - effective system matrix
%FSOL - Solution vector
%SOL - time hiostory solution of selected nodes
%gcoord - Coordinate values of each node
%nodes -nodal conectivity of each element
%index - a vector containing system dofs associated with each element
%bcdof-  a vector containing dofs associated with boundary conditions
%bcval - a vector containing boundary condition values associated with  the dofs in bcdof 


%input data for control parameters
Lx=50;                               %Length of the domain along X-axis
By=6;                                %Breadth of the domain along Y-axis
Nx=10;                               %Number of Elements along X-axis
Ny=2;                                %Number of Elements along Y-axis
NNel = 4;                            %number of nodes per element
Ndof = 1;                            %number of dofs per node
delta_t=0.1;                         %time step for transient analysis
in_time= 0;                          %initial time
end_time = 10;                       %termination time
Ntime = fix((end_time-in_time)/delta_t);    %number of time increment
a=10.4167;                              %thermal coefficient
b=25;
c=50;


%Function to generate the nodal connectivity,number of ellements and nodal
%coordinates
[coordinates,Nodes,Nnode,Nel] = RectangularMesh(Lx,By,Nx,Ny); 

sdof=Nnode*Ndof;                     %total system dofs

 %input data for boundary conditions
 bcdof(1)=1;                             %first node is constrained
 bcval(1)=b;                            %whose described value is 25
 bcdof(2)=11;                             %11th node is constrained
 bcval(2)=c;                            %whose described value is 40
 bcdof(3)=12;                             %12th node is constrained
 bcval(3)=b;                            %whose described value is 25
bcdof(4)=22;                             %22th node is constrained
 bcval(4)=c;                            %whose described value is 40
bcdof(5)=23;                             %23rd node is constrained
 bcval(5)=b;                            %whose described value is 25
bcdof(6)=33;                             %33rd node is constrained
 bcval(6)=c;                            %whose described value is 40

%initialization of matrices and vectors
FF=zeros(sdof,1);                        %initialization of system vector
FN=zeros(sdof,1);                        %initialization of effective system vector
fsol=zeros(sdof,1);                      %solution vector
sol=zeros(1,Ntime+1);                    %vector containing time history solution
KK =zeros(sdof,sdof);                    %initialization of system matrix
MM=zeros(sdof,sdof);                     %initialization of system matrix
KN =zeros(sdof,sdof);                    %effective system matrix
index =zeros(NNel*Ndof,1);               %initialization of index vector


%computation of element matrices,vectors and their assembly
for iel =1:Nel
    nd(1) = Nodes(iel,1);
    nd(2) =Nodes(iel,2);
    nd(3) =Nodes(iel,3);
    nd(4) =Nodes(iel,4);
    x1=coordinates(nd(1),1);
    y1=coordinates(nd(1),2);
    x2=coordinates(nd(2),1);
    y2=coordinates(nd(2),2);
    x3=coordinates(nd(3),1);
    y3=coordinates(nd(3),2);
    x4=coordinates(nd(4),1);
    y4=coordinates(nd(4),2);
  

    index =Systemdofs(nd,NNel,Ndof);
    K=Kelementmatrix(x1,x2,y1,y4);
    M=a*ElementMatrix(x1,x2,y1,y4);
    KK=Assembly1(KK,K,index);
    MM=Assembly2(MM,M,index);
    
end


%loop for time integration
for in = 1:sdof
    fsol(in)=25;
end
sol(1)=fsol(31);
KN=2*MM+delta_t*KK;

for it=1:Ntime
    FN=delta_t*FF+((2*MM)-(delta_t*KK))*fsol;
    [KN,FN]=ApplyBC(KN,FN,bcdof,bcval);
    fsol=KN\FN;
    sol(it+1)=fsol(31);
end
