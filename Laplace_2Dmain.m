%FEM code to solve the transient 2D Laplace equation
%10.4167*dT/dt = (d^2)T/dx^2 +(d^2)T/dy^2
%boundary conditions
%T(0,y,t)=25, T(50,y,t)=40, dT/dy (x,0,t)=0 ,dT/dy (x,6,t)=50(T-Tf)
%Initial condition: T(x,y,0)=25 over the domain.
%using bilinear rectangular elements and forward difference method
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
%K1- element matrix due to Cauchy-type flux
%F1 - element vector due to flux boundary condition
%index1 - index for nodal dofs with flux



%input data for control parameters
Lx=50;                              %Length of the domain along X-axis
By=6;                                %Breadth of the domain along Y-axis
Nx=10;                               %Number of Elements along X-axis
Ny=2;                                %Number of Elements along Y-axis
NNel = 4;                            %number of nodes per element
ndof = 1;                            %number of dofs per node
Nnode = 33;                          %total number of nodes in the system
sdof=Nnode*ndof;                     %total system dofs
delta_t=0.1;                         %time step for transient analysis
in_time= 0;                          %initial time
end_time = 1200;                       %termination time
Ntime = fix((end_time - in_time)/delta_t);    %number of time increment
a=10.4167;                              %thermal coefficient
nf=10;                                 %number of element boundaries with flux
nnels=2;                             %number of nodes per side of each element
b=298.15;                                %initial temperature
c=323.15;                                 %temperature at the end of the boundary
h=300;                                 %convection coefficient(air convection)
f=298.15;                                 %ambient temperature

%input data for nodal coordinate values

[coordinates,Nodes,Nnode,Nel] = RectangularMesh(Lx,By,Nx,Ny);

 %input data for boundary conditions
 bcdof(1)=1;                             %first node is constrained
 bcval(1)=c;                            %whose described value is 25
 bcdof(2)=11;                             %11th node is constrained
 bcval(2)=c;                            %whose described value is 40
 bcdof(3)=12;                             %12th node is constrained
 bcval(3)=c;                            %whose described value is 25
bcdof(4)=22;                             %22th node is constrained
 bcval(4)=c;                            %whose described value is 40
bcdof(5)=23;                             %23rd node is constrained
 bcval(5)=c;                            %whose described value is 25
bcdof(6)=33;                             %33rd node is constrained
 bcval(6)=c;                            %whose described value is 40
 
 %input for flux boundary conditions
 %nflx(i,j) i- element no: j- two side nodes
 nflx(1,1)=23; nflx(1,2)=24;
 nflx(2,1)=24; nflx(2,2)=25;
 nflx(3,1)=25; nflx(3,2)=26;
 nflx(4,1)=26; nflx(4,2)=27;
 nflx(5,1)=27; nflx(5,2)=28;
 nflx(6,1)=28; nflx(6,2)=29;
 nflx(7,1)=29; nflx(7,2)=30;
 nflx(8,1)=30; nflx(8,2)=31;
 nflx(9,1)=31; nflx(9,2)=32;
 nflx(10,1)=32; nflx(10,2)=33;

%initialization of matrices and vectors
FF=zeros(sdof,1);                        %initialization of system vector
FN=zeros(sdof,1);                        %initialization of effective system vector
fsol=zeros(sdof,1);                      %solution vector
sol=zeros(1,Ntime+1);                    %vector containing time history solution
KK =zeros(sdof,sdof);                    %initialization of system matrix
MM=zeros(sdof,sdof);                     %initialization of system matrix
KN =zeros(sdof,sdof);                    %effective system matrix
index =zeros(NNel*ndof,1);               %initialization of index vector
F1=zeros(nnels*ndof,1);                  %element flux vector
K1=zeros(nnels*ndof,nnels*ndof);         %flux matrix
index1=zeros(nnels*ndof,1);              %flux index vector


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
    xlen=x2-x1;
    ylen=y4-y1;
    
    index =Systemdofs(nd,NNel,ndof);
    K=Kmatrix(xlen,ylen);
    M=a*ElementMatrix(x1,x2,y1,y4);
    KK=Assembly1(KK,K,index);
    MM=Assembly2(MM,M,index);
    
end

%additional computation due to flux boundary condition
for ifx=1:nf
    nds(1)=nflx(ifx,1);
    nds(2)=nflx(ifx,2);
    x1=coordinates(nds(1),1);
    y1=coordinates(nds(1),2);
    x2=coordinates(nds(2),1);
    y2=coordinates(nds(2),2);
    eleng=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    
    index1=Fluxelementdof(nds,nnels,ndof);
    K1=h*CauchyFlux(eleng);
    F1=h*f*fef1l(eleng);
    [KK,FF]=Assembly3(KK,FF,K1,F1,index1);
    
end


%loop for time integration
for in = 1:sdof
    fsol(in)=298.15;
end
sol(1)=fsol(10);
KN=MM+delta_t*KK;

for it=1:Ntime
    FN=delta_t*FF+MM*fsol;
    [KN,FN]=ApplyBC(KN,FN,bcdof,bcval);
    fsol=KN\FN;
    sol(it+1)=fsol(10);
end

%plot the solution at node 31
time =0:delta_t:Ntime*delta_t;
plot(time,sol);
xlabel('Time (s)');
ylabel('Temperature (K)');
title('Transient Behaviour at Node 10');
    
