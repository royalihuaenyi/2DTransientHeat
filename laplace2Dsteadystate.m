%FEM code to solve the steadystate 2D Laplace equation
%f = -k[(d^2)U/dx^2 +(d^2)U/dy^2]
%boundary conditions
%U(0,y)=10, kdU/dy (x,0)=0 ,kdU/dy (2,y)=0, -kdU/dy (x,2)=10000
%Variable descriptions
%K=element matrix for spatial term (U,xx +U,yy)
%F - element forcing vector
%KK -system matrix of K
%FF- system forcing vector
%coordinates - Coordinate values of each node
%Nodes -nodal conectivity of each element
%index - a vector containing system dofs associated with each element
%bcdof-  a vector containing dofs associated with boundary conditions
%bcval - a vector containing boundary condition values associated with  the dofs in bcdof 


%input data for control parameters
Lx=2;                                %Length of the domain along X-axis
By=2;                                %Breadth of the domain along Y-axis
Nx=4;                                %Number of Elements along X-axis
Ny=4;                                %Number of Elements along Y-axis
NNel = 4;                            %number of nodes per element
Ndof = 1;                            %number of dofs per node
k=30;                                %thermal conductivity
bc=10;                               %boundary condition at x=0
f=1000;                              %forcing term
nf=4;                                %number of element boundaries with flux
NNels=2;                             %number of nodes per side of each element

%Function to generate the nodal connectivity,number of ellements and nodal
%coordinates
[coordinates,Nodes,Nnode,Nel] = RectangularMesh(Lx,By,Nx,Ny); 

sdof=Nnode*Ndof;                     %total system dofs

 %input data for boundary conditions
 bcdof(1)=1;                             %first node is constrained
 bcval(1)=bc;                            %whose described value is 10
 bcdof(2)=6;                             %6th node is constrained
 bcval(2)=bc;                            %whose described value is 10
 bcdof(3)=11;                             %11th node is constrained
 bcval(3)=bc;                            %whose described value is 10
 bcdof(4)=16;                             %16th node is constrained
 bcval(4)=bc;                            %whose described value is 10
 bcdof(5)=21;                             %21st node is constrained
 bcval(5)=bc;                            %whose described value is 10
 
 %input for flux boundary conditions
 nflux(1,1)=21; nflux(1,2)=22;
 nflux(2,1)=22; nflux(2,2)=23;
 nflux(3,1)=23; nflux(3,2)=24;
 nflux(4,1)=24; nflux(4,2)=25;


%initialization of matrices and vectors
FF=zeros(sdof,1);                        %initialization of system vector
KK =zeros(sdof,sdof);                    %initialization of system matrix
index =zeros(NNel*Ndof,1);               %initialization of index vector
f1=zeros(NNels*Ndof,1);                  %element flux vector
k1=zeros(NNels*Ndof,NNels*Ndof);         %flux matrix
index1=zeros(NNels*Ndof,1);              %flux index vector

%computation of element matrices,vectors and their assembly
for iel =1:Nel
    nd(1) = Nodes(iel,1);
    nd(2) = Nodes(iel,2);
    nd(3) = Nodes(iel,3);
    nd(4) = Nodes(iel,4);
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
    [F,FF]=Forcingvector(f,k);
    KK=Assembly1(KK,K,index);
    
end
    

% %additional computation due to flux boundary condition
for ifx=1:nf
    nds(1)=nflux(ifx,1);
    nds(2)=nflux(ifx,2);
   
    index1=Fluxelementdof(nds,NNels,Ndof);
    K1=FluxElementMatrix();
    F1=FluxElementForcingVector(k);
   [KK,FF]=Assembly3(KK,FF,K1,F1,index1);
end



 %applying boundary conditions
    [KK,FF]=ApplyBC(KK,FF,bcdof,bcval);
%  %solving the matrix equation
    fsol=KK\FF;
     