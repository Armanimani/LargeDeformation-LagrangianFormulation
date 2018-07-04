%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this program will calculate the HW 1 with aid of Updated Lagrangian     %
% Description and Full Newton-Raphson method                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear all;
format long;
%.......................Reading data from input file.......................
fid=fopen('12 elements.msh');
fgetl(fid);
fgetl(fid);
Nodes=fscanf(fid,'%g',4*[1 inf])';
fgetl(fid);
fgetl(fid);
fgetl(fid);
Elements=fscanf(fid,'%g',5*[1 inf])';
fgetl(fid);
fgetl(fid);
fgetl(fid);
Restrains=fscanf(fid,'%g',3*[1 inf])';
fgetl(fid);
fgetl(fid);
fgetl(fid);
Load_nodes=fscanf(fid,'%g',2*[1 inf])';
fgetl(fid);
fgetl(fid);
fgetl(fid);
q_increment=str2double(fgetl(fid));
%......................Loading gauss points/weights data...................
load Q_gauss_points;
load Q_gauss_weights;
%.....................some basic defination................................
number_nodes=size(Nodes,1);
number_elements=size(Elements,1);
number_restrains=size(Restrains,1);
number_load_nodes=size(Load_nodes,1);
sigma_global=zeros(4*number_elements,3);

fid1=fopen('output.plt','wt');

E=100000;
v=0.25;
D=E/(1+v)/(1-2*v)*[1-v v 0; v 1-v 0; 0 0 0.5*(1-2*v)];
t=1;

U_global=zeros(2*number_nodes,1);
U=zeros(2*number_nodes,1);
initial_coordinate=zeros(2*number_nodes,1);
total_load=zeros(2*number_nodes,1);
for i=1:number_elements
    output2(i,:)=Elements(i,2:5);
end
output2=output2';
%.......Creating initial coordinate matrix based on Nodes data.............
for i=1:number_nodes
    node=Nodes(i,1);
    initial_coordinate(2*node-1,1)=Nodes(node,2);
    initial_coordinate(2*node,1)=Nodes(node,3);
end

current_coordinate=initial_coordinate;
step_coordinate=initial_coordinate;

Index=zeros(2*number_nodes,1);
k_global=zeros(2*number_nodes,2*number_nodes);
load_index=zeros(2*number_nodes,1);
Incremental_loading=zeros(2*number_nodes,1);
%.......Creating Index matrix to dife the Restrained DOFs..................
for i=1:number_restrains;
    node=Restrains(i,1);
    if (Restrains(i,2)==1)
        Index(2*node-1,1)=1;
    end
    if (Restrains(i,3)==1)
        Index(2*node,1)=1;
    end
end
%.......Calculating the incremental nodal force matrix based on q..........
for i=1:number_load_nodes
    node1=Load_nodes(i,1);
    node2=Load_nodes(i,2);
    x1=Nodes(node1,2);
    x2=Nodes(node2,2);
    ds=abs(x1-x2);
    fy=-q_increment*ds*0.5;
    Incremental_loading(2*node1,1)=Incremental_loading(2*node1,1)+fy;
    Incremental_loading(2*node2,1)=Incremental_loading(2*node2,1)+fy;
end
%.........Finding the desired nodes to apply condition.....................
for i=1:number_nodes
    if (and(Nodes(i,2)==0,Nodes(i,3)==0))
        node1=i;
    elseif (and(Nodes(i,2)==120,Nodes(i,3)==0))
        node2=i;
    end
end

load_step_counter=0;
%.............Starting the load incremental loop ..........................
while (abs((current_coordinate(2*node2-1,1)-current_coordinate(2*node1-1,1)))>22)
    U=zeros(2*number_nodes,1);
    load_step_counter=load_step_counter+1
    k_global=zeros(2*number_nodes,2*number_nodes);
    step_coordinate=current_coordinate;
    %......Loop over elements to generate element_coordinate matrix........
    for el=1:number_elements
        for i=1:4
            node=Elements(el,i+1);
            initial_element_coordinate(i,:)=Nodes(node,1:3);
            current_element_coordinate(i,1)=node;
            current_element_coordinate(i,2)=current_coordinate(2*node-1,1);
            current_element_coordinate(i,3)=current_coordinate(2*node,1);
            step_element_coordinate(i,1)=node;
            step_element_coordinate(i,2)=step_coordinate(2*node-1,1);
            step_element_coordinate(i,3)=step_coordinate(2*node,1);
        end
        %........Finding appropirate node displacements of the element.....
        for i=1:4
            node=Elements(el,i+1);
            U_temp(2*i-1:2*i,1)=U(2*node-1:2*node,1);
        end
        %........Calculating the total stiffness matrix....................
        for i=1:4
            ni=Elements(el,i+1);
            for j=1:4
                nj=Elements(el,j+1);
                k_local=zeros(2,2);
                for gp=1:4
                    xi=Q_gauss_points(1,gp);
                    yi=Q_gauss_points(2,gp);
                    wi=Q_gauss_weights(1,gp);
                    [ N0,DNX0,DNY0,J0] = Shape_functions( xi,yi,step_element_coordinate );
                    [ N,DNX,DNY,J] = Shape_functions( xi,yi,current_element_coordinate );
                    S(1,gp)=sigma_global(4*el-4+gp,1);
                    S(2,gp)=sigma_global(4*el-4+gp,2);
                    S(3,gp)=sigma_global(4*el-4+gp,3);
                    [ Bi,Gi,A_tetai,B_bari ] = Matrices( DNX0,DNY0,i,U_temp );
                    [ Bj,Gj,A_tetaj,B_barj ] = Matrices( DNX0,DNY0,j,U_temp );
                    M_s=[S(1,gp) 0 S(3,gp) 0; 0 S(1,gp) 0 S(3,gp) ; S(3,gp) 0 S(2,gp) 0 ;0 S(3,gp) 0 S(2,gp)];
                    k_local=k_local+B_bari'*D*B_barj*wi*t*det(J0)+Gi'*M_s*Gj*wi*t*det(J0);
                end
                %...........Assembeling k_global matrix....................
                k_global(2*ni-1:2*ni,2*nj-1:2*nj)=k_global(2*ni-1:2*ni,2*nj-1:2*nj)+k_local;
            end
        end
    end
    output_cond(load_step_counter,1)=cond(k_global);
    psy1=-Incremental_loading;
    total_load=total_load+Incremental_loading;
    %...............Removing the coulumns and rows of restrained DOfs......
    ii=0;
    for i=1:2*number_nodes
        if (Index(i,1)==1)
            k_global(i-ii,:)=[];
            k_global(:,i-ii)=[];
            psy1(i-ii,:)=[];
            ii=ii+1;
        end
    end
    k_inverse=k_global^(-1);
    psy2=psy1;
    itteration_counter=1;
    %..................Calculating du based on k_global................
    while (norm(psy2)>0.00001);
        k_global=zeros(2*number_nodes,2*number_nodes);
        for el=1:number_elements
            for i=1:4
                node=Elements(el,i+1);
                initial_element_coordinate(i,:)=Nodes(node,1:3);
                current_element_coordinate(i,1)=node;
                current_element_coordinate(i,2)=current_coordinate(2*node-1,1);
                current_element_coordinate(i,3)=current_coordinate(2*node,1);
                step_element_coordinate(i,1)=node;
                step_element_coordinate(i,2)=step_coordinate(2*node-1,1);
                step_element_coordinate(i,3)=step_coordinate(2*node,1);
            end
            %........Finding appropirate node displacements of the element.....
            for i=1:4
                node=Elements(el,i+1);
                U_temp(2*i-1:2*i,1)=U(2*node-1:2*node,1);
            end
            %........Calculating the total stiffness matrix....................
            for i=1:4
                ni=Elements(el,i+1);
                for j=1:4
                    nj=Elements(el,j+1);
                    k_local=zeros(2,2);
                    for gp=1:4
                        xi=Q_gauss_points(1,gp);
                        yi=Q_gauss_points(2,gp);
                        wi=Q_gauss_weights(1,gp);
                        [ N0,DNX0,DNY0,J0] = Shape_functions( xi,yi,step_element_coordinate );
                        [ N,DNX,DNY,J] = Shape_functions( xi,yi,current_element_coordinate );
                        S(1,gp)=sigma_global(4*el-4+gp,1);
                        S(2,gp)=sigma_global(4*el-4+gp,2);
                        S(3,gp)=sigma_global(4*el-4+gp,3);
                        [ Bi,Gi,A_tetai,B_bari ] = Matrices( DNX0,DNY0,i,U_temp );
                        [ Bj,Gj,A_tetaj,B_barj ] = Matrices( DNX0,DNY0,j,U_temp );
                        M_s=[S(1,gp) 0 S(3,gp) 0; 0 S(1,gp) 0 S(3,gp) ; S(3,gp) 0 S(2,gp) 0 ;0 S(3,gp) 0 S(2,gp)];
                        k_local=k_local+B_bari'*D*B_barj*wi*t*det(J0)+Gi'*M_s*Gj*wi*t*det(J0);
                    end
                    %...........Assembeling k_global matrix....................
                    k_global(2*ni-1:2*ni,2*nj-1:2*nj)=k_global(2*ni-1:2*ni,2*nj-1:2*nj)+k_local;
                end
            end
        end
        ii=0;
        for i=1:2*number_nodes
            if (Index(i,1)==1)
                k_global(i-ii,:)=[];
                k_global(:,i-ii)=[];
                ii=ii+1;
            end
        end
        k_inverse=k_global^(-1);
        du=-k_inverse*psy1;
        ii=1;
        for i=1:2*number_nodes
            if (Index(i,1)==1)
                du_g(i,1)=0;
            else
                du_g(i,1)=du(ii,1);
                ii=ii+1;
            end
        end
        U=U+du_g;
        psy2=zeros(2*number_nodes,1);
        current_coordinate=current_coordinate+du_g;
        %............calculating psy function..............................
        for el=1:number_elements
            for i=1:4
                node=Elements(el,i+1);
                U_temp(2*i-1:2*i,1)=U(2*node-1:2*node,1);
            end
            for i=1:4
                node=Elements(el,i+1);
                initial_element_coordinate(i,:)=Nodes(node,1:3);
                current_element_coordinate(i,1)=node;
                current_element_coordinate(i,2)=current_coordinate(2*node-1,1);
                current_element_coordinate(i,3)=current_coordinate(2*node,1);
                step_element_coordinate(i,1)=node;
                step_element_coordinate(i,2)=step_coordinate(2*node-1,1);
                step_element_coordinate(i,3)=step_coordinate(2*node,1);
            end
            for i=1:4
                ni=Elements(el,i+1);
                psy2_local=zeros(2,1);
                for gp=1:4
                    xi=Q_gauss_points(1,gp);
                    yi=Q_gauss_points(2,gp);
                    wi=Q_gauss_weights(1,gp);
                    [ N0,DNX0,DNY0,J0] = Shape_functions( xi,yi,step_element_coordinate );
                    [ N,DNX,DNY,J] = Shape_functions( xi,yi,current_element_coordinate );
                    F=J'*(J0')^(-1);
                    dummy=0.5*(F'*F-eye(2));
                    strain_vector(1,1)=dummy(1,1);
                    strain_vector(2,1)=dummy(2,2);
                    strain_vector(3,1)=2*dummy(1,2);
                    S(:,gp)=D*strain_vector+sigma_global(4*el-4+gp,:)';
                    sigma_global_temp(4*el-4+gp,:)=S(:,gp)';
                    [ Bi,Gi,A_tetai,B_bari ] = Matrices( DNX0,DNY0,i,U_temp );
                    psy2_local=psy2_local+B_bari'*S(:,gp)*t*wi*det(J0);
                end
                %..................Assembeling psy2........................
                psy2(2*ni-1:2*ni,1)=psy2(2*ni-1:2*ni,1)+psy2_local;
            end
        end
        psy2=psy2-total_load;
        ii=0;
        %..............Removing restrained DOFs from psy2..................
        for i=1:2*number_nodes
            if (Index(i,1)==1)
                psy2(i-ii,:)=[];
                ii=ii+1;
            end
        end
        psy1=psy2;
        itteration_counter=itteration_counter+1;
    end
    sigma_global=sigma_global_temp;
    output_deflection(load_step_counter,1)=U_global(26,1);
    output_itteration_counter(load_step_counter,1)=itteration_counter;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % note that the following part will generate tec plot file format for %
    % output and it can be removed for faster running                     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %................Generating output file for tecplot....................
    clear output1;
    if (or(mod(load_step_counter,10)==0,load_step_counter==1))
        for i=1:number_nodes
            output1(i,1)=step_coordinate(2*i-1,1);
            output1(i,2)=step_coordinate(2*i,1);
            output1(i,6)=current_coordinate(2*i-1,1)-step_coordinate(2*i-1,1);
            output1(i,7)=current_coordinate(2*i,1)-step_coordinate(2*i,1);
            avg_number=0;
            Stress_temp=zeros(1,3);
            for el=1:number_elements
                for j=1:4
                    if (Elements(el,j+1)==i)
                        x_node=Nodes(i,2);
                        y_node=Nodes(i,3);
                        x_max=max(Nodes(Elements(el,2),2),Nodes(Elements(el,4),2));
                        x_min=min(Nodes(Elements(el,2),2),Nodes(Elements(el,4),2));
                        y_max=max(Nodes(Elements(el,2),3),Nodes(Elements(el,4),3));
                        y_min=min(Nodes(Elements(el,2),3),Nodes(Elements(el,4),3));
                        dx=abs(x_max-x_min);
                        dy=abs(y_max-y_min);
                        ds_max=dx;
                        x_c=x_min+dx/2;
                        y_c=y_min+dx/2;
                        for gp=1:4
                            xi=Q_gauss_points(1,gp);
                            yi=Q_gauss_points(2,gp);
                            x_gp=xi*dx/2+x_c;
                            y_gp=yi*dy/2+y_c;
                            ds=sqrt((x_gp-x_node)^2+(y_gp-y_node)^2);
                            if (ds<=ds_max)
                                ds_max=ds;
                                stress_temp=sigma_global(4*el-4+gp,1:3);
                            end
                        end
                        avg_number=avg_number+1;
                        Stress_temp=Stress_temp+stress_temp;
                    end
                end
            end
            output1(i,3:5)=Stress_temp'/avg_number;
        end
        output1=output1';
        fprintf(fid1,'SOLUTIONTIME=%d \n',load_step_counter);
        fprintf(fid1,'VARIABLES = "X" "Y" "Sx" "Sy" "Txy" "u"  "v"\n');
        fprintf(fid1,'ZONE N=  %d, E=  %d, ZONETYPE=FEQuadrilateral, DATAPACKING=POINT\n',number_nodes,number_elements);
        fprintf(fid1,'\n');
        fprintf(fid1,'%15.10f  %15.10f  %15.10f  %15.10f  %15.10f  %15.10f  %15.10f  \n',output1);
        fprintf(fid1,'\n');
        fprintf(fid1,'%d  %d  %d  %d  \n',output2);
    end
    U_global=U_global+U;
end