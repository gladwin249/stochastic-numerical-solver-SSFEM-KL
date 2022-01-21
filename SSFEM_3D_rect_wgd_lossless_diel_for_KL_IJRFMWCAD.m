%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NAME OF FILE: SSFEM_3D_rect_wgd_lossless_diel_for_KL_IJRFMWCAD.m
%
% PURPOSE: SSFEM for 3D dielectric wavegudie with epsilon_r
% variations in the dielectric.
%
% Written by Gladwin Jos
% 
% Date : 29/12/2021: Applied SSFEM to rectangular
% waveguide with dielectric block with KL
%                    for IJRFMWCAD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference: 
% Modelling of Spatial Permittivity Variations using Karhunen
% Loève Expansion for Stochastic Electromagnetic Problems

clc
clear all;  % Clear all variables
close all;  % Close all graphics 
tic;
% Loading the mesh file for the example in Figure 3
 mesh_file_name='3d_dielectric_brick_rectangular_waveguide_2_1_24_NE_coarser_Ne_2057.msh';
[node_list one_node_point two_node_line three_node_triangle ...
    four_node_tetrahedron ]=parse_gmsh_file_3D(mesh_file_name);
 
N_SSFEM=10000;   

node_list=node_list/100; % in meters
node_list=roundn(node_list,-6); % rounding the values to 10^-6
x=node_list(:,1);
y=node_list(:,2);
z=node_list(:,3);
x_max=max(x);
x_min=min(x);
y_max=max(y);
y_min=min(y);
z_min=min(z);
z_max=max(z);
if z_max == 0 
    temp=z_max;
    z_max=z_min;
    z_min=temp;
end   
input_len=z_min;
output_len=z_max;

% defining dielectric dimension inside the 
% rectangular waveguide (Ref: Figure 3)
b_value=1;
diel_len=0.8;
err=0.00001*b_value;
X_left=(0.556*b_value)/100-err;
X_right=(1.444*b_value)/100+err;  
Y_bottom=(0*b_value)/100;
Y_top=(0.399*b_value)/100;
Z_left=((2.4*b_value-diel_len*b_value)/2)/100;
Z_right= Z_left+(diel_len*b_value)/100;

% Generating the features of the tetrahedron finite element
% Ref: [30], [31] of the reference paper
Ne=length(four_node_tetrahedron) 
tet_edge=[1 2;1 3;1 4;2 3;4 2;3 4];
local_edge_array= edge_array_local_3D(four_node_tetrahedron,Ne,tet_edge);
[global_edge_array_final local_edge_array_elm_no] = edge_array_global_3D_unique(local_edge_array);
edge_element=element_edge_array_3D_reshape(Ne,local_edge_array_elm_no);
sign_edge_element=obtain_sign_edge_element(four_node_tetrahedron,Ne);
local_face_array= face_array_local_3D(four_node_tetrahedron,Ne);
[global_face_array local_face_array_elm_no global_face_array_element_no]= ...
    face_array_global_3D_unique( local_face_array);
face_element_array=element_face_array_3D_reshape(Ne,local_face_array_elm_no);
[face_surface face_surface_element]= ...
    surface_faces_3D(global_face_array,face_element_array,global_face_array_element_no,Ne);
edge_surface=surface_edges_3D_unique(face_surface);
edge_surface_element=edge_surface_element(edge_surface,global_edge_array_final);
[input_face_surface input_face_surface_element_no output_face_surface output_face_surface_element_no pec_surface]= ...
    input_output_face_surface(face_surface,node_list,face_surface_element,input_len,output_len);
input_normal_vector= ...
    find_normal_vector(input_face_surface,node_list,input_face_surface_element_no,four_node_tetrahedron,sign_edge_element);
output_normal_vector= ...
    find_normal_vector(output_face_surface,node_list,output_face_surface_element_no,four_node_tetrahedron,sign_edge_element);
No_of_edges=length(global_edge_array_final); 
[pec_edge_surface  pec_edge_element_no]=get_pec_edges(pec_surface,global_edge_array_final);
edge_element_no_input=get_input_output_edge_element(input_face_surface,global_edge_array_final);
edge_element_no_output=get_input_output_edge_element(output_face_surface,global_edge_array_final);
n=four_node_tetrahedron;
% generating the edge sign array    
abs_edge=abs(edge_element); 

% dimensions
a=x_max-x_min;
b=y_max-y_min;
E0=1;
width=a;
height=b;
freq=10*10^9;  % frequency of operation
mu_0=4*pi*10^(-7); 
epsilon_0=8.854187817*10^(-12);
k0=2*pi*freq*sqrt(mu_0*epsilon_0); 
if k0^2 >=(pi/a)^2
    k_z_10=sqrt(k0^2-(pi/a)^2);
else
    k_z_10=sqrt((pi/a)^2-k0^2);
end
gamma_input=k_z_10*1i;
mu_r=1+0*1i; 
 
% eigenvalues and eigenvectors obtained as in Appendix B
x_node_1=[X_left X_right];
x_node_2=[Y_bottom Y_top];
x_node_3=[Z_left Z_right];
corr_length_1=0.1*(x_node_1(2)-x_node_1(1));
corr_length_2=0.1*(x_node_2(2)-x_node_2(1));
corr_length_3=0.1*(x_node_3(2)-x_node_3(1)); 
corr_factor_1=corr_length_1*(2/(x_node_1(2)-x_node_1(1)));
corr_factor_2=corr_length_2*(2/(x_node_2(2)-x_node_2(1)));
corr_factor_3=corr_length_3*(2/(x_node_3(2)-x_node_3(1)));
No_of_terms_KL=15;
epsilon_r_dash_dielectric=4; % mean value of permittivity
epsilon_r_air=1;  
std_dev_1=0.1*epsilon_r_dash_dielectric;
std_dev_2=0.1*epsilon_r_dash_dielectric;
std_dev_3=0.1*epsilon_r_dash_dielectric;
[eig_value_3D,eig_vector_3D]= ...
    get_eigen_value_eigen_vector_3D_analytical ...
    (std_dev_1,std_dev_2,std_dev_3,x_node_1,x_node_2,x_node_3,corr_factor_1,corr_factor_2,corr_factor_3,No_of_terms_KL,x_node_2(2)-x_node_2(1));

% Mean value of permittivity for the tetrahedran elements
epsilon_r_final=zeros(Ne,1);
dielectric_count=0;
for e=1:Ne
     node1=n(e,1);
     node2=n(e,2);
     node3=n(e,3);
     node4=n(e,4);
     if x(node1)>=X_left && x(node1)<=X_right && ...
            x(node2)>=X_left && x(node2)<=X_right && ...  
            x(node3)>=X_left && x(node3)<=X_right && ...
            x(node4)>=X_left && x(node4)<=X_right && ...
            y(node1)>=Y_bottom && y(node1)<=Y_top && ...
            y(node2)>=Y_bottom && y(node2)<=Y_top && ...  
            y(node3)>=Y_bottom && y(node3)<=Y_top && ...
            y(node4)>=Y_bottom && y(node4)<=Y_top && ...
            z(node1)>=Z_left && z(node1)<=Z_right && ...
            z(node2)>=Z_left && z(node2)<=Z_right && ...  
            z(node3)>=Z_left && z(node3)<=Z_right && ...
            z(node4)>=Z_left && z(node4)<=Z_right 
                  epsilon_r_final(e,1)= epsilon_r_dash_dielectric; 
                  dielectric_count=dielectric_count+1; 
     else
              epsilon_r_final(e,1)=epsilon_r_air;   
     end     
end
% Generating Boundary matrices at the input and output port
% Ref: Appendix C of the reference paper
[B1ij,B2ij]=obtain_boundary_matrix_input_output_port(four_node_tetrahedron,node_list,edge_element, ...
    sign_edge_element,input_face_surface_element_no,output_face_surface_element_no,No_of_edges, ...
    output_face_surface,input_face_surface,input_normal_vector,output_normal_vector);
f_mean=obtain_excitation_matrix(No_of_edges,Ne,input_face_surface,edge_element, ...
    input_face_surface_element_no,four_node_tetrahedron,sign_edge_element, ...
    input_normal_vector,node_list,k_z_10,E0,width);
% Generating the mean stiffness matrix in equation (13) of the reference paper
% Ref: Appendix C of the reference paper
[K_mean]=obtain_K_mean_3D_wgd_SSFEM ...
    (mu_r,tet_edge,No_of_edges,n,node_list,epsilon_r_final, ...
    edge_element,sign_edge_element,B1ij,B2ij,k0,gamma_input);

% Applying PEC boundary conditions
K_mean(abs(pec_edge_element_no),:)=[];
K_mean(:,abs(pec_edge_element_no))=[];
f_mean(abs(pec_edge_element_no),:)=[]; 
 
No_of_rv=No_of_terms_KL;
xi_value=zeros(N_SSFEM,No_of_terms_KL);
for ii=1:No_of_rv
    seed_value=ii;
    xi_value(:,ii) = random_number_seed_gaussian(N_SSFEM,seed_value);   
end
% random field generation as in equation (4) of the reference paper
eig_vector_centroid=zeros(Ne,No_of_terms_KL);
for ii=1:No_of_terms_KL
    for e=1:Ne
         node1=n(e,1);
         node2=n(e,2);
         node3=n(e,3);
         node4=n(e,4);
         if x(node1)>=X_left && x(node1)<=X_right && ...
            x(node2)>=X_left && x(node2)<=X_right && ...  
            x(node3)>=X_left && x(node3)<=X_right && ...
            x(node4)>=X_left && x(node4)<=X_right && ...
            y(node1)>=Y_bottom && y(node1)<=Y_top && ...
            y(node2)>=Y_bottom && y(node2)<=Y_top && ...  
            y(node3)>=Y_bottom && y(node3)<=Y_top && ...
            y(node4)>=Y_bottom && y(node4)<=Y_top && ...
            z(node1)>=Z_left && z(node1)<=Z_right && ...
            z(node2)>=Z_left && z(node2)<=Z_right && ...  
            z(node3)>=Z_left && z(node3)<=Z_right && ...
            z(node4)>=Z_left && z(node4)<=Z_right
                  centroid_x=(x(node1)+x(node2)+x(node3)+x(node4))/4;
                  centroid_y=(y(node1)+y(node2)+y(node3)+y(node4))/4;
                  centroid_z=(z(node1)+z(node2)+z(node3)+z(node4))/4;
                  eig_vector_centroid(e,ii)= ...
               eig_vector_3D{ii}(centroid_x,centroid_y,centroid_z);
         else
                  eig_vector_centroid(e,ii)=0;
         end     
    end   
    random_field(:,ii) = sqrt(eig_value_3D(ii))*eig_vector_centroid(:,ii);   
end 
beta=zeros(Ne,No_of_terms_KL);
for ii=1:No_of_terms_KL 
    K_i(:,:,ii)=zeros(No_of_edges,No_of_edges);
    beta(:,ii)=-k0^2*random_field(:,ii);
end

% Generating the stochastic matrix as in equation (17)
for e=1:Ne 
    e
    l1=sqrt((x(n(e,tet_edge(1,1)))-x(n(e,tet_edge(1,2))))^2+ ...
        (y(n(e,tet_edge(1,1)))-y(n(e,tet_edge(1,2))))^2 +...
        (z(n(e,tet_edge(1,1)))-z(n(e,tet_edge(1,2))))^2);
    l2=sqrt((x(n(e,tet_edge(2,1)))-x(n(e,tet_edge(2,2))))^2+ ...
        (y(n(e,tet_edge(2,1)))-y(n(e,tet_edge(2,2))))^2 +...
        (z(n(e,tet_edge(2,1)))-z(n(e,tet_edge(2,2))))^2);
    l3=sqrt((x(n(e,tet_edge(3,1)))-x(n(e,tet_edge(3,2))))^2+ ...
        (y(n(e,tet_edge(3,1)))-y(n(e,tet_edge(3,2))))^2 +...
        (z(n(e,tet_edge(3,1)))-z(n(e,tet_edge(3,2))))^2);
    l4=sqrt((x(n(e,tet_edge(4,1)))-x(n(e,tet_edge(4,2))))^2+ ...
        (y(n(e,tet_edge(4,1)))-y(n(e,tet_edge(4,2))))^2 +...
        (z(n(e,tet_edge(4,1)))-z(n(e,tet_edge(4,2))))^2);
    l5=sqrt((x(n(e,tet_edge(5,1)))-x(n(e,tet_edge(5,2))))^2+ ...
        (y(n(e,tet_edge(5,1)))-y(n(e,tet_edge(5,2))))^2 +...
        (z(n(e,tet_edge(5,1)))-z(n(e,tet_edge(5,2))))^2);
    l6=sqrt((x(n(e,tet_edge(6,1)))-x(n(e,tet_edge(6,2))))^2+ ...
        (y(n(e,tet_edge(6,1)))-y(n(e,tet_edge(6,2))))^2 +...
        (z(n(e,tet_edge(6,1)))-z(n(e,tet_edge(6,2))))^2);
    le=[l1 l2 l3 l4 l5 l6];
    A=[1 1 1 1; x(n(e,1))  x(n(e,2)) x(n(e,3)) x(n(e,4)); ...
        y(n(e,1)) y(n(e,2)) y(n(e,3)) y(n(e,4)); ...
        z(n(e,1)) z(n(e,2)) z(n(e,3)) z(n(e,4))];
    Ve=1/6*(det(A));
    b_temp=A;
    b_temp(2,:)=[];
    b1=b_temp;b2=b_temp;b3=b_temp;b4=b_temp;
    b1(:,1)=[];
    b2(:,2)=[];
    b3(:,3)=[];
    b4(:,4)=[];
    b=[-det(b1) det(b2) -det(b3) det(b4)];
    c_temp=A;
    c_temp(3,:)=[];
    c1=c_temp;c2=c_temp;c3=c_temp;c4=c_temp;
    c1(:,1)=[];
    c2(:,2)=[];
    c3(:,3)=[];
    c4(:,4)=[];
    c=[det(c1) -det(c2) det(c3) -det(c4)];
    d_temp=A;
    d_temp(4,:)=[];
    d1=d_temp;d2=d_temp;d3=d_temp;d4=d_temp;
    d1(:,1)=[];
    d2(:,2)=[];
    d3(:,3)=[];
    d4(:,4)=[];
    d=[-det(d1) det(d2) -det(d3) det(d4)];
    Feij=zeros(6,6);
    for ii=1:4
        for jj=1:4
            f(ii,jj)=b(ii)*b(jj)+c(ii)*c(jj)+d(ii)*d(jj);
        end
    end
    Fe11=2*(f(2,2)-f(1,2)+f(1,1));
    Fe12=2*f(2,3)-f(2,1)-f(1,3)+f(1,1);
    Fe13=2*f(2,4)-f(2,1)-f(1,4)+f(1,1);
    Fe14=f(2,3)-f(2,2)-2*f(1,3)+f(1,2);
    Fe15=f(2,2)-f(2,4)-f(1,2)+2*f(1,4);
    Fe16=f(2,4)-f(2,3)-f(1,4)+f(1,3);
    Fe21=2*f(3,2)-f(3,1)-f(1,2)+f(1,1);
    Fe22=2*(f(3,3)-f(1,3)+f(1,1));
    Fe23=2*f(3,4)-f(3,1)-f(1,4)+f(1,1);
    Fe24=f(3,3)-f(3,2)-f(1,3)+2*f(1,2);
    Fe25=f(3,2)-f(3,4)-f(1,2)+f(1,4);
    Fe26=f(3,4)-f(3,3)-2*f(1,4)+f(1,3);
    Fe31=2*f(4,2)-f(4,1)-f(1,2)+f(1,1);
    Fe32=2*f(4,3)-f(4,1)-f(1,3)+f(1,1);
    Fe33=2*(f(4,4)-f(1,4)+f(1,1));
    Fe34=f(4,3)-f(4,2)-f(1,3)+f(1,2);
    Fe35=f(4,2)-f(4,4)-2*f(1,2)+f(1,4);
    Fe36=f(4,4)-f(4,3)-f(1,4)+2*f(1,3);
    Fe41=f(3,2)-2*f(3,1)-f(2,2)+f(2,1); 
    Fe42=f(3,3)-f(3,1)-f(2,3)+2*f(2,1); 
    Fe43=f(3,4)-f(3,1)-f(2,4)+f(2,1); 
    Fe44=2*(f(3,3)-f(2,3)+f(2,2)); 
    Fe45=f(3,2)-2*f(3,4)-f(2,2)+f(2,4);
    Fe46=f(3,4)-f(3,3)-2*f(2,4)+f(2,3);
    Fe51=f(2,2)-f(2,1)-f(4,2)+2*f(4,1); 
    Fe52=f(2,3)-f(2,1)-f(4,3)+f(4,1);  
    Fe53=f(2,4)-2*f(2,1)-f(4,4)+f(4,1); 
    Fe54=f(2,3)-f(2,2)-2*f(4,3)+f(4,2); 
    Fe55=2*(f(2,2)-f(2,4)+f(4,4));
    Fe56=f(2,4)-2*f(2,3)-f(4,4)+f(4,3);
    Fe61=f(4,2)-f(4,1)-f(3,2)+f(3,1);
    Fe62=f(4,3)-2*f(4,1)-f(3,3)+f(3,1);
    Fe63=f(4,4)-f(4,1)-f(3,4)+2*f(3,1); 
    Fe64=f(4,3)-2*f(4,2)-f(3,3)+f(3,2); 
    Fe65=f(4,2)-f(4,4)-2*f(3,2)+f(3,4); 
    Fe66=2*(f(4,4)-f(3,4)+f(3,3)); 
    
    Fe=[  Fe11 Fe12 Fe13 Fe14 Fe15 Fe16;
          Fe21 Fe22 Fe23 Fe24 Fe25 Fe26; 
          Fe31 Fe32 Fe33 Fe34 Fe35 Fe36; 
          Fe41 Fe42 Fe43 Fe44 Fe45 Fe46; 
          Fe51 Fe52 Fe53 Fe54 Fe55 Fe56;
          Fe61 Fe62 Fe63 Fe64 Fe65 Fe66;
        ];
     for ii=1:6
        for jj=1:6
            Feij(ii,jj)=((le(ii)*le(jj))/(720*Ve))*Fe(ii,jj);
        end
     end
    for kk=1:No_of_terms_KL
        K_matrix=K_i(:,:,kk);
        Ke_temp_matrix=Feij*beta(e,kk);
        for ii=1:6
          for jj=1:6 
                  K_matrix(abs_edge(e,ii),abs_edge(e,jj))=K_matrix(abs_edge(e,ii),abs_edge(e,jj))+ ...
                       (sign_edge_element(e,ii)*sign_edge_element(e,jj)* Ke_temp_matrix(ii,jj));
          end
        end
        K_i(:,:,kk)=K_matrix;
    end
    clear K_matrix;
end
  
% Applying Boundary condition for K1 and K2  
% Applying PEC boundary conditions 
No_of_edges_final=No_of_edges-length(pec_edge_element_no);
for ii=1:No_of_terms_KL
   K_i_bc(:,:,ii)=zeros(No_of_edges_final,No_of_edges_final); 
end
for ii=1:No_of_terms_KL
    K_matrix=K_i(:,:,ii);
    K_matrix(abs(pec_edge_element_no) ,:)=[]; % Deleting corresponding row 
    K_matrix(:,abs(pec_edge_element_no))=[]; % Deleting corresponding column
    K_i_bc(:,:,ii)=K_matrix;
end
clear K_matrix
clear K_i

% generating SSFEM matrix as in equation (24) of the reference paper
poly_order=1; % change accordingly 
no_of_hermite_poly=factorial(poly_order+No_of_rv)/ ...
    (factorial(poly_order)*factorial(No_of_rv));
Psi_matrix = get_orthogonal_polynomial_regularized_modified(no_of_hermite_poly,N_SSFEM,xi_value);
No_of_edges_final=No_of_edges-length(pec_edge_element_no);
K_stochastic_matrix =sparse(No_of_edges_final*no_of_hermite_poly,No_of_edges_final*no_of_hermite_poly);
f_stochastic_matrix=zeros(No_of_edges_final*no_of_hermite_poly,1);

 
Psi_j_square=zeros(no_of_hermite_poly,no_of_hermite_poly);
for ii=1:no_of_hermite_poly
    for jj=1:no_of_hermite_poly
        Psi_j_square(ii,jj)=get_Psi_j_square(Psi_matrix(:,ii),Psi_matrix(:,jj)); 
    end
end

for kk=1:no_of_hermite_poly
    kk
    for jj=1:no_of_hermite_poly 
        jj
            K_temp=K_mean*Psi_j_square(jj,kk); 
            Psi_ijk_value=zeros(No_of_rv,1);
            for ii=1:No_of_rv
                Psi_ijk_value(ii,1)= ...
                    get_Psi_i_Psi_j_Psi_k(Psi_matrix(:,ii+1),Psi_matrix(:,jj),Psi_matrix(:,kk));
            end 
            for ii=1:No_of_rv
                K_temp=K_temp+ K_i_bc(:,:,ii)*Psi_ijk_value(ii,1);
            end
            K_stochastic_matrix((kk-1)*No_of_edges_final+1:kk*No_of_edges_final,...
            (jj-1)*No_of_edges_final +1:jj*No_of_edges_final)=K_temp;       
    end
end
clear K_i_bc; % to save memory
f_stochastic_matrix(1:No_of_edges_final)=f_mean;
u_stochastic=K_stochastic_matrix\f_stochastic_matrix;
% Adding the unknown quantity obtained with the known quantity
u_stochastic_final=zeros(No_of_edges*no_of_hermite_poly,1);
for kk=1:no_of_hermite_poly 
    kk
    count=1;
    edge_count=1;
    u_stochastic_temp=zeros(No_of_edges,1);
    u_final=u_stochastic((kk-1)*No_of_edges_final+1:kk*No_of_edges_final);
    for ii=1:No_of_edges
        if edge_count<=length(pec_edge_element_no)
            if(ii==abs(pec_edge_element_no(edge_count)))
                u_stochastic_temp(ii)=0; 
                edge_count=edge_count+1;
            else
                u_stochastic_temp(ii)=u_final(count);
                count=count+1; 
            end 
        else
            u_stochastic_temp(ii)=u_final(count);
            count=count+1;        
        end
    end
    u_stochastic_final((kk-1)*No_of_edges +1 : kk*No_of_edges)=u_stochastic_temp;
end


bc_ABC_left_edge_no=unique(abs(edge_element_no_input)); 
bc_ABC_right_edge_no=unique(abs(edge_element_no_output)); 
left_Ez_stochastic=zeros(no_of_hermite_poly,1);
right_Ez_stochastic=zeros(no_of_hermite_poly,1);
for ii=1:no_of_hermite_poly
    Ey_final_final=u_stochastic_final((ii-1)*No_of_edges+1:ii*No_of_edges); 
    gamma_port1_temp=0;
    tau_port2_temp=0; 
    for e=1:Ne
        if ismember(e,input_face_surface_element_no) 
            [match_status face_no]=ismember(e,input_face_surface_element_no);
            face_surface_temp=input_face_surface(face_no,:);
            tetrahedron_element=four_node_tetrahedron(e,:);
            sign_edge_element_tetrahedron=sign_edge_element(e,:);  
            normal_vector=input_normal_vector(face_no,:); 
            edge_element_no_temp=abs(edge_element_no_input(face_no,:));
            gamma_port1=evaluate_surface_integral_gamma(face_surface_temp,input_normal_vector,node_list,sign_edge_element_tetrahedron,tetrahedron_element,width,edge_element_no_temp,Ey_final_final);
            gamma_port1=(gamma_port1*(2*exp(-1i*k_z_10*input_len))/(width*height*E0));
            gamma_port1_temp= gamma_port1_temp+gamma_port1;
        end
        if ismember(e,output_face_surface_element_no)
            [match_status face_no]=ismember(e,output_face_surface_element_no);
            face_surface_temp=output_face_surface(face_no,:);
            tetrahedron_element=four_node_tetrahedron(e,:);
            sign_edge_element_tetrahedron=sign_edge_element(e,:);  
            normal_vector=output_normal_vector(face_no,:); 
            edge_element_no_temp=abs(edge_element_no_output(face_no,:));
            tau_port2=evaluate_surface_integral_gamma(face_surface_temp,input_normal_vector,node_list,sign_edge_element_tetrahedron,tetrahedron_element,width,edge_element_no_temp,Ey_final_final);
            tau_port2=(tau_port2*(2*exp( 1i*k_z_10*output_len))/(width*height*E0));
            tau_port2_temp=tau_port2_temp+tau_port2;
        end
    end
    left_Ez_stochastic(ii)=gamma_port1_temp;
    right_Ez_stochastic(ii)=tau_port2_temp;
end
u_SSFEM_ref=zeros(N_SSFEM,1);
u_SSFEM_trans=zeros(N_SSFEM,1);
for jj=1:no_of_hermite_poly 
    u_SSFEM_ref=u_SSFEM_ref+left_Ez_stochastic(jj)*Psi_matrix(:,jj);
    u_SSFEM_trans=u_SSFEM_trans+right_Ez_stochastic(jj)*Psi_matrix(:,jj);
end 
gamma_fem_SSFEM=abs(u_SSFEM_ref-exp(-2*1i*k_z_10*input_len));
trans_fem_SSFEM=abs(u_SSFEM_trans);  
time_elapsed_SSFEM=toc;
if poly_order==2
    save('SSFEM_3D_multilayer_dielectric_block_10GHz_poly_order_2_IJRFMWCAD.mat');
end
if poly_order==1
    save('SSFEM_3D_multilayer_dielectric_block_10GHz_poly_order_1_IJRFMWCAD.mat');
end

% Plots
clear all;
close all;
poly_order=1;
if poly_order==2
    load('SSFEM_3D_multilayer_dielectric_block_10GHz_poly_order_2_IJRFMWCAD.mat');
end
if poly_order==1
    load('SSFEM_3D_multilayer_dielectric_block_10GHz_poly_order_1_IJRFMWCAD.mat');
end
[prob_dist_vector_fem_gamma_SSFEM, set_of_points_fem_gamma_SSFEM ]=ksdensity(gamma_fem_SSFEM);  
[prob_dist_vector_fem_trans_SSFEM ,set_of_points_fem_trans_SSFEM ]=ksdensity(trans_fem_SSFEM); 
figure(1)
colorstring = 'krbgmwy';
d_1(1)=plot(set_of_points_fem_gamma_SSFEM ,prob_dist_vector_fem_gamma_SSFEM,'LineWidth',1.5,'Color',colorstring(3));
xlabel('Magnitude of Reflection Coefficient  , $|\Gamma| $ ','Interpreter','latex');
ylabel('Probability density function  ','Interpreter','latex'); 
hold on;
figure(2)
colorstring = 'krbgmwy';
t_1(1)=plot(set_of_points_fem_trans_SSFEM,prob_dist_vector_fem_trans_SSFEM,'LineWidth',1.5,'Color',colorstring(3));
xlabel('Magnitude of Tranmission Coefficient  , $|T| $ ','Interpreter','latex');
ylabel('Probability density function  ','Interpreter','latex'); 
hold on;

% MC
% MCFEM_dielectric_loaded_rect_wgde_3D_4_IJRFMWCAD.m
load('MC_rect_2cm_1cm_24cm_for_IJRFMWCAD.mat');  
gamma_abs_MC=abs(gamma_port1_final);
tau_abs_MC=abs(tau_port2_final);
[prob_dist_vector_fem_MC, set_of_points_fem_MC ]=ksdensity(gamma_abs_MC);
figure(1)  
d_1(2)=plot(set_of_points_fem_MC,prob_dist_vector_fem_MC,'--','LineWidth',1.5,'Color',colorstring(2));
hold off;
legend([d_1(1) d_1(2) ],'SSFEM',' MCFEM - 10000','Interpreter','latex');
[prob_dist_vector_fem_trans_coeff_MC, set_of_points_fem_trans_coeff_MC ]=ksdensity(tau_abs_MC);
figure(2) 
t_1(2)=plot(set_of_points_fem_trans_coeff_MC,prob_dist_vector_fem_trans_coeff_MC,'--', 'LineWidth',1.5,'Color',colorstring(2));
hold off;
legend([t_1(1) t_1(2) ],'SSFEM',' MCFEM - 10000','Interpreter','latex');

figure(3)
subplot(1,2,1)
spy(K_mean)
title('FEM matrix')
subplot(1,2,2)
spy(K_stochastic_matrix);
title('SSFEM matrix')
hold on;
colorstring = 'krbgmwy';
for ii=1:no_of_hermite_poly
    x_cord=[0 No_of_edges_final*no_of_hermite_poly];
    y_cord=[No_of_edges_final*ii No_of_edges_final*ii];
    line(x_cord,y_cord,'LineWidth',1,'Color',colorstring(1));
    hold on;
end
for ii=1:no_of_hermite_poly
    y_cord=[0 No_of_edges_final*no_of_hermite_poly];
    x_cord=[No_of_edges_final*ii No_of_edges_final*ii];
    line(x_cord,y_cord,'LineWidth',1,'Color',colorstring(1));
    hold on;
end
