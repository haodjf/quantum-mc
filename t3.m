close all; clear; clc;
v_z=290;       %速度
d3=0.0003;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:2000000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7303')

close all; clear; clc;
v_z=290;       %速度
d3=0.0005;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:1500000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7305')

close all; clear; clc;
v_z=290;       %速度
d3=0.0008;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:1000000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7308')

close all; clear; clc;
v_z=290;       %速度
d3=0.001;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:1000000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7310')

close all; clear; clc;
v_z=290;       %速度
d3=0.0013;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:1000000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7313')

close all; clear; clc;
v_z=290;       %速度
d3=0.0015;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:1000000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7315')

close all; clear; clc;
v_z=290;       %速度
d3=0.0018;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:1000000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7318')

close all; clear; clc;
v_z=290;       %速度
d3=0.002;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:1000000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7320')

close all; clear; clc;
v_z=290;       %速度
d3=0.0023;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:1000000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7323')

close all; clear; clc;
v_z=290;       %速度
d3=0.0025;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:1000000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7325')

close all; clear; clc;
v_z=290;       %速度
d3=0.0028;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:1000000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7328')

close all; clear; clc;
v_z=290;       %速度
d3=0.003;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:1000000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7330')

close all; clear; clc;
v_z=290;       %速度
d3=0.0035;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:500000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7335')

close all; clear; clc;
v_z=290;       %速度
d3=0.004;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:500000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7340')

close all; clear; clc;
v_z=290;       %速度
d3=0.0045;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:500000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7345')

close all; clear; clc;
v_z=290;       %速度
d3=0.005;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:500000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7350')

close all; clear; clc;
v_z=290;       %速度
d3=0.006;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:500000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7360')

close all; clear; clc;
v_z=290;       %速度
d3=0.007;     %slit3宽度
hbar=1.0546*10^(-34);
tau=97.88*10^(-9);%寿命
gamma=1/tau; %w
c_vac=2.99792458*10^8;
nu_0=276736.600*10^9;
lambda_0=c_vac/nu_0;
k=2*pi/lambda_0;
m=4.0026*1.6606*10^(-27); %mass of helium-4
%m=4*1.672*10^(-27);    %mass of helium-4
epsilon_0=8.854*10^(-12);
vrec=hbar*k/m;
p=10;
i1=-p;
i2=p;
num1=i2-i1+1;      %态的数目
w0=1*10^(-3);      %束腰直径
cg=1/sqrt(3);%Clebsch-Gordan系数，2^3S_1->2^3P_0
spec_x=[-1:0.1:1];
dt=1/10*gamma^(-1); %as SMALL  as possible!，这里蒙卡时间间隔为十分之一自然线宽寿命
det_num=length(spec_x); %number of detunings                      
d1=0.0003;
d2=0.0003;               %狭缝宽度
diagm0=zeros(num1,num1);

%k1=-k
matrix_jump_m1=[diagm0,diagm0;...
    diag(ones(1,num1-1),-1),diagm0];

matrix_jump_m1=sparse(matrix_jump_m1);

%k1=+k
matrix_jump_p1=[diagm0,diagm0;...
    diag(ones(1,num1-1),+1),diagm0];

matrix_jump_p1=sparse(matrix_jump_p1);

%k1=0k
matrix_jump_0=[diagm0,diagm0;...
    diag(ones(1,num1)),diagm0];

matrix_jump_0=sparse(matrix_jump_0);         %不同自发辐射的跳跃矩阵
w=(-1+sqrt(3)*i)/2; %三次函数根的系数


time=w0/v_z; %acutal_time = 46.7*gamma^(-1)
step=floor(time/dt); %total time，天花板函数得到蒙卡次数
t1=0.48/v_z;            %slit1到slit2
t2=0.67/v_z;            %slit2到probe
t3=1.53/v_z;            %probe之后
asym=0;                 %频率不对称性
alpha=-0.0000;    

c=zeros(det_num*6,1);  %探测到的原子数目
spec_i=[0.15 0.3 0.45 0.6 0.75 0.9]; %功率
l=0.0000; %slit3中心位置

for l_det=1:1:6     %选择功率
i_isat=spec_i(l_det); %饱和度
omega=(1/sqrt(2))*cg*gamma*sqrt(i_isat/2);%拉比频率
for j_det=1:det_num  %再选择激光的失谐
count=0; %原子计数
parfor n=1:1:500000
y0=0;  %slit2处初始位置
v1=normrnd(0,0.8*v_z/2190); %slit2处速度分布
x_initial=y0+v1*t2;            %原子在probe处的横向位置 
  
        detuning=spec_x(j_det);   % (-922*v1/gamma); 不同的detuning

                no_jump_step=90; %为简化计算，90次更新一下作用矩阵
                no_jump_num=floor(step/no_jump_step);%step=floor(time/dt),time是渡越时间，dt是相邻两次模拟的间隔
               % no_jump_list=zeros(2*num1, 2*num1*no_jump_num);
                jump_num=0; %自发辐射数目
                jump_num1=0; %自发辐射数目，发生跳跃后更新作用矩阵
                vec0=zeros(2*num1,1);  
                vec0((3*num1+1)/2,1)=1;  %原子在不同态的概率分布
                vec2=zeros(num1,1);      %自发辐射后的概率分布
               j02_stop=step;            %蒙卡次数
                    j_nj=1;

                   random1=rand(step,1);%直接生成time/dt这么多个0~1上的随机数...
                    random2=rand(step,1);
                    random31=(8/3)*rand(step,1)-4/3;
                    random32=random31/2;
                    random4=2*pi*rand(step,1);
                    for j02=1:j02_stop         %蒙卡次数
                        
                        ep1=random1(j02);
                        ep2=random2(j02);       
                        ep4=random4(j02);

if mod(j02,no_jump_step)==1%简化Nojump计算
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %作用矩阵

 j_nj=j_nj+1;

end

if jump_num1<jump_num

jump_num1=jump_num1+1;
  
  delta=detuning*gamma;
  
  block_diag=diag(ones(1,num1));
  
  block_zeros=zeros(num1,num1);
  coe_relax=-hbar*(1i)*(gamma/2);
  
  H_relax= [block_diag*coe_relax, block_zeros;
                     block_zeros, block_zeros];

  H_unit=diag(ones(1,2*num1));
  
  block_coe_g=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m));%+ 1/2*hbar*delta);
  
  block_coe_e=diag((hbar*k*(i1:i2)+m*(v1+alpha*v_z)).^2/(2*m)- hbar*delta);

  coe_off_m1=hbar*omega/2*exp(-(j02*dt-time/2)^2/((time/2)^2));
                                     %sigma = time/2
  coe_off_p1=coe_off_m1;

  block_off=coe_off_m1*diag(ones(1,num1-1),+1); % m1, -hk, p+1
  block_offf=coe_off_p1*diag(ones(1,num1-1),-1);    % p1, +hk, p-1

  H_s= [block_coe_e, block_offf;
        block_off, block_coe_g];

  H_tot=H_s + H_relax;

  matrix_tot_1=H_unit - (1i)*(dt/hbar)*H_tot;

  %matrix_tot_2=H_unit - (1i)*(dt/hbar)*H_tot - 1/2*(dt/hbar)^2*H_tot*H_tot;

 nojump=sparse(matrix_tot_1);  %更新作用矩阵



end
                        
                        
                        
                        dp=gamma*dt*norm(vec0(1:num1))^2; %jump probability
 if ep1>dp %no jump

                           vec_evo=nojump*vec0; %evolve with H_tot 
                            %simply normalzation
                            vec1=vec_evo/norm(vec_evo); %归一化
                            
                        else %quantum jump
                           if ep2>2/3 
                            break
                           elseif ep2<1/3 %jump to |g>
                                ep31=acos(nthroot((3*random31(j02)/2+sqrt(9*random31(j02)*random31(j02)/4+1)),3)+nthroot((3*random31(j02)/2-sqrt(9*random31(j02)*random31(j02)/4+1)),3));
                                vec1=matrix_jump_0*vec0;
                                vec1=vec1/norm(vec1);
                                v1=v1+vrec*cos(ep4)*sin(ep31); %自发辐射导致横向速度改变
                                jump_num=jump_num+1;
                               
                                
                            else %jump to |0>
                               if (random32(j02)<0)
                        ep32=acos(-norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        else
                        ep32=acos(norm(w^2*(-3*random32(j02)/2+sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)+w*(-3*random32(j02)/2-sqrt(9*random32(j02)*random32(j02)/4-1))^(1/3)));
                        end
                                vec1=matrix_jump_0*vec0;
                                vec2=vec1(1+num1:2*num1);
                                vec2=vec2/norm(vec2);
                                v1=v1+vrec*cos(ep4)*sin(ep32); %自发辐射导致横向速度改变
                                
                             break
                                
                                
                            end
                            
                        end
    vec0=vec1;
                        
                    end %end of j02, step

probability=real(abs(vec2.'.^2)); %跳跃到不同态的概率
y=zeros(num1,1);
for j=1:num1

y(j,1)=x_initial+(v1+(j-(num1+1)/2)*vrec)*t3; %跳跃到不同态对应的slit3位置

if(y(j,1)<l+d3/2&&y(j,1)>l-d3/2)
count=count+probability(1,j);   %探测到的计数
end


end




end

c(j_det+det_num*(l_det-1),1)=count;  %不同功率，不同失谐对应的计数


end
end
save('bgl7370')
