%% HW 4. Problem 3

%% (A)
load('granule.mat')
Ri = 100;
Rm = 10000;
Cm = 1;
Er = 0;

indexcompartment        = granule(:,1);
typecompartment         = granule(:,2);
parentcompartmentmatrix = granule(:,7);

somaindex = indexcompartment(typecompartment==1);

%Take the subbranch compartment out and connect them to their parent branch
maincablelogic = false(size(granule,1),1);
for i = 2:size(granule)
    maincablelogic(i) = parentcompartmentmatrix(i)+1 ~= indexcompartment(i);
end

subbranchindex=indexcompartment(maincablelogic);
subbranchindex=[subbranchindex;size(granule,1)];

%The subbranch info in cell 
cellsize = size(subbranchindex)-1;
subbranch = cell(cellsize);
ss = size(subbranchindex,1);
for i = 1:ss
    if i ~= ss
    first = subbranchindex(i);
    last  = subbranchindex(i+1)-1;
    parent = parentcompartmentmatrix(first);
    subbranch{i} = granule([parent,first:last],:);
    end
end

subbranch = subbranch';

% The main branch's info
removesublogic =true(52,1);
removesublogic(subbranchindex(1):subbranchindex(end),1)=false;
maincable = granule(removesublogic,:);

wholecell = cell(ss,1);
wholecell{1}=maincable;
wholecell(2:size(subbranch,1)+1)=subbranch;

figure; grid on; hold on
set(gca, 'color', [0 0 0])
xlabel('x coordinate'); ylabel('y coordinate');zlabel('z coordinate')
title('The Dendrite Tree in 3D')
for i = 1:ss
    plot3(wholecell{i}(:,3),wholecell{i}(:,4),wholecell{i}(:,5),'g')   
    if any(wholecell{i}(:,2))==1
    somalogic = wholecell{i}(:,2)==1;
    plot3(wholecell{i}(somalogic,3),wholecell{i}(somalogic,4),wholecell{i}(somalogic,5),'r*')   
    end    
end

%% (B)
l = NaN(52,1);
l(1) = 0.0005;
for i=2:size(granule,1)
    f=granule(i,7);
    l(i,1)=(((granule(i,3)-granule(f,3))^2)+((granule(i,4)-granule(f,4))^2)...
        +((granule(i,5)-granule(f,5))^2))^.5;
end
sc = sqrt(granule(:,6)*Rm/(Ri*2));
satisfy = l<=0.1*sc; 
if satisfy == true(52,1)
disp('All compartments satisfy the constaint.');
end

%% (C)
matrix_radius_compartment = granule(:,6);


gi = pi .* matrix_radius_compartment.^2 ./ l ./ Ri;
cj = 2*pi.*matrix_radius_compartment.* l .*Cm;
gjm = 2*pi.*matrix_radius_compartment.* l ./ Rm;

Iapp = 1e-9;

n=size(granule,1);
A = zeros(n,n);
B = A;
v = zeros(n,1);
u = v;
u(35,1)=Iapp;

prebranch = granule(subbranchindex,7);
prebranch(end)=[];
prebranch = sort(prebranch);

endbranch = NaN(ss,1);
endbranch(1:end-1,1) = subbranchindex(1:end-1,1)-1;
endbranch(end)=granule(end,1);

daughterbranch = NaN(ss-1,3) ;
daughterbranch(:,1)=prebranch(:);
daughterbranch(:,2)=prebranch(:)+1;
for i = 2: size(wholecell,1)
    subparent = wholecell{i,1}(1,1);
    [~,loc]=ismember(subparent,daughterbranch(:,1));
    daughterbranch(loc,3)= wholecell{i,1}(2,1);
    
end
daughterbranch = daughterbranch(:,2:3);

tf = ismember(1,prebranch);
tf = double(tf);

h = 0;
for i = 1:n
    B(i,i)=(1/cj(i));
    f=granule(i,7);
    if i == tf
       h = h+1;
       A(i,i)=-(gjm(i)+gi(i+1)+gi(daughterbranch(h,2)))/cj(i);
       A(i,i+1)=gi(i+1)/cj(i);
       A(i,daughterbranch(h,2))=gi(daughterbranch(h,2))/cj(i);
    elseif i == 1 %Initial mode
        A(i,i)=-(gjm(i)+gi(i+1))/cj(i);
        A(i,i+1)=gi(i+1)/cj(i);
    elseif ismember(i,prebranch) %Prebranch mode
        h=h+1;
        A(i,f)= gi(i)/cj(i);
        A(i,i) =-(gi(i)+gjm(i)+gi(i+1)+gi(daughterbranch(h,2)))/cj(i);
        A(i,i+1)= gi(i+1)/cj(i);
        A(i,daughterbranch(h,2)) = gi(daughterbranch(h,2))/cj(i);
    elseif ismember(i,endbranch)% End of the subbranch
        A(i,f)= gi(i)/cj(i);
        A(i,i) =-(gi(i)+gjm(i))/cj(i);
    else 
        A(i,f)= gi(i)/cj(i);
        A(i,i)=-(gi(i)+gjm(i)+gi(i+1))/cj(i);
        A(i,i+1)= gi(i+1)/cj(i);    
    end
end

fulllength =zeros(length(endbranch),10000);
for i =1:length(endbranch)
    h =1;
    n = granule(endbranch(i),7);
    fulllength(i,h)=endbranch(i);
    while n>=0
        h=h+1;
        fulllength(i,h)=n;
        n=granule(n,7);
    end
end

v = -inv(A)*B*u;

fulllength1=fulllength(1,:);
fulllength1=fulllength1(fulllength1~=0);
fulllength1=flip(fulllength1);

l1 = zeros(length(fulllength1),1);
v1 = zeros(length(fulllength1),1);
for i = 1:length(fulllength1)
    v1(i,1)=v(fulllength1(i));
    for j=i:length(fulllength1)
       l1(j,1)=l1(j,1)+l(fulllength1(i));
    end
end


fulllength2=fulllength(2,:);
fulllength2=fulllength2(fulllength2~=0);
fulllength2=flip(fulllength2);

l2 = zeros(length(fulllength2),1);
v2 = zeros(length(fulllength2),1);
for i = 1:length(fulllength2)
    v2(i,1)=v(fulllength2(i));
    for j=i:length(fulllength2)
       l2(j,1)=l2(j,1)+l(fulllength2(i));
    end
end


fulllength3=fulllength(3,:);
fulllength3=fulllength3(fulllength3~=0);
fulllength3=flip(fulllength3);

l3 = zeros(length(fulllength3),1);
v3 = zeros(length(fulllength3),1);
for i = 1:length(fulllength3)
    v3(i,1)=v(fulllength3(i));
    for j=i:length(fulllength3)
       l3(j,1)=l3(j,1)+l(fulllength3(i));
    end
end

fulllength4=fulllength(4,:);
fulllength4=fulllength4(fulllength4~=0);
fulllength4=flip(fulllength4);

l4 = zeros(length(fulllength4),1);
v4 = zeros(length(fulllength4),1);
for i = 1:length(fulllength4)
    v4(i,1)=v(fulllength4(i));
    for j=i:length(fulllength4)
       l4(j,1)=l4(j,1)+l(fulllength4(i));
    end
end

fulllength5=fulllength(5,:);
fulllength5=fulllength5(fulllength5~=0);
fulllength5=flip(fulllength5);

l5 = zeros(length(fulllength5),1);
v5 = zeros(length(fulllength5),1);
for i = 1:length(fulllength5)
    v5(i,1)=v(fulllength5(i));
    for j=i:length(fulllength5)
       l5(j,1)=l5(j,1)+l(fulllength5(i));
    end
end

fulllength6=fulllength(6,:);
fulllength6=fulllength6(fulllength6~=0);
fulllength6=flip(fulllength6);

l6 = zeros(length(fulllength6),1);
v6 = zeros(length(fulllength6),1);
for i = 1:length(fulllength6)
    v6(i,1)=v(fulllength6(i));
    for j=i:length(fulllength6)
       l6(j,1)=l6(j,1)+l(fulllength6(i));
    end
end

figure; clf; hold on
plot(l1,v1,'g');
plot(l2,v2,'y');
plot(l3,v3,'p');
plot(l4,v4,'c');
plot(l5,v5,'r');
plot(l6,v6,'b');
xlabel('Dimensionless Length'); ylabel('V [mV]');
title('Steady State voltage across all branches vs. dimensionless distance from the soma');


%% (D)

dvdt =@(V) A*V  + B*u;
t_span =linspace(0,5e4,10000);
V0 = zeros(1,52);
tc = Rm * Cm; 
[t,V]=ode23(@(t,V) dvdt(V),t_span,V0);

figure; hold on
for i = 1:52
plot(t/tc,V(:,i))
end
xlabel('T'); ylabel('V(X) [mV]');
title('Voltage V(T) over time')

