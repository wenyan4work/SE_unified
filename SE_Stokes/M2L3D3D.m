scaleDownCheck=1.05;
scaleDownEquiv=2.95;
p=16;
filename='M2L3Dp16';
meshDownCheck=cubeMesh(p,[-(scaleDownCheck-1)/2,-(scaleDownCheck-1)/2,-(scaleDownCheck-1)/2],scaleDownCheck,0);
meshDownEquiv=cubeMesh(p,[-(scaleDownEquiv-1)/2,-(scaleDownEquiv-1)/2,-(scaleDownEquiv-1)/2],scaleDownEquiv,0);

Mmesh=6*(p-1)^2+2;
A=zeros(Mmesh*3);
M2L=zeros(Mmesh*3);

for i=1:Mmesh
    for j=1:Mmesh
        target=meshDownCheck(i,:);
        source=meshDownEquiv(j,:);
        rvec=target-source;
        rnorm=norm(rvec);
        A(3*i-3+1:3*i,3*j-3+1:3*j)=eye(3)/rnorm + transpose(rvec)*rvec/rnorm^3;
    end
end

[U,S,V] = svd(A);
[ms,ns]=size(S);
for ii=1:min([ms,ns])
    if(S(ii,ii) > 1e-15*S(1,1))
        S(ii,ii)=1/S(ii,ii);
    else
        S(ii,ii)=0;
    end
    
end
b=ones(Mmesh*3,1);
x=(V)*((S*U')*b);
max(A*x-b) % check backward stable

SE_opt.M=[64,64,64];

for i=1:Mmesh
    bright=zeros(3*Mmesh,3);
    i
    xpos=meshDownCheck;
    xext=meshDownCheck;
    xext(i,:)=[0,0,0]; % not overlap with j
    % fx
    [ucheck,~]=stokeslet_total_ext(xext,xpos(i,:),[1,0,0],15,SE_opt); % other
    bright(:,1)=reshape(ucheck',[3*Mmesh,1]);
    uself = stokeslet_total(1,xpos(i,:),[1,0,0],15,SE_opt); % self
    bright(3*i-3+1:3*i,1)=uself';
        % fy
    [ucheck,~]=stokeslet_total_ext(xext,xpos(i,:),[0,1,0],15,SE_opt); % other
    bright(:,2)=reshape(ucheck',[3*Mmesh,1]);
    uself = stokeslet_total(1,xpos(i,:),[0,1,0],15,SE_opt); % self
    bright(3*i-3+1:3*i,2)=uself';
        % fz
    [ucheck,~]=stokeslet_total_ext(xext,xpos(i,:),[0,0,1],15,SE_opt); % other
    bright(:,3)=reshape(ucheck',[3*Mmesh,1]);
    uself = stokeslet_total(1,xpos(i,:),[0,0,1],15,SE_opt); % self
    bright(3*i-3+1:3*i,3)=uself';


    M2L(:,3*i-3+1:3*i-3+3)= (V)*((S*U')*bright);
end

dlmwrite(filename, M2L, 'delimiter', '\n', 'precision', '%.18e');