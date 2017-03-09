function [u,ufft,ureal] = stokeslet_total( eval_idx,x,f,xi,opt )
%STOKESLET_TOTAL Summary of this function goes here
%   Detailed explanation goes here

xorigin=x;

% periodicity of source
[M,~]=size(x);
for t=1:M
    x(t,:)=x(t,:)-floor(x(t,:)./[opt.box(1),opt.box(2),opt.box(3)]).*[opt.box(1),opt.box(2),opt.box(3)];
end

% suggest for cubic box: xi=15
% opt.M=[128,128,128]
% opt.box=[1,1,1]
% opt.P=25
ufft=SE_Stokes(eval_idx,x,f,xi,opt);

ureal=zeros(numel(eval_idx),3);


% build a extended point cloud
xext=x;
for a=-1:1
    for b=-1:1
        for c=-1:1
            if(a==0)&&(b==0)&&(c==0)
                continue
            end
            xext=[xext;x+[a*opt.box(1),b*opt.box(2),c*opt.box(3)]];
        end
    end
end

pcloud=pointCloud(xext);


parfor i=1:numel(eval_idx)
    idx=eval_idx(i);
    target=x(idx,:);
    ureal(i,:)=-f(idx,:)*4*xi/sqrt(pi);

    
    % find neighbor
    nbId=findNeighborsInRadius(pcloud,target,0.5);
    [M,~]=size(nbId);
    for j=1:M
        source=xext(nbId(j),:);
        rvec=(target-source);
        if (norm(rvec)<1e-13)
            continue;
        end
        ureal(i,:)=ureal(i,:)+f(idx,:)*AEW(xi,rvec);
    end
end 

u=ufft+ureal;


unb=zeros(size(u));

[M,~]=size(xorigin);

parfor i=1:M
   target = xorigin(i,:);
   unb(i,:)=[0,0,0];
   for j=1:M
       source=xorigin(j,:);
       rvec=(target-source);

       for a=-2:2
            for b=-2:2
                for c=-2:2
                    rvec=target-source-[a*opt.box(1),b*opt.box(2),c*opt.box(3)];
                    if(norm(rvec)<1e-13)
                        continue
                    end
                    unb(i,:)=unb(i,:)+f(j,:)*(GStokes(rvec));
                end
            end
        end
   end
    
end

u=u-unb;

end

function A = AEW(xi,rvec)
% Hasimoto
    r=norm(rvec);
    A = 2*(xi*exp(-(xi^2)*(r^2))/(sqrt(pi)*r^2)+erfc(xi*r)/(2*r^3))*(r*r*eye(3)+transpose(rvec)*rvec) - 4*xi/sqrt(pi)*exp(-(xi^2)*(r^2))*eye(3) ;
end

function G=GStokes(rvec)
    r=norm(rvec);
    G=eye(3)/r + rvec'*rvec/r^3;
end
