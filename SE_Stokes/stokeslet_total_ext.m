
function [u,ufft,ureal] = stokeslet_total_ext( xext,x,f,xi,opt )
%STOKESLET_TOTAL Summary of this function goes here
%   Detailed explanation goes here

xextorigin=xext;
xorigin=x;

[T,~]=size(xext);
for t=1:T % periodicity of target
    xext(t,:)=xext(t,:)-floor(xext(t,:)./[opt.box(1),opt.box(2),opt.box(3)]).*[opt.box(1),opt.box(2),opt.box(3)];
end

% periodicity of source
[M,~]=size(x);
for t=1:M
    x(t,:)=x(t,:)-floor(x(t,:)./[opt.box(1),opt.box(2),opt.box(3)]).*[opt.box(1),opt.box(2),opt.box(3)];
end

% suggest for cubic box: xi=15
% opt.M=[128,128,128]
% opt.box=[1,1,1]
% opt.P=25
ufft=SE_Stokes_ext(xext,x,f,xi,opt);
[numext,~]=size(xext);
ureal=zeros(numext,3);

for i=1:numext
    target=xext(i,:);
    ureal(i,:)=[0,0,0]; % no self term for external points
    % find neighbor
    [M,~]=size(x);
    for j=1:M
        source=x(M,:);
        rvec=(target-source);
        if (norm(rvec)<1e-13)
            continue;
        end
        for a=-1:1
            for b=-1:1
                for c=-1:1
                    rvec=target-source-[a*opt.box(1),b*opt.box(2),c*opt.box(3)];
                    ureal(i,:)=ureal(i,:)+f(j,:)*(AEW(xi,rvec));
                end
            end
        end
        
    end
end 

u=ufft+ureal;

% unb=zeros(size(u));
% [numext,~]=size(xext);
% parfor i=1:numext
%    target = xextorigin(i,:);
%    unb(i,:)=[0,0,0];
%    [M,~]=size(xorigin);
%    for j=1:M
%        source=xorigin(j,:);
%        rvec=(target-source);
%        if(norm(rvec)<1e-13)
%            continue
%        end
%        for a=-2:2
%             for b=-2:2
%                 for c=-2:2
%                     rvec=target-source-[a*opt.box(1),b*opt.box(2),c*opt.box(3)];
%                     unb(i,:)=unb(i,:)+f(j,:)*(GStokes(rvec));
%                 end
%             end
%         end
%    end
%     
% end
% 
% u=u-unb;

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