% 
% [xc,wc]=clencurt(40);
% 
% xc=(xc+1)*0.5;
% wc=wc*0.5;

% uniform trapezoidal 

np=40;
xc=[0:1/np:1];
wc=ones(size(xc))*(1/np);
wc(1)=wc(1)*0.5;
wc(np+1)=wc(np+1)*0.5;


N=numel(xc);
xpos=zeros(N*N*N,3);
force=zeros(N*N*N,3);


force(1,:)=[N,0,0];

for i=1:N
    for j=1:N
        for k=1:N
            xpos((i-1)*N*N+(j-1)*N+k,:)=[xc(i),xc(j),xc(k)];
        end
    end
end
evalid=[1:1:N*N*N];
SE_opt.M=[128,128,128];
u=stokeslet_total_ext(xpos,[0.3123,0.4416,0.5219],[1,0,0],15,SE_opt);

unet=[0,0,0];

for i=1:N
    for j=1:N
        for k=1:N
            unet=unet+wc(i)*wc(j)*wc(k)*u((i-1)*N*N+(j-1)*N+k,:);
        end
    end
end

unet