function mesh = cubeMesh(p,center,alpha,depth)
    nPoint=6*(p-1)^2+2
    
    mesh=[-1,-1,-1];

    for i=0:(p-1)-1
        for j=0:(p-1)-1
            mesh=[mesh;[-1,(2.0 * (i + 1) - p + 1) / (p - 1),(2.0 * j - p + 1) / (p - 1)]];
        end
    end
    
    for i=0:(p-1)-1
        for j=0:(p-1)-1
            mesh=[mesh;[ (2.0 * i - p + 1) / (p - 1),-1,(2.0 * (j + 1) - p + 1) / (p - 1)]];
        end
    end
    
    for i=0:(p-1)-1
        for j=0:(p-1)-1
            mesh=[mesh;[(2.0 * (i + 1) - p + 1) / (p - 1),(2.0 * j - p + 1) / (p - 1),-1]];
        end
    end
    
    for i=0:((nPoint/2))-1
        mesh=[mesh;-mesh(i+1,:)];
    end
            
    r = 0.5*power(0.5,depth);
    b=alpha*r;
    for i=0:(nPoint-1)
        mesh(i+1,:) = (mesh(i+1,:)+1)*b+center;
    end
end


% 
%     nPoint=6*(p-1)**2+2
%     mesh=[]
%     mesh.append(np.array([-1,-1,-1]))
%     for i in range(p-1):
%         for j in range(p-1):
%             mesh.append(np.array([-1,(2.0 * (i + 1) - p + 1) / (p - 1),(2.0 * j - p + 1) / (p - 1)]))
%     for i in range(p-1):
%         for j in range(p-1):
%             mesh.append(np.array([ (2.0 * i - p + 1) / (p - 1),-1,(2.0 * (j + 1) - p + 1) / (p - 1)]))
%     for i in range(p-1):
%         for j in range(p-1):
%             mesh.append(np.array([(2.0 * (i + 1) - p + 1) / (p - 1),(2.0 * j - p + 1) / (p - 1),-1]))
%     for i in range(int(nPoint/2)):
%         mesh.append(-mesh[i])
%     r = 0.5*pow(0.5,depth)
%     b=alpha*r
%     for i in range(nPoint):
%         mesh[i] = (mesh[i]+1)*b+center