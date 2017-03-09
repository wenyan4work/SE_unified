p=16
nameIn='M2L3Dp'+str(p)
nameOut='M2L3D3Dp'+str(p)

with open(nameIn) as fin:
    with open(nameOut,'w') as fout:
        row=0
        col=0 
        for line in fin:
            fout.write(str(row)+' '+str(col)+' '+line)
            col=col+1
            if(col==3*(6*(p-1)**2+2) ):
                col=0
                row=row+1


