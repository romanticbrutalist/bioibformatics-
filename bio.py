#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 15:03:37 2019

@author: Kaan Hamilton
NOTE!!!
The program assumes that the short reads text file is in the 
same directory with no initial fasta declaration, the same goes for 
the reference gene. The short reads and reference gene are in the 
same directory as the source code saved as sr.txt and rs.txt
respectively
"""
import numpy as np
import random
import matplotlib.pyplot as plt
def match_score(a,b):
    if(a==b):
        if(a=='N'):
            return 0
        else:
            return 1
    else:
        return -3
    

lines = [line.rstrip('\n') for line in open('sr.txt')]#reads Short Reads file
sr=lines
lines = [line.rstrip('\n') for line in open('rs.txt')]#reads Reference Gene file
rs=lines
rs="".join(rs)
s=sr[1113344]


def nw(s,rs):
    n=len(s)+1
    m=len(rs)+1
    f=np.zeros((n,m),dtype=int)
    p = np.zeros((n,m),dtype = str)
    gp=-2


    for i in range(n):
        f[i][0]=i*gp
        
        for i in range(m):
            f[0][i]=i*gp

    for i in range(1,n):
        for j in range(1,m):
            f[i][j]=max(f[i-1][j-1]+match_score(s[i-1],rs[j-1]),f[i-1][j]+gp,f[i][j-1]+gp)
        
#print(f)
#np.savetxt('test.out', f, delimiter='  ')
        

    aa=""
    ab=""
    i=n-1
    j=m-1

    while(i>0 or j>0):
            if(i>0 and j>0 and f[i][j]==f[i-1][j-1]+match_score(s[i-1],rs[j-1])):
                aa=s[i-1]+aa
                ab=rs[j-1]+ab
                p[i][j]='d'
                i-=1
                j-=1
            elif(i>0 and f[i][j]==f[i-1][j]+gp):
                aa=s[i-1]+aa
                ab="-"+ab
                i-=1
                p[i][j]='h'
            else:
                aa="-"+aa
                ab=rs[j-1]+ab
                j-=1
                p[i][j]='v'
    
    return(aa,ab,f)


ar=np.zeros(4197,dtype=int)
for l in range(5000):
    s=sr[random.choice(range(len(sr)))]           
    mx=0
    mxi=0
    xi=0

    for i in range(0,len(rs),10+len(s)):
        cs=0
        w, e, r = nw(s,rs[i:i+10+len(s)])
        #print(w+'\n'+e)
        # print("\n ########### \n")
        fst=0
        lst=0
        for j in range(len(w)):
            if w[j]!='-':
                fst=j
                break
    
        for k in range(len(w)):
            if w[k]!='-':
                lst=k
        for u in range(fst,lst+1):
            if(w[u]=='N'):
                cs+=0
            elif(w[u]!='-'):
                cs+=1
            elif(w[u]=='-'):
                cs-=2
        if(cs>mx):
            mx=cs
            mxi=xi
        
        #print("\n #####  "+str(i)+" | "+str(xi*len(s))+" | "+str(np.amax(r))+" | "+str(cs)+" | "+str(mxi)+" | "+str(mx) +"  ###### ")
        xi+=1
    if(mx>30):
        for o in range(mxi*len(s),(mxi+1)*len(s)+10):
            ar[o]+=1
   #print('!!!'+str(fst)+"|"+str(lst)+'!!!\n')
    
plt.plot(range(1,4198),ar)
plt.show()
np.savetxt('bio.txt', ar, delimiter=',')
    

            

