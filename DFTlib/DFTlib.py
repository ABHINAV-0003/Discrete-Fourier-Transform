# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 23:12:30 2017

@author: LALIT ARORA
"""

from sympy import *
from pandas import *
import cmath
import math

def power(a,p):
    temp=pow(a,p)
    realtemp=temp.real
    imgtemp=temp.imag
    realtemp=math.ceil(realtemp*1000)/1000
    imgtemp=math.ceil(imgtemp*1000)/1000
    if (realtemp==0 or realtemp==-0) and imgtemp!=0:
        return (1j*imgtemp)
    elif realtemp!=0 and (imgtemp==1j*0 or imgtemp==-1j*0):
        return (realtemp)
    else:
        temp=realtemp+1j*imgtemp
        return (temp)
    
    
def wp(N):
    
    x=math.ceil(cos(pi*(2/N))*100)/100
    y=math.ceil(sin(pi*(2/N))*100)/100
    if x==0 and y!=0:
        return (-1j*y)
    elif x!=0 and y==0:
        return x
    else:
        return (x-1j*y)
        
def cwp(N):
    x=math.ceil(cos(pi*(2/N))*100)/100
    y=math.ceil(sin(pi*(2/N))*100)/100
    if x==0 and y!=0:
        return (1j*y)
    elif x!=0 and y==0:
        return x
    else:
        return (x+1j*y)
    
def TwiddleMatrix(N,st):
    finalmat=[]
    k=0
    if st=="TM":
        for i in range(N):
            temp=[]
            for j in range(N):
                temp.append(power(wp(N),(k*j)))
            k=k+1
            finalmat.append(temp)
        return finalmat
    elif st=="CTM":
        for i in range(N):
            temp=[]
            for j in range(N):
                temp.append(power(cwp(N),(k*j)))
            k=k+1
            finalmat.append(temp)
        return finalmat
    else:
        return("Error Occured!..")
# TM is for user
def TM(N,st):
    matrix=TwiddleMatrix(N,st)
    if matrix=="Error Occured!..":
        return "Error Occured!.."
    else:
        return DataFrame(matrix)

def DFT(seq):
    N=len(seq)
    matrix=TwiddleMatrix(N,"TM")
    multiply=[]
    for j in range(N):
        temp=0
        for i in range(N):
            temp=temp+(seq[i]*matrix[j][i])
        multiply.append(temp)
    return DataFrame(multiply)

def DFTcoder(seq):
    N=len(seq)
    matrix=TwiddleMatrix(N,"TM")
    multiply=[]
    for j in range(N):
        temp=0
        for i in range(N):
            temp=temp+(seq[i]*matrix[j][i])
        multiply.append(temp)
    return (multiply)

def IDFTcoder(seq):
    N=len(seq)
    matrix=TwiddleMatrix(N,"CTM")
    multiply=[]
    for j in range(N):
        temp=0
        for i in range(N):
            temp=temp+(seq[i]*matrix[j][i])
        multiply.append(temp/N)
    return (multiply)

def IDFT(seq):
    N=len(seq)
    matrix=TwiddleMatrix(N,"CTM")
    multiply=[]
    for j in range(N):
        temp=0
        for i in range(N):
            temp=temp+(seq[i]*matrix[j][i])
        multiply.append(temp/N)
    return DataFrame(multiply)

def circonv(seq1,seq2):
    N1=len(seq1)
    N2=len(seq2)
    if N1>N2:
        temp=N1-N2
        for i in range(temp):
            seq2.append(0)
    elif N2>N1:
        temp=N2-N1
        for i in range(temp):
            seq1.append(0)
    
    
    dftseq1=DFTcoder(seq1)
    dftseq2=DFTcoder(seq2)
    idftmatrix=[]
    for i in range(len(dftseq1)):
        idftmatrix.append(dftseq1[i]*dftseq2[i])
    
    retmatrix=IDFTcoder(idftmatrix)
    return DataFrame(retmatrix)

def lcrotate(l, n):
    return l[n:] + l[:n]

def rcrotate(l, n):
    return l[-n:] + l[:-n]

def linconv(seq1,seq2):
    N1=len(seq1)
    N2=len(seq2)
    L=N1+N2-1
    for i in range(N2-1):
        seq1.append(0)
    for i in range(N1-1):
        seq2.append(0)
    matrix=[]
    temp=[]
    temp.append(seq1[0])
    
    revseq1=seq1[::-1]
    for i in range(L-1):
        temp.append(revseq1[i])
    matrix.append(temp)
    for i in range(1,L):
        matrix.append(rcrotate(matrix[i-1],1))
    finalmat=[]
    for i in range(L):
        temp=0
        for j in range(L):
            temp=temp+matrix[i][j]*seq2[j]
        finalmat.append(temp)
        
    return DataFrame(finalmat)


def blocks(seq1,seq2):
    N2=len(seq2)
    N1=len(seq1)
    no=math.ceil(float(N1/N2))
    exlength=((N2))*(no+1)
    matrix=[]
    k=0
    for i in range(exlength-N1):
        seq1.append(0)
    
    while(k<exlength):
        temp=[]
        for j in range((2*N2)-1):
            if(1-N2+k<0):
                temp.append(0)
            else:
                temp.append(seq1[1-N2+k])
            k=k+1
            if j==((2*N2)-2):
                k=k-(N2-1)
        
        matrix.append(temp)
    return matrix
    
        
def convos(matrix,seq):
    L=len(matrix[0])
    anothermatrix=[]
    for i in range(len(matrix)):
        temp=[]
        sendmatrix=[]
        temp.append(matrix[i][0])
        revseq1=matrix[i][::-1]
        for i in range(L-1):
            temp.append(revseq1[i])
        sendmatrix.append(temp)
        for i in range(1,L):
            sendmatrix.append(rcrotate(sendmatrix[i-1],1))
        finalmat=[]
        for i in range(L):
            temp=0
            for j in range(L):
                temp=temp+sendmatrix[i][j]*seq[j]
            finalmat.append(temp)
        anothermatrix.append(finalmat)
    return anothermatrix


def osconv(seq1,seq2):
    N2=len(seq2)
    N1=len(seq1)
    if N1>N2:
        matrix=blocks(seq1,seq2)
        print("Blocks")
        print(DataFrame(matrix))
        length=len(matrix[0])
        diff=length-N2
        for j in range(diff):
            seq2.append(0)
        print(seq2)
        omat=convos(matrix,seq2)
        print("Convolved vectors before overlaping")
        print(DataFrame(omat))
        return(overlap(omat,N2))
    else:
        matrix=blocks(seq2,seq1)
        print("Blocks")
        print(DataFrame(matrix))
        length=len(matrix[0])
        diff=length-N1
        for j in range(diff):
            seq1.append(0)
        print(seq1)
        omat=convos(matrix,seq1)
        print("Convolved vectors before overlaping")
        print(DataFrame(omat))
        print("Overlaped and convolved output")
        return (overlap(omat,N1))
        
    

def overlap(mat,N):
    linmat=[]
    L=len(mat[0])
    for i in range(len(mat)):
        for j in range(N-1,L):
            linmat.append(mat[i][j])
    return linmat 



def blocks_add(seq1,seq2):
    N1=len(seq1)
    N2=len(seq2)
    no=math.ceil(float(N1/N2))
    exlength=no*N2
    for i in range(exlength-N1):
        seq1.append(0)
    matrix=[]
    k=0
    for i in range(no):
        temp=[]
        for j in range(N2):
            temp.append(seq1[k])
            k=k+1
        for t in range(N2-1):
            temp.append(0)
        matrix.append(temp)
    return matrix

def convos_add(matrix,seq):
    L=len(matrix[0])
    anothermatrix=[]
    for i in range(len(matrix)):
        temp=[]
        sendmatrix=[]
        temp.append(matrix[i][0])
        revseq1=matrix[i][::-1]
        for i in range(L-1):
            temp.append(revseq1[i])
        sendmatrix.append(temp)
        for i in range(1,L):
            sendmatrix.append(rcrotate(sendmatrix[i-1],1))
        finalmat=[]
        for i in range(L):
            temp=0
            for j in range(L):
                temp=temp+sendmatrix[i][j]*seq[j]
            finalmat.append(temp)
        anothermatrix.append(finalmat)
    return anothermatrix

def oaconv(seq1,seq2):
    N1=len(seq1)
    N2=len(seq2)
    if N1>N2:
        matrix=blocks_add(seq1,seq2)
        print("Blocks")
        print(DataFrame(matrix))
        length=len(matrix[0])
        for j in range(length-N2):
            seq2.append(0)
        mat=convos_add(matrix,seq2)
        print("After Convolution")
        print(DataFrame(mat))
        print("Added and Convolved..")
        return add(mat,N2)
    else:
        matrix=blocks_add(seq2,seq1)
        print("Blocks")
        print(DataFrame(matrix))
        length=len(matrix[0])
        for j in range(length-N1):
            seq1.append(0)
        mat=convos_add(matrix,seq1)
        print("After Convolution")
        print(DataFrame(mat))
        print("Added and Convolved..")
        return add(mat,N1)
        
def add(mat,N2):
    L=len(mat)
    linmat=[]
    L1=len(mat[0])
    for i in range(L1-N2+1):
        linmat.append(mat[0][i])
        
    for i in range(1,L):
        for j in range(len(mat[i])):
            if j<(L1-N2):
                linmat.append(mat[i][j]+mat[i-1][j+L1+1-N2])
            if j<L1-(N2-1) and j>(L1-(N2+1)):
                linmat.append(mat[i][j])
            
    for t in range(L1-(N2-1),L1):
        linmat.append(mat[i][t])
    
    return linmat  

if __name__ == '__main__':
    # write your code here
    
    
    

