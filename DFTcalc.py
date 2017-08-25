# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 20:20:59 2017

@author: LALIT ARORA

DFT Calculator for N bit finite sequence 
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

print("Enter the value of N  ")
N=int(input())
finalmat=[]
k=0
for i in range(N):
 
    temp=[]
    for j in range(N):
        temp.append(power(wp(N),(k*j)))
    k=k+1
    finalmat.append(temp)

def takeNbitseq():
    print("Enter N bit sequence......")
    ls = input().strip().split(',')
    for i in range(len(ls)):
        ls[i]=int(ls[i])
    return ls

Nbitseq=takeNbitseq()

print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
print()
print("No of discrete bits in input:  ",N)
print()
print("Twiddle Matrix is :-  ")
print (DataFrame(finalmat))

multiply=[]
for j in range(N):
    temp=0
    for i in range(N):
        temp=temp+(Nbitseq[i]*finalmat[j][i])
    multiply.append(temp)
print()
print ("Input sequence is :",Nbitseq)
print()
print("DFT of the input sequence is :- ")
print(DataFrame(multiply))
print()
print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")