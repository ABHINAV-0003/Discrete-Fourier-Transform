# Discrete-Fourier-Transform
## **Python script for calculating DFT of N bit finite sequence.**

This script will help you to calculate Discrete Fourier Transform of N bit finite Sequence .

You have to enter N - Number of bits in sequence
Enter the sequence of N bits seperated by commas ','. (Assumed : First element is at origin.)

Script will generate Twiddle Matrix and DFT for input N bit sequence.

**Sample Input:**

Enter the value of N  
4
Enter N bit sequence......
-1,1,-1,1

**Sample Output:**

No of discrete bits in input:   4

Twiddle Matrix is :-  
     0        1    2        3
0  1.0   (1+0j)  1.0   (1+0j)
1  1.0  (-0-1j) -1.0       1j
2  1.0  (-1+0j)  1.0  (-1+0j)
3  1.0       1j -1.0  (-0-1j)

Input sequence is : [-1, 1, -1, 1]

DFT of the input sequence is :- 
         0
0       0j
1       0j
2  (-4+0j)
3       0j

Requirements:
1. Python 3.0 or more
2. Sympy, panda and cmath libraries installed.
