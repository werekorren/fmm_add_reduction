# fmm_add_reduction

This repo contains c code to minimize the number of linear operators in fast multiplication alrorithms.
Input - a fast multiplication algorithm of Strassen-type - is read from a file.
Addition reduction is then performed to minimize the number of additions.
The available strategies for addition reduction are:
 - brute force
 - greedy vanilla
 - greedy potential

Output is a full specification of the addition-reduced algorithm.
The algorithm is verified (symbolically/algebraically) to output the correct result (the matrix product).
As optional output, it is possible to print both the input and output algorithms in LaTeX format.

For additional details, please refer to the paper.


# Sample output
Sample output when run with default parameters on the Laderman algorithm:

<pre>
Using file algorithms/other/Laderman-333-23-98.txt

Naive algorithm in compact form:  
(k, l, m, q) = (3, 3, 3, 23)

A =  
  1  1  0 -1  0  1 -1 -1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  
  1  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  1  0  0  0  0  
  1  0  0  0  0  0  0  0  0  1  0 -1  1  1  0 -1  1  0  0  0  0  0  0  
 -1 -1  0  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  
 -1  0  1  1  1  0  0  0  0 -1  0  0  0  0  0  1  0  1  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  1 -1  1  0  1  0  0  0  
  0  0  0  0  0  0  1  1  1 -1  0  0  0  0  0  0  0  0  0  0  0  1  0  
 -1  0  0  0  0  0  1  0  1 -1  1  1  0  0  1  0  0  0  0  0  0  0  0  
 -1  0  0  0  0  0  0  0  0  0  0  1 -1  0  1  0  0  0  0  0  0  0  1  

B =  
  0  0 -1  1 -1  1  1  0 -1  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  
  0 -1  1 -1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  
  0  0  0  0  0  0 -1  1  1  0  1  0  0  0  0  0  0  0  0  0  1  0  0  
  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0  
  1  1 -1  1  0  0  0  0  0  0 -1  1  1  0  0  0  0  0  0  0  0  0  0  
  0  0 -1  0  0  0  1 -1  0  1 -1  0  0  0  0  1  1  0  0  0  0  0  0  
  0  0 -1  0  0  0  0  0  0  0 -1  1  0  1 -1  1  0 -1  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0  0  1 -1 -1  0  1  0  0  0  0  1  0  0  0  
  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0 -1 -1  1  0  0  0  0  1

M =  
  0  1  0  0  0  0  0  0  0  
  0  0  0  1  1  0  0  0  0  
  0  0  0  1  0  0  0  0  0  
  0  1  0  1  1  0  0  0  0  
  0  1  0  0  1  0  0  0  0  
  1  1  1  1  1  0  1  0  1  
  0  0  1  0  0  0  1  0  1  
  0  0  0  0  0  0  1  0  1  
  0  0  1  0  0  0  0  0  1  
  0  0  1  0  0  0  0  0  0  
  0  0  0  0  0  0  1  0  0  
  0  1  0  0  0  0  1  1  0  
  0  0  0  0  0  0  1  1  0  
  1  1  1  1  0  1  1  1  0  
  0  1  0  0  0  0  0  1  0  
  0  0  1  1  0  1  0  0  0  
  0  0  0  1  0  1  0  0  0  
  0  0  1  0  0  1  0  0  0  
  1  0  0  0  0  0  0  0  0  
  0  0  0  0  1  0  0  0  0  
  0  0  0  0  0  1  0  0  0  
  0  0  0  0  0  0  0  1  0  
  0  0  0  0  0  0  0  0  1

Algorithm correctness verification before reduction: ok

A:  
For alpha = 0.000000 = k2/k1 = [approx] =    0/    1 = 0.000000:   28 ->   18  
For alpha = 0.100000 = k2/k1 = [approx] =    1/   10 = 0.100000:   28 ->   16  
For alpha = 0.200000 = k2/k1 = [approx] =    1/    5 = 0.200000:   28 ->   16  
For alpha = 0.300000 = k2/k1 = [approx] =    3/   10 = 0.300000:   28 ->   16  
For alpha = 0.400000 = k2/k1 = [approx] =    2/    5 = 0.400000:   28 ->   16  
For alpha = 0.500000 = k2/k1 = [approx] =    1/    2 = 0.500000:   28 ->   16  
Best greedy potential params found for A are alpha = 0.100000, or ( k1, k2) = (  10,   1), reducing from  28 to  16 additions

B:  
For alpha = 0.000000 = k2/k1 = [approx] =    0/    1 = 0.000000:   28 ->   18  
For alpha = 0.100000 = k2/k1 = [approx] =    1/   10 = 0.100000:   28 ->   16  
For alpha = 0.200000 = k2/k1 = [approx] =    1/    5 = 0.200000:   28 ->   16  
For alpha = 0.300000 = k2/k1 = [approx] =    3/   10 = 0.300000:   28 ->   16  
For alpha = 0.400000 = k2/k1 = [approx] =    2/    5 = 0.400000:   28 ->   16  
For alpha = 0.500000 = k2/k1 = [approx] =    1/    2 = 0.500000:   28 ->   16  
Best greedy potential params found for B are alpha = 0.100000, or ( k1, k2) = (  10,   1), reducing from  28 to  16 additions

C:  
For alpha = 0.000000 = k2/k1 = [approx] =    0/    1 = 0.000000:   42 ->   34  
For alpha = 0.100000 = k2/k1 = [approx] =    1/   10 = 0.100000:   42 ->   34  
For alpha = 0.200000 = k2/k1 = [approx] =    1/    5 = 0.200000:   42 ->   30  
For alpha = 0.300000 = k2/k1 = [approx] =    3/   10 = 0.300000:   42 ->   30  
For alpha = 0.400000 = k2/k1 = [approx] =    2/    5 = 0.400000:   42 ->   31  
For alpha = 0.500000 = k2/k1 = [approx] =    1/    2 = 0.500000:   42 ->   32  
Best greedy potential params found for C are alpha = 0.200000, or ( k1, k2) = (   5,   1), reducing from  42 to  30 additions

Reduced algorithm in compact form:  
(k, l, m, q) = (3, 3, 3, 23)

A =  
  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  |  1  1  0  0  0  0  0  0  
  1  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  1  0  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  |  0  0  1  1  0  0  0  0  
  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  | -1  0  0  0  0  0  0  0  
  0  0  1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  |  0  0  0  0  1  1  0  0  
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  1  0  0  0  |  0  0 -1  0  0  0  0  0  
  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  1  0  |  0 -1  0  0  0  0  0  0  
  0  0  0  0  0  0  0  0  1  0  1  0  0  0  1  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0  1  1  
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  1  |  0  0  0 -1  0  0  0  0  
------------------------------------------------------------------------------------------------  
  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  |  0  0  0  0 -1  0  0  0  
  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0 -1  0  
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  |  0  0  0  0  0 -1  0  0  
  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0  0 -1  
 -1  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  1  0  0  0  0  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  1  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0  0  0  
 -1  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0  0  0  

B =  
  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  |  1  1  0  0  0  0  0  0  
  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  | -1  0  0  0  0  0  0  0  
  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  |  0 -1  0  0  0  0  0  0  
  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0  |  0  0  0  0  0  0  0  0  
  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  |  0  0  1  0  1  0  0  0  
  0  0  0  0  0  0  0 -1  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  |  0  0  0  1  0  1  0  0  
  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1  0  0 -1  0  0  0  0  0  |  0  0  0  0  0  0  1  1  
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  1  0  0  0  |  0  0 -1  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  1  |  0  0  0 -1  0  0  0  0  
------------------------------------------------------------------------------------------------  
  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  |  0  0  0  0  1  0  0  0  
  0  0  0  0  0  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  1  0  0  
  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0  1  0  
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  |  0  0  0  0  0  0  0  1  
  0  0 -1  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  1  0  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0  0 -1  1  0  0  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0 -1  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  |  0  0  0  0  0  0  0  0  

M =  
  0  1  0  0  0  0  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0  |  0  0  0  0  1  0  0  0  
  0  0  0  1  0  0  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0  |  1  0  0  0  0  0  0  0  
  0  1  0  0  1  0  0  0  0  |  0  0  0  0  0  0  0  0  
  1  0  0  0  0  0  0  0  0  |  1  1  0  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0  |  0  1  0  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  1  0  0  
  0  0  1  0  0  0  0  0  1  |  0  0  0  0  0  0  0  0  
  0  0  1  0  0  0  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  1  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0  |  0  0  1  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0  1  0  
  1  0  0  0  0  0  0  0  0  |  0  0  1  1  0  0  0  0  
  0  1  0  0  0  0  0  1  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  0  0  0  |  0  0  0  1  0  0  0  0  
  0  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0  0  1  
  0  0  1  0  0  1  0  0  0  |  0  0  0  0  0  0  0  0  
  1  0  0  0  0  0  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  1  0  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  1  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  0  1  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  0  0  1  |  0  0  0  0  0  0  0  0  
------------------------------------------------------  
  0  1  0  0  0  0  0  0  0  |  0  0  0  0  1  0  0  0  
  0  0  1  0  0  0  0  0  0  |  0  0  0  0  0  1  0  0  
  0  1  0  0  0  0  0  0  0  |  0  0  0  0  0  0  1  0  
  0  0  1  0  0  0  0  0  0  |  0  0  0  0  0  0  0  1  
  0  0  0  1  1  0  0  0  0  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  1  0  1  |  0  0  0  0  0  0  0  0  
  0  0  0  0  0  0  1  1  0  |  0  0  0  0  0  0  0  0  
  0  0  0  1  0  1  0  0  0  |  0  0  0  0  0  0  0  0  

Algorithm correctness verification after reduction: ok

Algorithm uses  28 +  28 +  42 =  98 additions [naive]  
Algorithm uses  16 +  16 +  30 =  62 additions after reduction with alpha = ( 0.100000, 0.100000, 0.200000) [greedy potential]  
               ---------------------  
Total savings:  12 +  12 +  12 =  36 additions (36.73%)  

Reduction runtime: 0.040000 sec


Process returned 0 (0x0)   execution time : 0.182 s
Press any key to continue.
</pre>
