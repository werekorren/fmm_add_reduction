# Introduction

This repo contains c code to minimize the number of linear operators in fast multiplication algorithms.
Input - a fast multiplication algorithm of Strassen-type - is read from a file.
Addition reduction is then performed to minimize the number of additions.
The available strategies for addition reduction are:
 - brute force
 - greedy vanilla
 - greedy potential

Output is a full specification of the addition-reduced algorithm.
The algorithm is verified (symbolically/algebraically) to output the correct result (the matrix product).
As optional output, it is possible to print both the input and output algorithms in LaTeX format.

The fmm prefix in fmm_add_reduction is for fast matrix multiplication.

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

In the compact form of the naive algorithm printed above, the $\mathbb{A}$, $\mathbb{B}$ and $\mathbb{M}$ correspond to the matrices $\mathcal{A}$, $\mathcal{B}$ and $\mathcal{C}$ in the targeted computation $\mathcal{A}\mathcal{B}=\mathcal{C}$.

The dimensions of $\mathbb{A}$ in the printout are $kl\times{}q$.
Each row represents an entry in the original matrix $\mathcal{A}$ (of dimension $k\times{}l$, so entries $\mathcal{A}_0, \mathcal{A}_1, \mathcal{A}_2, \ldots$), and each column represents an $\mathbb{M}_i, with 0\leq{}i<q,$ expression in the fast matrix multiplication equation system.
A nono-empty entry at $(i, j)$ thus indicates that $\mathcal{A}_i$ is utilized in the expression for $\mathcal{M}_j$.

Correspondingly, the dimensions of $\mathbb{B}$ in the printout are $lm\times{}q$, with each row representing an entry in the original matrix $\mathcal{B}$ (of dimension $l\times{}m$), and each column represents an $\mathbb{M}_i, with 0\leq{}i<q,$ expression in the fast matrix multiplication equation system.

For $\mathbb{M}$, the dimensions in the printout are $q\times{}km$, and each row represents an $\mathbb{M}_i$, with $0\leq{}i<q$, expression in the fast matrix multiplication equation system, with each column representing an entry in $\mathcal{C}$.

The equation system printouts after the reduction are larger, showing an additional space ($t$-space) delimited by vertical and horizontal bars.
The $t$-space details all the variable substitutions that have been introduced to minimize the total number of additions in the system.
For example, looking at the first (index zero) column in the $t$-space of $\mathbb{A}$, it specifies the variable substitution as $u_0=\mathcal{A}_0-\mathcal{A}_3$

# Command-line parameters
<pre>
> fmm_add_reduction.exe [file name] [reduction method] [optional greedy potential parameters] [verbosity level] [latex]
</pre>

### file name
Name of file containing the fast matrix multiplication algorithm that you want to reduce.
Three different file formats are supported.
Sample input files are provided in the $\texttt{algorithms}$ folder.

### reduction method
Three different reduction methods are supported:
 - brute force (use 'bf' or 'brute-force' or 'bruteforce')
 - greedy vanilla (use 'v' or 'gv' or 'vanilla' or "greedy vanilla")
 - greedy potential (use 'p'  or 'gp' or 'potential' or "greedy potential")

### optional greedy potential parameters
Optional greedy potential parameters can be given in three ways.
- Omitting them.
This option uses the default greedy potential parameters, which is to use alpha values from alpha_start=0.0 to alpha_end=0.5 in steps of 0.1 (alpha_steps=5) -- the same for all matrices $\mathcal{A}$, $\mathcal{B}$ and $\mathcal{C}$.
- Specifying 1 (one) set of alpha values (alpha_start, alpha_end, alpha_step).
This option will use the specified alpha values, but will use the same paramater set for each of $\mathcal{A}$, $\mathcal{B}$ and $\mathcal{C}$.
- Specifying 3 (three) consecutive sets of alpha values.
This option lets you specify separate alpha values for each of $\mathcal{A}$, $\mathcal{B}$ and $\mathcal{C}$.

### verbosity level
A number from 0 (silent) and up. A higher number indicates more intermediate output during the search process.
For brute force search, verbosity > 1 prints search statistics at corresponding max depth of the search tree.
For greedy vanilla and greedy potential, intermediate results are printed.

### latex
Add 'latex' option to get LateX printouts of both the initial and the reduced equation systems.

### Examples (Windows console)
Single quotes are not written out, while double quotes are (to allow spaces).</br>
Any of the following commands will brute-force the original Laderman equation system.
<pre>
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt b
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt bf
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt brute-force
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt "brute force"
</pre>
Any of the following commands will run Greedy Vanilla on the original Laderman equation system.
<pre>
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt v
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt gv
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt vanilla
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt "greedy vanilla"
</pre>
Any of the following commands will run Greedy Potential on the original Laderman equation system with the default alpha values as parameters.
<pre>
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt p
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt gp
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt potential
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt "greedy potential"
</pre>
You can specify one set of alpha values to be used with Greedy Potential for all $\mathcal{A}$, $\mathcal{B}$ and $\mathcal{C}$ as follows.
<pre>
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt p 0 0.2 4
</pre>
In this case alpha values $\set{ 0, 0.05, 0.1, 0.15, 0.2 }$ will be used for all $\mathcal{A}$, $\mathcal{B}$ and $\mathcal{C}$.</br>
You can also specify three sets of alpha values to be used with Greedy Potential for each of $\mathcal{A}$, $\mathcal{B}$ and $\mathcal{C}$, respectively as follows.
<pre>
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt p 0 0.2 4 0 0.05 5 0.25 0.37 12
</pre>
In this case alpha values $\set{ 0, 0.05, 0.1, 0.15, 0.2 }$ are used for $\mathcal{A}$, $\set{ 0, 0.01, 0.02, 0.03, 0.04, 0.05 }$ for $\mathcal{B}$ and $\set{ 0.25, 0.26, ..., 0.36, 0.37 }$ for $\mathcal{C}$.</br>
And finally, adding verbosity level 2 and LaTeX printouts can be done as follows.
<pre>
> fmm_add_reduction.exe algorithms\other\Laderman-333-23-98.txt p 0 0.2 4 0 0.05 5 0.25 0.37 12 2 latex
</pre>
