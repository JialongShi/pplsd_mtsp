# pplsd_mtsp
# PPLS/D algorithm for solving the multiobjective TSP problem
# Jialong Shi (jialong.shi@xjtu.edu.cn)



## Introduction

The code is distributed for research use. The author reserves all rights to the code.

In this program, PPLS/D (Parallel Pareto Local Search based on Decomposition) is implemented to optimize the multiobjective Traveling Salesman Problem (mTSP). It returns an approximation of the Pareto optimal solution set of the input mtsp instance. 

PPLS/D is a parallel algorithm with multiple parallel processes. Hence this program requires the MPI tool, such as MPICH or OpenMPI. In this document we use MPICH as example to show how to compile and use this program.

** We also provide a sequential version of this program, which does not require MPI and hence is more friendly for debugging. Please refer to     https://github.com/JialongShi/pplsd_mtsp_sequential

Relevant literature:

[1] Shi J, Zhang Q, Sun J. PPLS/D: Parallel Pareto local search based on decomposition[J]. IEEE transactions on cybernetics, 2020, 50(3): 1060-1071.

[2] Shi J, Zhang Q, Derbel B, et al. Using parallel strategies to speed up Pareto local search[C] Asia-Pacific Conference on Simulated Evolution and Learning. 2017: 62-74.



## The mTSP problem

In the mTSP, G = (V, E) is a fully connected graph where V is its node set with n nodes and E is the edge set. Each edge in E is corresponded to m different costs. Here we denote the m costs of the edge connecting node i and node j as {c_{i,j,1},c_{i,j,2},...,c_{i,j,m}} and all costs are larger than 0. A feasible solution x = (x(1),x(2),...x(n)) is a permutation of the n nodes which represents a Hamilton cycle passing through every node in V exactly once. The mTSP problem can be formalized as follows

     minimize f_k(x) = c_{x(1),x(2),k} + c_{x(2),x(3),k} + ... + c_{x(n-1),x(n),k} + c_{x(n),x(1),k},    k  = 1,...,m.
	 
where f_k(x) is the k-th objective function. In this program, the neighborhood move in the mTSP is based on the 2-Opt move, in which two edges in the current solution are replaced by two other edges.


## File list (the following files should be in the same directory)

- Source code: main.cpp  problem.h  problem.cpp  solution.h  solution.cpp  archive.h  archive.cpp  misc.h  misc.cpp  pplsd.h  pplsd.cpp

- Problem file examples: example_mtsp_m2_n50.tsp example_mtsp_m3_n30.tsp

- Weight vector file examples: example_wv_m2_6procs.txt example_wv_m3_15procs.txt

- cmake file: CMakeLists.txt



## Requires

As mentioned earlier, this program requires MPI tool like MPICH. To install MPICH and verify the installation:
```
sudo apt-get install mpich
which mpic++
```

CMakeLists.txt has been provided for compiling the program by cmake, so better make sure that cmake has been installed. However, cmake is not compulsory, you can also compile the program without cmake. We will show the commands in the next section.



## How to compile it

You can compile the project by running

```
cd <directory-of-this-program>
cmake CMakeLists.txt
make
```

Then an executable named ‘pplsd’ will appear in the directory. 

!!! If the cmake command or the make command is failed, you can directly compile the program by the following command. 

```
mpic++  -o  pplsd  -O2  main.cpp  problem.cpp  solution.cpp  archive.cpp  misc.cpp  pplsd.cpp
```



## How to use it

You can run the program by the 'mpiexec' command of MPICH:

```
mpiexec -n <process_num> ./pplsd  <problem_filename>  <weight_vector_filename>   <max_runtime>
```

** Describes **
- <process_num>: This argument indicates the parallel process number of PPLS/D. You better make sure that this number is not larger than your computer's core number.

- <problem_filename>: This argument is the input mtsp problem filename. We have given two example problem files in the package. The format of the problem file is

         NAME: <problem_name>
         TYPE: mTSP
         COMMENT: <comment>
         OBJECTIVE_NUM: <M>
         DIMENSION: <N>
		 EDGE_WEIGHT_TYPE : EUC_2D
         NODE_COORD_SECTION
		 <node_index> <coord_x_1> <coord_y_1> <coord_x_2> <coord_y_2> ... <coord_x_M> <coord_y_M>
         ...
		 EOF
		 
We can see that in the mTSP file each node has M different (x,y) coordinate, which will generate M different costs for each edge.

-  <weight_vector_filename>: This argument is the file describes the weight vectors of different PPLS/D processes. The format is 

     <process_ID>  <weight_value_1>  <weight_value_2> ... <weight_value_M>

Note here that the process ID must start from 0, the weight vector dimension M should be equal to the objective number of the input problem, and the weight vector number (i.e. the line number) should be equal to the process number <process_num>. We have given two example problem files in the package.

- <max_runtime>: This argument is the runtime budget. The units are seconds. Note here that PPLS/D may converge and finish before the time budget.

** Example **

```
mpiexec -n 6 ./pplsd  example_mtsp_m2_n50.tsp  example_wv_m2_6procs.txt  120
```

The above command runs a 6-process PPLS/D on a 2-objective mtsp. The runtime budget is 120s.



## Outputs

After finish, the program will create a file named ‘result_final.txt’ in the same directory. It lists the solutions in the final archive of PPLS/D. The first line is the solution number in the final archive. The second line are the objective number M and the dimension N of the problem. From the third line to the end, each line indicates a solution. In each line, the first M values is the objective function vector of the solution and the last N values are the variable values of the solution. The format can be summarize as

     <solution_num>
     <obj_num_M>  <dim_N>
     <obj_value_1> <obj_value_2> ... <obj_value_M> <var_1> <var_2> ... <var_N>
     ... ...

The program also will output the sub-archive of each process, which are listed in the files ‘result_process0.txt’, ‘results_process_1.txt’, ...


## Visualization

For 2-objective or 3-objective problems, the following MATLAB script can help to visualize the result. Remember to modify the values of ‘m’, ‘procs_num’ and ‘dir’ before using the script.

```
clear all
close all
clc

m = 3;
procs_num = 15;
dir = '.';

figure()
for procs_index = 0:procs_num-1
    filename = sprintf('%s\\result_process%d.txt',dir,procs_index);
    res = read_result_file(filename);
    if m == 2
        scatter(res(:,1),res(:,2));
        hold on
    elseif m == 3
        scatter3(res(:,1),res(:,2),res(:,3));
        hold on
    end
end
hold off


function result  = read_result_file( filename )

    result = [];
    FID = fopen(filename, 'r');
    if FID == -1 
        disp(['ERROR: cannot open file',filename]);
    else
        archiveSize = fscanf(FID, '%d', 1);
        if archiveSize > 0
            m = fscanf(FID, '%d',1);
            n = fscanf(FID, '%d',1);
            result = zeros(archiveSize,m+n);
            for line = 1:archiveSize
                fit = fscanf(FID, '%lf',m);
                sol = fscanf(FID, '%d',n);
                result(line,:) = [fit; sol]';
            end
        end
        fclose(FID);
    end

end
```


