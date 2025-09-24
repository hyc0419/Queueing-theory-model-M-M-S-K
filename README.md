# Queueing-theory-model-M-M-S-K
M/M/s/k Queue Simulation
#Overview
This C program simulates an M/M/s/k queueing system, where:

M: Poisson arrival process (exponential inter-arrival times)
M: Exponential service times
s: Number of servers (defined as S)
k: Maximum system capacity, including those in service and in queue (defined as K)

The program generates random inter-arrival and service times, simulates packet arrivals and departures, and computes performance metrics such as average system time (W), average queue time (Wq), average number of customers in the system (L), and average number in queue (Lq). Both simulation-based and analytical results are calculated and output to the console and a file named mmsk.csv.
#Features

Random Number Generation: Uses a linear congruential generator (lcgrand) to produce pseudo-random numbers for generating exponential inter-arrival and service times.
Queueing Model: Simulates a multi-server queue with limited capacity (k), where packets are dropped if the system is full.
Performance Metrics:
L: Average number of packets in the system.
Lq: Average number of packets in the queue.
W: Average time a packet spends in the system.
Wq: Average time a packet spends in the queue.


Analytical Calculations: Computes theoretical values for L, Lq, W, and Wq using M/M/s/k queueing formulas.
Output: Results are printed to the console and saved to mmsk.csv for further analysis.

#Constants and Parameters

MODLUS, MULT1, MULT2: Parameters for the random number generator.
S: Number of servers (default: 2).
K: System capacity (default: 8).
lambda: Arrival rate (default: 25.0 packets per minute).
mu: Service rate (default: 40.0 packets per minute).
n: Number of packets to simulate (default: 1,000,000).

#Program Structure

Random Number Generator:
lcgrand(int stream): Generates a random number for a given stream.
lcgrandst(long zset, int stream): Sets the seed for a stream.
lcgrandgt(int stream): Retrieves the current seed for a stream.


Queueing Simulation:
struct node: Represents a packet with arrival time, service time, departure time, and a pointer to the next packet.
struct hd: Represents a server with usage status and departure time.
struct Queue: Represents the queue holding packets waiting for service.
poisson_intertime(double rate): Generates exponential inter-arrival or service times.
getarrivalpacketinfo(double lambda, double mu, int n): Generates a linked list of n packets with arrival and service times.
arrive: Handles packet arrivals, assigning them to free servers or the queue.
departure: Handles packet departures, updating server status and processing queued packets.


Analytical Functions:
factorial(int n): Computes factorial for probability calculations.
P0(float lambda, float mu, int s, int k): Computes the probability of zero customers in the system.
Pn(float lambda, float mu, int s, int k, int n): Computes the probability of n customers in the system.
Lq_math(float lambda, float mu, int s, int k): Computes the average number in queue (Lq).
L_math(float lambda, float mu, int s, int k): Computes the average number in the system (L).


Main Function:
Initializes servers and queue.
Simulates packet arrivals and departures.
Computes and outputs simulation-based and analytical metrics.



#Compilation and Execution

Requirements: A C compiler (e.g., gcc).
Compile:gcc -o mmsk mmsk.c -lm

The -lm flag links the math library for functions like log and pow.
Run:./mmsk


Output:
Console: Prints L, Lq, W, Wq (simulation and analytical results).
File: Appends results to mmsk.csv in the format:L_sim,Lq_sim,W_sim,Wq_sim,L_math,Lq_math,W_math,Wq_math





#Usage Notes

The program uses a fixed seed array (zrng) for reproducibility. Modify zrng or the random number generator functions for different random sequences.
The simulation assumes a stable system (i.e., lambda / (s * mu) < 1).
The output file mmsk.csv is opened in append mode (a+), so running the program multiple times will append new rows to the file.
To change the simulation parameters, modify the constants S, K, lambda, mu, or n in the code.

#Example Output
For the default parameters (lambda = 25.0, mu = 40.0, S = 2, K = 8, n = 1,000,000):
0.73512456, 0.11023456, 0.02940578, 0.00440988
0.73123456, 0.10876543, 0.02924938, 0.00435062


The first line shows simulation results (L, Lq, W, Wq).
The second line shows analytical results.
The same values are appended to mmsk.csv.

#Limitations

The random number generator is deterministic and uses a predefined seed array, which may limit randomness for large simulations.
The program assumes packets are dropped when the system is full (queue length + servers = K).
Floating-point precision may affect analytical calculations for very large or small values of lambda or mu.

Future Improvements

Add command-line arguments to customize lambda, mu, S, K, and n.
Implement a more robust random number generator.
Include error handling for file operations and invalid parameter inputs.
Add support for multiple runs to compute confidence intervals for simulation results.
