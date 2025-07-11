OFDM/TDM Communication System
Project Overview
This project implements the communication system described in the paper "Performance of OFDM/TDM with MMSE-FDE Using Pilot-assisted Channel Estimation." The objective is to recreate Figure 6 from the paper, which illustrates the average Bit Error Rate (BER) for various Eb/N0 values across five distinct channel estimation scenarios:

OFDM/TDM Pilot and Non-Ideal CE
OFDM/TDM Pilot and Ideal CE
OFDM with TDM Pilot and Non-Ideal CE
OFDM with TDM Pilot and Ideal CE
OFDM with FDM Pilot

The implementation is written in C, with separate files handling different cases.
File Descriptions

OFDM_TDM_case1and2.c: Implements cases 1 and 2. Set the case_type variable to 1 for OFDM/TDM with CE or 2 for OFDM/TDM with Ideal CE.
Conventional_OFDM_TDM_case3and4.c: Implements cases 3 and 4. Set the case_type variable to 3 for OFDM with TDM Pilot or 4 for OFDM with TDM Pilot with Ideal CE.
Conventional_OFDM_FDM_case5.c: Implements case 5 (OFDM with FDM Pilot).
ber_data.txt: Output file storing simulation results. Each line follows the format: case_type Eb/N0 BER, where case_type is an integer (1-5), Eb/N0 is in dB, and BER is the bit error rate.
ber_data_all.txt: Aggregates results from all cases in order for plotting.
Plotting_python.txt: Includes a Python script to generate BER vs. Eb/N0 graphs from ber_data_all.txt.

Usage Instructions

Compile the C files: Use a C compiler (e.g., GCC) to compile each .c file.gcc -o case1and2 OFDM_TDM_case1and2.c
gcc -o case3and4 Conventional_OFDM_TDM_case3and4.c
gcc -o case5 Conventional_OFDM_FDM_case5.c


Run a simulation: Execute the desired program, setting case_type as needed. For example:
For case 1: In OFDM_TDM_case1and2.c, set case_type = 1 and run:./case1and2




Collect results: Simulation outputs are appended to ber_data.txt.
Repeat for all cases: Adjust case_type and run each program accordingly.
Aggregate data: Copy the contents of ber_data.txt into ber_data_all.txt in case order (1 through 5).
Generate graphs: Run the Python script in Plotting_python.txt to produce graphs matching Figure 6.

Known Issues

Case 1 (OFDM/TDM with CE): The BER curve does not align with Figure 6 from the paper, despite adhering to the described channel estimation method.
