import numpy as np
import matplotlib.pyplot as plt

# Load data from the text file
data = np.loadtxt('ber_data_all.txt')

# Extract columns
case_types = data[:, 0].astype(int)  # Case type (1 to 5)
eb_n0_db = data[:, 1]                # Eb/N0 in dB
avg_ber = data[:, 2]                 # Average BER

# Create a figure
plt.figure(figsize=(10, 6))

# Plot each case
labels = [
    'OFDM/TDM with CE',
    'OFDM/TDM with Ideal CE',
    'OFDM with TDM Pilot',
    'OFDM with TDM Pilot with Ideal CE',
    'OFDM with FDM Pilot'
]
markers = ['o-', 's-', '^-', 'd-', '*-']  # Different markers for each case

for case in range(1, 6):
    mask = (case_types == case)
    plt.semilogy(eb_n0_db[mask], avg_ber[mask], markers[case-1], label=labels[case-1])

# Customize the plot
plt.xlabel('Eb/N0 (dB)')
plt.ylabel('Average BER')
plt.title('Average BER vs. Average Eb/N0 (dB) - Figure 6')
plt.grid(True, which="both", ls="--")
plt.legend()
plt.savefig('ber_plot.png')  # Save the plot
plt.show()