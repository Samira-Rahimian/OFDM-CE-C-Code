#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

//This code implements Cases 3 and 4: 
// 3- OFDM with FDM Pilot
// 4- OFDM with TDM Pilot with Ideal CE
// Please choose your case_type!

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define N 256    // Subcarriers per OFDM symbol
#define CP 16    // Cyclic prefix length
#define K 16     // OFDM symbols per frame
#define L 8      // Channel paths
#define BETA_DB 3.0 // HPA saturation level in dB


// Function prototypes
double complex generate_qpsk();
double complex generate_noise(double sigma);
void generate_channel(double complex* h);
void generate_chu_sequence(double complex* P, int size);
void ifft(double complex* X, double complex* x, int size);
void fft(double complex* x, double complex* X, int size);
void add_cyclic_prefix(double complex* x, double complex* x_cp, int size, int cp_len);
void remove_cyclic_prefix(double complex* r_cp, double complex* r, int size, int cp_len);
void convolve_with_channel(double complex* x_cp, double complex* h, double complex* r_cp, int len, int num_paths);


// Generate QPSK symbol: (±1 ± j) / sqrt(2)
double complex generate_qpsk() {
    double re = (rand() % 2 ? 1 : -1) / sqrt(2);
    double im = (rand() % 2 ? 1 : -1) / sqrt(2);
    return re + im * I;
}

// Generate complex Gaussian noise
double complex generate_noise(double sigma) {
    double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
    double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
    double re = sigma * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
    double im = sigma * sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
    return re + im * I;
}

// Generate L-path Rayleigh fading channel
void generate_channel(double complex* h) {
    double sigma = sqrt(0.5 / L);
    for (int l = 0; l < L; l++) {
        double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
        double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
        double z0 = sigma * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
        double z1 = sigma * sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
        h[l] = z0 + I * z1;
    }
}

// Generate Chu sequence for pilot
void generate_chu_sequence(double complex* P, int size) {
    for (int m = 0; m < size; m++) {
        double angle = M_PI * m * m / size;
        P[m] = cos(angle) + I * sin(angle); // e^(j * pi * m^2 / N)
    }
}

// IFFT: Scaled by 1/sqrt(size)
void ifft(double complex* X, double complex* x, int size) {
    for (int n = 0; n < size; n++) {
        x[n] = 0;
        for (int k = 0; k < size; k++) {
            double angle = 2 * M_PI * k * n / size;
            x[n] += X[k] * (cos(angle) + I * sin(angle));
        }
        x[n] /= sqrt(size);
    }
}

// FFT: Scaled by 1/sqrt(size)
void fft(double complex* x, double complex* X, int size) {
    for (int k = 0; k < size; k++) {
        X[k] = 0;
        for (int n = 0; n < size; n++) {
            double angle = -2 * M_PI * k * n / size;
            X[k] += x[n] * (cos(angle) + I * sin(angle));
        }
        X[k] /= sqrt(size);
    }
}

// Add cyclic prefix
void add_cyclic_prefix(double complex* x, double complex* x_cp, int size, int cp_len) {
    for (int t = 0; t < cp_len; t++) {
        x_cp[t] = x[size - cp_len + t];
    }
    for (int t = 0; t < size; t++) {
        x_cp[cp_len + t] = x[t];
    }
}
// Apply HPA
void apply_hpa(double complex* signal, int len, double A) {
    for (int t = 0; t < len; t++) {
        double mag = cabs(signal[t]);
        if (mag > A) {
            signal[t] = A * signal[t] / mag;
        }
    }
}

// Remove cyclic prefix
void remove_cyclic_prefix(double complex* r_cp, double complex* r, int size, int cp_len) {
    for (int t = 0; t < size; t++) {
        r[t] = r_cp[cp_len + t];
    }
}

// Convolve with channel
void convolve_with_channel(double complex* x_cp, double complex* h, double complex* r_cp, int len, int num_paths) {
    for (int t = 0; t < len; t++) {
        r_cp[t] = 0;
        for (int l = 0; l < num_paths && t - l >= 0; l++) {
            r_cp[t] += h[l] * x_cp[t - l];
        }
    }
}

// Compute true channel frequency response H[k]
void compute_true_channel(double complex* h, double complex* H_true, int size, int num_paths) {
    for (int k = 0; k < size; k++) {
        H_true[k] = 0;
        for (int l = 0; l < num_paths; l++) {
            double angle = -2 * M_PI * k * l / size;
            H_true[k] += h[l] * (cos(angle) + I * sin(angle));
        }
    }
}
// MMSE FDE
void mmse_fde(double complex* R, double complex* Y_eq, double complex* H_Eq, double sigma2_Eq, double alpha, double beta) {
    for (int m = 0; m < N; m++) {
        double complex h = H_Eq[m];
        double h_mag2 = creal(h) * creal(h) + cimag(h) * cimag(h);
        double denom = h_mag2 * (1 - exp(-beta * beta)) + 2 * sigma2_Eq;
        double complex w = (alpha * conj(h)) / denom;
        Y_eq[m] = R[m] * w;
    }
}

int main() {
    srand(time(NULL));
    // Case selection (3 to 4)
    int case_type = 3; 

    // HPA parameters
    double beta_dB = BETA_DB;
    double A = sqrt(pow(10, beta_dB / 10.0));
    double beta = pow(10, beta_dB / 20.0);
    double erfc_beta = erfc(beta);
    double alpha = 1 - exp(-beta* beta) + (sqrt(M_PI) * beta / 2.0) * erfc_beta;


    double complex P[N];
    generate_chu_sequence(P, N); // Generate Chu sequence for pilot

    FILE* fp = fopen("ber_data.txt", "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    for (double Eb_N0_dB = 5; Eb_N0_dB <= 30; Eb_N0_dB += 5) {
        double Eb_N0 = pow(10, Eb_N0_dB / 10.0);
        double Es = 1.0; // Symbol energy (QPSK has unit energy)
        double Eb = Es / 2.0; // 2 bits per symbol
        double N0 = Eb / Eb_N0;
        double sigma2 = N0 / 2.0; // Noise variance per real/imag component
        double sigma = sqrt(sigma2);
        double total_ber = 0.0;
        int num_frames = 1000;

        for (int frame_idx = 0; frame_idx < num_frames; frame_idx++) {
            // Generate channel
            double complex h[L];
            generate_channel(h);

            // Pilot symbol (k=0)
            double complex X_pilot[N];
            for (int m = 0; m < N; m++) {
                X_pilot[m] = P[m];
            }
            double complex x_pilot[N];
            ifft(X_pilot, x_pilot, N);
            double complex x_pilot_cp[N + CP];
            add_cyclic_prefix(x_pilot, x_pilot_cp, N, CP);
            apply_hpa(x_pilot_cp, N + CP, A);
            double complex r_pilot_cp[N + CP];
            convolve_with_channel(x_pilot_cp, h, r_pilot_cp, N + CP, L);
            for (int t = 0; t < N + CP; t++) {
                r_pilot_cp[t] += generate_noise(sigma);
            }
            double complex r_pilot[N];
            remove_cyclic_prefix(r_pilot_cp, r_pilot, N, CP);
            double complex R_pilot[N];
            fft(r_pilot, R_pilot, N);
            double complex H_est[N];
            double complex N_e[N];
            double sum_power = 0.0;
            double sigma2_e;
            for (int m = 0; m < N; m++) {
                H_est[m] = R_pilot[m] / P[m];
                // Noise power estimation: N_e(m) = R_pilot(m) - H_est(m) * P(m)
                N_e[m] = R_pilot[m] - H_est[m] * P[m];
                sum_power += creal(N_e[m]) * creal(N_e[m]) + cimag(N_e[m]) * cimag(N_e[m]);
            }
            sigma2_e = sum_power / N;
            double complex H_true[N];
            compute_true_channel(h, H_true, N, L);
            // Data symbols (k=1 to K-1)
            int bit_errors = 0;
            for (int k = 1; k < K; k++) {
                double complex X[N];
                for (int m = 0; m < N; m++) {
                    X[m] = generate_qpsk();
                }
                double complex x[N];
                ifft(X, x, N);
                double complex x_cp[N + CP];
                add_cyclic_prefix(x, x_cp, N, CP);
                apply_hpa(x_cp, N + CP, A);
                double complex r_cp[N + CP];
                convolve_with_channel(x_cp, h, r_cp, N + CP, L);
                for (int t = 0; t < N + CP; t++) {
                    r_cp[t] += generate_noise(sigma);
                }
                double complex r[N];
                remove_cyclic_prefix(r_cp, r, N, CP);
                double complex R[N];
                fft(r, R, N);
                double complex Y_eq[N];
                double sigma2_Eq;
                double complex H_Eq[N];
                if (case_type == 3) {
                    for (int m = 0; m < N; m++) {
                        H_Eq[m] = H_est[m];
                    }
                    sigma2_Eq = sigma2_e;
                }
                else if (case_type == 4) {
                    for (int m = 0; m < N; m++) {
                        H_Eq[m] = H_true[m];
                    }
                    sigma2_Eq = sigma2;
                }
                mmse_fde(R, Y_eq, H_Eq, sigma2_Eq, alpha, beta);
                for (int m = 0; m < N; m++) {
                    double re = creal(Y_eq[m]);
                    double im = cimag(Y_eq[m]);
                    double complex X_decided = (re > 0 ? 1 : -1) / sqrt(2) + I * (im > 0 ? 1 : -1) / sqrt(2);
                    if (creal(X[m]) * creal(X_decided) < 0) bit_errors++;
                    if (cimag(X[m]) * cimag(X_decided) < 0) bit_errors++;
                }
            }
            double number_of_bits = 2.0 * (K - 1) * N;
            total_ber += (double)bit_errors / number_of_bits;
        }
        double avg_ber = total_ber / num_frames;
        printf("Eb/N0 = %.1f dB, BER = %e\n", Eb_N0_dB, avg_ber);
        fprintf(fp, "%d %.1f %e\n", case_type, Eb_N0_dB, avg_ber);
    }
    fclose(fp);
    return 0;
}