#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>

//This code implements Cases 1 and 2: 
// 1- OFDM/TDM 
// 2- OFDM/TDM with Ideal CE
// Please choose your case_type accordingly!
//This code considers correlation in channel generation.

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define Nc 256    // Frame length
#define Nm 16     // Subcarriers per slot (IFFT size)
#define K 16      // Number of slots
#define GI 16     // Guard interval
#define L 8       // Number of channel paths
#define BETA_DB 3.0 // HPA saturation level in dB
#define P 1.0     // Transmit power coefficient
#define GAMMA 0.05

// Function prototypes
double complex generate_qpsk();
void ifft(double complex* X, double complex* x, int N);
void fft(double complex* x, double complex* X, int N);
void ifft_interp(double complex* X, double complex* x, int N);
void fft_interp(double complex* x, double complex* X, int N, int M);
void generate_ofdm_tdm_frame(double complex frame[Nc], double complex X[K][Nm], double complex* pilot_symbols);
void add_cyclic_prefix(double complex frame[Nc], double complex tx_signal[Nc + GI]);
void apply_hpa(double complex* tx_signal, int len, double A);
double complex generate_noise(double sigma);
void generate_channel(double complex* h, double rho, double complex* h_prev);
void apply_channel(double complex* tx_signal, double complex* r, double complex* h, double sigma);
void receive_frame(double complex* r, double complex* R);
void channel_estimation(double complex* r, double complex* H_est, double complex* pilot_symbols, double* sigma2_e);
void mmse_fde(double complex* R, double complex* Y_eq, double complex* H_est, double sigma2, double alpha, double beta);
void recover_time_domain(double complex* Y_eq, double complex* time_signal, int N);
void demodulate_slot(double complex* time_signal, double complex* X_hat, int slot);
void compute_true_channel(double complex* h, double complex* H_true, int size, int num_paths);


// Generate QPSK symbol: (±1 ± j) / sqrt(2)
double complex generate_qpsk() {
    double re = (rand() % 2 ? 1 : -1) / sqrt(2);
    double im = (rand() % 2 ? 1 : -1) / sqrt(2);
    return re + im * I;
}

// Scaled IFFT: x[n] = (1/sqrt(N)) * sum_{k=0}^{N-1} X[k] * exp(j*2*pi*k*n / N)
void ifft(double complex* X, double complex* x, int N) {
    for (int n = 0; n < N; n++) {
        x[n] = 0;
        for (int k = 0; k < N; k++) {
            double angle = 2 * M_PI * k * n / N;
            x[n] += X[k] * (cos(angle) + I * sin(angle));
        }
        x[n] /= sqrt(N);
    }
}

// Scaled FFT: X[k] = (1/sqrt(N)) * sum_{n=0}^{N-1} x[n] * exp(-j*2*pi*k*n / N)
void fft(double complex* x, double complex* X, int N) {
    for (int k = 0; k < N; k++) {
        X[k] = 0;
        for (int n = 0; n < N; n++) {
            double angle = -2 * M_PI * k * n / N;
            X[k] += x[n] * (cos(angle) + I * sin(angle));
        }
        X[k] /= sqrt(N);
    }
}

// IFFT with 1/N scaling for interpolation
void ifft_interp(double complex* X, double complex* x, int N) {
    for (int n = 0; n < N; n++) {
        x[n] = 0;
        for (int k = 0; k < N; k++) {
            double angle = 2 * M_PI * k * n / N;
            x[n] += X[k] * (cos(angle) + I * sin(angle));
        }
        x[n] /= N;
    }
}

// FFT with 1/N scaling for interpolation
void fft_interp(double complex* x, double complex* X, int N, int M) {
    for (int k = 0; k < N; k++) {
        X[k] = 0;
        for (int t = 0; t < M; t++) {
            double angle = -2 * M_PI * k * t / N;
            X[k] += x[t] * (cos(angle) + I * sin(angle));
        }
        X[k] /= N;
    }
}

// Generate OFDM/TDM frame with Chu pilot in slot K-1
void generate_ofdm_tdm_frame(double complex frame[Nc], double complex X[K][Nm], double complex* pilot_symbols) {
    double complex s_slot[Nm];
    for (int k = 0; k < K; k++) {
        for (int m = 0; m < Nm; m++) {
            if (k == K - 1) {
                X[k][m] = pilot_symbols[m];
            }
            else {
                X[k][m] = generate_qpsk();
            }
        }
        ifft(X[k], s_slot, Nm);
        for (int t = 0; t < Nm; t++) {
            frame[k * Nm + t] = s_slot[t];
        }
    }
}

// Add Nm-sample cyclic prefix
void add_cyclic_prefix(double complex frame[Nc], double complex tx_signal[Nc + GI]) {
    for (int t = 0; t < GI; t++) {
        tx_signal[t] = frame[Nc - GI + t];
    }
    for (int t = 0; t < Nc; t++) {
        tx_signal[GI + t] = frame[t];
    }
}

// Apply HPA soft limiter with saturation amplitude A
void apply_hpa(double complex* tx_signal, int len, double A) {
    for (int t = 0; t < len; t++) {
        double mag = cabs(tx_signal[t]);
        if (mag > A) {
            tx_signal[t] = A * tx_signal[t] / mag;
        }
    }
}

// Generate complex Gaussian noise
double complex generate_noise(double sigma) {
    double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
    double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
    double re = sigma * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
    double im = sigma * sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
    return re + im * I;
}

/* Generate L - path Rayleigh fading channel without correlation
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
*/
// Generate L-path Rayleigh fading channel with correlation and slow fading
void generate_channel(double complex* h, double rho, double complex* h_prev) {
    double sig = sqrt(0.5 / L);
    for (int l = 0; l < L; l++) {
        double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
        double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
        double w_re = sig * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
        double w_im = sig * sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
        double complex w = w_re + I * w_im;
        h[l] = rho * h_prev[l] + sqrt(1 - rho * rho) * w;
        h_prev[l] = h[l];
    }
}


// Apply channel and noise
void apply_channel(double complex* tx_signal, double complex* r, double complex* h, double sigma) {
    int len_tx = Nc + GI;
    for (int t = 0; t < len_tx; t++) {
        r[t] = 0;
        for (int l = 0; l < L && t - l >= 0; l++) {
            r[t] +=  h[l] * tx_signal[t - l];
        }
        r[t] += generate_noise(sigma);
    }
}

// Remove GI and apply FFT
void receive_frame(double complex* r, double complex* R) {
    double complex r_frame[Nc];
    for (int t = 0; t < Nc; t++) {
        r_frame[t] = r[GI + t];
    }
    fft(r_frame, R, Nc);
}


// Channel estimation with noise power estimation
void channel_estimation(double complex* r, double complex* H_est, double complex* pilot_symbols, double* sigma2_e) {
    //Time Domain Filtering on GI
    static double complex prev_r[Nm] = { 0 };
    static int first_call = 1;
    double gamma = GAMMA;
    double complex r_g[Nm];

    if (first_call) {
        for (int t = 0; t < Nm; t++) {
            r_g[t] = r[t];
            prev_r[t] = r_g[t];
        }
        first_call = 0;
    }
    else {
        for (int t = 0; t < Nm; t++) {
            r_g[t] = gamma * r[t] + (1 - gamma) * prev_r[t];
            prev_r[t] = r_g[t];
        }
    }

    // Compute R_g(q)
    double complex R_g[Nm];
    for (int q = 0; q < Nm; q++) {
        R_g[q] = 0;
        for (int t = 0; t < Nm; t++) {
            double angle = -2 * M_PI * q * t / Nm;
            R_g[q] += r_g[t] * (cos(angle) + I * sin(angle));
        }
    }

    // Compute H_g(q) = R_g(q) / P(q)
    double complex H_g[Nm];
    for (int q = 0; q < Nm; q++) {
        H_g[q] = R_g[q] / pilot_symbols[q];
    }

    // Compute h(t) using ifft_interp
    double complex h[Nm];
    ifft_interp(H_g, h, Nm);

    // Compute H_e(n) using fft_interp
    fft_interp(h, H_est, Nc, Nm);

    // Noise power estimation: N_e(q) = R_g(q) - H_e(q) * P(q)
    double complex N_e[Nm];
    double sum_power = 0.0;
    for (int q = 0; q < Nm; q++) {
        N_e[q] = R_g[q] - H_g[q] * pilot_symbols[q];
        sum_power += creal(N_e[q]) * creal(N_e[q]) + cimag(N_e[q]) * cimag(N_e[q]);
    }
    *sigma2_e = sum_power / Nm;
}

// MMSE Frequency-Domain Equalization
void mmse_fde(double complex* R, double complex* Y_eq, double complex* H_est, double sigma2, double alpha, double beta) {
    double beta_linear = pow(10, beta / 20.0);
    for (int m = 0; m < Nc; m++) {
        double complex h = H_est[m];
        double h_mag2 = creal(h) * creal(h) + cimag(h) * cimag(h);
        double denom = h_mag2 * (1 - exp(-beta_linear * beta_linear)) + 2 * sigma2;
        double complex w = (alpha * conj(h)) / denom;
        Y_eq[m] = R[m] * w;
    }
}

// Recover time-domain signal using IFFT
void recover_time_domain(double complex* Y_eq, double complex* time_signal, int N) {
    ifft(Y_eq, time_signal, N);
}

// Demodulate each slot using FFT
void demodulate_slot(double complex* time_signal, double complex* X_hat, int slot) {
    double complex slot_time[Nm];
    for (int t = 0; t < Nm; t++) {
        slot_time[t] = time_signal[slot * Nm + t];
    }
    fft(slot_time, X_hat, Nm);
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
int main() {
    srand(time(NULL));
    // Case selection (1 to 2)
    int case_type = 1;

    double complex X[K][Nm];
    double complex pilot_symbols[Nm];
    double complex frame[Nc];
    double complex h[L];
    double complex h_prev[L];
    double complex tx_signal[Nc + GI];
    double complex r[Nc + GI];
    double complex H_est[Nc];
    double sigma2_e;
    double complex R[Nc];
    double complex Y_eq[Nc];
    double complex time_signal[Nc];
    double complex X_hat[Nm];
    double complex H_true[Nc];
    double complex H_Eq[Nc];
    double complex sigma2_Eq;




     // Generate Chu pilot sequence
    for (int m = 0; m < Nm; m++) {
        double angle = M_PI * m * m / Nc;
        pilot_symbols[m] = cos(angle) + I * sin(angle);
    }

    // Compute correlation coefficient rho
    //double rho = j0(2 * M_PI * FD_TS);
    double rho = 0.9999999013;
    // Initialize first channel
    double sig = sqrt(0.5 / L);
    for (int l = 0; l < L; l++) {
        double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
        double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
        double z0 = sig * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2);
        double z1 = sig * sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
        h_prev[l] = z0 + I * z1;
    }



    // HPA parameters
    double beta_dB = BETA_DB;
    double A = sqrt(pow(10, beta_dB / 10.0));
    double beta_linear = pow(10, beta_dB / 20.0);
    double erfc_beta = erfc(beta_linear);
    double alpha = 1 - exp(-beta_linear * beta_linear) + (sqrt(M_PI) * beta_linear / 2.0) * erfc_beta;

    // Open file for writing
    FILE* fp = fopen("ber_data.txt", "w");
    if (fp == NULL) {
        printf("Error opening file!\n");
        return 1;
    }

    // Eb/N0 sweep
    for (double Eb_N0_dB = 5; Eb_N0_dB <= 30; Eb_N0_dB += 5) {
        double Eb_N0 = pow(10, Eb_N0_dB / 10.0);
        double number_of_bits = 2.0 * (K - 1) * Nm;
        double E_frame_avg = (Nc + GI) * P;
        double Eb = E_frame_avg / number_of_bits;
        double N0 = Eb / Eb_N0;
        double sigma2 = N0 / 2.0;
        double sigma = sqrt(sigma2);
        double total_ber = 0.0;
        int num_frames = 1000;

        for (int frame_idx = 0; frame_idx < num_frames; frame_idx++) {
            generate_ofdm_tdm_frame(frame, X, pilot_symbols);
            add_cyclic_prefix(frame, tx_signal);
            apply_hpa(tx_signal, Nc + GI, A);
            for (int t = 0; t < Nc + GI; t++) {
                tx_signal[t] *= sqrt(P);
            }
            //generate_channel(h);
            generate_channel(h, rho, h_prev);
            apply_channel(tx_signal, r, h, sigma);
            receive_frame(r, R);
            channel_estimation(r, H_est, pilot_symbols, &sigma2_e);
            compute_true_channel(h, H_true, Nc, L);
           
            if (case_type == 1) {
                for (int m = 0; m < Nc; m++) {
                    H_Eq[m] = H_est[m];
                }
                sigma2_Eq = sigma2_e;
            }
            else if (case_type == 2) {
                for (int m = 0; m < Nc; m++) {
                    H_Eq[m] = H_true[m];
                }
                sigma2_Eq = sigma2;
            }
            mmse_fde(R, Y_eq, H_Eq, sigma2_Eq, alpha, beta_linear); // Use estimated noise power
            recover_time_domain(Y_eq, time_signal, Nc);

            int bit_errors = 0;
            for (int k = 0; k < K - 1; k++) {
                demodulate_slot(time_signal, X_hat, k);
                for (int m = 0; m < Nm; m++) {
                    double complex Y = X_hat[m];
                    double re = creal(Y);
                    double im = cimag(Y);
                    double complex X_decided = (re > 0 ? 1 : -1) / sqrt(2) + I * (im > 0 ? 1 : -1) / sqrt(2);
                    double complex X_trans = X[k][m];
                    if (creal(X_trans) * creal(X_decided) < 0) bit_errors++;
                    if (cimag(X_trans) * cimag(X_decided) < 0) bit_errors++;
                }
            }
            total_ber += (double)bit_errors / number_of_bits;
        }
        double avg_ber = total_ber / num_frames;
        printf("Eb/N0 = %.1f dB, BER = %e\n", Eb_N0_dB, avg_ber);
        fprintf(fp, "%d %.1f %e\n", case_type, Eb_N0_dB, avg_ber);
    }

    // Close the file
    fclose(fp);
    return 0;
}