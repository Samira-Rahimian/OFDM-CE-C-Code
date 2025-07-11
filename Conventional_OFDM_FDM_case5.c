#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>



//This code implements Cases 5: 
// 5- OFDM with TDM Pilot
//No case selection is required.

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

//Nm= Nc/PILOT_SPACING=16 // num of pilot subcarriers
#define Nc 256     // Subcarriers per OFDM symbol
#define K 16      // Number of OFDM symbols
#define GI 16     // Guard interval per symbol
#define L 8       // Channel paths
#define BETA_DB 3.0 // HPA saturation in dB
#define P 1.0     // Transmit power
#define PILOT_SPACING 16 // Pilot every 16th subcarrier

// QPSK symbol generation
double complex generate_qpsk() {
    double re = (rand() % 2 ? 1 : -1) / sqrt(2);
    double im = (rand() % 2 ? 1 : -1) / sqrt(2);
    return re + im * I;
}

// Scaled IFFT
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

// Scaled FFT
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


// Generate OFDM symbol with FDM pilots
void generate_ofdm_symbol(double complex X[Nc], double complex s[Nc + GI]) {
    for (int m = 0; m < Nc; m++) {
        if (m % PILOT_SPACING == 0) {
            double angle = M_PI * m * m / Nc;
            X[m] = cos(angle) + I * sin(angle); // Pilot
        }
        else {
            X[m] = generate_qpsk(); // Data
        }
    }
    double complex s_no_cp[Nc];
    ifft(X, s_no_cp, Nc);
    for (int t = 0; t < GI; t++) {
        s[t] = s_no_cp[Nc - GI + t];
    }
    for (int t = 0; t < Nc; t++) {
        s[GI + t] = s_no_cp[t];
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

// Generate noise
double complex generate_noise(double sigma) {
    double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
    double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
    return sigma * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2) +
        I * sigma * sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
}

// Generate channel
void generate_channel(double complex* h) {
    double sigma = sqrt(0.5 / L);
    for (int l = 0; l < L; l++) {
        double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
        double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
        h[l] = sigma * sqrt(-2 * log(u1)) * cos(2 * M_PI * u2) +
            I * sigma * sqrt(-2 * log(u1)) * sin(2 * M_PI * u2);
    }
}

// Apply channel
void apply_channel(double complex* tx_signal, double complex* r, double complex* h, double sigma, int len) {
    for (int t = 0; t < len; t++) {
        r[t] = 0;
        for (int l = 0; l < L && t - l >= 0; l++) {
            r[t] += sqrt(2) * h[l] * tx_signal[t - l];
        }
        r[t] += generate_noise(sigma);
    }
}

//Channel estimation with high-resolution interpolation
void channel_estimation_fdm(double complex* R, double complex* H_est, int pilot_spacing, double* sigma2_e) {
    int Nm = Nc / pilot_spacing; // Nm = 256 / 16 = 16
    double complex H_pilot[Nm];

    // Compute channel estimates at pilot positions
    for (int q = 0; q < Nm; q++) {
        int m = q * pilot_spacing;
        double angle = M_PI * m * m / Nc;
        double complex P_m = cos(angle) + I * sin(angle);
        H_pilot[q] = R[m] / P_m;
    }

    // Compute channel impulse response h(t) using Nm-point IFFT
    double complex h[Nm];
    ifft_interp(H_pilot, h, Nm);
    // Compute interpolated channel estimates H_est[n] for all Nc subcarriers
    fft_interp(h, H_est, Nc, Nm);
    // Noise power estimation: N_e(q) = R(q) - H_est(q) * P(q)
    double complex N_e[Nm];
    double sum_power = 0.0;
    for (int q = 0; q < Nm; q++) {
        int n = q * pilot_spacing; // Map pilot index to full frequency index
        double angle = M_PI * n * n / Nc;
        double complex P_n = cos(angle) + I * sin(angle);
        N_e[q] = R[n] - H_est[n] * P_n;
        sum_power += creal(N_e[q]) * creal(N_e[q]) + cimag(N_e[q]) * cimag(N_e[q]);
    }
    *sigma2_e = sum_power / Nm;
}


// MMSE FDE
void mmse_fde(double complex* R, double complex* Y_eq, double complex* H_est, double sigma2, double alpha, double beta) {
    for (int m = 0; m < Nc; m++) {
        double complex h = H_est[m];
        double h_mag2 = creal(h) * creal(h) + cimag(h) * cimag(h);
        double denom = h_mag2 * (1 - exp(-beta * beta)) + 2 * sigma2;
        double complex w = (alpha * conj(h)) / denom;
        Y_eq[m] = R[m] * w;
    }
}

int main() {
    srand(time(NULL));

    int case_type = 5;

    int frame_len = K * (Nc + GI);
    double complex X[K][Nc];
    double complex tx_signal[frame_len];
    double complex r[frame_len];
    double complex h[L];
    double complex R[Nc];
    double complex H_est[Nc];
    double complex Y_eq[Nc];
    double sigma2_e;

    double beta_dB = BETA_DB;
    double A = sqrt(pow(10, beta_dB / 10.0));
    double beta = pow(10, beta_dB / 20.0);
    double erfc_beta = erfc(beta);
    double alpha = 1 - exp(-beta * beta) + (sqrt(M_PI) * beta / 2.0) * erfc_beta;

    FILE* fp = fopen("ber_data.txt", "w");
    for (double Eb_N0_dB = 5; Eb_N0_dB <= 30; Eb_N0_dB += 5) {
        double Eb_N0 = pow(10, Eb_N0_dB / 10.0);
        double num_data_subs = Nc - (Nc / PILOT_SPACING);
        double num_bits = 2.0 * K * num_data_subs;
        double E_frame_avg = K * (Nc + GI) * P;
        double Eb = E_frame_avg / num_bits;
        double N0 = Eb / Eb_N0;
        double sigma2 = N0 / 2.0;
        double sigma = sqrt(sigma2);
        double total_ber = 0.0;
        int num_frames = 1000;

        for (int frame_idx = 0; frame_idx < num_frames; frame_idx++) {
            // Transmit
            for (int k = 0; k < K; k++) {
                double complex s[Nc + GI];
                generate_ofdm_symbol(X[k], s);
                apply_hpa(s, Nc + GI, A);
                for (int t = 0; t < Nc + GI; t++) {
                    tx_signal[k * (Nc + GI) + t] = s[t] * sqrt(P);
                }
            }
            generate_channel(h);
            apply_channel(tx_signal, r, h, sigma, frame_len);

            // Receive and process per symbol
            int bit_errors = 0;
            for (int k = 0; k < K; k++) {
                double complex r_symbol[Nc];
                for (int t = 0; t < Nc; t++) {
                    r_symbol[t] = r[k * (Nc + GI) + GI + t];
                }
                fft(r_symbol, R, Nc);
                channel_estimation_fdm(R, H_est, PILOT_SPACING, &sigma2_e);
                mmse_fde(R, Y_eq, H_est, sigma2_e, alpha, beta);
                for (int m = 0; m < Nc; m++) {
                    if (m % PILOT_SPACING != 0) {
                        double re = creal(Y_eq[m]);
                        double im = cimag(Y_eq[m]);
                        double complex X_decided = (re > 0 ? 1 : -1) / sqrt(2) +
                            I * (im > 0 ? 1 : -1) / sqrt(2);
                        double complex X_trans = X[k][m];
                        if (creal(X_trans) * creal(X_decided) < 0) bit_errors++;
                        if (cimag(X_trans) * cimag(X_decided) < 0) bit_errors++;
                    }
                }
            }
            total_ber += (double)bit_errors / num_bits;
        }
        double avg_ber = total_ber / num_frames;
        printf("Eb/N0 = %.1f dB, BER = %e\n", Eb_N0_dB, avg_ber);
        fprintf(fp, "%d %.1f %e\n", case_type, Eb_N0_dB, avg_ber);
    }
    fclose(fp);
    return 0;
}