#include <iostream>
#include <vector>
#include <string>
#include <gmpxx.h>
#include <chrono>

#include "./lib/Struct.h"
#include "./lib/RLWE.h"
#include "./lib/Group.h"
#include "./lib/Authentic.h"
#include "./lib/Control.h"

using namespace std;

// ==================== Hyper Parameter ==================== //
const int poly_degree = (int)powl(2, 14);
const int plain_bits = 42;
const int cipher_bits = 256;
const int group_bits = 3072;

const double r_scale = 0.0001;
const double s_scale = 0.001;
// ========================================================= //


int main()
{   
    // ===================== Check Runtime ==================== //
    auto stc = chrono::high_resolution_clock::now();
    auto edc = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::nanoseconds>(edc - stc);
    double run_time = duration.count() / 1000000;
    // ======================================================== //

    
    stc = chrono::high_resolution_clock::now();
    // ====================== Parameter ====================== //
    mpz_class plain_mod;
    prime_handler::find_ntt_prime(poly_degree, plain_bits, plain_mod);
    mpz_class psi_p;
    prime_handler::find_ntt_root(poly_degree, plain_mod, psi_p);

    mpz_class cipher_mod;
    prime_handler::find_ntt_prime(poly_degree, cipher_bits, cipher_mod);
    mpz_class psi_c;
    prime_handler::find_ntt_root(poly_degree, cipher_mod, psi_c);

    mpz_class group_mod, group_gen;
    prime_handler::find_schnorr_prime(cipher_mod, group_bits,  group_mod);
    prime_handler::find_schnorr_gen(cipher_mod, group_mod, group_gen);
    // ======================================================== //
    edc = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::nanoseconds>(edc - stc);
    run_time = duration.count() / 1000000;
    cout << "Parameter setting time: " << run_time << "ms" << endl;
    // cout << plain_mod << endl;
    // cout << cipher_mod << endl;
    // cout << group_mod << endl;
    // cout << group_gen << endl;



    stc = chrono::high_resolution_clock::now();
    // ================= Control (ARX) matrix ================== //
    // Original(real number field) ARX matrix representation
    vector<double> P1({-0.3844, 6.5970}); // H*(F-R*H)^3*G
    vector<double> P2({1.9504, -32.9119}); // H*(F-R*H)^2*G
    vector<double> P3({-2.8663, 50.1128}); // H*(F-R*H)^1*G
    vector<double> P4({1.3015, -23.9663}); // H*(F-R*H)^0*G
    vector<double> Q1({-0.0945}); // H*(F-R*H)^3*R
    vector<double> Q2({0.5509}); // H*(F-R*H)^2*R
    vector<double> Q3({-1.3776}); // H*(F-R*H)^1*R
    vector<double> Q4({1.9065}); // H*(F-R*H)^0*R

    // Change to integer field
    vector<int64_t> P1i(2); P1i[0] = (int64_t)(P1[0] / s_scale); P1i[1] = (int64_t)(P1[1] / s_scale);
    vector<int64_t> P2i(2); P2i[0] = (int64_t)(P2[0] / s_scale); P2i[1] = (int64_t)(P2[1] / s_scale);
    vector<int64_t> P3i(2); P3i[0] = (int64_t)(P3[0] / s_scale); P3i[1] = (int64_t)(P3[1] / s_scale);
    vector<int64_t> P4i(2); P4i[0] = (int64_t)(P4[0] / s_scale); P4i[1] = (int64_t)(P4[1] / s_scale);
    vector<int64_t> Q1i(1); Q1i[0] = (int64_t)(Q1[0] / s_scale);
    vector<int64_t> Q2i(1); Q2i[0] = (int64_t)(Q2[0] / s_scale);
    vector<int64_t> Q3i(1); Q3i[0] = (int64_t)(Q3[0] / s_scale);
    vector<int64_t> Q4i(1); Q4i[0] = (int64_t)(Q4[0] / s_scale);

    // Integer matrix
    vector<vector<int64_t>> P({P1i, P2i, P3i, P4i}); 
    vector<vector<int64_t>> Q({Q1i, Q2i, Q3i, Q4i}); 
    // ========================================================= //
    edc = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::nanoseconds>(edc - stc);
    run_time = duration.count() / 1000000;
    cout << "Control (ARX) matrix setting time: " << run_time << "ms" << endl;


    stc = chrono::high_resolution_clock::now();
    // ============ Encrypted (ARX) controller ready =========== //
    arx* arx_ctrl = new arx();

    // secret key generation
    poly* sk = new poly(poly_degree);
    random_handler::secret_key(sk);

    // P and Q matrix encryption and link pointer in arx class inner variable
    vector<int64_t> pod_matrix(poly_degree, 0LL);
    for(int i = 0; i < 4; i++)
    {
        pod_matrix[0] = P[i][0];
        pod_matrix[1] = P[i][1];

        poly* pre_packed_data = new poly(poly_degree);
        batch_encoder::encode(pod_matrix, plain_mod, psi_p, pre_packed_data);

        poly* pre_plaintext = new poly(poly_degree);
        poly_handler::pack_2_plain(pre_packed_data, pre_plaintext);

        cipher* ciphertext = new cipher(pre_plaintext->ring_dim, plain_mod, cipher_mod, psi_p, psi_c);
        crypto_handler::encrypt(pre_plaintext, sk, ciphertext);

        arx_ctrl->P_y[i] = ciphertext;

        delete pre_packed_data;
        delete pre_plaintext;
    }

    pod_matrix[1] = 0;
    for(int i = 0; i < 4; i++)
    {
        pod_matrix[0] = Q[i][0];

        poly* pre_packed_data = new poly(poly_degree);
        batch_encoder::encode(pod_matrix, plain_mod, psi_p, pre_packed_data);

        poly* pre_plaintext = new poly(poly_degree);
        poly_handler::pack_2_plain(pre_packed_data, pre_plaintext);

        cipher* ciphertext = new cipher(pre_plaintext->ring_dim, plain_mod, cipher_mod, psi_p, psi_c);
        crypto_handler::encrypt(pre_plaintext, sk, ciphertext);

        arx_ctrl->Q_u[i] = ciphertext;

        delete pre_packed_data;
        delete pre_plaintext;
    }

    pod_matrix[0] = 0;

    poly* pre_packed_data = new poly(poly_degree);
    batch_encoder::encode(pod_matrix, plain_mod, psi_p, pre_packed_data);

    poly* pre_plaintext = new poly(poly_degree);
    poly_handler::pack_2_plain(pre_packed_data, pre_plaintext);

    cipher* ciphertext = new cipher(pre_plaintext->ring_dim, plain_mod, cipher_mod, psi_p, psi_c);
    crypto_handler::encrypt(pre_plaintext, sk, ciphertext);

    arx_ctrl->temp = ciphertext;
    arx_ctrl->zero_set();

    delete pre_packed_data;
    delete pre_plaintext;
    // ========================================================= //
    edc = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::nanoseconds>(edc - stc);
    run_time = duration.count() / 1000000;
    cout << "controller encryption time: " << run_time << "ms" << endl;



    stc = chrono::high_resolution_clock::now();
    // ==================== Authenticator test ==================== //
    authentic* auth = new authentic(poly_degree, cipher_mod, group_mod, group_gen);
    auth->make_ekf(arx_ctrl->P_y, arx_ctrl->Q_u);
    // ============================================================ //
    edc = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::nanoseconds>(edc - stc);
    run_time = duration.count() / 1000000;
    cout << "Authenticator test time: " << run_time << "ms" << endl;



    cout << "!! Simulation test start !!" << endl;
    // ===================== Simulation test ====================== //
    // make plant class
    vector<double> x_init({0, 0, 0.1, 0});
    plant* plt = new plant(x_init);

    // variables for signal encryption
    poly* plaintext = new poly(poly_degree);
    poly* packed_data = new poly(poly_degree);
    vector<int64_t> result_data(poly_degree, 0LL);
    int64_t temp_u;
    vector<double> real_u(1);

    // variables for singal transaction
    vector<int64_t> plant_output(poly_degree, 0LL);
    vector<int64_t> control_input(poly_degree, 0LL);
    cipher* plt_out = new cipher(plaintext->ring_dim, plain_mod, cipher_mod, psi_p, psi_c);
    cipher* ctrl_in = new cipher(plaintext->ring_dim, plain_mod, cipher_mod, psi_p, psi_c);

    // ready to authentic
    mpz_class previous_pf;
    std::vector<poly*> initial_crypto_state_stack(16);
    for(int i = 0; i < 8; i++)
    {
        if(i < 4)
        {
            initial_crypto_state_stack[2 * i] = new poly(poly_degree);
            initial_crypto_state_stack[2 * i]->coeff = arx_ctrl->mem_y_pre[i]->ciphertext[0]->coeff;
            initial_crypto_state_stack[2 * i + 1]= new poly(poly_degree);
            initial_crypto_state_stack[2 * i + 1]->coeff = arx_ctrl->mem_y_pre[i]->ciphertext[1]->coeff;
        }
        else
        {
            initial_crypto_state_stack[2 * i] = new poly(poly_degree);
            initial_crypto_state_stack[2 * i]->coeff = arx_ctrl->mem_u_pre[i - 4]->ciphertext[0]->coeff;
            initial_crypto_state_stack[2 * i + 1]= new poly(poly_degree);
            initial_crypto_state_stack[2 * i + 1]->coeff = arx_ctrl->mem_u_pre[i - 4]->ciphertext[1]->coeff;
        }
    }
    poly* initial_crypto_state = poly_handler::poly_recur_concat(initial_crypto_state_stack);
    group_handler::group_dot(auth->g_r_1, initial_crypto_state, previous_pf);
    for(auto& factor : initial_crypto_state_stack)
    {
        if(factor != nullptr)
        {
            delete factor;
            factor = nullptr;
        }
    }
    initial_crypto_state_stack.clear();
    delete initial_crypto_state;

    // authentic pass check
    bool pass = false;

    int iter = 20;

    // CSV save
    FILE* ps = fopen("date(y1_y2_u_vc_enct_vct).csv", "w");
    if(ps == nullptr)
    {
        cout << "Error: Unable to open file" << endl;
        exit(1);
    }
    fprintf(ps, "y1,y2,u,vc,enct,vct\n");

    auto enc_stc = std::chrono::high_resolution_clock::now();
    auto enc_edc = std::chrono::high_resolution_clock::now();
    auto enc_duration = chrono::duration_cast<chrono::nanoseconds>(edc - stc);
    double enc_run_time = duration.count() / 1000000;

    auto vc_stc = std::chrono::high_resolution_clock::now();
    auto vc_edc = std::chrono::high_resolution_clock::now();
    auto vc_duration = chrono::duration_cast<chrono::nanoseconds>(edc - stc);
    double vc_run_time = duration.count() / 1000000;

    for(int i = 0; i < iter; i++)
    {
        // iteration time check
        stc = std::chrono::high_resolution_clock::now();

        // get plant output and calculation control input
        plt->output();
        arx_ctrl->calc();

        // get control input
        crypto_handler::decrypt(arx_ctrl->calc_res, sk, plaintext);
        batch_encoder::decode(plaintext, plain_mod, psi_p, result_data);
        temp_u = result_data[0] + result_data[1];

        real_u[0] = (double)temp_u * r_scale * s_scale;

        enc_stc = std::chrono::high_resolution_clock::now();
        // plant output encryption
        plant_output[0] = (int64_t)(plt->y[0] / r_scale);
        plant_output[1] = (int64_t)(plt->y[1] / r_scale);
        batch_encoder::encode(plant_output, plain_mod, psi_p, plaintext);
        crypto_handler::encrypt(plaintext, sk, plt_out);

        // re-encryption
        control_input[0] = (int64_t)(real_u[0] / r_scale);   
        batch_encoder::encode(control_input, plain_mod, psi_p, plaintext);
        crypto_handler::encrypt(plaintext, sk, ctrl_in);
        enc_edc = std::chrono::high_resolution_clock::now();
        enc_duration = chrono::duration_cast<chrono::nanoseconds>(enc_edc - enc_stc);
        enc_run_time = enc_duration.count() / 1000000;

        // plant state update and controller memory update
        plt->state_update(real_u);
        arx_ctrl->mem_update(plt_out, ctrl_in);
        if(i >= 50 && i < 70)
        {
            arx_ctrl->mem_y_new[3]->ciphertext[0]->coeff[1] += 1;
        }

        vc_stc = std::chrono::high_resolution_clock::now();
        // Authentication
        auth->generate_proof(arx_ctrl->mem_y_new, arx_ctrl->mem_u_new, arx_ctrl->mem_y_pre, arx_ctrl->mem_u_pre);
        pass = auth->verifying_proof(plt_out, ctrl_in, arx_ctrl->calc_res, previous_pf);
        vc_edc = std::chrono::high_resolution_clock::now();
        vc_duration = chrono::duration_cast<chrono::nanoseconds>(vc_edc - vc_stc);
        vc_run_time = vc_duration.count() / 1000000;

        // save data
        fprintf(ps, "%lf,%lf,%lf,%d,%lf,%lf\n",plt->y[0], plt->y[1], real_u[0], (int)pass, enc_run_time, vc_run_time);

        // debug print
        edc = std::chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::nanoseconds>(edc - stc);
        run_time = duration.count() / 1000000;
        cout << "iter [" << i+1 << "] | time: " << run_time << "ms" << " | u: " << real_u[0] << " | y: " << plt->y[0] << ", " << plt->y[1] << endl;
        cout << "pass: " << pass << " | enc run time: " << enc_run_time << "ms | vc run time: " << vc_run_time << "ms" << endl;
    }
    // ============================================================ //



    // ====================== Free variables ====================== //
    delete sk;

    delete arx_ctrl;

    delete auth;

    delete plt;
    delete plaintext;
    delete packed_data;
    delete plt_out;
    delete ctrl_in;

    fclose(ps);
    // ============================================================ //

    return 0;
}
