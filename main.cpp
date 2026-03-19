#include <iostream>
#include <vector>
#include <string>
#include <gmpxx.h>
#include <chrono>

#include "./lib/Struct.h"
#include "./lib/RLWE.h"
#include "./lib/Group.h"
#include "./lib/Authentic.h"
// #include "./lib/Control.h"

using namespace std;

// ==================== Hyper Parameter ==================== //
const int poly_degree = (int)powl(2, 14);
const int plain_bits = 42;
const int cipher_bits = 256;
const int group_bits = 3072;
// ========================================================= //

int main()
{
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



    // ====================== Encryption ====================== //
    vector<int64_t> pod_matrix(poly_degree, 2LL);
    poly* packed_data = new poly(poly_degree);
    poly* plaintext = new poly(poly_degree);
    cipher* ciphertext = new cipher(poly_degree, plain_mod, cipher_mod, psi_p, psi_c);
    poly* sk = new poly(poly_degree);
    random_handler::secret_key(sk);
    
    batch_encoder::encode(pod_matrix, plain_mod, psi_p, packed_data);
    poly_handler::pack_2_plain(packed_data, plaintext);
    crypto_handler::encrypt(plaintext, sk, ciphertext);
    // ======================================================== //


    
    // ==================== Representation ==================== //
    vr_cipher* vrciphertext = new vr_cipher(poly_degree, plain_mod, cipher_mod);
    format_transform_handler::cipher_2_vr_cipher(ciphertext, vrciphertext);

    mr_cipher* mrciphertext = new mr_cipher(poly_degree, plain_mod, cipher_mod);
    format_transform_handler::cipher_2_mr_cipher(ciphertext, mrciphertext);

    vr_cipher* res_vrciphertext = new vr_cipher(poly_degree, plain_mod, cipher_mod);
    crypto_handler::eval_mvr_mul(mrciphertext, vrciphertext, res_vrciphertext);
    // ======================================================== //



    // ==================== On group structure ==================== //
    poly* concat_vrciphertext = new poly(2 * poly_degree);
    poly_handler::poly_concat(vrciphertext->ciphertext_1, vrciphertext->ciphertext_2, concat_vrciphertext);

    gvec* gociphertext = new gvec(concat_vrciphertext->ring_dim, group_mod, group_gen);
    group_handler::poly_2_gvec(concat_vrciphertext, gociphertext);
    // ============================================================ //



    // ======================== De-flag ======================= //
    cipher* re_cipher = new cipher(poly_degree, plain_mod, cipher_mod, psi_p, psi_c);
    format_transform_handler::vr_cipher_2_cipher(res_vrciphertext, re_cipher);
    // ======================================================== //



    // ==================== Authenticator test ==================== //
    authentic* auth = new authentic(poly_degree, cipher_mod, group_mod, group_gen);
    auth->make_ekf();
    // ============================================================ //



    // ======================== Decrypt ======================= //
    poly* res_plain = new poly(poly_degree);
    poly* res_pack = new poly(poly_degree);
    vector<int64_t> res_matrix(poly_degree, 0LL);

    crypto_handler::decrypt(re_cipher, sk, res_plain);
    poly_handler::plain_2_pack(res_plain, res_pack);
    batch_encoder::decode(res_pack, plain_mod, psi_p, res_matrix);
    // ======================================================== //

    cout << res_matrix[0] << endl;

    delete packed_data;
    delete plaintext;
    delete ciphertext;
    delete sk;

    delete vrciphertext;
    delete mrciphertext;
    delete res_vrciphertext;

    delete concat_vrciphertext;
    delete gociphertext;

    delete auth;
    
    delete re_cipher;

    delete res_plain;
    delete res_pack;

    return 0;
}

// Previous Description
// int main()
// {
//     // ================================================================================ //
//     //                               Setting for Control                                //
//     // ================================================================================ //

//     vector<double> P1({-0.3844, 6.5970}); // H*(F-R*H)^3*G
//     vector<double> P2({1.9504, -32.9119}); // H*(F-R*H)^2*G
//     vector<double> P3({-2.8663, 50.1128}); // H*(F-R*H)^1*G
//     vector<double> P4({1.3015, -23.9663}); // H*(F-R*H)^0*G

//     vector<double> Q1({-0.0945}); // H*(F-R*H)^3*R
//     vector<double> Q2({0.5509}); // H*(F-R*H)^2*R
//     vector<double> Q3({-1.3776}); // H*(F-R*H)^1*R
//     vector<double> Q4({1.9065}); // H*(F-R*H)^0*R

//     double r_scale = 0.0001;
//     double s_scale = 0.001;

//     vector<int64_t> P1i(2); P1i[0] = (int64_t)(P1[0] / s_scale); P1i[1] = (int64_t)(P1[1] / s_scale);
//     vector<int64_t> P2i(2); P2i[0] = (int64_t)(P2[0] / s_scale); P2i[1] = (int64_t)(P2[1] / s_scale);
//     vector<int64_t> P3i(2); P3i[0] = (int64_t)(P3[0] / s_scale); P3i[1] = (int64_t)(P3[1] / s_scale);
//     vector<int64_t> P4i(2); P4i[0] = (int64_t)(P4[0] / s_scale); P4i[1] = (int64_t)(P4[1] / s_scale);

//     vector<int64_t> Q1i(1); Q1i[0] = (int64_t)(Q1[0] / s_scale);
//     vector<int64_t> Q2i(1); Q2i[0] = (int64_t)(Q2[0] / s_scale);
//     vector<int64_t> Q3i(1); Q3i[0] = (int64_t)(Q3[0] / s_scale);
//     vector<int64_t> Q4i(1); Q4i[0] = (int64_t)(Q4[0] / s_scale);

//     vector<vector<int64_t>> P({P1i, P2i, P3i, P4i}); 
//     vector<vector<int64_t>> Q({Q1i, Q2i, Q3i, Q4i}); 

//     vector<double> x_init({0, 0, 0.1, 0});

//     // ================================================================================ //
//     //                               Setting for Encrypt                                //
//     // ================================================================================ //
//     // set poly degree and plaintext modulus bits size
//     int poly_degree = (int)powl(2, 12);
//     int plain_bits = 42;
//     int cipher_bits = 256;
//     int group_bits = 3072;

//     mpz_class plain_mod;
//     prime_handler::find_ntt_prime(poly_degree, plain_bits, plain_mod);
//     mpz_class psi_p;
//     prime_handler::find_ntt_root(poly_degree, plain_mod, psi_p);

//     mpz_class cipher_mod;
//     prime_handler::find_ntt_prime(poly_degree, cipher_bits, cipher_mod);
//     mpz_class psi_c;
//     prime_handler::find_ntt_root(poly_degree, cipher_mod, psi_c);

//     mpz_class group_mod, group_gen;
//     prime_handler::find_schnorr_prime(cipher_mod, group_bits,  group_mod);
//     prime_handler::find_schnorr_gen(cipher_mod, group_mod, group_gen);

//     cout << "ℹ️ Plaintext modulus: \n" << plain_mod << endl << endl;
//     cout << "ℹ️ Ciphertext modulus: \n" << cipher_mod << endl << endl;
//     cout << "ℹ️ Schnorr Group safety modulus: \n" << group_mod << endl << endl;
//     cout << "ℹ️ Schnorr Group safety generator: \n" << group_gen << endl << endl;

//     // ================================================================================ //
//     //                               Simulation setting                                 //
//     // ================================================================================ //

//     // time trip set
//     auto start = std::chrono::steady_clock::now();
//     auto end = std::chrono::steady_clock::now();
//     std::chrono::duration<double, std::milli> ms_double = end - start;

//     // set cipher secret key
//     poly* sk = new poly(poly_degree);
//     random_handler::secret_key(sk);

//     // set verifiable key 
//     poly* r_0 = new poly(24 * poly_degree);
//     random_handler::not_zero_public_key(cipher_mod, r_0);
//     poly* r_1 = new poly(24 * poly_degree);
//     random_handler::not_zero_public_key(cipher_mod, r_1);
//     poly* s = new poly(3 * poly_degree);
//     random_handler::not_zero_public_key(cipher_mod, s);
    
//     // set controller on ciphertext space
//     arx* ctrl = new arx();
//     vc* ek = new vc();

//     // set plant on plaintext(Real) space
//     plant* plt = new plant(x_init);

//     // set cipher initializer
//     start = std::chrono::steady_clock::now();

//     vector<int64_t> pod_matrix(poly_degree, 0LL);
//     for(int i = 0; i < 4; i++)
//     {
//         pod_matrix[0] = P[i][0];
//         pod_matrix[1] = P[i][1];

//         poly* pre_packed_data = new poly(poly_degree);
//         batch_encoder::encode(pod_matrix, plain_mod, psi_p, pre_packed_data);

//         poly* pre_plaintext = new poly(poly_degree);
//         poly_handler::pack_2_plain(pre_packed_data, pre_plaintext);

//         cipher* ciphertext = new cipher(pre_plaintext->ring_dim, plain_mod, cipher_mod, psi_p, psi_c);
//         crypto_handler::encrypt(pre_plaintext, sk, ciphertext);

//         ctrl->P_y[i] = ciphertext;

//         delete pre_packed_data;
//         delete pre_plaintext;
//     }

//     pod_matrix[1] = 0;
//     for(int i = 0; i < 4; i++)
//     {
//         pod_matrix[0] = Q[i][0];

//         poly* pre_packed_data = new poly(poly_degree);
//         batch_encoder::encode(pod_matrix, plain_mod, psi_p, pre_packed_data);

//         poly* pre_plaintext = new poly(poly_degree);
//         poly_handler::pack_2_plain(pre_packed_data, pre_plaintext);

//         cipher* ciphertext = new cipher(pre_plaintext->ring_dim, plain_mod, cipher_mod, psi_p, psi_c);
//         crypto_handler::encrypt(pre_plaintext, sk, ciphertext);

//         ctrl->Q_u[i] = ciphertext;

//         delete pre_packed_data;
//         delete pre_plaintext;
//     }

//     pod_matrix[0] = 0;
//     poly* pre_packed_data = new poly(poly_degree);
//     batch_encoder::encode(pod_matrix, plain_mod, psi_p, pre_packed_data);

//     poly* pre_plaintext = new poly(poly_degree);
//     poly_handler::pack_2_plain(pre_packed_data, pre_plaintext);

//     cipher* ciphertext = new cipher(pre_plaintext->ring_dim, plain_mod, cipher_mod, psi_p, psi_c);
//     crypto_handler::encrypt(pre_plaintext, sk, ciphertext);

//     ctrl->temp = ciphertext;

//     delete pre_packed_data;
//     delete pre_plaintext;

//     ctrl->zero_set();

//     end = std::chrono::steady_clock::now();
//     ms_double = end - start;

//     cout << "✅ Nominal (encrypted) controller set done | 🕒 " <<  ms_double.count() << "ms" << endl;

//     // set controller side value link to verifiable computation object
//     start = std::chrono::steady_clock::now();

//     // ek->arx_coe_set(ctrl);
//     ek->link_memory(ctrl);
//     ek->set_ekf(ctrl, r_0, r_1, s);

//     end = std::chrono::steady_clock::now();
//     ms_double = end - start;

//     cout << "EKF gen-time | 🕒 " <<  ms_double.count() << "ms" << endl;

//     ek->up_2_g(group_mod, group_gen);

//     end = std::chrono::steady_clock::now();
//     ms_double = end - start;

//     cout << "✅ Verifiable computation set done | 🕒 " <<  ms_double.count() << "ms" << endl << endl;

//     // ================================================================================ //
//     //                                     Simulation                                   //
//     // ================================================================================ //

//     // set maximum iter
//     int iter = 2;

//     // set needs variable
//     poly* plaintext = new poly(poly_degree);
//     poly* packed_data = new poly(poly_degree);
//     vector<int64_t> result_data(poly_degree, 0LL);
//     int64_t temp_u;
//     vector<double> real_u(1);

//     vector<int64_t> plant_output(poly_degree, 0LL);
//     vector<int64_t> control_input(poly_degree, 0LL);
//     cipher* plt_out = new cipher(plaintext->ring_dim, plain_mod, cipher_mod, psi_p, psi_c);
//     cipher* ctrl_in = new cipher(plaintext->ring_dim, plain_mod, cipher_mod, psi_p, psi_c);


//     cout << "ℹ️ Iteration Start" << endl;

//     // iteration start
//     for(int i = 0; i < iter; i++)
//     {
//         // iteration time check
//         start = std::chrono::steady_clock::now();

//         // get plant output and calculation control input
//         plt->output();
//         ctrl->calc();

//         // get control input
//         crypto_handler::decrypt(ctrl->calc_res, sk, plaintext);
//         batch_encoder::decode(plaintext, plain_mod, psi_p, result_data);
//         temp_u = result_data[0] + result_data[1];

//         real_u[0] = (double)temp_u * r_scale * s_scale;

//         // plant output encryption
//         plant_output[0] = (int64_t)(plt->y[0] / r_scale);
//         plant_output[1] = (int64_t)(plt->y[1] / r_scale);
//         batch_encoder::encode(plant_output, plain_mod, psi_p, plaintext);
//         crypto_handler::encrypt(plaintext, sk, plt_out);

//         // re-encryption
//         control_input[0] = (int64_t)(real_u[0] / r_scale);   
//         batch_encoder::encode(control_input, plain_mod, psi_p, plaintext);
//         crypto_handler::encrypt(plaintext, sk, ctrl_in);

//         // plant state update and controller memory update
//         plt->state_update(real_u);
//         ctrl->mem_update(plt_out, ctrl_in);

//         // debug print
//         end = std::chrono::steady_clock::now();
//         ms_double = end - start;
//         cout << "🔄️ iter : "<< i + 1 << " 🕒 : " << ms_double.count() <<"ms " << "🖥️ u: " << real_u[0] << " y: " << plt->y[0] << ", " << plt->y[1] << endl;
//     }

//     // memory free
//     delete sk;
//     delete ctrl;
//     delete ek;
//     delete plt;
    
//     delete plaintext;
//     delete packed_data;

//     delete plt_out;
//     delete ctrl_in;
    
//     return 0;
// }