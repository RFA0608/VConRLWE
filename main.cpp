#include <iostream>
#include <vector>
#include <string>
#include <gmpxx.h>
#include <chrono>

#include "RLWE.h"

using namespace std;

int main()
{
    // ============================================================================= //
    // Setting for Encryption
    // ============================================================================= //
    // set poly degree and plaintext modulus bits size
    int poly_degree = (int)powl(2, 13);
    int plain_bits = 40;
    int cipher_bits = 210;

    mpz_class plain_mod;
    prime_handler::find_ntt_prime(poly_degree, plain_bits, plain_mod);
    mpz_class g_p;
    prime_handler::find_primitive_root(plain_mod, g_p);

    mpz_class cipher_mod;
    prime_handler::find_ntt_prime(poly_degree, cipher_bits, cipher_mod);
    mpz_class g_c;
    prime_handler::find_primitive_root(cipher_mod, g_c);

    // ============================================================================= //
    // Test
    // ============================================================================= //
    vector<int64_t> pod_matrix(poly_degree, 0LL);
    pod_matrix[0] = 1;
    pod_matrix[1] = -2;
    pod_matrix[2] = 4;
    pod_matrix[3] = -8;
    pod_matrix[4] = 11;

    poly* packed_data = new poly(poly_degree), * result_data = new poly(poly_degree);
    batch_encoder::encode(pod_matrix, plain_mod, g_p, packed_data);
    
    poly* p1 = new poly(poly_degree);
    poly* p2 = new poly(poly_degree);
    poly* p3 = new poly(poly_degree);

    poly_handler::plain_2_poly(packed_data, p1);
    poly_handler::plain_2_poly(packed_data, p2);

    cipher* c1 = new cipher(p1->ring_dim, plain_mod, cipher_mod, g_p, g_c);
    cipher* c2 = new cipher(p2->ring_dim, plain_mod, cipher_mod, g_p, g_c);
    cipher* c3 = new cipher(p2->ring_dim, plain_mod, cipher_mod, g_p, g_c);
    cipher* c4 = new cipher(p2->ring_dim, plain_mod, cipher_mod, g_p, g_c);

    poly* sk = new poly(poly_degree);
    random_handler::secret_key(sk);

    auto start = std::chrono::high_resolution_clock::now();

    crypto_handler::encrypt(p1, sk, c1);
    crypto_handler::encrypt(p2, sk, c2);
    crypto_handler::eval_mul(c1, c2, c3);
    crypto_handler::eval_add(c3, c3, c4);

    crypto_handler::decrypt(c4, sk, p3);
    poly_handler::poly_2_plain(p3, result_data);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> ms_double = end - start;
    std::cout << "runtime: " << ms_double.count() << " ms" << std::endl;
    
    vector<int64_t> repod_matrix(poly_degree, 0LL);
    batch_encoder::decode(result_data, plain_mod, g_p, repod_matrix);

    for(int i = 0; i < 5; i++)
    {
        cout << repod_matrix[i] << endl;
    }

    delete packed_data;
    delete result_data;
    delete p1;
    delete p2;
    delete p3;
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete sk;
    return 0;
}