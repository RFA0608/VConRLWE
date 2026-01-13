#include <iostream>
#include <vector>
#include <string>
#include <gmpxx.h>

#include "seal/seal.h"
#include "RLWE.h"

using namespace std;
using namespace seal;

int main()
{
    // ============================================================================= //
    // Setting for Encryption
    // ============================================================================= //
    // set poly degree and plaintext modulus bits size
    int poly_degree = (int)powl(2, 14);
    int plain_modulus = 42;

    // make context from parameters
    EncryptionParameters parms(scheme_type::bgv);
    size_t poly_modulus_degree = (size_t)poly_degree;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, plain_modulus));
    SEALContext context(parms);

    // set batch encoder 
    BatchEncoder batch_encoder(context);
    size_t slot_count = batch_encoder.slot_count();

    // ============================================================================= //
    // Test
    // ============================================================================= //
    vector<int64_t> pod_matrix(slot_count, 0LL);
    pod_matrix[0] = 1;
    pod_matrix[1] = -2;
    pod_matrix[2] = 4;
    pod_matrix[3] = -8;
    pod_matrix[4] = 11;

    Plaintext packed_data, result_data(poly_degree);
    batch_encoder.encode(pod_matrix, packed_data);
    
    poly* p1 = new poly(poly_degree);
    poly* p2 = new poly(poly_degree);
    poly* p3 = new poly(poly_degree);

    mpz_class plain_mod(to_string(parms.plain_modulus().value()));
    mpz_class g_p;
    prime_handler::find_primitive_root(plain_mod, g_p);

    mpz_class cipher_mod;
    prime_handler::find_ntt_prime(poly_degree, 210, cipher_mod);
    mpz_class g_c;
    prime_handler::find_primitive_root(cipher_mod, g_c);

    poly_handler::plain_2_poly(packed_data, p1);
    poly_handler::plain_2_poly(packed_data, p2);

    cipher* c1 = new cipher(p1->ring_dim, plain_mod, cipher_mod, g_p, g_c);
    cipher* c2 = new cipher(p2->ring_dim, plain_mod, cipher_mod, g_p, g_c);
    cipher* c3 = new cipher(p2->ring_dim, plain_mod, cipher_mod, g_p, g_c);
    cipher* c4 = new cipher(p2->ring_dim, plain_mod, cipher_mod, g_p, g_c);

    poly* sk = new poly(poly_degree);
    random_handler::secret_key(sk);

    crypto_handler::encrypt(p1, sk, c1);
    crypto_handler::encrypt(p2, sk, c2);
    crypto_handler::eval_mul(c1, c2, c3);
    crypto_handler::eval_add(c3, c3, c4);

    crypto_handler::decrypt(c4, sk, p3);

    poly_handler::poly_2_plain(p3, result_data);
    
    vector<int64_t> repod_matrix(slot_count, 0LL);
    batch_encoder.decode(result_data, repod_matrix);

    for(int i = 0; i < 5; i++)
    {
        cout << repod_matrix[i] << endl;
    }

    delete p1;
    delete p2;
    delete p3;
    return 0;
}