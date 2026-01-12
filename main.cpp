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

    mpz_t plain_mod;
    mpz_init_set_str(plain_mod, to_string(parms.plain_modulus().value()).c_str(), 10);

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

    poly_set(packed_data, p1);

    mpz_t* cipher_mod = find_ntt_prime(p1->size, 110);
    mpz_out_str(stdout, 10, *cipher_mod);
    cout << endl;
    cout <<"is prime : " << is_prime(*cipher_mod) << endl;

    cipher* c = new cipher(p1->size, plain_mod, *cipher_mod);
    poly* sk = secret_key(p1->size);

    encrypt(p1, sk, c);
    decrypt(c, sk, p2);
    poly_get(p2, result_data);
    // poly_set(packed_data, p2);

    // // poly_add(p1, p2, p3);
    // poly_mul(p1, p2, p3);
    // poly_mod(p3, mod, p3);

    // mpz_t rop;
    // mpz_init(rop);
    // find_primitive_root(mod, rop);
    // int r = is_prime(mod);


    // poly_mul_ntt(p1, p2, mod, rop, p3);
    // poly_get(p3, result_data);

    vector<int64_t> repod_matrix(slot_count, 0LL);
    batch_encoder.decode(result_data, repod_matrix);

    for(int i = 0; i < 5; i++)
    {
        cout << repod_matrix[i] << endl;
    }

    delete p1;
    delete p2;
    delete p3;
    delete c;
    delete sk;

    mpz_clears(plain_mod, *cipher_mod, NULL);
    delete cipher_mod;

    return 0;
}