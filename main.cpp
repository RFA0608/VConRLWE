#include <iostream>
#include <string>
#include <gmpxx.h>

#include "seal/seal.h"
#include <vector>

#include "Struct.h"

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

    mpz_t mod;
    mpz_init_set_str(mod, to_string(parms.plain_modulus().value()).c_str(), 10);

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
    poly_set(packed_data, p2);

    // poly_add(p1, p2, p3);
    poly_mul(p1, p2, p3);
    poly_mod(p3, mod, p3);

    mpz_t rop;
    mpz_init(rop);
    find_primitive_root(mod, rop);

    poly_mul_ntt(p1, p2, mod, rop, p3);

    poly_get(p3, result_data);

    vector<int64_t> repod_matrix(slot_count, 0LL);
    batch_encoder.decode(result_data, repod_matrix);

    for(int i = 0; i < 5; i++)
    {
        cout << repod_matrix[i] << endl;
    }

    p1->~poly();
    p2->~poly();
    p3->~poly();

    return 0;
}