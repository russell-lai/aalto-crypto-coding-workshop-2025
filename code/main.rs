struct PKEParams {
    secpar: u32,    // security parameter
    f:u32,          // conductor
    phi: u32,       // degree
    n: u32,         // module rank
    q: u32,         // ciphertext modulus
    ell: u32,       // plaintext dimension
    p: u32,         // plaintext modulus
    b: u32,         // parameter for centred binomial distribution
}

fn main() {
    println!("Hello, World!");
}