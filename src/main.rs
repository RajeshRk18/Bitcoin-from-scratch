mod ecc;
mod finitefield;
mod traits;

use ecc::*;
use finitefield::*;

use num_bigint::{BigUint, BigInt};
use num_traits::FromPrimitive;
use rand::Rng;

fn y(point: BigInt) -> (BigInt, BigInt) {
    let m_prime: BigInt = BigInt::parse_bytes(
        "115792089237316195423570985008687907852837564279074904382605163141518161494337".as_bytes(),
        10,
    )
    .unwrap();

    let x = point.clone();
    let y = (point.pow(3) + BigInt::from_u8(7).unwrap()).sqrt() % m_prime;

    (x, y)
}

fn main() {
    let _bitcoin_curve = Secp256k1::new(7, 1);

    let _ff = FiniteField::new(BigInt::from_u32(11).unwrap());

    //let generators = ff.find_generators();

    /*let mut range = rand::thread_rng();
    let points: Vec<(BigUint, BigUint)> = (0..100000)
        .into_iter()
        .map(|_| {
            let x = BigUint::from_u128(range.gen::<u128>()).unwrap();
            if let Some((x, y)) = in_curve(x.clone()) {
                (x, y)
            } else {
                (
                    BigUint::from_u128(0).unwrap(),
                    BigUint::from_u128(0).unwrap(),
                )
            }
        })
        .collect();

    let filtered: Vec<&(BigUint, BigUint)> = points
        .iter()
        .filter(|(x, y)| {
            *x != BigUint::from_u128(0).unwrap() && *y != BigUint::from_u128(0).unwrap()
        })
        .collect();
    println!("{:?}", filtered);*/

    //println!("{:?}", ff.order());
    //println!("{:?}", generators);
    //println!("{:?}", ff);
    let (x, y) = y(BigInt::parse_bytes("55066263022277343669578718895168534326250603453777594175500187360389116729240".as_bytes(), 10).unwrap());

    println!("x: {:?}\ny: {:?}", x, y);

}