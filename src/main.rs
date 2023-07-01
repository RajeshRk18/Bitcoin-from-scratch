mod ecc;
mod fieldelement;
mod finitefield;
mod sec;
mod traits;
mod utils;




fn main() {
    /*let bitcoin_curve = Secp256k1::new(7, 1);

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

    let res = bitcoin_curve.mul(BigInt::from_u8(1).unwrap(), (BigInt::parse_bytes("55066263022277343669578718895168534326250603453777594175500187360389116729240".as_bytes(), 10).unwrap(), BigInt::parse_bytes("32670510020758816978083085130507043184471273380659243275938904335757337482424".as_bytes(), 10).unwrap()));
    println!("{:?}", res);*/
}
