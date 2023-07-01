use crate::fieldelement::*;
use num_bigint::{BigInt, ToBigInt};
use std::str::FromStr;
use std::thread::sleep;
use std::time::Duration;

pub fn field_element<T: ToBigInt>(value: T) -> FieldElement {
    FieldElement::new(value.to_bigint().unwrap())
}

pub fn felt_from_str(value: &str) -> FieldElement {
    FieldElement::new(BigInt::parse_bytes(value.as_bytes(), 10).unwrap())
}

// As we work under prime field, there always exists multiplicative inverse.
pub fn ext_euclid(mut a: FieldElement, mut b: FieldElement) -> FieldElement {
    if gcd(a.clone(), b.clone()) != 1 {
        return field_element(0);
    } else {
        let mut s = field_element(1);
        let mut s_old = field_element(0);
        let mut t = field_element(0);
        let mut t_old = field_element(1);

        while b != 0 {
            let q = a;
            a = b.clone();
            b = &a % &b;
            let temp = s.clone();
            s = s_old - &q * &s;
            s_old = temp;
            let temp1 = t.clone();
            t = t_old - q * t;
            t_old = temp1;
        }

        if s_old < field_element(0) {
            s_old += b;
            return s_old;
        } else {
            return s_old;
        }
    }
}

pub fn modulus(a: FieldElement, p: FieldElement) -> FieldElement {
    let mut r = a.clone();
    if r < p {
        r
    } else {
        while r > p {
            r = &r - &p;
        }
        r
    }
}

pub fn gcd(mut a: FieldElement, mut b: FieldElement) -> FieldElement {
    while b != 0 {
        let t = b.clone();
        b = &a % &b;
        a = t;
    }
    return a;
}

pub fn new_bigint<T: ToBigInt>(num: T) -> BigInt
where
    BigInt: From<T>,
{
    num.into()
}

pub fn to_bits(num: BigInt) -> Vec<u8> {
    let mut bits = vec![0; 256];
    let mut k = num.clone();

    for i in (0..256).rev() {
        let bit = k.clone() & new_bigint(1);
        k /= new_bigint(2);

        bits[i] = u8::from_str(&bit.to_string()).unwrap();

        sleep(Duration::from_nanos(150));
    }
    bits
}
