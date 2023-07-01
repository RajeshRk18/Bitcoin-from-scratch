#![allow(unused_imports, unused_must_use, dead_code)]

use crate::fieldelement::*;
use crate::finitefield::*;
use crate::traits::EllipticCurve;
use crate::utils::*;

use lazy_static::lazy_static;
use num_bigint::{BigInt, ToBigInt};
use num_traits::FromPrimitive;
use std::str::FromStr;
use thiserror::Error;

lazy_static! {
    pub static ref ORDER: FieldElement = FieldElement::new(
        BigInt::parse_bytes(
            "115792089237316195423570985008687907852837564279074904382605163141518161494337"
                .as_bytes(),
            10
        )
        .unwrap()
    );
    pub static ref PRIME: FieldElement = FieldElement::new(
        BigInt::parse_bytes(
            "115792089237316195423570985008687907853269984665640564039457584007908834671663"
                .as_bytes(),
            10
        )
        .unwrap()
    );
    pub static ref G_X: FieldElement = FieldElement::new(
        BigInt::parse_bytes(
            "55066263022277343669578718895168534326250603453777594175500187360389116729240"
                .as_bytes(),
            10
        )
        .unwrap()
    );
    pub static ref G_Y: FieldElement = FieldElement::new(
        BigInt::parse_bytes(
            "32670510020758816978083085130507043184471273380659243275938904335757337482424"
                .as_bytes(),
            10
        )
        .unwrap()
    );
}

#[derive(Debug, PartialEq, Eq, PartialOrd)]
pub struct Secp256k1 {
    pub degree0_coeff: FieldElement,
    pub degree3_coeff: FieldElement,
}

impl Secp256k1 {
    pub fn new(degree0_coeff: i32, degree3_coeff: i32) -> Self {
        assert!(Self::is_indiscriminant(0, degree0_coeff));
        Self {
            degree0_coeff: FieldElement::new(BigInt::from_i32(degree0_coeff).unwrap()),
            degree3_coeff: FieldElement::new(BigInt::from_i32(degree3_coeff).unwrap()),
        }
    }

    #[inline(always)]
    pub fn order(&self) -> FieldElement {
        PRIME.clone()
    }

    pub fn generator() -> <Self as EllipticCurve>::Affine {
        (G_X.clone(), G_Y.clone())
    }

    fn zero() -> FieldElement {
        field_element(0)
    }

    fn one() -> FieldElement {
        field_element(1)
    }

    fn is_indiscriminant(a: i32, b: i32) -> bool {
        let op = 4 * a.pow(3) + 27 * b.pow(2);
        if op != 0 {
            return true;
        }
        false
    }

    pub fn mul(
        &self,
        scalar: BigInt,
        point: <Self as EllipticCurve>::Affine,
    ) -> <Self as EllipticCurve>::Affine {
        let to_proj = self.to_projective(point);

        let proj_result = self.scalar_mul(scalar, to_proj);
        let proj_result_to_affine = self.to_affine(proj_result);

        proj_result_to_affine
    }
}

impl EllipticCurve for Secp256k1 {
    type Affine = (FieldElement, FieldElement);
    type Projective = (FieldElement, FieldElement, FieldElement);

    #[inline(always)]
    fn is_on_curve(&self, point: Self::Affine) -> bool {
        let lhs = point.1.pow(2) % PRIME.clone();
        let rhs = (self.degree3_coeff.clone() * point.0.pow(3) + self.degree0_coeff.clone())
            % PRIME.clone();
        if lhs != rhs {
            return false;
        }
        true
    }

    #[inline(always)]
    fn to_projective(&self, point: Self::Affine) -> Self::Projective {
        (point.0, point.1, Self::one())
    }

    #[inline(always)]
    fn to_affine(&self, point: Self::Projective) -> Self::Affine {
        let x = point.0.clone();
        let y = point.1.clone();
        let z = point.2.clone();

        if z == Self::zero() {
            return (Self::zero(), Self::zero());
        } else {
            let z_inv = ext_euclid(z, PRIME.clone());
            let z_inv_squared = z_inv.clone() * z_inv.clone() % PRIME.clone();
            let z_inv_cubed = z_inv_squared.clone() * z_inv.clone() % PRIME.clone();

            let x = x * z_inv_squared % PRIME.clone();
            let y = y * z_inv_cubed % PRIME.clone();

            return (x, y);
        }
    }

    #[inline(always)]
    fn scalar_mul(&self, scalar: BigInt, point: Self::Projective) -> Self::Projective {
        assert!(scalar != new_bigint(0) && FieldElement::from(scalar.clone()) < ORDER.clone());
        let mut result = (Self::zero(), Self::zero(), Self::zero());
        let mut addend = point.clone();
        let mut scalar = scalar;

        while <BigInt as Into<FieldElement>>::into(scalar.clone()) >= Self::zero() {
            if scalar.clone() & new_bigint(1) == new_bigint(1) {
                result = self.point_add(result.clone(), addend.clone());
            }
            scalar >>= 1;
            dbg!(&scalar);
            dbg!(&result);
            addend = self.point_double(addend.clone());
        }
        result
    }

    #[inline(always)]
    fn point_add(&self, point1: Self::Projective, point2: Self::Projective) -> Self::Projective {
        let (x1, y1, z1) = point1.clone();
        let (x2, y2, z2) = point2;

        let z1_squared = z1.pow(2).clone() % PRIME.clone();
        let z2_squared = z2.pow(2).clone() % PRIME.clone();

        let u1 = x1 * &z2_squared % PRIME.clone();
        let u2 = x2 * &z1_squared % PRIME.clone();

        let s1 = y1 * &(&z2_squared * &z2) % PRIME.clone();
        let s2 = y2 * &(&z1_squared * &z1) % PRIME.clone();

        if u1 == u2 {
            if s1 != s2 {
                return (Self::zero(), Self::zero(), Self::zero());
            } else {
                return self.point_double(point1);
            }
        }

        let h = &u2 - &u1;
        let i = h.pow(3) % PRIME.clone();
        let j = (h * &i) % PRIME.clone();
        let r = (&s2 - &s1) % PRIME.clone();
        let v = (u1 * &i) % PRIME.clone();

        let x3 = (r.clone() * &r - j.clone() - field_element(2) * &v) % PRIME.clone();
        let y3 = (r * &(&v - &x3) - s1 * &j) % PRIME.clone();
        let z3 = ((z1 + z2).pow(2) - z1_squared - z2_squared) % PRIME.clone();

        (x3, y3, z3)
    }

    #[inline(always)]
    fn point_double(&self, point: Self::Projective) -> Self::Projective {
        dbg!(point.1.clone());
        if &point.1 == &Self::zero() {
            return (Self::zero(), Self::zero(), Self::zero());
        }
        let (x, y, z) = point;

        let a = (&x * &x) % PRIME.clone();
        let b = (&y * &y) % PRIME.clone();
        let c = (&b * &b) % PRIME.clone();

        let d = (field_element(2) * &(&x + &b) * (&x + &b) - a.clone() - c.clone()) % PRIME.clone();
        let e = (field_element(3) * &a) % PRIME.clone();
        let f = (&e * &e) % PRIME.clone();

        let x3 = (&f - &(field_element(2) * &d)) % PRIME.clone();
        let y3 = (e * &(&d - &x3) - field_element(8) * &c) % PRIME.clone();
        let z3 = (field_element(2) * y * z) % PRIME.clone();

        dbg!(&x3, &y3, &z3);

        (x3, y3, z3)
    }
}

#[derive(Debug, Error)]
#[error("The point {{(0, 1):?}} is not on the curve!")]
pub struct PointError;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fieldelement::*;
    use num_bigint::BigInt;
    use num_traits::FromPrimitive;
    pub(crate) fn gcd(mut a: i128, mut b: i128) -> i128 {
        while b != 0 {
            let t = b;
            b = a % b;
            a = t;
        }
        return a;
    }
    pub(crate) fn ext_euclid(mut a: i128, mut b: i128) -> i128 {
        if gcd(a, b) != 1 {
            return 0;
        } else {
            let mut s = 1;
            let mut s_old = 0;
            let mut t = 0;
            let mut t_old = 1;

            while b != 0 {
                let q = a;
                a = b;
                b = a % b;
                let temp = s;
                s = s_old - q * s;
                s_old = temp;
                let temp1 = t;
                t = t_old - q * t;
                t_old = temp1;
            }

            if s_old < 0 {
                s_old += b;
                return s_old;
            } else {
                return s_old;
            }
        }
    }

    #[test]
    fn scalar_test() {
        let curve = Secp256k1::new(7, 1);
        let point = (field_element(5), field_element(7));
        let scalar = new_bigint(2);
        let result = curve.mul(scalar, point);
        assert_eq!(result, (field_element(2), field_element(5)));
    }

    #[test]
    fn modinv() {
        assert_eq!(ext_euclid(3, 17), 1);
    }

    fn bigint(x: &str) -> BigInt {
        BigInt::parse_bytes(x.as_bytes(), 10).expect("cannot create bigint from the gn input")
    }

    #[test]
    fn test_vectors() {
        let curve = Secp256k1::new(7, 1);
        let point = (G_X.clone(), G_Y.clone());

        assert_eq!(
            curve.mul(bigint("2"), point),
            (
                FieldElement::new(bigint(
                    "89565891926547004231252920425935692360644145829622209833684329913297188986597"
                )),
                FieldElement::new(bigint(
                    "12158399299693830322967808612713398636155367887041628176798871954788371653930"
                ))
            )
        );
    }

    #[test]
    fn to_bits_repr() {
        use rand::{thread_rng, Rng};
        use std::time::{Duration, Instant};
        fn to_bits(k: BigInt) -> Vec<u8> {
            let mut k = k;
            let mut bits = Vec::with_capacity(256);
            let mut count = 0;
            let ti = Instant::now();
            while k > new_bigint(0) {
                bits.push(u8::from_str(&(k.clone() % new_bigint(2)).to_string()).unwrap());
                k /= new_bigint(2);
                count += 1;
            }
            dbg!(ti.elapsed());
            let bits = bits.into_iter().rev().collect::<Vec<u8>>();
            bits
        }
        assert_eq!(to_bits(new_bigint(1)), vec![1]);
        assert_eq!(to_bits(new_bigint(40)), vec![1, 0, 1, 0, 0, 0]);
        assert_eq!(to_bits(new_bigint(4)), vec![1, 0, 0]);
        // https://www.rapidtables.com/convert/number/decimal-to-binary.html?x=72488970228380509287422715226575535698893157 [represented in least significant bit]
        let bit = vec![
            1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1,
            0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1,
            0, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0,
            0, 1, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1,
            1,
        ];
        let bits = bit.into_iter().rev().collect::<Vec<u8>>();
        assert_eq!(
            to_bits(bigint("72488970228380509287422715226575535698893157")),
            bits
        );
        assert_eq!(
            to_bits(bigint(
                "3987468754290932415981468458752646776237488329483678564598965"
            )),
            vec![
                1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1,
                0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1,
                0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0,
                1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0,
                1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1,
                0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0,
                0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0,
                1, 1, 0, 1, 0, 1
            ]
        );
        assert_eq!(to_bits(new_bigint(5)), vec![1, 0, 1]);
        assert_eq!(
            to_bits(new_bigint(1000)),
            vec![1, 1, 1, 1, 1, 0, 1, 0, 0, 0]
        );
    }

    // Just tried out to implement constant time algorithm for converting decimal into binary reps. It is not technically constant time but the runtime bound is between 25 ms and 35 ms.
    // It drastically improves in consistency over non constant time algorithm
    #[test]
    fn constant_time_bits_conversion() {
        use std::thread::sleep;
        use std::time::{Duration, Instant};
        fn to_bits(k: BigInt) -> Vec<u8> {
            let mut bits = vec![0; 256];
            let ti = Instant::now();
            let mut k = k.clone();

            for i in (0..256).rev() {
                let bit = k.clone() & new_bigint(1);
                k >>= 1;

                bits[i] = u8::from_str(&bit.to_string()).unwrap();

                sleep(Duration::from_nanos(150));
            }

            dbg!(ti.elapsed());
            bits
        }
        debug_assert_eq!(
            to_bits(new_bigint(1)),
            vec![
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1
            ]
        );
        assert_eq!(
            to_bits(new_bigint(40)),
            vec![
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                1, 0, 0, 0
            ]
        );
        assert_eq!(
            to_bits(new_bigint(4)),
            vec![
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 1, 0, 0
            ]
        );
        assert_eq!(
            to_bits(bigint("72488970228380509287422715226575535698893157")),
            vec![
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1,
                0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1,
                0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0,
                1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1,
                0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0,
                0, 1, 0, 1
            ]
        );
        assert_eq!(
            to_bits(bigint(
                "3987468754290932415981468458752646776237488329483678564598965"
            )),
            vec![
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
                0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0,
                1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1,
                1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1, 1,
                0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 0,
                1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 0, 1, 0, 1,
                1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0,
                1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1,
                0, 1, 0, 1
            ]
        );
        assert_eq!(
            to_bits(new_bigint(5)),
            vec![
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 1, 0, 1
            ]
        );
        assert_eq!(
            to_bits(new_bigint(1000)),
            vec![
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0,
                1, 0, 0, 0
            ]
        );
    }
}
/*let mut bytes = k.to_bytes_be().1;
let mut binary = Vec::with_capacity(256);

let t = Instant::now();
while !bytes.is_empty() {
    let byte = bytes.pop().unwrap();
    let mut inner_bits = Vec::with_capacity(8);
    for i in (0..8).rev() {
        let bit = (byte >> i) & 1 == 1;
        let bit = u8::from(bit);
        binary.push(bit);
    }
}
dbg!(t.elapsed());
binary*/
