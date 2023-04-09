#![allow(unused_imports, unused_must_use, dead_code)]

use crate::traits::EllipticCurve;
use crate::finitefield::*;

use lazy_static::lazy_static;
use num_bigint::{BigInt, ToBigInt};
use num_traits::FromPrimitive;
use thiserror::Error;

lazy_static! {
    pub static ref ORDER: BigInt = BigInt::parse_bytes(
        "115792089237316195423570985008687907852837564279074904382605163141518161494337".as_bytes(),
        10
    )
    .unwrap();

    pub static ref PRIME: BigInt = BigInt::parse_bytes(
        "115792089237316195423570985008687907853269984665640564039457584007908834671663".as_bytes(),
        10
    )
    .unwrap();

    pub static ref G_X: FieldElement = FieldElement(BigInt::parse_bytes(
        "55066263022277343669578718895168534326250603453777594175500187360389116729240".as_bytes(),
        10
    )
    .unwrap());
    pub static ref G_Y: FieldElement = FieldElement(BigInt::parse_bytes(
        "32670510020758816978083085130507043184471273380659243275938904335757337482424".as_bytes(),
        10
    )
    .unwrap());
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct Secp256k1 {
    pub degree0_coeff: BigInt,
    pub degree3_coeff: BigInt,
}

impl Secp256k1 {
    pub fn new(degree0_coeff: i32, degree3_coeff: i32) -> Self {
        assert!(Self::is_indiscriminant(0, degree0_coeff));
        Self {
            degree0_coeff: BigInt::from_i32(degree0_coeff).unwrap(),
            degree3_coeff: BigInt::from_i32(degree3_coeff).unwrap(),
        }
    }

    #[inline(always)]
    pub fn order(&self) -> BigInt {
        PRIME.to_bigint().unwrap()
    }

    pub fn generator() -> (FieldElement, FieldElement) {
        (G_X.clone(), G_Y.clone())
    }

    fn zero() -> BigInt {
        BigInt::from_u8(0).unwrap()
    }

    fn one() -> BigInt {
        BigInt::from_u8(1).unwrap()
    }

    fn is_indiscriminant(a: i32, b: i32) -> bool {
        let op = 4 * a.pow(3) + 27 * b.pow(2);
        if op != 0 {
            return true;
        }
        false
    }

    fn is_infinity(&self, point: (BigInt, BigInt)) -> bool {
        if point.0 == Self::zero() && point.1 == Self::zero() {
            return true;
        }
        false
    }

    pub fn point_at_infinity() -> (BigInt, BigInt) {
        (BigInt::from_u8(0).unwrap(), BigInt::from_u8(0).unwrap())
    }

    fn ext_euclid(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
        if b == 0i8.to_bigint().unwrap() {
            return (a.clone(), 1i8.to_bigint().unwrap(), 0i8.to_bigint().unwrap());
        } else {
            let (gcd, x1, y1) = Self::ext_euclid(b.clone(), a.clone() % b.clone());
            let x = y1.clone();
            let y = x1.clone() - (a / b) * y1.clone();
            return (gcd, x, y);
        }
    }

    fn mod_inv(&self, a: BigInt, b: BigInt) -> Option<BigInt> {
        let (gcd, x, _) = Self::ext_euclid(a, b.clone());

        if gcd != 1.to_bigint().unwrap() {
            return None;
        } else {
            return Some((x % b.clone() + b.clone()) % b.clone());
        }
    }

}

impl EllipticCurve for Secp256k1 {
    fn check_point_not_at_infinity(
        &self,
        point1: (BigInt, BigInt),
        point2: (BigInt, BigInt),
    ) -> bool {
        if point1.0 == point2.0 {
            return false;
        }
        true
    }

    #[inline(always)]
    fn is_on_curve(&self, point: (BigInt, BigInt)) -> bool {
        let lhs = point.1.pow(2) % PRIME.to_bigint().unwrap();
        dbg!(lhs.clone());
        let rhs = (self.degree3_coeff.clone() * point.0.pow(3) + self.degree0_coeff.clone()) % PRIME.to_bigint().unwrap();
        dbg!(rhs.clone());
        if lhs != rhs {
            return false;
        }
        true
    }

    #[inline(always)]
    fn point_add(
        &self,
        point1: (BigInt, BigInt),
        point2: (BigInt, BigInt),
    ) -> Result<(BigInt, BigInt), PointError> {
        // Returns error if the point is not on the curve
        /*if !self.is_infinity(point1.clone()) {
            if !self.is_on_curve(point1.clone()) {
                return Err(PointError(point1.0.clone(), point1.1.clone()));
            }
        }

        if !self.is_infinity(point2.clone()) {
            if !self.is_on_curve(point2.clone()) {
                return Err(PointError(point2.0.clone(), point2.1.clone()));
            }
        }*/

        // checks if both points are at infinity
        //assert!(point1.0.clone() != Self::zero() && point2.0.clone() != Self::zero() && point1.1.clone() != Self::zero());

        // Returns the other point if the point is at infinity (0, 0)
        if self.is_infinity(point1.clone()) {
            return Ok(point2);
        }
        if self.is_infinity(point2.clone())  {
            return Ok(point1);
        }

        let (x1, y1) = point1;
        let (x2, y2) = point2;

        if (x1 == x2) && (y1 == -y2.clone()) {
            return Ok(Self::point_at_infinity());
        } 
        let slope;
        if (&x1 == &x2) && (&y1 == &y2) {
            slope = (BigInt::from_u8(3).unwrap() * x1.pow(2) / BigInt::from_u8(2).unwrap() * &y1) % self.order();
        } else {
            slope = ((&y2 - &y1) / (&x2 - &x1)) % self.order();
        }

        let x3 = (slope.pow(2) - &x1 - &x2) % self.order();
        let y3 = (slope * (&x1 - &x3) - &y1) % self.order();

        Ok((x3, y3))
        /*elif x1 == x2 && y1 == y2 {

        }
        if x1 != x2 {
            let slope = (&y2 - &y1) / (&x2 - &x1);
            let x3 = slope.pow(2) - &x1 - &x2;
            let y3 = slope * (&x1 - &x3) - &y1;
            return Ok((x3, y3));
        } else if (x1 == x2) && (y1 == -y2.clone()) {
            return Ok((BigInt::from(1 as usize), BigInt::from(0 as usize)));
        } else if ((x1 == x2) && (y1 == y2)) || (x1 == x2) {
            let slope: BigInt = BigInt::from_u32(3).unwrap() * x1.clone().pow(2)
                / BigInt::from_u32(2).unwrap()
                * y1.clone();
            let x3: BigInt = slope.clone().pow(2) - BigInt::from_u32(2).unwrap() * x1.clone();
            let y3: BigInt = slope * (x1.clone() - x3.clone()) - y1.clone();
            return Ok((x3, y3));
        } else {
            return Ok((BigInt::from(1 as usize), BigInt::from(0 as usize)));
        }*/
    }

    #[inline(always)]
    fn point_double(&self, point: (BigInt, BigInt)) -> (BigInt, BigInt) {
        if self.is_infinity(point.clone()) {
            return Self::point_at_infinity();
        } else {
            let slope = ((BigInt::from_u8(3).unwrap() * point.0.clone().pow(2)) / BigInt::from_u8(2).unwrap() * point.1.clone()) % self.order();
            let x3 = (slope.pow(2) - BigInt::from_u8(2).unwrap() * point.0.clone()) % self.order();
            let y3 = (slope * (x3.clone() - point.0) - point.1.clone()) % self.order();
            return (x3, y3);
        }
    }

    #[inline(always)]
    fn scalar_mul(&self, mut scalar: BigInt, point: (BigInt, BigInt)) -> Result<(BigInt, BigInt), PointError> {
        if &scalar == &Self::zero() {
            return Ok(Self::point_at_infinity());
        }
        let mut q = point.clone();

        let mut r = Self::point_at_infinity();

        /*for bit in scalar.to_str_radix(2).chars().rev().skip(1) {
            q = self.point_double(q.clone());
    
            if bit == '1' {
                r = self.point_add(r.clone(), q.clone()).unwrap();
            } 
        }*/
        while scalar > BigInt::from_u8(0).unwrap() {
            if scalar.clone() % 2 == Self::one() {
                r = self.point_add(r.clone(), q.clone())?;
            } else {
                scalar /= BigInt::from_u32(2).unwrap();
                q = self.point_double(q.clone());
            }
        }

        Ok(r)
    }
}

#[derive(Debug, Error, PartialEq)]
#[error("The point {{(0, 1):?}} is not on the curve!")]
pub struct PointError(BigInt, BigInt);

#[cfg(test)]
mod tests {
    use super::*;
    use num_traits::FromPrimitive;
    // Mock curve?
    struct TestCurve {
        b: u32,
        order: i128,
    }

    impl TestCurve {
        fn is_on_curve(&self, point: (BigInt, BigInt)) -> bool {
            let lhs = point.1.pow(2) % BigInt::from_i128(self.order).unwrap();
            let rhs = (point.0.pow(3) + self.b) % BigInt::from_i128(self.order).unwrap();
            if lhs != rhs {
                return false;
            }
            true
        }
    }


    pub(crate) fn gcd(mut a: i128, mut b: i128) -> i128 {
        while b!= 0 {
            let t = b;
            b = a % b;
            a = t;
        }
        return a;
    }
    pub(crate) fn ext_euclid(mut a: i128,mut b: i128) -> i128 {
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
    fn modinv() {
        assert_eq!(ext_euclid(3, 17), 1);
    }

    #[test]
    fn test_vectors() {
        fn new_bigint(x: &str) -> BigInt {
            BigInt::parse_bytes(x.as_bytes(), 10).expect("cannot create bigint from the gn input")
        }

        fn new_u32(x: u32) -> BigInt {
            BigInt::from_u32(x).expect("cannot create bigint from the gn input") 
        }

        let curve = Secp256k1::new(7, 1);

        let point = (new_bigint("112711660439710606056748659173929673102114977341539408544630613555209775888121"), new_bigint("25583027980570883691656905877401976406448868254816295069919888960541586679410"));

        assert_eq!(curve.scalar_mul(new_u32(3), point), Ok((new_u32(0), new_u32(0))));
    }
}
