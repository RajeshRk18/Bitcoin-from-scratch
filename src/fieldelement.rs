#![allow(unused_variables, dead_code)]

use crate::ecc::PRIME;
use crate::utils::*;
use num_bigint::BigInt;
use std::cmp::Ordering;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Rem, Sub};

#[derive(Debug, Clone)]
pub struct FieldElement(BigInt);

impl FieldElement {
    pub fn new(value: BigInt) -> Self {
        FieldElement(value)
    }
    pub fn pow(&self, exp: u32) -> FieldElement {
        self.0.pow(exp).into()
    }

    pub fn gcd(&mut self) -> FieldElement {
        let mut p = field_element(15);
        while p != field_element(0) {
            dbg!(&p);
            let temp = p.clone();
            p = modulus(self.clone(), p.clone());
            *self = temp;
        }
        return self.clone();
    }
}

impl Add<FieldElement> for FieldElement {
    type Output = FieldElement;
    fn add(self, rhs: Self) -> Self::Output {
        ((self.0 + rhs.0) % PRIME.0.clone()).into()
    }
}

impl<'a> Add<&'a FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn add(self, rhs: &'a FieldElement) -> Self::Output {
        ((self.0.clone() + rhs.0.clone()) % PRIME.0.clone()).into()
    }
}
impl<'a> Add<&'a FieldElement> for FieldElement {
    type Output = FieldElement;
    fn add(self, rhs: &'a FieldElement) -> Self::Output {
        ((self.0 + rhs.0.clone()) % PRIME.0.clone()).into()
    }
}

impl Sub<FieldElement> for FieldElement {
    type Output = FieldElement;
    fn sub(self, rhs: Self) -> Self::Output {
        let result = (self.0.clone() + (PRIME.0.clone() - rhs.0.clone())) % PRIME.0.clone();
        result.into()
    }
}

impl<'a> Sub<&'a FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn sub(self, rhs: &'a FieldElement) -> Self::Output {
        let result = (self.0.clone() + (PRIME.0.clone() - rhs.0.clone())) % PRIME.0.clone();
        result.into()
    }
}

impl<'a> Sub<&'a FieldElement> for FieldElement {
    type Output = FieldElement;
    fn sub(self, rhs: &'a FieldElement) -> Self::Output {
        let result = (self.0.clone() + (PRIME.0.clone() - rhs.0.clone())) % PRIME.0.clone();
        result.into()
    }
}

impl Mul<FieldElement> for FieldElement {
    type Output = FieldElement;
    fn mul(self, rhs: Self) -> Self::Output {
        ((self.0.clone() * rhs.0.clone()) % PRIME.0.clone()).into()
    }
}

impl<'a> Mul<&'a FieldElement> for FieldElement {
    type Output = FieldElement;
    fn mul(self, rhs: &'a FieldElement) -> Self::Output {
        ((self.0.clone() * rhs.0.clone()) % PRIME.0.clone()).into()
    }
}

impl<'a> Mul<&'a FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn mul(self, rhs: &'a FieldElement) -> Self::Output {
        ((self.0.clone() * rhs.0.clone()) % PRIME.0.clone()).into()
    }
}

impl Div<FieldElement> for FieldElement {
    type Output = FieldElement;
    fn div(self, rhs: Self) -> Self::Output {
        self.0.div(rhs.0).into()
    }
}

impl<'a> Div<&'a FieldElement> for FieldElement {
    type Output = FieldElement;
    fn div(self, rhs: &'a FieldElement) -> Self::Output {
        self.0.clone().div(rhs.0.clone()).into()
    }
}
impl<'a> Div<&'a FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn div(self, rhs: &'a FieldElement) -> Self::Output {
        self.0.clone().div(rhs.0.clone()).into()
    }
}

impl Neg for FieldElement {
    type Output = FieldElement;
    fn neg(self) -> Self::Output {
        //(modulo - value) % modulo
        let reduced = (PRIME.0.clone() - self.0.clone()) % PRIME.0.clone();
        FieldElement(reduced)
    }
}

impl<'a> Neg for &'a FieldElement {
    type Output = FieldElement;
    fn neg(self) -> Self::Output {
        let reduced = (PRIME.0.clone() - self.0.clone()) % PRIME.0.clone();
        FieldElement(reduced)
    }
}

impl MulAssign for FieldElement {
    fn mul_assign(&mut self, rhs: Self) {
        self.0 *= rhs.0;
    }
}

impl AddAssign for FieldElement {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone().add(rhs);
    }
}

impl Rem for FieldElement {
    type Output = FieldElement;
    fn rem(self, rhs: Self) -> Self::Output {
        (self.0 % rhs.0).into()
    }
}

impl<'a> Rem<&'a FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn rem(self, rhs: &'a FieldElement) -> Self::Output {
        (self.0.clone() % rhs.0.clone()).into()
    }
}

impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }

    fn ne(&self, other: &Self) -> bool {
        &self.0 != &other.0
    }
}

impl PartialEq<i32> for FieldElement {
    fn eq(&self, other: &i32) -> bool {
        self == &FieldElement(BigInt::from(other.clone()))
    }

    fn ne(&self, other: &i32) -> bool {
        self != &FieldElement(BigInt::from(other.clone()))
    }
}

impl PartialOrd for FieldElement {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.0.cmp(&other.0))
    }
    fn ge(&self, other: &Self) -> bool {
        &self.0 >= &other.0
    }

    fn le(&self, other: &Self) -> bool {
        &self.0 <= &other.0
    }

    fn gt(&self, other: &Self) -> bool {
        &self.0 > &other.0
    }

    fn lt(&self, other: &Self) -> bool {
        &self.0 < &other.0
    }
}

impl Eq for FieldElement {}

impl From<BigInt> for FieldElement {
    fn from(value: BigInt) -> Self {
        FieldElement::new(value)
    }
}

impl<'a> From<&'a BigInt> for FieldElement {
    fn from(value: &'a BigInt) -> Self {
        FieldElement::new(value.clone())
    }
}

#[cfg(test)]
mod field_element {
    use super::*;

    #[test]
    fn test_neg() {
        let neg = -field_element(34);
        assert_eq!(
            neg,
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671629"
            )
        );
        assert_eq!(
            -field_element(456),
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671207"
            )
        );
    }

    #[test]
    fn neg_op() {
        let lhs = field_element(900);
        let rhs = -field_element(34);
        assert_eq!((lhs + rhs) % PRIME.clone(), field_element(866));
    }

    #[test]
    fn sub() {
        let lhs = field_element(45);
        let rhs = field_element(34);
        assert_eq!(lhs - rhs, field_element(11));
    }

    #[test]
    fn big_neg() {
        let lhs = -felt_from_str("115792089237316195423570985008687907853269984665640564039457584007908834671629869898598");
        assert_eq!(lhs, felt_from_str("33130101402"));
    }

    #[test]
    fn test_exp() {
        let lhs = field_element(45) - field_element(30) * field_element(3);

        debug_assert_eq!(
            lhs,
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671618"
            )
        )
    }

    #[test]
    fn sub_exp() {
        let lhs = field_element(-45);
        let rhs = felt_from_str(
            "115792089237316195423570985008687907853269984665640564039457584007908834671618",
        );

        debug_assert_eq!(lhs, rhs);
    }
}
