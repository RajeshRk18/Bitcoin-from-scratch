#![allow(unused_variables, dead_code)]

use crate::ecc::PRIME;
use crate::utils::*;
use num_bigint::{BigUint, BigInt};
use std::cmp::Ordering;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Rem, Sub};

#[derive(Debug, Clone)]
pub struct FieldElement(BigUint);

impl FieldElement {
    pub fn new(value: BigUint) -> Self {
        FieldElement(value)
    }
    pub fn pow(&self, exp: u32) -> FieldElement {
        self.0.pow(exp).into()
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
        let result: BigUint =
            (self.0.clone() + (PRIME.0.clone() - rhs.0.clone())) % PRIME.0.clone();
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

impl PartialEq<u32> for FieldElement {
    fn eq(&self, other: &u32) -> bool {
        self == &FieldElement(BigUint::from(other.clone()))
    }

    fn ne(&self, other: &u32) -> bool {
        self != &FieldElement(BigUint::from(other.clone()))
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

impl From<BigUint> for FieldElement {
    fn from(value: BigUint) -> Self {
        FieldElement::new(value)
    }
}

impl<'a> From<&'a BigUint> for FieldElement {
    fn from(value: &'a BigUint) -> Self {
        FieldElement::new(value.clone())
    }
}

impl From<BigInt> for FieldElement {
    fn from(value: BigInt) -> Self {
        FieldElement::new(value.to_biguint().unwrap())
    }
}

impl<'a> From<&'a BigInt> for FieldElement {
    fn from(value: &'a BigInt) -> Self {
        FieldElement::new(value.to_biguint().unwrap().clone())
    }
}

#[cfg(test)]
mod ff_test {
    use super::*;

    #[test]
    fn test_neg1() {
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
    fn test_neg2() {
        debug_assert_eq!(
            -field_element(45),
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671618",
            )
        );

        debug_assert_eq!(
            -felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671618"
            ),
            field_element(45)
        )
    }

    #[test]
    fn test_neg_op() {
        let lhs = field_element(900);
        let rhs = -field_element(34);
        assert_eq!((lhs + rhs) % PRIME.clone(), field_element(866));
    }

    #[test]
    fn test_sub() {
        let lhs = field_element(45);
        let rhs = field_element(34);
        assert_eq!(lhs - rhs, field_element(11));
    }

    #[test]
    fn test_big_neg1() {
        let lhs = -felt_from_str("115792089237316195423570985008687907853269984665640564039457584007875704570261");
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
    fn test_two_sub() {
        debug_assert_eq!(
            -field_element(45) - field_element(3),
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671615"
            )
        )
    }

    #[test]
    fn test_mulassign() {
        let mut a = field_element(5);
        let b = field_element(34);
        a *= b;
        assert_eq!(a, field_element(170));
    }

    #[test]
    fn test_partial_eq() {
        let a = field_element(6763567);
        let b = felt_from_str("8753976666665397589478555");
        let c = field_element(34);
        let d = felt_from_str("6763567");

        assert!(a == d);
        assert!(c < d);
        assert!(b > a);
        assert!(c <= a);
    }
}
