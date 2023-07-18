use crate::ecc::Secp256k1;
use crate::traits::Secp256k1Curve;
use crate::utils::ext_euclid;
use num_bigint::{BigInt, BigUint};
use std::cmp::Ordering;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Rem, Sub};

#[derive(Debug, Clone)]
pub struct FieldElement(BigUint);

impl FieldElement {
    pub fn new(value: BigUint) -> Self {
        FieldElement(value)
    }

    pub fn floor_div(a: &Self, b: &Self) -> Self {
        (&a.0 / &b.0).into()
    }
}

//      ---ADD---
impl<'a> Add<&'a FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn add(self, rhs: &'a FieldElement) -> Self::Output {

        let add: FieldElement = (&self.0 + &rhs.0).into();
        add % Secp256k1::scalar_field_order()
    }
}
impl Add<FieldElement> for FieldElement {
    type Output = FieldElement;
    fn add(self, rhs: Self) -> Self::Output {
        &self + &rhs
    }
}
impl<'a> Add<&'a FieldElement> for FieldElement {
    type Output = FieldElement;
    fn add(self, rhs: &'a FieldElement) -> Self::Output {
        &self + rhs
    }
}

//      ---SUB---
impl<'a> Sub<&'a FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn sub(self, rhs: &'a FieldElement) -> Self::Output {
        let prime = Secp256k1::scalar_field_order();
        ((&self.0 + (prime.0.clone() - &rhs.0)) % prime.0.clone()).into()
    }
}
impl Sub<FieldElement> for FieldElement {
    type Output = FieldElement;
    fn sub(self, rhs: Self) -> Self::Output {
        &self - &rhs
    }
}
impl<'a> Sub<&'a FieldElement> for FieldElement {
    type Output = FieldElement;
    fn sub(self, rhs: &'a FieldElement) -> Self::Output {
        &self - rhs
    }
}

//      ---MUL---
impl<'a> Mul<&'a FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn mul(self, rhs: &'a FieldElement) -> Self::Output {
        let mul: FieldElement = (&self.0 * &rhs.0).into(); 
        mul % Secp256k1::scalar_field_order()
    }
}
impl Mul<FieldElement> for FieldElement {
    type Output = FieldElement;
    fn mul(self, rhs: Self) -> Self::Output {
        &self * &rhs
    }
}
impl<'a> Mul<&'a FieldElement> for FieldElement {
    type Output = FieldElement;
    fn mul(self, rhs: &'a FieldElement) -> Self::Output {
        &self * rhs
    }
}

//      ---MUL INVERSE---
impl<'a> Div<&'a FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn div(self, rhs: &'a FieldElement) -> Self::Output {
        self * &ext_euclid(rhs, &Secp256k1::scalar_field_order())
    }
}
impl Div<FieldElement> for FieldElement {
    type Output = FieldElement;
    fn div(self, rhs: Self) -> Self::Output {
        &self / &rhs
    }
}
impl<'a> Div<&'a FieldElement> for FieldElement {
    type Output = FieldElement;
    fn div(self, rhs: &'a FieldElement) -> Self::Output {
        &self / rhs
    }
}

//      ---NEGATE---
impl<'a> Neg for &'a FieldElement {
    type Output = FieldElement;
    fn neg(self) -> Self::Output {
        //(modulo - value) % modulo
        Secp256k1::scalar_field_order() - self
    }
}
impl Neg for FieldElement {
    type Output = FieldElement;
    fn neg(self) -> Self::Output {
        -&self
    }
}

//      ---MULASSIGN---
impl MulAssign for FieldElement {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone() * rhs;
    }
}

//      ---ADDASSIGN---
impl AddAssign for FieldElement {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}

//      ---MODULO--- 
impl Rem for FieldElement {
    type Output = FieldElement;
    fn rem(self, rhs: Self) -> Self::Output {
        &self % &rhs
    }
}
impl<'a> Rem<&'a FieldElement> for &'a FieldElement {
    type Output = FieldElement;
    fn rem(self, rhs: &'a FieldElement) -> Self::Output {
        (&self.0 % &rhs.0).into()
    }
}

//      ---EQUAL & NOT-EQUAL---
impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }

    fn ne(&self, other: &Self) -> bool {
        self.0 != other.0
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

//      ---GREATER-THAN, LESS-THAN, GREATER/EQUAL, LESSER/EQUAL---
impl PartialOrd for FieldElement {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.0.cmp(&other.0))
    }
    fn ge(&self, other: &Self) -> bool {
        self.0 >= other.0
    }

    fn le(&self, other: &Self) -> bool {
        self.0 <= other.0
    }

    fn gt(&self, other: &Self) -> bool {
        self.0 > other.0
    }

    fn lt(&self, other: &Self) -> bool {
        self.0 < other.0
    }
}

impl Eq for FieldElement {}

//      ---TYPE CONVERSIONS---
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
    use crate::utils::*;

    #[test]
    fn test_neg1() {
        let neg = -felt_from_uint(34);
        assert_eq!(
            neg,
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671629"
            )
        );
        assert_eq!(
            -felt_from_uint(456),
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671207"
            )
        );
    }

    #[test]
    fn test_neg2() {
        debug_assert_eq!(
            -felt_from_uint(45),
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671618",
            )
        );

        debug_assert_eq!(
            -felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671618"
            ),
            felt_from_uint(45)
        )
    }

    #[test]
    fn test_neg_op() {
        let pos = felt_from_uint(900);
        let neg = -felt_from_uint(34);
        assert_eq!(pos + neg, felt_from_uint(866));
    }

    #[test]
    fn test_sub() {
        let lhs = felt_from_uint(45);
        let rhs = felt_from_uint(34);
        assert_eq!(lhs - rhs, felt_from_uint(11));
    }

    #[test]
    fn test_big_neg1() {
        let lhs = -felt_from_str(
            "115792089237316195423570985008687907853269984665640564039457584007875704570261",
        );
        assert_eq!(lhs, felt_from_str("33130101402"));
    }

    #[test]
    fn test_exp() {
        let lhs = felt_from_uint(45) - felt_from_uint(30) * felt_from_uint(3);

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
            -felt_from_uint(45) - felt_from_uint(3),
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671615"
            )
        )
    }

    #[test]
    fn test_mulassign() {
        let mut a = felt_from_uint(5);
        let b = felt_from_uint(34);
        a *= b;
        assert_eq!(a, felt_from_uint(170));
    }

    #[test]
    fn test_partial_eq() {
        let a = felt_from_uint(6763567);
        let b = felt_from_str("8753976666665397589478555");
        let c = felt_from_uint(34);
        let d = felt_from_str("6763567");

        assert!(a == d);
        assert!(c < d);
        assert!(b > a);
        assert!(c <= a);
    }

    #[test]
    fn test_exp2() {
        let exp = felt_from_uint(25) * felt_from_uint(5) - felt_from_uint(32) + felt_from_uint(21);
        debug_assert_eq!(exp, felt_from_uint(114));
    }

    #[test]
    fn test_exp3() {
        let exp = -felt_from_uint(45) * felt_from_uint(2) + felt_from_uint(67);
        debug_assert_eq!(
            exp,
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671640"
            )
        );
    }

    #[test]
    fn test_result_one() {
        let a = felt_from_str("7137376184023522026654217343440040540351594155614695530712440441403759244341");
        let b = felt_from_str("18233444644265725414720095600783733811551140822142748479967321742263849032310");

        debug_assert_eq!(a * b, felt_from_uint(3));
    }

    #[test]
    fn test_div() {
        let a = felt_from_str("18233444644265725414720095600783733811551140822142748479967321742263849032310");
        let b = felt_from_str("12158399299693830322967808612713398636155367887041628176798871954788371653930");

        debug_assert_eq!(&a / &b , felt_from_str("3581453595367386406056576018017432656403623753535430210497256105630064126918"));
    }
}
