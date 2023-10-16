use crate::ecc::Secp256k1;
use crate::traits::{Secp256k1Curve, ZeroIt};
use crate::utils::*;
use num_bigint::{BigInt, BigUint, ToBigInt, RandomBits};
use num_traits::FromPrimitive;
use rand::distributions::Distribution;
use rand_chacha::rand_core::OsRng;
use std::cmp::Ordering;
use std::ops::{Add, AddAssign, Div, Mul, MulAssign, Neg, Rem, Sub, SubAssign};

#[derive(Debug, Clone)]
pub struct FieldElement(BigUint);

impl FieldElement {
    pub fn new(value: BigUint) -> Self {
        FieldElement(value)
    }

    pub fn floor_div(a: &Self, b: &Self) -> Self {
        (&a.0 / &b.0).into()
    }

    pub fn from_bytes(bytes: &[u8]) -> FieldElement {
        FieldElement(BigUint::from_bytes_le(bytes))
    }

    pub fn gen_random() -> Self {
        let mut rng = OsRng {};
        let big_rand = RandomBits::new(255);

        let value = big_rand.sample(&mut rng);
        FieldElement(value)
    }

    pub fn one() -> Self {
        felt_from_uint(1)
    }

    pub fn zero() -> Self {
        felt_from_uint(0)
    }
}

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

impl MulAssign for FieldElement {
    fn mul_assign(&mut self, rhs: Self) {
        *self = self.clone() * rhs;
    }
}

impl AddAssign for FieldElement {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}

impl SubAssign for FieldElement {
    fn sub_assign(&mut self, rhs: Self) {
        *self = self.clone() - rhs
    }
}

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


// constant time comparisions to be implemented from https://gist.github.com/sneves/10845247
impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }

    fn ne(&self, other: &Self) -> bool {
        !self.eq(other)
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

impl Default for FieldElement {
    fn default() -> Self {
        Self::zero()
    }
}

//      ---ZEROIZE---
impl ZeroIt for FieldElement {
    // This will be used to clear any secret data (eg., private key) from memory after being used.

    // NOTE: This still wont securely zero out as the BigInt in num_bigint crate implements Copy trait
    // Value will be copied when moved and thus value at original memory will not be cleared.

    // TODO: To use crypto-bigint crate which provides Boxed limbs and thus pointer can be cloned when needed
    // to be moved.
    fn zeroize(&mut self) {
        use std::ptr::write_volatile;
        use std::sync::atomic;

        unsafe {
            write_volatile(&mut *self, Self::default());
        }
        atomic::compiler_fence(atomic::Ordering::SeqCst);
    }
}

#[cfg(test)]
mod ff_test {
    use crate::{ecc::Secp256k1, traits::Secp256k1Curve, utils::*, fieldelement::FieldElement};

    #[test]
    fn test_negation() {
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
        assert_eq!(
            -felt_from_uint(45),
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671618",
            )
        );

        assert_eq!(
            -felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671618"
            ),
            felt_from_uint(45)
        );
        let lhs = -felt_from_str(
            "115792089237316195423570985008687907853269984665640564039457584007875704570261",
        );
        assert_eq!(lhs, felt_from_str("33130101402"));
    }

    #[test]
    fn test_sub() {
        let lhs = felt_from_uint(45);
        let rhs = felt_from_uint(34);
        assert_eq!(lhs - rhs, felt_from_uint(11));
        debug_assert_eq!(
            -felt_from_uint(45) - felt_from_uint(3),
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671615"
            )
        );
        debug_assert_eq!(
            felt_from_uint(275386243) - felt_from_uint(275386243),
            felt_from_uint(0)
        );
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
    fn test_mul() {
        let mut a = felt_from_uint(5);
        let b = felt_from_uint(34);
        a *= b;
        assert_eq!(a, felt_from_uint(170));

        let a = felt_from_str(
            "7137376184023522026654217343440040540351594155614695530712440441403759244341",
        );
        let b = felt_from_str(
            "18233444644265725414720095600783733811551140822142748479967321742263849032310",
        );

        debug_assert_eq!(a * b, felt_from_uint(1));
        let x = felt_from_str("76312198732308146610779301658447197634344129726058341717334777383093884774635");
        debug_assert_eq!(&x * &x, felt_from_str("0"));
    }

    #[test]
    fn test_edge_cases() {
        // Test zero
        let zero = felt_from_uint(0);
        let non_zero = felt_from_uint(1729);

        assert_eq!(&zero + &zero, zero);
        assert_eq!(&zero - &zero, zero);
        assert_eq!(&zero * &non_zero, zero);
        assert_eq!(&non_zero * &zero, zero);
        assert_eq!(&zero / &non_zero, zero);

        // Test one
        let one = felt_from_uint(1);

        assert_eq!(&non_zero * &one, non_zero);
        assert_eq!(&non_zero / &one, non_zero);

        // Test scalar field order
        let scalar_order = Secp256k1::scalar_field_order();

        assert_eq!(&scalar_order + &zero, zero);
        assert_eq!(&scalar_order - &scalar_order, zero);

        // Test large values
        let large_value = felt_from_str("1234567890123456789012345678901234567890");
        let another_large_value = felt_from_str("9876543210987654321098765432109876543210");

        let expected_product = felt_from_str(
            "35093743783979003143549847474448534855363696641363699189191142970683623002285",
        );
        assert_eq!(&large_value * &another_large_value, expected_product);
    }

    #[test]
    fn test_division() {
        // Test cases with known inverses
        let a = FieldElement::gen_random();
        let b = FieldElement::gen_random();
        let c = FieldElement::gen_random();
        let d = FieldElement::gen_random();
        let e = FieldElement::gen_random();
        let f = FieldElement::gen_random();
        assert_ne!(&d / &e, d);
        assert_ne!(&e / &d, e);

        assert_ne!(&e / &f, e);
        assert_ne!(&f / &e, f);

        assert_ne!(&d / &f, d);
        assert_ne!(&f / &d, f);

        // Test division by one and division by self
        let one = felt_from_uint(1);

        assert_eq!(&a / &one, a);
        assert_eq!(&b / &one, b);
        assert_eq!(&c / &one, c);
        assert_eq!(&d / &one, d);
        assert_eq!(&e / &one, e);
        assert_eq!(&f / &one, f);

        assert_eq!(&a / &a, one);
        assert_eq!(&b / &b, one);
        assert_eq!(&c / &c, one);
        assert_eq!(&d / &d, one);
        assert_eq!(&e / &e, one);
        assert_eq!(&f / &f, one);

        let a = felt_from_str(
            "18233444644265725414720095600783733811551140822142748479967321742263849032310",
        );
        let b = felt_from_str(
            "12158399299693830322967808612713398636155367887041628176798871954788371653930",
        );

        assert_eq!(
            &a / &b,
            felt_from_str(
                "3581453595367386406056576018017432656403623753535430210497256105630064126918"
            )
        );
    }

    #[test]
    fn test_exps() {
        let a = FieldElement::gen_random();
        let b = felt_from_str(
            "57896044618658097711785492504343953926634992332820282019728792003954417335808",
        );
        let c = felt_from_str("987654321098765432109876543210");
        let d = FieldElement::gen_random();
        let e = felt_from_str(
            "115792089237316195423570985008687907853269984665640564039457584007908834671657",
        );

        // Test expression: (a + b) * c - d / e
        let expected_result1 = ((&a + &b) * &c) - (&d / &e);
        let actual_result1 = (&a * &c + &b * &c) - &d / &e;
        assert_eq!(actual_result1, expected_result1);

        // Test expression: (c + d) / (a - b)
        let expected_result2 = (&c + &d) / (&a - &b);
        let actual_result2 = (&c / &(&a - &b)) + (&d / &(&a - &b));
        assert_eq!(actual_result2, expected_result2);

        // random exps
        let exp = -felt_from_uint(45) * felt_from_uint(2) + felt_from_uint(67);
        debug_assert_eq!(
            exp,
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671640"
            )
        );
        let exp = felt_from_uint(25) * felt_from_uint(5) - felt_from_uint(32) + felt_from_uint(21);
        debug_assert_eq!(exp, felt_from_uint(114));

        let lhs = felt_from_uint(45) - felt_from_uint(30) * felt_from_uint(3);

        debug_assert_eq!(
            lhs,
            felt_from_str(
                "115792089237316195423570985008687907853269984665640564039457584007908834671618"
            )
        )
    }
}
