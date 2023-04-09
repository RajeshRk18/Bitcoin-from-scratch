#![allow(unused_variables, dead_code, unused_imports)]
use crate::ecc::PRIME;

use num_bigint::{BigInt, ToBigInt};
use num_traits::{WrappingAdd, FromPrimitive, ToPrimitive};
use std::collections::HashSet;
use std::ops::{Add, Sub, Mul, Div, MulAssign, Rem};

#[derive(Debug, Clone)]
pub struct FieldElement(pub BigInt);

#[derive(Debug, Clone)]
pub struct FiniteField {
    pub elements: Vec<FieldElement>,
    pub order: BigInt
}

impl FiniteField {
    pub fn new(order: BigInt) -> Self {
        let elements = {
            let mut field = Vec::new();
            let mut count = Self::zero();
            while count < order {
                field.push(FieldElement(count.clone()));
                count += Self::one();
            }
            field
        };

        Self { elements, order }
    }

    pub fn wrapping_add(&self, a: FieldElement, b: FieldElement) -> FieldElement {
        ((a.0 + b.0) % self.order.clone() ).into()
    }

    fn new_bigint(a: &str) -> BigInt {
        BigInt::parse_bytes(a.as_bytes(), 10).expect("Error parsing the value into BigInt")
    }

    pub fn wrapping_mul(&self, a: FieldElement, b: FieldElement) -> FieldElement {
        ((a.0 * b.0) % self.order.clone()).into()
    }

    fn divmod(&self, a: FieldElement) -> FieldElement {
        (a.0 % self.order.clone()).into()
    }

    fn zero() -> BigInt {
        BigInt::from_u8(0).unwrap()
    }

    fn one() -> BigInt {
        BigInt::from_u8(1).unwrap()
    }

    pub fn find_generators(&self) -> Vec<FieldElement> {
        let mut generators = Vec::<FieldElement>::new();

        for gen in self.elements.iter().skip(1) {
            dbg!(&gen);
            let mut gen_elements = HashSet::<BigInt>::new();
            let mut power = Self::zero();
            while power < self.order.clone()  {
                if power == Self::zero() {
                    gen_elements.insert(Self::one().clone());
                    continue;
                }
                let multiply = self.square_and_multiply(gen.clone(), FieldElement(power.clone()));
                dbg!(&multiply);
                if !gen_elements.insert(multiply) {
                    break;
                }
                power += Self::one();
            }

            if BigInt::from_usize(gen_elements.len()).unwrap() == self.order.clone() - Self::one() {
                generators.push(gen.clone());
            } else {
                continue;
            }
        }
        generators
    }

    fn square_and_multiply(&self, base: FieldElement, mut power: FieldElement) -> BigInt {
        let mut result = FieldElement(BigInt::from_u8(1).unwrap());
        let mut base = self.divmod(base);

        while power.0 > Self::zero() {
            if result.0.clone() % Self::new_bigint("2") == Self::one() {
                result = self.wrapping_mul(result.clone(), base.clone());
            }
            base = self.wrapping_mul(base.clone(), base.clone());
            power = (power.0 / Self::new_bigint("2")).into();
        }

        result.0
    }
}

impl Add for FieldElement {
    type Output = FieldElement;
    fn add(self, rhs: Self) -> Self::Output {
        ((self.0 + rhs.0) % PRIME.to_bigint().unwrap()).into()
    }
}

impl Add for &FieldElement {
    type Output = FieldElement;
    fn add(self, rhs: Self) -> Self::Output {
        ((self.0.clone() + rhs.0.clone()) % PRIME.to_bigint().unwrap()).into()
    }
}

impl Sub for FieldElement {
    type Output = FieldElement;
    fn sub(self, rhs: Self) -> Self::Output {
        ((self.0 - rhs.0) % PRIME.to_bigint().unwrap()).into()
    }
}

impl Sub for &FieldElement {
    type Output = FieldElement;
    fn sub(self, rhs: Self) -> Self::Output {
        ((self.0.clone() - rhs.0.clone()) % PRIME.to_bigint().unwrap()).into()
    }
}

impl Mul for FieldElement {
    type Output = FieldElement;
    fn mul(self, rhs: Self) -> Self::Output {
        ((self.0 * rhs.0) % PRIME.to_bigint().unwrap()).into()
    }
}

impl Mul for &FieldElement {
    type Output = FieldElement;
    fn mul(self, rhs: Self) -> Self::Output {
        ((self.0.clone() * rhs.0.clone()) % PRIME.to_bigint().unwrap()).into()
    }
}

impl Div for FieldElement {
    type Output = FieldElement;
    fn div(self, rhs: Self) -> Self::Output {
        let div = self.0.div(rhs.0);
        (div % PRIME.to_bigint().unwrap()).into()
    }
}

impl Div for &FieldElement {
    type Output = FieldElement;
    fn div(self, rhs: Self) -> Self::Output {
        let div = self.0.clone().div(rhs.0.clone());
        (div % PRIME.to_bigint().unwrap()).into()
    }
}

impl MulAssign for FieldElement {
    fn mul_assign(&mut self, rhs: Self) {
        self.0 *= rhs.0;
    }
}

impl PartialEq for FieldElement {
    fn eq(&self, other: &Self) -> bool {
        self.0 == other.0
    }
}

impl From<BigInt> for FieldElement {
    fn from(value: BigInt) -> Self {
        FieldElement(value)
    }
}

impl Rem for FieldElement {
    type Output = FieldElement;
    fn rem(self, rhs: Self) -> Self::Output {
        (self.0 % rhs.0).into()
    }
}

impl From<FieldElement> for u32 {
    fn from(value: FieldElement) -> Self {
        value.0.to_u32().unwrap()
    }
}

pub fn add<T>(a: T, b: T) -> T
where
    T: WrappingAdd,
{
    a.wrapping_add(&b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;
    proptest! {
        #[test]
        fn test_add(a: i64, b: i64) {
        assert_eq!(add(a, b), a.wrapping_add(b));
        }
    }

    #[test]
    fn field_element() {
        let big_int = BigInt::parse_bytes(
            "115792089237316195423570985008687907852837564279074904382605163141518161494337"
                .as_bytes(),
            10,
        )
        .unwrap();
        println!("{:?}", big_int.bits());
    }

    #[test]
    fn generator() {
        let prime = FiniteField::new(BigInt::from_u32(23).unwrap());
        let generators: Vec<u32> = prime.find_generators().into_iter().map(|elem| u32::from(elem)).collect();

        assert_eq!(vec![2, 3, 4, 5, 6, 8, 9, 12, 13, 16, 18], generators);
    }
}
