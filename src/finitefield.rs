/*#![allow(unused_variables, dead_code, unused_imports)]
use crate::ecc::PRIME;
use crate::fieldelement::FieldElement;

use num_bigint::{BigInt, ToBigInt};
use num_traits::{WrappingAdd, FromPrimitive, ToPrimitive};
use std::collections::HashSet;
use std::ops::{Add, Sub, Mul, Div, MulAssign, Rem};

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
        (a + b) % self.order.clone().into()
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

impl From<FieldElement> for u32 {
    fn from(value: FieldElement) -> Self {
        value.0.to_u32().unwrap()
    }
}*/

/*#[inline(always)]
fn point_add(
    &self,
    point1: (FieldElement, FieldElement),
    point2: (FieldElement, FieldElement),
) -> Result<(FieldElement, FieldElement), PointError> {

    if self.is_infinity(point1.clone()) {
        return Ok(point2);
    }
    if self.is_infinity(point2.clone())  {
        return Ok(point1);
    }

    let (x1, y1) = point1;
    let (x2, y2) = point2;

    if (x1, y1) == (x2, -y2) {
        return Ok(Self::point_at_infinity());
    }
    let slope;
    if (&x1 == &x2) && (&y1 == &y2) {
        slope = (field_element(3) * x1.pow(2) / field_element(2) * y1.clone()) % self.order();
    } else {
        slope = ((&y2 - &y1) / (&x2 - &x1)) % self.order();
    }

    let x3 = slope.pow(2).into() - &x1 - &x2;
    let y3 = slope.into() * (&x3 - &x1) + &y1;

    Ok((x3, y3))
}

#[inline(always)]
fn point_double(&self, point: (FieldElement, FieldElement)) -> (FieldElement, FieldElement) {
    if self.is_infinity(point.clone()) {
        return Self::point_at_infinity();
    } else {
        let slope = (field_element(3).unwrap() * point.0.clone().pow(2) / field_element(2).unwrap()) * point.1.clone();
        let x3 = slope.pow(2) - field_element(2) * point.0.clone();
        let y3 = slope * (x3.clone() - point.0) - point.1.clone();
        return (x3, y3);
    }
}

#[inline(always)]
fn scalar_mul(&self, mut scalar: FieldElement, point: (FieldElement, FieldElement)) -> Result<(FieldElement, FieldElement), PointError> {
    if &scalar == &Self::zero() {
        return Ok(Self::point_at_infinity());
    }

    let mut q = point.clone();

    let mut r = Self::point_at_infinity();
    let count = 0;
    while scalar > FieldElement::from_u8(0).unwrap() {
        dbg!(count);
        if scalar.clone() & Self::one() == Self::one() {
            r = self.point_add(r.clone(), q.clone())?;
        }
        scalar >>= 1;
        q = self.point_double(q.clone());
    }

    Ok(r)
}

fn point_add(
        &self,
        point1: (FieldElement, FieldElement),
        point2: (FieldElement, FieldElement),
    ) -> Result<(FieldElement, FieldElement), PointError> {
        if self.is_infinity(point1.clone()) {
            return Ok(point2);
        }
        if self.is_infinity(point2.clone()) {
            return Ok(point1);
        }

        let (x1, y1) = point1;
        let (x2, y2) = point2;
        if x1 == x2 && y1 != y2 {
            return Self::point_at_infinity();
        } else {

        }
        let slope = (&y2 - &y1) * self.mod_inv(&x2 - &x1, self.order()).unwrap() % self.order();
        let x3 = (slope.pow(2) - &x1 - &x2) % self.order();
        let y3 = (slope * (&x1 - &x3) - &y1) % self.order();
        Ok((x3, y3))
}
fn point_double(&self, point: (FieldElement, FieldElement)) -> (FieldElement, FieldElement) {
    let (x, y) = point.clone();
    if self.is_infinity(point) || y == Self::zero() {
        return (Self::zero(), Self::zero());
    }
    let slope = (FieldElement::from_u8(3).unwrap() * x.pow(2) / FieldElement::from_u8(2).unwrap() * &y) % self.order();
    let x3 = (slope.pow(2) - FieldElement::from_u8(2).unwrap() * &x) % self.order();
    let y3 = (slope * (&x - &x3) - &y) % self.order();
    (x3, y3)
}

fn scalar_mul(&self, scalar: FieldElement, point: (FieldElement, FieldElement)) -> Result<(FieldElement, FieldElement), PointError> {
    let mut result = Self::point_at_infinity();
    let mut addend = point.clone();
    let mut scalar = scalar.clone();
    while scalar != 0.to_bigint().unwrap() {
        if scalar.clone() % 2.to_bigint().unwrap() == 1.to_bigint().unwrap() {
            result = self.point_add(result.clone(), addend.clone())?;
        }
        addend = self.point_double(addend.clone());
        scalar = scalar >> 1;
    }
    Ok(result)
}*/

/*//write ext_euclid function
fn ext_euclid(a: FieldElement, b: FieldElement) -> (FieldElement, FieldElement, FieldElement) {
    if b == 0.to_bigint().unwrap().into() {
        return (a.clone(), 1.to_bigint().unwrap().into(), 0.to_bigint().unwrap().into());
    }
    let (gcd, x, y) = Self::ext_euclid(b.clone(), a.clone() % b.clone());
    (gcd, y.clone(), x - (a / b) * y)
}

//write mod_inv function
fn mod_inv(&self, a: FieldElement, b: FieldElement) -> Option<FieldElement> {
    let (gcd, x, _) = Self::ext_euclid(a, b.clone());

    if gcd != 1.to_bigint().unwrap().into() {
        return None;
    } else {
        return Some((x % b.clone() + b.clone()) % b.clone());
    }
}*/

/*impl TestCurve {
     fn ext_euclid(a: FieldElement, b: FieldElement) -> Projective {
         if b == 0i8.to_bigint().unwrap() {
             return (a.clone(), 1i8.to_bigint().unwrap(), 0i8.to_bigint().unwrap());
         } else {
             let (gcd, x1, y1) = Self::ext_euclid(b.clone(), a.clone() % b.clone());
             let x = y1.clone();
             let y = x1.clone() - (a / b) * y1.clone();
             return (gcd, x, y);
         }
     }

     fn mod_inv(&self, a: FieldElement, b: FieldElement) -> Option<FieldElement> {
         let (gcd, x, _) = Self::ext_euclid(a, b.clone());

         if gcd != 1.to_bigint().unwrap() {
             return None;
         } else {
             return Some((x % b.clone() + b.clone()) % b.clone());
         }
     }
 }

 impl EllipticCurve for TestCurve {
     fn is_on_curve(&self, point: (FieldElement, FieldElement)) -> bool {
         let lhs = point.1.pow(2) % FieldElement::from_i128(self.order).unwrap();
         let rhs = (point.0.pow(3) + self.b) % FieldElement::from_i128(self.order).unwrap();
         if lhs != rhs {
             return false;
         }
         true
     }

     fn check_point_not_at_infinity(
         &self,
         point1: (FieldElement, FieldElement),
         point2: (FieldElement, FieldElement),
     ) -> bool {
         if point1.0 == point2.0 {
             return false;
         }
         true
     }

     fn point_add(
         &self,
         point1: (FieldElement, FieldElement),
         point2: (FieldElement, FieldElement),
     ) -> Result<(FieldElement, FieldElement), PointError> {
         if !self.is_on_curve(point1.clone()) {
             return Err(PointError(point2.0.clone(), point2.1.clone()));
         }
         if !self.is_on_curve(point2.clone()) {
             return Err(PointError(point2.0.clone(), point2.1.clone()));
         }
         let (x1, y1) = point1;
         let (x2, y2) = point2;

         if x1 != x2 {
             let slope = (y2.clone() - y1.clone()) / (x2.clone() - x1.clone());
             let x3 =
                 (slope.pow(2) - x1.clone() - x2.clone()) % FieldElement::from(self.order as isize);
             let y3 = (slope * (x1.clone() - x3.clone()) - y1.clone())
                 % FieldElement::from(self.order as isize);
             return Ok((x3, y3));
         } else if ((x1 == x2) && (y1 == -y2.clone())) || ((x1 == x2) && (y1 != y2)) {
             return Ok((FieldElement::from(1 as usize), FieldElement::from(0 as usize)));
         } else if ((x1 == x2) && (y1 == y2)) || ((y1 == y2) && (x1 != x2)) {
             let slope: FieldElement = (FieldElement::from_u32(3).unwrap() * x1.clone().pow(2)
                 / FieldElement::from_u32(2).unwrap()
                 * y1.clone())
                 % FieldElement::from(self.order as isize);
             let x3: FieldElement = (slope.clone().pow(2) - FieldElement::from_u32(2).unwrap() * x1.clone())
                 % FieldElement::from(self.order as isize);
             let y3: FieldElement = (slope * (x1.clone() - x3.clone()) - y1.clone())
                 % FieldElement::from(self.order as isize);
             return Ok((x3, y3));
         } else {
             return Ok((FieldElement::from(1 as usize), FieldElement::from(0 as usize)));
         }
     }

     fn point_double(&self, point: (FieldElement, FieldElement)) -> (FieldElement, FieldElement) {
         if point.1 == Self::zero() {
             return Self::point_at_infinity();
         } else {
             let slope = (FieldElement::from_u8(3).unwrap() * point.0.clone() + self.degree0_coeff.clone()) / FieldElement::from_u8(2).unwrap() * point.1.clone();
             let x3 = slope.pow(2) - FieldElement::from_u8(2).unwrap() * point.0.clone();
             let y3 = point.1 + slope * (x3.clone() - point.0);
             return (x3, y3);
         }
     }

     fn scalar_mul(&self, mut scalar: FieldElement, point: (FieldElement, FieldElement)) -> Result<(FieldElement, FieldElement), PointError> {
         assert!(scalar.clone() != Self::zero());
         let mut q = (FieldElement::from_u8(0).unwrap(), FieldElement::from_u8(0).unwrap());

         let mut r = point.clone();

         while scalar.clone() > FieldElement::from_u8(0).unwrap() {
             if scalar.clone() % 2 == Self::one() {
                 q = self.point_add(q.clone(), r.clone())?;
             } else {
                 scalar = scalar.clone() >> 1;
                 r = self.point_double(r.clone());
             }
         }

         Ok(q)
     }

 }

 #[test]
 fn in_curve() {
     let test_curve = TestCurve { b: 7, order: 223 };
     assert!(test_curve.is_on_curve((
         FieldElement::from_u32(192).unwrap(),
         FieldElement::from_u32(105).unwrap()
     )));
     assert!(!test_curve.is_on_curve((
         FieldElement::from_u32(200).unwrap(),
         FieldElement::from_u32(119).unwrap()
     )));
 }

 #[test]
 fn add() {
     let test_curve = TestCurve { b: 7, order: 223 };
     debug_assert_eq!(
         test_curve.point_add(
             (
                 FieldElement::from_u32(192).unwrap(),
                 FieldElement::from_u32(105).unwrap()
             ),
             (FieldElement::from_u32(17).unwrap(), FieldElement::from_u32(56).unwrap())
         ),
         Ok((
             FieldElement::from_i32(-209).unwrap(),
             FieldElement::from_i32(-105).unwrap()
         ))
     );
 }

/*#[test]
 fn mul() {
     let test_curve = TestCurve { b: 7, order: 223 };
     debug_assert_eq!(
         test_curve.scalar_mul(
             3,
             (
                 FieldElement::from_u32(192).unwrap(),
                 FieldElement::from_u32(105).unwrap()
             )
         ),
         (
             FieldElement::from_i32(18).unwrap(),
             FieldElement::from_i32(189).unwrap()
         )
     );
 }*/*/

/*for bit in scalar.to_str_radix(2).chars().rev().skip(1) {
    q = self.point_double(q.clone());

    if bit == '1' {
        r = self.point_add(r.clone(), q.clone()).unwrap();
    }
}*/

/*
#[inline(always)]
fn point_add(
    &self,
    point1: (FieldElement, FieldElement),
    point2: (FieldElement, FieldElement),
) -> Result<(FieldElement, FieldElement), PointError> {

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
        slope = (FieldElement::from_u8(3).unwrap() * x1.pow(2) / FieldElement::from_u8(2).unwrap() * &y1) % self.order();
    } else {
        slope = ((&y2 - &y1) / (&x2 - &x1)) % self.order();
    }

    let x3 = (slope.pow(2) - &x1 - &x2) % self.order();
    let y3 = (slope * (&x1 - &x3) - &y1) % self.order();

    Ok((x3, y3))
}
#[inline(always)]
fn point_double(&self, point: (FieldElement, FieldElement)) -> (FieldElement, FieldElement) {
    if self.is_infinity(point.clone()) {
        return Self::point_at_infinity();
    } else {
        let slope = ((FieldElement::from_u8(3).unwrap() * point.0.clone().pow(2)) / FieldElement::from_u8(2).unwrap() * point.1.clone()) % self.order();
        let x3 = (slope.pow(2) - FieldElement::from_u8(2).unwrap() * point.0.clone()) % self.order();
        let y3 = (slope * (x3.clone() - point.0) - point.1.clone()) % self.order();
        return (x3, y3);
    }
}
    #[inline(always)]
fn scalar_mul(&self, mut scalar: FieldElement, point: (FieldElement, FieldElement)) -> Result<(FieldElement, FieldElement), PointError> {
    if &scalar == &Self::zero() {
        return Ok(Self::point_at_infinity());
    }
    let mut q = point.clone();

    let mut r = Self::point_at_infinity();
    let count = 0;
    while scalar > FieldElement::from_u8(0).unwrap() {
        dbg!(count);
        if scalar.clone() & Self::one() == Self::one() {
            r = self.point_add(r.clone(), q.clone())?;
        } else {
            scalar >>= 1;
            q = self.point_double(q.clone());
        }
    }

    Ok(r)
}
 */
