use crate::fieldelement::*;
use crate::traits::EllipticCurve;
use crate::utils::*;
use crate::error::PointError;

use lazy_static::lazy_static;
use num_bigint::{BigUint, BigInt};
use num_traits::FromPrimitive;

lazy_static! {
    pub static ref ORDER: FieldElement = FieldElement::new(
        BigUint::parse_bytes(
            "115792089237316195423570985008687907852837564279074904382605163141518161494337"
                .as_bytes(),
            10
        )
        .unwrap()
    );
    pub static ref PRIME: FieldElement = FieldElement::new(
        BigUint::parse_bytes(
            "115792089237316195423570985008687907853269984665640564039457584007908834671663"
                .as_bytes(),
            10
        )
        .unwrap()
    );
    pub static ref G_X: FieldElement = FieldElement::new(
        BigUint::parse_bytes(
            "55066263022277343669578718895168534326250603453777594175500187360389116729240"
                .as_bytes(),
            10
        )
        .unwrap()
    );
    pub static ref G_Y: FieldElement = FieldElement::new(
        BigUint::parse_bytes(
            "32670510020758816978083085130507043184471273380659243275938904335757337482424"
                .as_bytes(),
            10
        )
        .unwrap()
    );
}

#[derive(Debug, PartialEq, Eq, PartialOrd)]
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
        scalar: BigUint,
        point: <Self as EllipticCurve>::Affine,
    ) -> Result<<Self as EllipticCurve>::Affine, PointError> {
        if !self.is_on_curve(point.clone()) {
            return Err(PointError(point));
        }
        let to_proj = self.to_projective(point);

        let proj_result = self.scalar_mul(scalar, to_proj);
        let proj_result_to_affine = self.to_affine(proj_result);

        Ok(proj_result_to_affine)
    }

    fn is_negative(point1: <Self as EllipticCurve>::Projective, point2: <Self as EllipticCurve>::Projective) -> bool {
        if (&point1.0 == &point2.0) && (&point1.2 == &point2.2) && (point1.1.clone() == -point2.1.clone()) {
            return true;
        }
        false
    }
}

impl EllipticCurve for Secp256k1 {
    type Affine = (FieldElement, FieldElement);
    type Projective = (FieldElement, FieldElement, FieldElement);
    type Scalar = BigUint;

    #[inline(always)]
    fn is_on_curve(&self, point: Self::Affine) -> bool {
        let lhs = point.1.pow(2) % PRIME.clone();
        let rhs: FieldElement = (FieldElement::new(self.degree3_coeff.clone().to_biguint().unwrap()) * point.0.pow(3) + FieldElement::new(self.degree0_coeff.clone().to_biguint().unwrap()))
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
            let z_inv_squared = &z_inv * &z_inv;
            let z_inv_cubed = &z_inv_squared * &z_inv;

            let x = x * z_inv_squared;
            let y = y * z_inv_cubed;

            return (x, y);
        }
    }

    #[inline(always)]
    fn scalar_mul(&self, scalar: BigUint, point: Self::Projective) -> Self::Projective {

        assert!(
            scalar != new_BigUint(0 as u32) && FieldElement::from(scalar.clone()) < ORDER.clone()
        );

        let mut result = (Self::zero(), Self::zero(), Self::zero());
        let mut addend = point.clone();
        let scalar_bits = to_bits(scalar);
        let mut dummy_point = (Self::zero(), Self::zero(), Self::zero());
        for bit in scalar_bits.clone().into_iter().rev() {
            addend = self.point_double(addend.clone());
            dbg!(&addend);
            if bit == 1u8 {
                result = self.point_add(result.clone(), addend.clone());
                dbg!(&result);
            
            } else {
                // Dummy op to make it constant time.
                dummy_point = self.point_add(dummy_point.clone(), addend.clone());
            }
        }
        result
    }

    #[inline(always)]
    fn point_add(&self, point1: Self::Projective, point2: Self::Projective) -> Self::Projective {
        let (x1, y1, z1) = point1.clone();
        let (x2, y2, z2) = point2.clone();

        if Self::is_negative(point1.clone(), point2) {
            return self.point_double(point1);
        }

        let t = &z1 * &z1;
        let t2 = &z2 * &z2;

        let u1 = x1 * &t2;
        let u2 = x2 * &t;
        let u3 = &u1 - &u2;
        let u4 = &y1 * &z2.pow(3);
        let u5 = &y2 * &z1.pow(3);
        let u6 = &u4 - &u5;

        if &u3 == &Self::zero() {
            if &u6 == &Self::zero() {
                return self.point_double(point1.clone());
            } else {
                // Dummy op to make it constant time.
                let dummy_op = self.point_double(point1);
                return (Self::zero(), Self::one(), Self::zero());
            }
        }

        let u7 = &u1 + &u2;
        let u8 = &u4 + &u5;
        let z3 = &z1 * &z2 * &u3;
        let x3 = u6.pow(2) - &u7 * &u3.pow(2);
        let u9 = &u7 * &u3.pow(2) - field_element(2) * &x3;
        let y3 = (&u9 * &u6 - &u8 * &u3.pow(3)) / field_element(2);

        (x3, y3, z3)
    }

    #[inline(always)]
    fn point_double(&self, point: Self::Projective) -> Self::Projective {
        if &point.2 == &Self::zero() {
            return (Self::zero(), Self::one(), Self::zero());
        }
        let (x, y, z) = point;

        let u1 = field_element(3) * &x.pow(3);
        let z3 = field_element(2) * &y * &z;
        let u2 = field_element(4) * &x * &y.pow(2);
        let x3 = &u1.pow(2) - &(field_element(2) * &u2);
        let u3 = field_element(8) * &y.pow(4);
        let y3 = u1 * (u2 - &x3) - u3;

        (x3, y3, z3)
    }
}

#[cfg(test)]
mod ecc_tests {
    use super::*;

    #[test]
    fn scalar_test() {
        let curve = Secp256k1::new(7, 1);
        let point = (field_element(5), field_element(7));
        let scalar = new_BigUint(2 as u32);
        let result = curve.mul(scalar, point).unwrap();
        assert_eq!(result, (field_element(2), field_element(5)));
    }

    #[test]
    fn test_vectors() {
        let curve = Secp256k1::new(7, 1);
        let point = (G_X.clone(), G_Y.clone());

        assert_eq!(
            curve.mul(BigUint("2"), point).unwrap(),
            (
                FieldElement::new(BigUint(
                    "89565891926547004231252920425935692360644145829622209833684329913297188986597"
                )),
                FieldElement::new(BigUint(
                    "12158399299693830322967808612713398636155367887041628176798871954788371653930"
                ))
            )
        );
    }

    #[test]
    fn run_test() {
        let curve = Secp256k1::new(7, 1);
        let point = (G_X.clone(), G_Y.clone());
        assert_eq!(
            curve.mul(BigUint("2"), point).unwrap(),
            (
                FieldElement::new(BigUint(
                    "89565891926547004231252920425935692360644145829622209833684329913297188986597"
                )),
                FieldElement::new(BigUint(
                    "12158399299693830322967808612713398636155367887041628176798871954788371653930"
                ))
            )
        );
    }
}
