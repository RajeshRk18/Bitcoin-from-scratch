use crate::fieldelement::FieldElement;
use num_bigint::BigInt;

pub trait EllipticCurve {
    type Affine;
    type Projective;
    fn is_on_curve(&self, point: Self::Affine) -> bool;

    fn to_projective(&self, point: Self::Affine) -> Self::Projective;

    fn point_add(&self, point1: Self::Projective, point2: Self::Projective) -> Self::Projective;

    fn point_double(&self, point: Self::Projective) -> Self::Projective;

    fn scalar_mul(&self, scalar: BigInt, point: Self::Projective) -> Self::Projective;

    fn to_affine(&self, point: Self::Projective) -> Self::Affine;
}

pub trait FiniteField {
    fn new(order: FieldElement) -> Self;

    fn order(&self) -> FieldElement;

    fn find_generators(&self) -> Vec<FieldElement>;

    fn is_generator(&self, generator: FieldElement) -> bool;

    fn is_in_field(&self, element: FieldElement) -> bool;

    fn add(&self, element1: FieldElement, element2: FieldElement) -> FieldElement;

    fn sub(&self, element1: FieldElement, element2: FieldElement) -> FieldElement;

    fn mul(&self, element1: FieldElement, element2: FieldElement) -> FieldElement;

    fn div(&self, element1: FieldElement, element2: FieldElement) -> FieldElement;

    fn pow(&self, element: FieldElement, exponent: FieldElement) -> FieldElement;

    fn sqrt(&self, element: FieldElement) -> Option<FieldElement>;

    fn inv(&self, element: FieldElement) -> Option<FieldElement>;
}
pub trait SecSerde {
    fn serialize(x: FieldElement, y: FieldElement) -> Vec<u8>;
    fn deserialize(bytes: &[u8]) -> Result<(FieldElement, FieldElement), String>;
}

pub trait Compression {
    fn compress(point: (FieldElement, FieldElement)) -> Vec<u8>;
    fn decompress(bytes: &[u8]) -> Result<(FieldElement, FieldElement), String>;
}
pub trait KeyPair {
    fn new() -> Self;
    fn private_key(&self) -> FieldElement;
    fn public_key(&self) -> (FieldElement, FieldElement);
}
