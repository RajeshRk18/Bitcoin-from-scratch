use crate::fieldelement::FieldElement;

pub trait IsEllipticCurve {
    type Affine;
    type Projective;
    type Scalar;

    fn is_valid(&self, point: &Self::Affine) -> bool;

    fn is_indiscriminant() -> bool;

    fn is_negative(point1: &Self::Projective, point2: &Self::Projective) -> bool;

    fn is_identity(p: &Self::Projective) -> bool;

    fn to_projective(&self, point: Self::Affine) -> Self::Projective;

    fn point_add(&self, point1: Self::Projective, point2: Self::Projective) -> Self::Projective;

    fn point_double(&self, point: Self::Projective) -> Self::Projective;

    fn scalar_gen_mul(&self, scalar: Self::Scalar) -> Self::Projective;

    fn to_affine(&self, point: Self::Projective) -> Self::Affine;
}
pub trait Secp256k1Curve: IsEllipticCurve {
    fn a() -> FieldElement; // 0

    fn b() -> FieldElement;

    fn scalar_field_order() -> FieldElement;

    fn base_field_order() -> FieldElement;

    fn one() -> FieldElement;

    fn zero() -> FieldElement;

    fn generator() -> Self::Affine;

    fn identity() -> Self::Projective;
}

pub trait SecSerde {
    fn serialize(x: FieldElement, y: FieldElement) -> Vec<u8>;
    fn deserialize(bytes: &[u8]) -> Result<(FieldElement, FieldElement), String>;
}

pub trait Compression: IsEllipticCurve {
    fn compress(point: Self::Affine) -> Vec<u8>;
    fn decompress(bytes: &[u8]) -> Result<Self::Affine, String>;
}

pub trait KeyPair {
    fn new() -> Self;
    fn private_key(&self) -> FieldElement;
    fn public_key(&self) -> (FieldElement, FieldElement);
}
