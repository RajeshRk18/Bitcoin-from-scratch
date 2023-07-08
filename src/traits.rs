use crate::fieldelement::FieldElement;

pub trait EllipticCurve {
    type Affine;
    type Projective;
    type Scalar;

    fn is_on_curve(&self, point: Self::Affine) -> bool;

    fn to_projective(&self, point: Self::Affine) -> Self::Projective;

    fn point_add(&self, point1: Self::Projective, point2: Self::Projective) -> Self::Projective;

    fn point_double(&self, point: Self::Projective) -> Self::Projective;

    fn scalar_mul(&self, scalar: Self::Scalar, point: Self::Projective) -> Self::Projective;

    fn to_affine(&self, point: Self::Projective) -> Self::Affine;
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
