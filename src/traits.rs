use crate::fieldelement::FieldElement;

pub trait IsEllipticCurve {
    type Affine;
    type Jacobian;
    type Scalar;

    fn is_valid(&self, point: &Self::Affine) -> bool;

    fn is_indiscriminant() -> bool;

    fn is_negative(point1: &Self::Jacobian, point2: &Self::Jacobian) -> bool;

    fn is_identity(p: &Self::Jacobian) -> bool;

    fn to_jacobian(&self, point: Self::Affine) -> Self::Jacobian;

    fn point_add(point1: &Self::Jacobian, point2: &Self::Jacobian) -> Self::Jacobian;

    fn point_double(point: &Self::Jacobian) -> Self::Jacobian;

    fn scalar_gen_mul(&self, scalar: Self::Scalar) -> Self::Jacobian;

    fn to_affine(&self, point: Self::Jacobian) -> Self::Affine;
}
pub trait Secp256k1Curve: IsEllipticCurve {
    fn a() -> FieldElement; // 0

    fn b() -> FieldElement;

    fn scalar_field_order() -> FieldElement;

    fn base_field_order() -> FieldElement;

    fn generator() -> Self::Affine;

    fn identity() -> Self::Jacobian;
}

pub trait Signable: IsEllipticCurve {
    fn sign<T: AsRef<[u8]>>(&self, data: T) -> (FieldElement, FieldElement);
    fn verify(&self, r: FieldElement, s: FieldElement) -> FieldElement;
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

// Implements zeroize from scratch without depending on zeroize crate
pub trait ZeroIt {
    fn zeroize(&mut self);
}