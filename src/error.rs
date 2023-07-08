use thiserror::Error;

use crate::fieldelement::FieldElement;

#[derive(Debug, Error)]
#[error("The point {{(0, 1):?}} is not on the curve!")]
pub struct PointError(pub (FieldElement, FieldElement));
