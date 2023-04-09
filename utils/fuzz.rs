#![no_main]
use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
#[derive(Arbitrary, Debug)]
struct Input {
 s: String,
}

fuzz_target!(|input: Input| {
 use fuzzme::parse_integer;
 parse_integer(&input.s);
});

use proptest::prelude::*;
proptest! {
 #[test]
 fn test_quicksort_proptest(
 vec in prop::collection::vec(prop::num::i64::ANY, 0..1000)
 ) {
    use quicksort_proptest::Quicksort;
    let mut vec_sorted = vec.clone();
    vec_sorted.sort();
    let mut vec_quicksorted = vec.clone();
    vec_quicksorted.quicksort();
    assert_eq!(vec_quicksorted, vec_sorted);
 }
}

   // square and multiply
        /*for bit in (0..64).rev() {
            if (power.clone().shr(FieldElement(bit))).bitand(FieldElement(1)) == FieldElement(1) {
                result = self.divmod(result.mul(base.clone()));
            }

            base.mul_assign(base.clone());
            base = self.divmod(base);
        }

        result.0*/
        /*let mut a = g;
            let mut b = FieldElement(1);

            while power.clone().gt(&FieldElement(0)) {
                if power.clone().rem(FieldElement(2)) == FieldElement(1) {
                    b.mul_assign(a.clone());
                    b = self.divmod(b);
                    a = self.divmod(a.clone().mul(a));
                    power = power.clone().div(FieldElement(2));

                } else {
                    a = self.divmod(a.clone().mul(a));
                    power.div_assign(FieldElement(2));
                }
            }

            b.0
        }*/