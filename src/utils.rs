//-------------------------------------------------------------------------
// @file utils.rs
//
// @date 01/06/21 11:38:58
// @author Martin Noblia
// @email mnoblia@disroot.org
//
// @brief
//
// @detail
//
//  Licence:
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
//--------------------------------------------------------------------------
use num::Float;

// NOTE(elsuizo:2020-06-02): the following function
// is a translation of the implementation that is here:
// https://floating-point-gui.de/errors/comparison/

/// a comparison function for floating point values
///
pub fn nearly_equal<T: Float>(a: T, b: T, epsilon: T) -> bool {
    let abs_a = a.abs();
    let abs_b = b.abs();
    let abs_diff = (a - b).abs();
    let zero = T::zero();
    // short-cut, handles infinity
    if a == b { true }

    else if a == zero || b == zero || (abs_a + abs_b < T::min_positive_value()) {
        // a or b is zero or both are extremely close to it
        // relative error is less meaningful here
        abs_diff < (epsilon * T::min_positive_value())
    } else {
      abs_diff / T::min(abs_a + abs_b, T::max_value()) < epsilon
    }
}

// TODO(elsuizo:2020-08-31): implement the Display trait for this type
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct MaxMin<T> {
    pub max: (T, usize),
    pub min: (T, usize),
}

/// generic function to fin min, max values and the position in a slice
pub fn find_max_min<T: std::cmp::PartialOrd + Copy>(slice: &[T]) -> MaxMin<T> {
    let mut max = &slice[0];
    let mut min = &slice[0];

    let mut max_pos: usize = 0;
    let mut min_pos: usize = 0;

    for index in 1..slice.len() {
        if slice[index] < *min { min = &slice[index]; min_pos = index;}
        if slice[index] > *max { max = &slice[index]; max_pos = index;}
    }

    MaxMin{max: (*max, max_pos), min: (*min, min_pos)}
}
