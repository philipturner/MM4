import XCTest
import MM4

// Complex functions for mutating state, which change the position and/or
// velocity. The tests cover several edge cases. For example, angular velocity
// should rotate with the atom positions when the object is rotated. Test some
// angles besides trivial 90 degree angles, which are only possible with the
// underlying quaternion and cross product math.
//
// Rigorously test the copy-on-write semantics. Ensure changes in one object
// don't cause changes to another copy created earlier, and vice versa. Such
// bugs would be extremely nasty and time-consuming to track down.
