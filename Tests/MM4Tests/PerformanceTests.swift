//
//  PerformanceTests.swift
//
//
//  Created by Philip Turner on 10/18/23.
//

import Foundation

// The initial tests should be for forcefield correctness. At another time, in
// the distant future, it will be optimized for performance, removing the O(n^2)
// scaling.
//
// TODO: Setup profiling tests. See whether multicore CPU could help for
// setting up quick energy minimizations of large nanosystems. The test suite
// should be compiled in release mode for performance tests, otherwise they
// will be disabled using '#if DEBUG'.
