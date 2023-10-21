//
//  MM4Error.swift
//
//
//  Created by Philip Turner on 10/20/23.
//

import Foundation

/// An error returned by the MM4 library.
///
/// All errors thrown in MM4 are of this type, although some may be subclasses.
public class MM4Error: Error {
  /// The description of the error.
  public let description: String
  
  /// For compatibility with the C API, store an explicit C string.
  internal let _cDescripton: UnsafeMutablePointer<CChar>
  
  /// Initialize the error object.
  public init(description: String) {
    self.description = description
    self._cDescripton = .allocate(capacity: description.count)
    strcpy(_cDescripton, description)
  }
  
  /// Deinitialize the error object.
  deinit {
    _cDescripton.deallocate()
  }
}
