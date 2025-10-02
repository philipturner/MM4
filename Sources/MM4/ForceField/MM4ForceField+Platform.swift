//
//  MM4ForceField+Platform.swift
//  MM4
//
//  Created by Philip Turner on 10/2/25.
//

import OpenMM

extension MM4ForceField {
  // Lazily initialized if users don't load OpenMM in client code and specify
  // a platform to prove they loaded OpenMM. Enables a more ergonomic API for
  // scripting workflows, where boilerplate for plugin loading no longer needs
  // to occur at the top of each script.
  //
  // It is a very easy hazard that OpenMM could fall back to the "Reference"
  // platform, without the user ever knowing. This could be especially
  // problematic for people without much experience using or benchmarking
  // OpenMM.
  nonisolated(unsafe)
  private static var defaultPlatform: OpenMM_Platform?
  
  private static func createDefaultPlatform() -> OpenMM_Platform {
    let pluginsDirectory = OpenMM_Platform.defaultPluginsDirectory
    guard let pluginsDirectory else {
      fatalError("Could not find the OpenMM plugins directory.")
    }
    
    #if os(macOS)
    let pluginFile = pluginsDirectory + "/" + "libOpenMMOpenCL.dylib"
    #elseif os(Windows)
    let pluginFile = pluginsDirectory + "/" + "OpenMMOpenCL.dll"
    #else
    #error("Linux is no longer supported.")
    #endif
    OpenMM_Platform.loadPluginLibrary(file: pluginFile)
    
    let platforms = OpenMM_Platform.platforms
    for platform in platforms {
      if platform.name == "OpenCL" {
        return platform
      }
    }
    
    fatalError("""
      Could not find the OpenCL platform.
      Plugins directory: \(pluginsDirectory)
      Platform count: \(platforms.count)
      Platforms: \(platforms.map(\.name))
      """)
  }
  
  // If the entered platform is valid, return it. Otherwise, return the
  // default platform.
  static func fallback(_ input: OpenMM_Platform?) -> OpenMM_Platform {
    if let input {
      return input
    }
    
    if let defaultPlatform {
      return defaultPlatform
    }
    
    let output = createDefaultPlatform()
    defaultPlatform = output
    return output
  }
}
