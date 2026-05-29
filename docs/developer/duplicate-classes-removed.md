# Duplicate Classes Between metapilot-core and metapilot-windows-gui

## Summary

`metapilot-windows-gui` depends on `metapilot-core` but contained duplicate copies of some classes with the same package and class name. Those duplicates were removed so the GUI uses the single definition from core.

Current version target: `1.0.0`.

## Removed from metapilot-windows-gui (use metapilot-core instead)

### Package `uk.ac.quadram.metapilot.quant.flashLFQ`
- `FlashLfqPsm.java`
- `FlashLfqPsmReader.java`
- `FlashLfqQuanPepReader.java`
- `FlashLfqQuanPeptide.java`
- `FlashLfqQuanProReader.java`
- `FlashLfqQuanProtein.java`
- `FlashLfqTask.java`
- `package-info.java` (optional; core has one)

### Package `uk.ac.quadram.metapilot.resources`
- `ResourceManager.java`

## Why this is safe

- `metapilot-windows-gui` has a Maven dependency on `metapilot-core`, so these types are available from the core JAR at compile and runtime.
- All GUI code that imported these classes (e.g. `DDAMagPfindTask`, `DDAMagSageTask`, `DDAIgcPFindTask`, `MetaDiaTask`, `MetaDiaTargetTask`) continues to work; the same fully qualified names now resolve to the core module.
- The removed files were either identical to core or effectively the same; no GUI-specific behavior was lost.

## If you need GUI-specific behavior

If the GUI ever needs a different implementation for one of these types, options are:
- Extend the core class in the GUI module and use the subclass where needed, or
- Add a GUI-specific class in a different package (e.g. `...gui.FlashLfqQuanPepReaderGui`).
