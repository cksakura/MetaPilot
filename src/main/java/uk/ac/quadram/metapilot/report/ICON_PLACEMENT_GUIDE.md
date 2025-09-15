# Icon Placement Guide for MetaPilot Reports

## Icon Files Location

Place your icon files in the `resources/` directory within the report folder:

```
report/
├── index.html
├── resources/
│   ├── favicon.ico          # Browser tab icon (16x16 or 32x32 pixels)
│   ├── logo_metapilot.png   # Main logo (120px width recommended)
│   └── logo_metapilot_icon.png  # Small icon (18px height recommended)
└── reports/
    └── rmd_html/
        └── [generated HTML reports]
```

## Icon File Requirements

### favicon.ico
- **Purpose**: Browser tab icon
- **Format**: ICO format
- **Size**: 16x16 or 32x32 pixels
- **Location**: `resources/favicon.ico`

### logo_metapilot.png
- **Purpose**: Main logo in the header
- **Format**: PNG format (transparent background recommended)
- **Size**: 120px width (height will scale proportionally)
- **Location**: `resources/logo_metapilot.png`

### logo_metapilot_icon.png
- **Purpose**: Small icon in the navigation bar
- **Format**: PNG format (transparent background recommended)
- **Size**: 18px height (width will scale proportionally)
- **Location**: `resources/logo_metapilot_icon.png`

## Customization

To use your own icons:

1. Replace the existing icon files in the `resources/` directory
2. Keep the same filenames to maintain compatibility
3. Ensure the images are optimized for web use (compressed but high quality)
4. Use transparent backgrounds for PNG files for better integration

## Color Scheme

The report now uses the Quadram Institute color scheme:
- **Primary**: Buddha Gold (#B5BD00)
- **Secondary**: Malibu (#8ED1FC)  
- **Text**: Black (#000000)

Your icons should complement these colors for the best visual experience.
