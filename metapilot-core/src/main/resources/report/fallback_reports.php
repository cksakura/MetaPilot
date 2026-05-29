<?php
// Fallback report listing script
// This can be used when the manifest system fails

$reportsDir = 'reports/';
$reports = [];

// Scan the reports directory for HTML files
if (is_dir($reportsDir)) {
    $files = scandir($reportsDir);
    foreach ($files as $file) {
        if (pathinfo($file, PATHINFO_EXTENSION) === 'html' && $file !== 'index.html') {
            $filePath = $reportsDir . $file;
            $fileTime = filemtime($filePath);
            
            $reports[] = [
                'title' => ucfirst(str_replace(['_', '-'], ' ', pathinfo($file, PATHINFO_FILENAME))),
                'file' => $file,
                'summary' => 'Generated report file',
                'created_at' => date('c', $fileTime),
                'tags' => ['generated']
            ];
        }
    }
}

// Sort by creation time (newest first)
usort($reports, function($a, $b) {
    return strtotime($b['created_at']) - strtotime($a['created_at']);
});

// Create manifest-like structure
$manifest = [
    'generated_at' => date('c'),
    'reports' => $reports
];

// Output as JSON
header('Content-Type: application/json');
echo json_encode($manifest, JSON_PRETTY_PRINT);
?>
