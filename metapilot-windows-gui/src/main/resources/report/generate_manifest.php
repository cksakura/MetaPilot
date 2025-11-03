<?php
header('Content-Type: application/json');
header('Cache-Control: no-cache, must-revalidate');

// Function to get report title from filename
function getReportTitle($filename) {
    $titles = [
        'report_ID_summary.html' => 'Identification Summary',
        'report_peptides_summary.html' => 'Peptides Summary', 
        'report_proteinGroups_summary.html' => 'Proteins Summary',
        'report_taxonomy_summary.html' => 'Taxonomy Summary',
        'report_function_summary.html' => 'Function Summary',
        'bubble_chart_report.html' => 'Bubble Chart Report',
        'Genome.html' => 'Genome Report',
        'functions.html' => 'Functions Report'
    ];
    return $titles[$filename] ?? str_replace(['_', '.html'], [' ', ''], $filename);
}

// Function to get report summary from filename
function getReportSummary($filename) {
    $summaries = [
        'report_ID_summary.html' => 'Summary of identification results including PSM counts and identification rates',
        'report_peptides_summary.html' => 'Detailed analysis of identified peptides with quantification data',
        'report_proteinGroups_summary.html' => 'Protein identification and quantification results',
        'report_taxonomy_summary.html' => 'Taxonomic analysis showing organism distribution',
        'report_function_summary.html' => 'Functional annotation of identified proteins',
        'bubble_chart_report.html' => 'Interactive hierarchical visualization of taxonomic data',
        'Genome.html' => 'Genome-level analysis and statistics',
        'functions.html' => 'Functional annotation results'
    ];
    return $summaries[$filename] ?? 'Generated report from MetaPilot analysis';
}

// Function to get report tags from filename
function getReportTags($filename) {
    $tags = [
        'report_ID_summary.html' => ['identification', 'summary'],
        'report_peptides_summary.html' => ['peptides', 'quantification'],
        'report_proteinGroups_summary.html' => ['proteins', 'quantification'],
        'report_taxonomy_summary.html' => ['taxonomy', 'organisms'],
        'report_function_summary.html' => ['functions', 'annotation'],
        'bubble_chart_report.html' => ['bubble-chart', 'taxonomy', 'visualization'],
        'Genome.html' => ['genome', 'taxonomy'],
        'functions.html' => ['functions', 'annotation']
    ];
    return $tags[$filename] ?? ['report'];
}

// Scan the reports directory
$reportsDir = __DIR__ . '/reports';
$reports = [];

if (is_dir($reportsDir)) {
    $files = scandir($reportsDir);
    
    foreach ($files as $file) {
        if (pathinfo($file, PATHINFO_EXTENSION) === 'html') {
            $filePath = $reportsDir . '/' . $file;
            $fileTime = filemtime($filePath);
            
            $reports[] = [
                'title' => getReportTitle($file),
                'file' => $file,
                'summary' => getReportSummary($file),
                'created_at' => date('c', $fileTime),
                'tags' => getReportTags($file)
            ];
        }
    }
}

// Create manifest
$manifest = [
    'generated_at' => date('c'),
    'version' => '1.0',
    'reports' => $reports
];

// Output JSON
echo json_encode($manifest, JSON_PRETTY_PRINT);
?>
