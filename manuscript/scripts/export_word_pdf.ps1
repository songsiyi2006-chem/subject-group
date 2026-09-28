param([string]$DocumentDirectory = (Join-Path $PSScriptRoot '../documents'),
      [string[]]$Names = @('manuscript_english', 'manuscript_chinese', 'supporting_information_english', 'supporting_information_chinese'))
# Uses an existing Microsoft Word installation. Opens only the four named files.
$ErrorActionPreference = 'Stop'
$documentRoot = (Resolve-Path -LiteralPath $DocumentDirectory).Path
$allowedNames = @('manuscript_english', 'manuscript_chinese', 'supporting_information_english', 'supporting_information_chinese')
foreach ($entry in $Names) { if ($entry -notin $allowedNames) { throw 'Unexpected document name' } }
$word = New-Object -ComObject Word.Application
$word.Visible = $false
$word.DisplayAlerts = 0
$records = @()
try {
    foreach ($name in $Names) {
        $inputPath = Join-Path $documentRoot ($name + '.docx')
        $outputPath = Join-Path $documentRoot ($name + '.pdf')
        $document = $word.Documents.Open($inputPath, $false, $true)
        try {
            $document.Repaginate()
            $document.Fields.Update() | Out-Null
            $document.ExportAsFixedFormat($outputPath, 17)
            $records += [pscustomobject]@{name=$name; pages=$document.ComputeStatistics(2); bytes=(Get-Item -LiteralPath $outputPath).Length}
        } finally {
            $document.Close(0)
            [Runtime.InteropServices.Marshal]::ReleaseComObject($document) | Out-Null
        }
    }
} finally {
    $word.Quit()
    [Runtime.InteropServices.Marshal]::ReleaseComObject($word) | Out-Null
}
$records | ConvertTo-Json
