param(
    [string]$Version = '2.1.0',
    [string]$ReleaseDate = '20260828',
    [string]$Profile = 'example_data'
)

$ErrorActionPreference = 'Stop'
$sourceDir = (Resolve-Path -LiteralPath $PSScriptRoot).Path
$releaseRoot = Join-Path (Split-Path $sourceDir -Parent) 'releases'
$packageName = "EELS-ZEST_v${Version}_${ReleaseDate}_${Profile}"
$stageDir = Join-Path $releaseRoot $packageName
$packageDir = Join-Path $stageDir 'EELS-ZEST'
$zipFile = Join-Path $releaseRoot ($packageName + '.zip')
$zipChecksumFile = $zipFile + '.sha256'

New-Item -ItemType Directory -Force -Path $releaseRoot | Out-Null
$releaseRootResolved = (Resolve-Path -LiteralPath $releaseRoot).Path
if (Test-Path -LiteralPath $stageDir) {
    $stageResolved = (Resolve-Path -LiteralPath $stageDir).Path
    $prefix = $releaseRootResolved.TrimEnd('\') + '\'
    if (-not $stageResolved.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to remove staging path outside release root: $stageResolved"
    }
    Remove-Item -LiteralPath $stageResolved -Recurse -Force
}
if (Test-Path -LiteralPath $zipFile) {
    $zipResolved = (Resolve-Path -LiteralPath $zipFile).Path
    $prefix = $releaseRootResolved.TrimEnd('\') + '\'
    if (-not $zipResolved.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to remove archive outside release root: $zipResolved"
    }
    Remove-Item -LiteralPath $zipResolved -Force
}
if (Test-Path -LiteralPath $zipChecksumFile) {
    $checksumResolved = (Resolve-Path -LiteralPath $zipChecksumFile).Path
    $prefix = $releaseRootResolved.TrimEnd('\') + '\'
    if (-not $checksumResolved.StartsWith($prefix, [StringComparison]::OrdinalIgnoreCase)) {
        throw "Refusing to remove checksum outside release root: $checksumResolved"
    }
    Remove-Item -LiteralPath $checksumResolved -Force
}

New-Item -ItemType Directory -Force -Path $packageDir | Out-Null

$rootFiles = @(
    'README.md', 'README_zh-CN.md',
    'VERSION.txt', 'VERSION_NOTES.md', 'NOTICE.md',
    'PACKAGE_CONTENTS.md', 'RELEASE_AUDIT.md',
    'LICENSE', 'CITATION.cff',
    '.gitignore', 'build_release_package.ps1',
    'index.html', 'logo.png', 'intmodes.mat', 'PK1.mat'
)
foreach ($name in $rootFiles) {
    $src = Join-Path $sourceDir $name
    if (-not (Test-Path -LiteralPath $src)) { throw "Missing release file: $name" }
    Copy-Item -LiteralPath $src -Destination $packageDir
}

# MATLAB runtime and launchers.
$runtimeFiles = @(
    'DM3Import.m', 'DM3Import0.m', 'EELS.m',
    'eels_fit_datacube.m', 'eels_fit_spectrum.m',
    'EELS_PearsonSparse_FullFit.m',
    'EELS_TwoStage_Toolbox.m', 'EELS_ZEST_App.m',
    'EELS_ZEST_DataOps.m', 'EELS_ZEST_Export.m',
    'EELS_ZEST_CutoffDialog.m', 'EELS_ZEST_ExportDialog.m',
    'EELS_ZEST_StackTool.m', 'EELS_ZEST_PeakFit.m', 'EELS_ZEST_PeakFitTool.m',
    'run_eels_twostage_toolbox.m', 'run_eels_zest.m'
)
foreach ($name in $runtimeFiles) {
    $src = Join-Path $sourceDir $name
    if (-not (Test-Path -LiteralPath $src)) { throw "Missing runtime file: $name" }
    Copy-Item -LiteralPath $src -Destination $packageDir
}

# Example scripts and bundled data.
$examplesOut = Join-Path $packageDir 'examples'
New-Item -ItemType Directory -Force -Path $examplesOut | Out-Null
$exampleFiles = @(
    'README.md',
    'zest_gui_launch.m',
    'zest_single_stage1_background_removal.m',
    'zest_single_stage2_thermometry.m',
    'zest_datacube_full_stage1.m',
    'zest_datacube_full_stage2_parallel.m',
    'zest_reuse_stage1_for_stage2.m',
    'zest_sparse_intensity_analysis_from_result.m',
    'zest_synthetic_smoke.m'
)
foreach ($name in $exampleFiles) {
    $src = Join-Path $sourceDir "examples\$name"
    if (-not (Test-Path -LiteralPath $src)) { throw "Missing public example: $name" }
    Copy-Item -LiteralPath $src -Destination $examplesOut
}
$exampleDataOut = Join-Path $examplesOut 'data'
New-Item -ItemType Directory -Force -Path $exampleDataOut | Out-Null
$exampleDataSource = Join-Path $sourceDir 'examples\data\example_data.mat'
if (-not (Test-Path -LiteralPath $exampleDataSource)) {
    throw "Missing example data: $exampleDataSource"
}
Copy-Item -LiteralPath $exampleDataSource -Destination $exampleDataOut

# Current documentation and the bundled example figures.
$docsOut = Join-Path $packageDir 'docs'
New-Item -ItemType Directory -Force -Path $docsOut | Out-Null
$docFiles = @(
    'README.md',
    'GUI_TOOLBOX_GUIDE.md',
    'LOSS_AND_ATTENTION_THEORY_v1.3.md',
    'theory.html',
    'EELS_ZEST_2p1_workflow.png', 'EELS_ZEST_2p1_workflow.svg'
)
foreach ($name in $docFiles) {
    $src = Join-Path $sourceDir "docs\$name"
    if (-not (Test-Path -LiteralPath $src)) { throw "Missing public document: $name" }
    Copy-Item -LiteralPath $src -Destination $docsOut
}
$exampleCase = Join-Path $sourceDir 'docs\example_data_figures'
if (-not (Test-Path -LiteralPath $exampleCase)) {
    throw "Missing example figure directory: $exampleCase"
}
Copy-Item -LiteralPath $exampleCase -Destination $docsOut -Recurse

# Keep internal numbered dataset labels out of the release.
$textExtensions = @('.m', '.md', '.html', '.txt', '.csv', '.ps1', '.py', '.json', '.xml', '.yml', '.yaml')
$disallowedPattern = '(?i)\bdatacube\d+\b'
$numberedDatasetHits = @()
Get-ChildItem -LiteralPath $packageDir -Recurse -File |
    Where-Object { $textExtensions -contains $_.Extension.ToLowerInvariant() } |
    ForEach-Object {
        $matches = Select-String -LiteralPath $_.FullName -Pattern $disallowedPattern -AllMatches
        foreach ($match in $matches) {
            $numberedDatasetHits += [pscustomobject]@{
                File = $match.Path.Substring($packageDir.Length + 1)
                Line = $match.LineNumber
                Text = $match.Line.Trim()
            }
        }
    }
if ($numberedDatasetHits.Count -gt 0) {
    $details = ($numberedDatasetHits | Format-Table -AutoSize | Out-String).Trim()
    throw "Numbered internal dataset label found in release:`n$details"
}

$numberedDatacubes = @(Get-ChildItem -LiteralPath $packageDir -Recurse -File |
    Where-Object { $_.Name -match '(?i)^datacube\d' })
if ($numberedDatacubes.Count -gt 0) {
    throw "Numbered internal dataset filename found in release: $($numberedDatacubes.FullName -join ', ')"
}

# Generate deterministic file inventory and hashes. Manifest and checksum files
# exclude themselves to avoid recursive hash changes.
$manifestRows = Get-ChildItem -LiteralPath $packageDir -Recurse -File |
    Sort-Object FullName |
    ForEach-Object {
        $relative = $_.FullName.Substring($packageDir.Length + 1).Replace('\', '/')
        [pscustomobject]@{
            Path = $relative
            Bytes = $_.Length
            SHA256 = (Get-FileHash -LiteralPath $_.FullName -Algorithm SHA256).Hash.ToLowerInvariant()
        }
    }
$manifestRows | Export-Csv -LiteralPath (Join-Path $packageDir 'RELEASE_MANIFEST.csv') -NoTypeInformation -Encoding utf8

$checksumLines = $manifestRows | ForEach-Object { "{0}  {1}" -f $_.SHA256, $_.Path }
[IO.File]::WriteAllLines((Join-Path $packageDir 'SHA256SUMS.txt'), $checksumLines, [Text.UTF8Encoding]::new($false))

Compress-Archive -LiteralPath $packageDir -DestinationPath $zipFile -CompressionLevel Optimal

$zipHash = (Get-FileHash -LiteralPath $zipFile -Algorithm SHA256).Hash.ToLowerInvariant()
[IO.File]::WriteAllText(
    $zipChecksumFile,
    ("{0}  {1}`n" -f $zipHash, (Split-Path $zipFile -Leaf)),
    [Text.UTF8Encoding]::new($false))
$fileCount = (Get-ChildItem -LiteralPath $packageDir -Recurse -File).Count
$bytes = (Get-ChildItem -LiteralPath $packageDir -Recurse -File | Measure-Object Length -Sum).Sum

[pscustomobject]@{
    PackageDirectory = $packageDir
    Archive = $zipFile
    ArchiveChecksum = $zipChecksumFile
    Files = $fileCount
    UncompressedMB = [math]::Round($bytes / 1MB, 2)
    ArchiveMB = [math]::Round((Get-Item -LiteralPath $zipFile).Length / 1MB, 2)
    ArchiveSHA256 = $zipHash
}
