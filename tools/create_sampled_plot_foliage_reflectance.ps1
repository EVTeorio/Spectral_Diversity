Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$outDir = Join-Path $root "reports/tables/pca_loading_spectral_regions"
$outPath = Join-Path $outDir "sampled_plot_foliage_reflectance_50m.csv"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

$wavelengths = @(398..998 | Where-Object { (($_ - 398) % 5) -eq 0 })
$manualExcluded = @("sub50_80", "sub50_79", "sub50_71", "sub50_70", "sub50_62", "sub50_53")
$rasterDir = Join-Path $root "Quad_Spectra/50m_smooth_5nm"
$files = @(Get-ChildItem -Path $rasterDir -File | Where-Object {
  $_.Extension -eq "" -and $manualExcluded -notcontains $_.Name
} | Sort-Object Name | Select-Object -First 12)

$samplePerRaster = 20
$maxCandidates = 1200
$samples = 680
$lines = 680
$pixels = $samples * $lines
$bands = 121
$bandBytes = $pixels * 4
$band563 = [int]((563 - 398) / 5)
$shadowThreshold = 0.0305476

$sums = New-Object double[] $bands
$count = 0
$rng = [System.Random]::new(42)

foreach ($file in $files) {
  $reader = [System.IO.BinaryReader]::new([System.IO.File]::Open($file.FullName, [System.IO.FileMode]::Open, [System.IO.FileAccess]::Read, [System.IO.FileShare]::Read))
  try {
    $accepted = 0
    $candidates = 0
    while ($accepted -lt $samplePerRaster -and $candidates -lt $maxCandidates) {
      $candidates++
      $idx = $rng.Next(0, $pixels)
      $reader.BaseStream.Seek(($band563 * $bandBytes) + ($idx * 4), [System.IO.SeekOrigin]::Begin) | Out-Null
      $v563 = [double]$reader.ReadSingle()
      if ([double]::IsNaN($v563) -or [double]::IsInfinity($v563) -or $v563 -le $shadowThreshold) {
        continue
      }

      $vals = New-Object double[] $bands
      $ok = $true
      for ($b = 0; $b -lt $bands; $b++) {
        $reader.BaseStream.Seek(($b * $bandBytes) + ($idx * 4), [System.IO.SeekOrigin]::Begin) | Out-Null
        $v = [double]$reader.ReadSingle()
        if ([double]::IsNaN($v) -or [double]::IsInfinity($v) -or $v -le 0) {
          $ok = $false
          break
        }
        $vals[$b] = $v
      }
      if (-not $ok) {
        continue
      }

      for ($b = 0; $b -lt $bands; $b++) {
        $sums[$b] += $vals[$b]
      }
      $accepted++
      $count++
    }
  } finally {
    $reader.Close()
  }
}

$rows = for ($b = 0; $b -lt $bands; $b++) {
  [pscustomobject]@{
    wavelength_nm = $wavelengths[$b]
    mean_reflectance = if ($count -gt 0) { $sums[$b] / $count } else { [double]::NaN }
    sampled_pixels = $count
    source = "50m_smooth_5nm_sunlit_sample"
  }
}

$rows | Export-Csv -Path $outPath -NoTypeInformation
Write-Output $outPath
