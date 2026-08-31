Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$tableDir = Join-Path $root "reports/tables/pca_loading_spectral_regions"
$outDir = Join-Path $root "Documents/Tables and Figures"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

Add-Type -AssemblyName System.Drawing

$fontFamily = "Arial"
$ink = [System.Drawing.Color]::FromArgb(22, 26, 32)
$muted = [System.Drawing.Color]::FromArgb(84, 91, 101)
$grid = [System.Drawing.Color]::FromArgb(218, 223, 230)
$pc1Color = [System.Drawing.Color]::FromArgb(24, 89, 153)
$pc2Color = [System.Drawing.Color]::FromArgb(185, 79, 43)
$shadePc1 = [System.Drawing.Color]::FromArgb(46, 24, 89, 153)
$shadePc2 = [System.Drawing.Color]::FromArgb(48, 185, 79, 43)
$associationMode = $env:PCA_BRIGHTNESS_ASSOCIATION_MODE
if ($associationMode -eq $null -or $associationMode -eq "") {
  $associationMode = "least"
}

function New-Font($size, $style = [System.Drawing.FontStyle]::Regular) {
  return [System.Drawing.Font]::new($fontFamily, [single]$size, $style)
}

function New-Brush($color) {
  return [System.Drawing.SolidBrush]::new($color)
}

function New-Pen($color, $width = 1) {
  return [System.Drawing.Pen]::new($color, [single]$width)
}

function Draw-Text($g, $text, $x, $y, $size = 22, $color = $ink, $style = [System.Drawing.FontStyle]::Regular) {
  $font = New-Font $size $style
  $brush = New-Brush $color
  $g.DrawString($text, $font, $brush, [single]$x, [single]$y)
  $font.Dispose()
  $brush.Dispose()
}

function Draw-Centered($g, $text, $x, $y, $w, $h, $size = 22, $color = $ink, $style = [System.Drawing.FontStyle]::Regular) {
  $font = New-Font $size $style
  $brush = New-Brush $color
  $fmt = [System.Drawing.StringFormat]::new()
  $fmt.Alignment = [System.Drawing.StringAlignment]::Center
  $fmt.LineAlignment = [System.Drawing.StringAlignment]::Center
  $g.DrawString($text, $font, $brush, [System.Drawing.RectangleF]::new($x, $y, $w, $h), $fmt)
  $fmt.Dispose()
  $font.Dispose()
  $brush.Dispose()
}

function Scale-X($x, $xmin, $xmax, $left, $right) {
  return $left + (([double]$x - [double]$xmin) / ([double]$xmax - [double]$xmin)) * ($right - $left)
}

function Scale-Y($y, $ymin, $ymax, $top, $bottom) {
  return $bottom - (([double]$y - [double]$ymin) / ([double]$ymax - [double]$ymin)) * ($bottom - $top)
}

function Get-ZScores($values) {
  $nums = @($values | ForEach-Object { [double]$_ })
  $mean = ($nums | Measure-Object -Average).Average
  $ss = 0.0
  foreach ($v in $nums) { $ss += ([double]$v - $mean) * ([double]$v - $mean) }
  $sd = [Math]::Sqrt($ss / ($nums.Count - 1))
  return @($nums | ForEach-Object { ([double]$_ - $mean) / $sd })
}

function Get-Wavelength-Color($nm) {
  $w = [double]$nm
  if ($w -lt 450) { return [System.Drawing.Color]::FromArgb(92, 64, 172) }
  if ($w -lt 495) { return [System.Drawing.Color]::FromArgb(34, 96, 190) }
  if ($w -lt 570) { return [System.Drawing.Color]::FromArgb(36, 145, 77) }
  if ($w -lt 590) { return [System.Drawing.Color]::FromArgb(218, 181, 47) }
  if ($w -lt 620) { return [System.Drawing.Color]::FromArgb(226, 134, 39) }
  if ($w -lt 750) { return [System.Drawing.Color]::FromArgb(185, 41, 38) }
  return [System.Drawing.Color]::FromArgb(112, 22, 22)
}

function Add-LoadingZ($rows, $reflectanceByWavelength) {
  $pc1Z = Get-ZScores @($rows | ForEach-Object { $_.pc1_loading })
  $pc2Z = Get-ZScores @($rows | ForEach-Object { $_.pc2_loading })
  for ($i = 0; $i -lt $rows.Count; $i++) {
    $nmKey = [string]([int][double]$rows[$i].wavelength_nm)
    $rows[$i] | Add-Member -NotePropertyName pc1_loading_z -NotePropertyValue $pc1Z[$i] -Force
    $rows[$i] | Add-Member -NotePropertyName pc2_loading_z_recomputed -NotePropertyValue $pc2Z[$i] -Force
    $rows[$i] | Add-Member -NotePropertyName foliage_reflectance -NotePropertyValue ([double]$reflectanceByWavelength[$nmKey]) -Force
  }
  return $rows
}

function Draw-Series($g, $rows, $field, $xmin, $xmax, $ymin, $ymax, $left, $top, $right, $bottom, $color, $width = 4) {
  $pen = New-Pen $color $width
  for ($i = 1; $i -lt $rows.Count; $i++) {
    $x1 = Scale-X ([double]$rows[$i - 1].wavelength_nm) $xmin $xmax $left $right
    $y1 = Scale-Y ([double]$rows[$i - 1].$field) $ymin $ymax $top $bottom
    $x2 = Scale-X ([double]$rows[$i].wavelength_nm) $xmin $xmax $left $right
    $y2 = Scale-Y ([double]$rows[$i].$field) $ymin $ymax $top $bottom
    $g.DrawLine($pen, [single]$x1, [single]$y1, [single]$x2, [single]$y2)
  }
  $pen.Dispose()
}

function Draw-Reflectance-Series($g, $rows, $xmin, $xmax, $ymin, $ymax, $left, $top, $right, $bottom) {
  $vals = @($rows | ForEach-Object { [double]$_.foliage_reflectance })
  $min = ($vals | Measure-Object -Minimum).Minimum
  $max = ($vals | Measure-Object -Maximum).Maximum
  $scaledValues = @($vals | ForEach-Object { $ymin + ((([double]$_ - $min) / ($max - $min)) * ($ymax - $ymin)) })
  for ($i = 1; $i -lt $rows.Count; $i++) {
    $nm = ([double]$rows[$i - 1].wavelength_nm + [double]$rows[$i].wavelength_nm) / 2.0
    $pen = New-Pen (Get-Wavelength-Color $nm) 8
    $x1 = Scale-X ([double]$rows[$i - 1].wavelength_nm) $xmin $xmax $left $right
    $x2 = Scale-X ([double]$rows[$i].wavelength_nm) $xmin $xmax $left $right
    $y1 = Scale-Y $scaledValues[$i - 1] $ymin $ymax $top $bottom
    $y2 = Scale-Y $scaledValues[$i] $ymin $ymax $top $bottom
    $g.DrawLine($pen, [single]$x1, [single]$y1, [single]$x2, [single]$y2)
    $pen.Dispose()
  }
}

function Draw-Shaded-Region($g, $startNm, $endNm, $color, $xmin, $xmax, $left, $top, $right, $bottom) {
  $brush = New-Brush $color
  $x1 = Scale-X $startNm $xmin $xmax $left $right
  $x2 = Scale-X $endNm $xmin $xmax $left $right
  $g.FillRectangle($brush, [single]$x1, [single]$top, [single]($x2 - $x1), [single]($bottom - $top))
  $brush.Dispose()
}

function Get-ContiguousRegions($wavelengths) {
  $sorted = @($wavelengths | Sort-Object)
  $regions = @()
  $start = $null
  $previous = $null
  foreach ($nm in $sorted) {
    if ($start -eq $null) {
      $start = [int]$nm
      $previous = [int]$nm
    } elseif ([int]$nm -eq ($previous + 5)) {
      $previous = [int]$nm
    } else {
      $regions += @{ Start = $start; End = $previous; Label = ("{0}-{1} nm" -f $start, $previous) }
      $start = [int]$nm
      $previous = [int]$nm
    }
  }
  if ($start -ne $null) {
    $regions += @{ Start = $start; End = $previous; Label = ("{0}-{1} nm" -f $start, $previous) }
  }
  return @($regions)
}

function Get-TopLoadingRegions($rows, $field, $count) {
  if ($count -le 0) {
    return @()
  }
  $items = @($rows | Sort-Object -Property @{ Expression = { [Math]::Abs([double]$_.$field) }; Descending = $true } |
    Select-Object -First $count |
    ForEach-Object { [int][double]$_.wavelength_nm })
  return @(Get-ContiguousRegions $items)
}

function Draw-Panel($g, $rows, $title, $pc1Region, $pc2Region, $associationLabel, $left, $top, $right, $bottom) {
  $xmin = 398.0
  $xmax = 998.0
  $allY = @()
  foreach ($row in $rows) {
    $allY += [double]$row.pc1_loading_z
    $allY += [double]$row.pc2_loading_z_recomputed
  }
  $absMax = [Math]::Ceiling((($allY | ForEach-Object { [Math]::Abs($_) }) | Measure-Object -Maximum).Maximum)
  $ymin = -1 * [Math]::Max(2, $absMax)
  $ymax = [Math]::Max(2, $absMax)

  foreach ($region in @($pc1Region)) {
    Draw-Shaded-Region $g $region.Start $region.End $shadePc1 $xmin $xmax $left $top $right $bottom
  }
  foreach ($region in @($pc2Region)) {
    Draw-Shaded-Region $g $region.Start $region.End $shadePc2 $xmin $xmax $left $top $right $bottom
  }

  $gridPen = New-Pen $grid 1
  $axisPen = New-Pen $ink 2
  $zeroPen = New-Pen ([System.Drawing.Color]::FromArgb(145, 88, 96, 105)) 2
  $zeroPen.DashStyle = [System.Drawing.Drawing2D.DashStyle]::Dash

  foreach ($xTick in @(400, 500, 600, 700, 800, 900, 1000)) {
    $x = Scale-X $xTick $xmin $xmax $left $right
    $g.DrawLine($gridPen, [single]$x, [single]$top, [single]$x, [single]$bottom)
    Draw-Centered $g ("{0}" -f $xTick) ($x - 35) ($bottom + 10) 70 26 16 $muted
  }
  for ($yTick = $ymin; $yTick -le $ymax; $yTick += 2) {
    $y = Scale-Y $yTick $ymin $ymax $top $bottom
    $g.DrawLine($gridPen, [single]$left, [single]$y, [single]$right, [single]$y)
    Draw-Centered $g ("{0:0}" -f $yTick) ($left - 76) ($y - 13) 62 26 16 $muted
  }
  $zeroY = Scale-Y 0 $ymin $ymax $top $bottom
  $g.DrawLine($zeroPen, [single]$left, [single]$zeroY, [single]$right, [single]$zeroY)

  Draw-Reflectance-Series $g $rows $xmin $xmax $ymin $ymax $left $top $right $bottom
  Draw-Series $g $rows "pc1_loading_z" $xmin $xmax $ymin $ymax $left $top $right $bottom $pc1Color 4.2
  Draw-Series $g $rows "pc2_loading_z_recomputed" $xmin $xmax $ymin $ymax $left $top $right $bottom $pc2Color 4.2

  $g.DrawLine($axisPen, [single]$left, [single]$bottom, [single]$right, [single]$bottom)
  $g.DrawLine($axisPen, [single]$left, [single]$top, [single]$left, [single]$bottom)

  Draw-Text $g $title ($left + 8) ($top - 48) 28 $ink ([System.Drawing.FontStyle]::Bold)
  $gridPen.Dispose()
  $axisPen.Dispose()
  $zeroPen.Dispose()
}

$loadings = @(Import-Csv -Path (Join-Path $tableDir "pca_pc1_pc2_loadings_by_wavelength.csv"))
$reflectanceRows = @(Import-Csv -Path (Join-Path $tableDir "sampled_plot_foliage_reflectance_50m.csv"))
$reflectanceByWavelength = @{}
foreach ($row in $reflectanceRows) {
  $reflectanceByWavelength[[string]([int][double]$row.wavelength_nm)] = [double]$row.mean_reflectance
}

$rawRows = Add-LoadingZ @($loadings | Where-Object { $_.pca_basis -eq "regular_PCA" } | Sort-Object {[double]$_.wavelength_nm}) $reflectanceByWavelength
$stdRows = Add-LoadingZ @($loadings | Where-Object { $_.pca_basis -eq "standardized_PCA" } | Sort-Object {[double]$_.wavelength_nm}) $reflectanceByWavelength

$brightnessCorrelationByPc = @{}
foreach ($path in @(
  (Join-Path $root "reports/tables/pc1_mean_reflectance_correlation/50m_regular_PCA_pc_mean_reflectance_correlation_summary.csv"),
  (Join-Path $root "reports/tables/pc1_mean_reflectance_correlation/50m_standardized_PCA_pc_mean_reflectance_correlation_summary.csv")
)) {
  if (Test-Path $path) {
    foreach ($row in @(Import-Csv -Path $path | Where-Object { $_.analysis_level -eq "quadrat_mean_level" })) {
      $brightnessCorrelationByPc[("{0}|{1}" -f $row.pca_basis, $row.pc_axis)] = [Math]::Abs([double]$row.pearson_r)
    }
  }
}

$width = 2450
$height = 1850
$bmp = [System.Drawing.Bitmap]::new($width, $height, [System.Drawing.Imaging.PixelFormat]::Format32bppArgb)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$g.Clear([System.Drawing.Color]::Transparent)
$g.SmoothingMode = [System.Drawing.Drawing2D.SmoothingMode]::AntiAlias
$g.TextRenderingHint = [System.Drawing.Text.TextRenderingHint]::ClearTypeGridFit

Draw-Text $g "PCA loading z-scores across wavelength" 145 80 42 $ink ([System.Drawing.FontStyle]::Bold)
Draw-Text $g "PC1 and PC2 loading structure with a scaled plot-level foliage reflectance reference" 145 140 24 $muted

$left = 250
$right = 2275
$top1 = 330
$bottom1 = 840
$top2 = 1120
$bottom2 = 1630

if ($associationMode -eq "none") {
  $rawPc1 = @()
  $rawPc2 = @()
  $stdPc1 = @()
  $stdPc2 = @()
  $associationLabel = ""
} elseif ($associationMode -eq "most") {
  $referenceBrightnessR = [double]$brightnessCorrelationByPc["regular_PCA|PC1"]
  $maxHighlightedBands = 25
  $rawPc1Count = [Math]::Max(1, [int][Math]::Round($maxHighlightedBands * ([double]$brightnessCorrelationByPc["regular_PCA|PC1"] / $referenceBrightnessR)))
  $rawPc2Count = [Math]::Max(1, [int][Math]::Round($maxHighlightedBands * ([double]$brightnessCorrelationByPc["regular_PCA|PC2"] / $referenceBrightnessR)))
  $stdPc1Count = [Math]::Max(1, [int][Math]::Round($maxHighlightedBands * ([double]$brightnessCorrelationByPc["standardized_PCA|PC1"] / $referenceBrightnessR)))
  $stdPc2Count = [Math]::Max(1, [int][Math]::Round($maxHighlightedBands * ([double]$brightnessCorrelationByPc["standardized_PCA|PC2"] / $referenceBrightnessR)))

  $rawPc1 = @(Get-TopLoadingRegions $rawRows "pc1_loading_z" $rawPc1Count)
  $rawPc2 = @(Get-TopLoadingRegions $rawRows "pc2_loading_z_recomputed" $rawPc2Count)
  $stdPc1 = @(Get-TopLoadingRegions $stdRows "pc1_loading_z" $stdPc1Count)
  $stdPc2 = @(Get-TopLoadingRegions $stdRows "pc2_loading_z_recomputed" $stdPc2Count)
  $associationLabel = "brightness-associated region"
} else {
  $rawPc1 = @(@{ Start = 933; End = 998; Label = "933-998 nm" })
  $rawPc2 = @(@{ Start = 643; End = 683; Label = "643-683 nm" })
  $stdPc1 = @(@{ Start = 843; End = 903; Label = "843-903 nm" })
  $stdPc2 = @(@{ Start = 753; End = 778; Label = "753-778 nm" })
  $associationLabel = "least associated with brightness"
}

Draw-Panel $g $rawRows "Raw PCA" $rawPc1 $rawPc2 $associationLabel $left $top1 $right $bottom1
Draw-Panel $g $stdRows "Vector-normalized PCA" $stdPc1 $stdPc2 $associationLabel $left $top2 $right $bottom2

Draw-Centered $g "Wavelength (nm)" 775 1740 900 44 27 $ink ([System.Drawing.FontStyle]::Bold)

$font = New-Font 27 ([System.Drawing.FontStyle]::Bold)
$brush = New-Brush $ink
$format = [System.Drawing.StringFormat]::new()
$format.Alignment = [System.Drawing.StringAlignment]::Center
$format.LineAlignment = [System.Drawing.StringAlignment]::Center
$format.FormatFlags = [System.Drawing.StringFormatFlags]::DirectionVertical
$g.DrawString("Loading z-score", $font, $brush, [System.Drawing.RectangleF]::new(55, 680, 54, 600), $format)
$format.Dispose()
$font.Dispose()
$brush.Dispose()

$legendX = 1410
$legendY = 95
$pc1Pen = New-Pen $pc1Color 5
$pc2Pen = New-Pen $pc2Color 5
$reflectPen = New-Pen ([System.Drawing.Color]::FromArgb(112, 22, 22)) 8
$pc1Brush = New-Brush $shadePc1
$pc2Brush = New-Brush $shadePc2

$g.DrawLine($pc1Pen, [single]$legendX, [single]$legendY, [single]($legendX + 90), [single]$legendY)
Draw-Text $g "PC1 loading z-score" ($legendX + 112) ($legendY - 20) 20 $ink
$g.DrawLine($pc2Pen, [single]$legendX, [single]($legendY + 45), [single]($legendX + 90), [single]($legendY + 45))
Draw-Text $g "PC2 loading z-score" ($legendX + 112) ($legendY + 25) 20 $ink
$g.DrawLine($reflectPen, [single]$legendX, [single]($legendY + 90), [single]($legendX + 90), [single]($legendY + 90))
Draw-Text $g "Plot foliage reflectance (scaled)" ($legendX + 112) ($legendY + 70) 20 $ink
if ($associationMode -ne "none") {
  $g.FillRectangle($pc1Brush, [single]$legendX, [single]($legendY + 132), 90, 24)
  Draw-Text $g ("PC1 " + $associationLabel) ($legendX + 112) ($legendY + 122) 20 $ink
  $g.FillRectangle($pc2Brush, [single]$legendX, [single]($legendY + 174), 90, 24)
  Draw-Text $g ("PC2 " + $associationLabel) ($legendX + 112) ($legendY + 164) 20 $ink
}

$pc1Pen.Dispose()
$pc2Pen.Dispose()
$reflectPen.Dispose()
$pc1Brush.Dispose()
$pc2Brush.Dispose()

if ($associationMode -eq "none") {
  $outPath = Join-Path $outDir "15_pca_loading_zscores_no_shaded_regions.png"
} elseif ($associationMode -eq "most") {
  $outPath = Join-Path $outDir "14_pca_loading_zscores_most_brightness_associated_regions.png"
} else {
  $outPath = Join-Path $outDir "13_pca_loading_departures_across_wavelength.png"
}
$bmp.Save($outPath, [System.Drawing.Imaging.ImageFormat]::Png)
$g.Dispose()
$bmp.Dispose()

Write-Output $outPath
